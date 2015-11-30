import logging
import itertools
try:
    logger = logging.getLogger("peak_grouping")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

from sqlalchemy.ext.baked import bakery
from sqlalchemy import func, bindparam, select

import numpy as np

from sklearn.linear_model import LogisticRegression

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, HypothesisSampleMatch, Decon2LSPeakGroup,
    PeakGroupDatabase, PeakGroupMatch, TempPeakGroupMatch, JointPeakGroupMatch,
    PeakGroupMatchToJointPeakGroupMatch, PeakGroupScoringModel)

from glycresoft_sqlalchemy.utils.collectiontools import flatten


from .common import (
    ppm_error, centroid_scan_error_regression, expanding_window,
    expected_a_peak_regression)

T_TempPeakGroupMatch = TempPeakGroupMatch.__table__
TPeakGroupMatch = PeakGroupMatch.__table__
T_JointPeakGroupMatch = JointPeakGroupMatch.__table__


query_oven = bakery()


ClassifierType = LogisticRegression


def _sklearn_classifier_to_scoring_model(classifier):
    return PeakGroupScoringModel.from_parameter_vector(classifier.coef_)


def _sklearn_classifier_from_scoring_model(model, *args, **kwargs):
    classifier = ClassifierType(*args, **kwargs)
    classifier.coef_ = model.to_parameter_vector()
    return classifier


def _group_unmatched_peak_groups_by_shifts(groups, mass_shift_map, grouping_error_tolerance=2e-5):

    mass_shift_range = []
    for mass_shift, count_range in mass_shift_map.items():
        shift = mass_shift.mass
        for shift_count in range(1, count_range + 1):
            mass_shift_range.append((shift * shift_count))

    try:
        lower_edge = min(mass_shift_range)
        upper_edge = max(mass_shift_range)
    except ValueError:
        lower_edge = 0
        upper_edge = 0

    while(len(groups) > 0):
        current_group = groups.pop(0)
        max_mass = (current_group.weighted_monoisotopic_mass + upper_edge) * (1 + grouping_error_tolerance)
        min_mass = (current_group.weighted_monoisotopic_mass + lower_edge) * (1 - grouping_error_tolerance)
        current_group_accumulator = []
        remaining_accumulator = []

        while(len(groups) > 0):
            group = groups.pop(0)
            if group.weighted_monoisotopic_mass > max_mass:
                remaining_accumulator.append(group)
                continue
            elif min_mass <= group.weighted_monoisotopic_mass <= max_mass:
                matched = False
                for mass_shift, count_range in mass_shift_map.items():
                    shift = mass_shift.mass
                    for shift_count in range(1, count_range + 1):
                        current_shift = (shift * shift_count)
                        error = ppm_error(
                            current_group.weighted_monoisotopic_mass + current_shift,
                            group.weighted_monoisotopic_mass)
                        if abs(error) <= grouping_error_tolerance:
                            matched = True
                            group.mass_shift = mass_shift
                            group.mass_shift_count = shift_count
                            current_group_accumulator.append(group)
                            break
                if not matched:
                    remaining_accumulator.append(group)
            elif group.weighted_monoisotopic_mass < min_mass:
                remaining_accumulator.append(group)
                break

        current_group_accumulator.append(current_group)
        yield (current_group_accumulator)
        remaining_accumulator.extend(groups)
        groups = remaining_accumulator


def _get_groups_by_composition_ids(session, hypothesis_sample_match_id):
    composition_id_q = session.query(PeakGroupMatch.theoretical_match_id).filter(
        PeakGroupMatch.hypothesis_sample_match_id).group_by(
        PeakGroupMatch.theoretical_match_id).order_by(
        PeakGroupMatch.weighted_monoisotopic_mass.desc())
    for _composition_id in composition_id_q:
        _composition_id = _composition_id[0]
        if _composition_id is None:
            continue
        yield session.query(PeakGroupMatch).filter(
            PeakGroupMatch.theoretical_match_id == _composition_id,
            PeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id).all()


def _merge_groups(group_matches, minimum_abundance_ratio=0.01):
    scan_count_total = 0
    min_scan = float("inf")
    max_scan = 0
    charge_states = set()
    scan_times = set()
    average_a_to_a_plus_2_ratio = 0
    total_volume = 0
    average_signal_to_noise = 0
    n = 0.
    n_modification_states = 0
    merged_peak_data = {
        "peak_ids": [],
        "intensities": [],
        "scan_times": []
    }
    maximum_volume = max(g.total_volume for g in group_matches)
    minimum_abundance = minimum_abundance_ratio * maximum_volume

    for peak_group in group_matches:
        if peak_group.total_volume < minimum_abundance:
            continue
        peak_data = peak_group.peak_data
        merged_peak_data['peak_ids'].extend(peak_data['peak_ids'])
        merged_peak_data['intensities'].extend(peak_data['intensities'])
        merged_peak_data['scan_times'].extend(peak_data['scan_times'])

        n_peaks = len(peak_data['peak_ids'])
        n += n_peaks
        scan_count_total += peak_group.scan_count
        total_volume += peak_group.total_volume
        n_modification_states += 1

        average_a_to_a_plus_2_ratio += peak_group.average_a_to_a_plus_2_ratio * n_peaks
        average_signal_to_noise += peak_group.average_signal_to_noise * n_peaks

        min_scan = min(min_scan, peak_group.first_scan_id)
        max_scan = max(max_scan, peak_group.last_scan_id)

        scan_times.update(peak_data['scan_times'])
        charge_states.update(peak_data.get("charge_states", ()))

    average_signal_to_noise /= n
    average_a_to_a_plus_2_ratio /= n
    if len(charge_states) != 0:
        charge_state_count = len(charge_states)
    else:
        charge_state_count = max(g.charge_state_count for g in group_matches)

    scan_times = sorted(scan_times)
    windows = expanding_window(scan_times)
    window_densities = []
    for window in windows:
        window_max_scan = window[-1]
        window_min_scan = window[0]
        window_scan_count = len(window)
        window_scan_density = window_scan_count / (
            float(window_max_scan - window_min_scan) + 15.) if window_scan_count > 1 else 0
        if window_scan_density != 0:
            window_densities.append(window_scan_density)
    if len(window_densities) != 0:
        scan_density = max(window_densities)  # sum(window_densities) / float(len(window_densities))
    else:
        scan_density = 0.

    instance_dict = {
        "first_scan_id": min_scan,
        "last_scan_id": max_scan,
        "scan_density": scan_density,
        "centroid_scan_estimate": sum(scan_times) / n,
        "average_a_to_a_plus_2_ratio": average_a_to_a_plus_2_ratio,
        "average_signal_to_noise": average_signal_to_noise,
        "charge_state_count": charge_state_count,
        "modification_state_count": n_modification_states,
        "total_volume": total_volume,
        "scan_count": scan_count_total,
        "peak_data": merged_peak_data,
        "fingerprint": ':'.join(map(str, scan_times)),
        "weighted_monoisotopic_mass": group_matches[0].weighted_monoisotopic_mass,
        "hypothesis_sample_match_id": group_matches[0].hypothesis_sample_match_id,
        "theoretical_match_id": group_matches[0].theoretical_match_id
    }
    return instance_dict, [p.id for p in group_matches]


def join_unmatched(session, hypothesis_sample_match_id, grouping_error_tolerance=2e-5, minimum_abundance_ratio=0.01):
    unmatched = session.query(PeakGroupMatch).filter(
        PeakGroupMatch.theoretical_match_id == None,
        PeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id).order_by(
        PeakGroupMatch.weighted_monoisotopic_mass.desc()).all()
    if len(unmatched) == 0:
        return

    mass_shift_map = unmatched[0].hypothesis_sample_match.parameters['mass_shift_map']
    conn = session.connection()
    for bunch in _group_unmatched_peak_groups_by_shifts(unmatched, mass_shift_map):
        group, member_ids = _merge_groups(bunch, minimum_abundance_ratio)
        group['matched'] = False
        group['theoretical_match_id'] = None
        joint_id = conn.execute(T_JointPeakGroupMatch.insert(), group).lastrowid
        conn.execute(PeakGroupMatchToJointPeakGroupMatch.insert(),
                     [{"peak_group_id": i, "joint_group_id": joint_id} for i in member_ids])


def join_matched(session, hypothesis_sample_match_id, minimum_abundance_ratio=0.01):
    conn = session.connection()
    for bunch in _get_groups_by_composition_ids(session, hypothesis_sample_match_id):
        group, member_ids = _merge_groups(bunch, minimum_abundance_ratio)
        group['matched'] = True
        group['theoretical_match_id'] = bunch[0].theoretical_match_id
        group['theoretical_match_type'] = bunch[0].theoretical_match_type
        joint_id = conn.execute(T_JointPeakGroupMatch.insert(), group).lastrowid
        conn.execute(PeakGroupMatchToJointPeakGroupMatch.insert(),
                     [{"peak_group_id": i, "joint_group_id": joint_id} for i in member_ids])


def estimate_trends(session, hypothesis_sample_match_id):
    '''
    After assigning peak group features, impute the global
    trend for peak and scan shapes
    '''
    logger.info("Estimating peak trends")

    conn = session.connect()

    cen_alpha, cen_beta = centroid_scan_error_regression(
        session, source_model=JointPeakGroupMatch, filter_fn=lambda q: q.filter(
            JointPeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id))

    expected_a_alpha, expected_a_beta = expected_a_peak_regression(
        session, source_model=JointPeakGroupMatch, filter_fn=lambda q: q.filter(
            JointPeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id))

    update_expr = T_JointPeakGroupMatch.update().values(
        centroid_scan_error=func.abs(
            T_JointPeakGroupMatch.c.centroid_scan_estimate - (
                cen_alpha + cen_beta * T_JointPeakGroupMatch.c.weighted_monoisotopic_mass)),
        a_peak_intensity_error=func.abs(
            T_JointPeakGroupMatch.c.average_a_to_a_plus_2_ratio - (
                expected_a_alpha + expected_a_beta * T_JointPeakGroupMatch.c.weighted_monoisotopic_mass))
            ).where(
        T_JointPeakGroupMatch.c.hypothesis_sample_match_id == hypothesis_sample_match_id)

    max_weight = conn.execute(select([func.max(T_JointPeakGroupMatch.c.weighted_monoisotopic_mass)])).scalar()
    slices = [0] + [max_weight * float(i)/10. for i in range(1, 11)]
    for i in range(1, len(slices)):
        lower = slices[i - 1]
        upper = slices[i]
        logger.info("Updating slice %f-%f", lower, upper)
        step = update_expr.where(
            T_JointPeakGroupMatch.c.weighted_monoisotopic_mass.between(
                lower, upper))
        conn.execute(step)
        session.commit()
        conn = session.connection()

    conn = session.connection()
    lower = slices[len(slices) - 1]
    step = update_expr.where(
        T_JointPeakGroupMatch.c.weighted_monoisotopic_mass >= lower)
    conn.execute(step)
    session.commit()


def logit(X):
    return np.log(X) - np.log(1 - X)


def inverse_logit(X):
    return np.exp(X) / (np.exp(X) + 1)


class PeakGroupMassShiftJoiningClassifier(PipelineModule):
    features = [
        T_JointPeakGroupMatch.c.charge_state_count,
        T_JointPeakGroupMatch.c.scan_density,
        T_JointPeakGroupMatch.c.modification_state_count,
        T_JointPeakGroupMatch.c.total_volume,
        T_JointPeakGroupMatch.c.a_peak_intensity_error,
        T_JointPeakGroupMatch.c.centroid_scan_error,
        T_JointPeakGroupMatch.c.scan_count,
        T_JointPeakGroupMatch.c.average_signal_to_noise
    ]

    label = [T_JointPeakGroupMatch.c.matched]

    ids = [T_JointPeakGroupMatch.c.id]

    classifier = None

    def __init__(
            self, database_path, observed_ions_path,
            sample_run_id=None, hypothesis_sample_match_id=None,
            search_type="TheoreticalGlycanComposition",
            match_tolerance=2e-5, minimum_abundance_ratio=0.01,
            use_legacy_coefficients=True,
            n_processes=4):
        self.manager = self.manager_type(database_path)
        self.lcms_database = self.manager_type(observed_ions_path)
        session = self.manager.session()
        self.sample_run_id = sample_run_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        hypothesis_sample_match = session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
        self.mass_shift_map = hypothesis_sample_match.parameters['mass_shift_map']
        self.match_tolerance = match_tolerance
        self.minimum_abundance_ratio = minimum_abundance_ratio
        self.use_legacy_coefficients = use_legacy_coefficients
        self.n_processes = n_processes

    def create_joins(self):
        session = self.manager.session()
        self.inform("Joining Unmatched Peak Groups")
        join_unmatched(session, self.hypothesis_sample_match_id,
                       self.match_tolerance, self.minimum_abundance_ratio)
        self.inform("Joining Matched Peak Groups")
        join_matched(session, self.hypothesis_sample_match_id, self.minimum_abundance_ratio)
        session.commit()
        session.close()

    def estimate_trends(self):
        '''
        After assigning peak group features, impute the global
        trend for peak and scan shapes
        '''
        session = self.manager.session()
        hypothesis_sample_match_id = self.hypothesis_sample_match_id
        logger.info("Estimating peak trends")

        conn = session.connection()

        cen_alpha, cen_beta = centroid_scan_error_regression(
            session, source_model=JointPeakGroupMatch, filter_fn=lambda q: q.filter(
                JointPeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id))

        expected_a_alpha, expected_a_beta = expected_a_peak_regression(
            session, source_model=JointPeakGroupMatch, filter_fn=lambda q: q.filter(
                JointPeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id))

        update_expr = T_JointPeakGroupMatch.update().values(
            centroid_scan_error=func.abs(
                T_JointPeakGroupMatch.c.centroid_scan_estimate - (
                    cen_alpha + cen_beta * T_JointPeakGroupMatch.c.weighted_monoisotopic_mass)),
            a_peak_intensity_error=func.abs(
                T_JointPeakGroupMatch.c.average_a_to_a_plus_2_ratio - (
                    expected_a_alpha + expected_a_beta * T_JointPeakGroupMatch.c.weighted_monoisotopic_mass))
                ).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == hypothesis_sample_match_id)

        max_weight = conn.execute(select([func.max(T_JointPeakGroupMatch.c.weighted_monoisotopic_mass)])).scalar()
        slices = [0] + [max_weight * float(i)/10. for i in range(1, 11)]
        for i in range(1, len(slices)):
            lower = slices[i - 1]
            upper = slices[i]
            logger.info("Updating slice %f-%f", lower, upper)
            step = update_expr.where(
                T_JointPeakGroupMatch.c.weighted_monoisotopic_mass.between(
                    lower, upper))
            conn.execute(step)
            session.commit()
            conn = session.connection()

        conn = session.connection()
        lower = slices[len(slices) - 1]
        step = update_expr.where(
            T_JointPeakGroupMatch.c.centroid_scan_error == None)
        conn.execute(step)
        session.commit()

    def construct_model_matrix(self):
        features = self.features
        ids = self.ids
        label = self.label

        data_model_session = self.manager.session()
        conn = data_model_session.connection()

        id_vec = flatten(conn.execute(select(ids).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall())
        feature_matrix = np.array(conn.execute(select(features).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall(),
            dtype=np.float64)
        label_vector = np.array(conn.execute(select(label).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall())

        return [label_vector, id_vec, feature_matrix]

    def patch_coefficient_sign(self):
        betas = self.classifier.coef_
        i_total_volume = [i for i, col in enumerate(self.features) if col.name == 'total_volume'][0]
        if betas[0][i_total_volume] < 0:
            betas[0][i_total_volume] *= -1.5
        self.classifier.coef_ = betas

    def legacy_model_coefficient(self):
        betas = [[0.65784006, 1.345376317, 0.219899787,
                  0.0000383503, -0.000349839, -0.001610525,
                  -0.000947516, 0.011453828]]
        betas = np.array(betas)
        self.classifier.coef_ = betas

    def fit_and_score(self):
        classifier = ClassifierType()
        label_vector, id_vec, feature_matrix = self.construct_model_matrix()
        classifier.fit(feature_matrix, label_vector.ravel())
        self.classifier = classifier
        if self.use_legacy_coefficients:
            self.legacy_model_coefficient()
        else:
            self.patch_coefficient_sign()
        scores = classifier.predict_proba(feature_matrix)[:, 1]
        update_data = [{"b_id": id_key, "ms1_score": float(score)} for id_key, score in itertools.izip(id_vec, scores)]
        stmt = T_JointPeakGroupMatch.update().where(
            T_JointPeakGroupMatch.c.id == bindparam("b_id")).values(
            ms1_score=bindparam("ms1_score"))
        data_model_session = self.manager.session()
        data_model_session.execute(stmt, update_data)
        data_model_session.commit()

    def coefficients(self):
        if self.classifier is None:
            return {}
        return {f.name: v for f, v in zip(self.features, self.classifier.coef_[0])}

    def transfer_peak_groups(self):
        data_model_session = self.manager.session()
        lcms_database_session = self.lcms_database.session()
        peak_group_labels = [
            'id',
            'sample_run_id',
            'charge_state_count',
            'scan_count',
            'first_scan_id',
            'last_scan_id',
            'scan_density',
            'weighted_monoisotopic_mass',
            'total_volume',
            'average_a_to_a_plus_2_ratio',
            'a_peak_intensity_error',
            'centroid_scan_estimate',
            'centroid_scan_error',
            'average_signal_to_noise',
            'matched',
            "peak_data"
        ]

        stmt = lcms_database_session.query(
                *[getattr(Decon2LSPeakGroup, label) for label in peak_group_labels]).filter(
                Decon2LSPeakGroup.sample_run_id == self.sample_run_id)

        id_stmt = data_model_session.query(
            PeakGroupMatch.peak_group_id).filter(
            PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            PeakGroupMatch.matched).selectable

        batch = lcms_database_session.connection().execute(stmt.selectable)

        conn = data_model_session.connection()

        # Move all Decon2LSPeakGroups, regardless of whether or not they matched to
        # the temporary table.
        while True:
            items = batch.fetchmany(10000)
            if len(items) == 0:
                break
            mapped_items = [dict(zip(peak_group_labels, row)) for row in items]
            conn.execute(T_TempPeakGroupMatch.insert(), mapped_items)
            data_model_session.commit()
            conn = data_model_session.connection()

        data_model_session.commit()
        conn = data_model_session.connection()

        id_stmt = data_model_session.query(
            PeakGroupMatch.peak_group_id).filter(
            PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            PeakGroupMatch.matched).selectable

        if not len(data_model_session.connection().execute(id_stmt).fetchmany(2)) == 2:
            raise ValueError("Hypothesis-Sample Match ID matches maps no PeakGroupMatches")

        # Use the presence of a PeakGroupMatch.matched == True row to indicate whether something was a match
        # and fill out the :attr:`TempPeakGroupMatch.matched` column
        update_stmt = T_TempPeakGroupMatch.update().where(
            T_TempPeakGroupMatch.c.id.in_(id_stmt)).values(matched=True)
        data_model_session.connection().execute(update_stmt)

        update_stmt = T_TempPeakGroupMatch.update().where(
            ~T_TempPeakGroupMatch.c.id.in_(id_stmt)).values(matched=False)
        data_model_session.connection().execute(update_stmt)

        # Copy the Decon2LSPeakGroups that did not match anything as PeakGroupMatches
        # with null theoretical group matches
        move_stmt = data_model_session.query(
            *[getattr(TempPeakGroupMatch, label) for label in peak_group_labels]).filter(
            ~TempPeakGroupMatch.matched).selectable

        def transform(row):
            out = dict(zip(peak_group_labels, row))
            peak_group_match_id = out.pop('id')
            out["theoretical_match_type"] = None
            out['matched'] = False
            out['hypothesis_sample_match_id'] = self.hypothesis_sample_match_id
            out['peak_group_id'] = peak_group_match_id
            assert out['weighted_monoisotopic_mass'] is not None
            return out

        conn = data_model_session.connection()
        batch = conn.execute(move_stmt)
        while True:
            items = batch.fetchmany(10000)
            if len(items) == 0:
                break
            conn.execute(TPeakGroupMatch.insert(), list(map(transform, items)))
        data_model_session.commit()

    def clear_peak_groups(self):
        """Delete all TempPeakGroupMatch rows.
        """
        data_model_session = self.manager.session()
        data_model_session.query(TempPeakGroupMatch).delete()
        data_model_session.commit()
        data_model_session.close()

    def run(self):
        self.clear_peak_groups()
        self.transfer_peak_groups()
        self.create_joins()
        self.estimate_trends()
        self.fit_and_score()


class PeakGroupClassification(PipelineModule):
    features = [
        T_TempPeakGroupMatch.c.charge_state_count,
        T_TempPeakGroupMatch.c.scan_density,
        T_TempPeakGroupMatch.c.scan_count,
        T_TempPeakGroupMatch.c.total_volume,
        T_TempPeakGroupMatch.c.a_peak_intensity_error,
        T_TempPeakGroupMatch.c.centroid_scan_error,
        T_TempPeakGroupMatch.c.average_signal_to_noise
    ]

    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=None, hypothesis_sample_match_id=None,
                 model_parameters=None):
        self.database_manager = self.manager_type(database_path)
        self.lcms_database = PeakGroupDatabase(observed_ions_path)
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.model_parameters = model_parameters
        self.classifier = None

    def transfer_peak_groups(self):
        """Copy Decon2LSPeakGroup entries from :attr:`observed_ions_manager`
        into :attr:`database_manager` across database boundaries so that they
        may be involved in the same table queries.

        Replicate the TempPeakGroupMatch rows for groups that did not match any
        database entries as full PeakGroupMatch rows with :attr:`PeakGroupMatch.matched` == `False`
        for post-processing.
        """
        self.clear_peak_groups()
        data_model_session = self.database_manager.session()
        lcms_database_session = self.lcms_database.session()
        peak_group_labels = [
            'id',
            'sample_run_id',
            'charge_state_count',
            'scan_count',
            'first_scan_id',
            'last_scan_id',
            'scan_density',
            'weighted_monoisotopic_mass',
            'total_volume',
            'average_a_to_a_plus_2_ratio',
            'a_peak_intensity_error',
            'centroid_scan_estimate',
            'centroid_scan_error',
            'average_signal_to_noise',
            'matched',
            "peak_data"
        ]

        stmt = lcms_database_session.query(
                *[getattr(Decon2LSPeakGroup, label) for label in peak_group_labels]).filter(
                Decon2LSPeakGroup.sample_run_id == self.sample_run_id)
        batch = lcms_database_session.connection().execute(stmt.selectable)

        conn = data_model_session.connection()

        # Move all Decon2LSPeakGroups, regardless of whether or not they matched to
        # the temporary table.
        while True:
            items = batch.fetchmany(10000)
            if len(items) == 0:
                break
            mapped_items = [dict(zip(peak_group_labels, row)) for row in items]
            conn.execute(T_TempPeakGroupMatch.insert(), mapped_items)
            data_model_session.commit()
            conn = data_model_session.connection()

        data_model_session.commit()
        conn = data_model_session.connection()

        id_stmt = data_model_session.query(
            PeakGroupMatch.peak_group_id).filter(
            PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            PeakGroupMatch.matched).selectable

        if not len(data_model_session.connection().execute(id_stmt).fetchmany(2)) == 2:
            raise ValueError("Hypothesis-Sample Match ID matches maps no PeakGroupMatches")

        # Use the presence of a PeakGroupMatch.matched == True row to indicate whether something was a match
        # and fill out the :attr:`TempPeakGroupMatch.matched` column
        update_stmt = T_TempPeakGroupMatch.update().where(
            T_TempPeakGroupMatch.c.id.in_(id_stmt)).values(matched=True)
        data_model_session.connection().execute(update_stmt)

        update_stmt = T_TempPeakGroupMatch.update().where(
            ~T_TempPeakGroupMatch.c.id.in_(id_stmt)).values(matched=False)
        data_model_session.connection().execute(update_stmt)

        # Copy the Decon2LSPeakGroups that did not match anything as PeakGroupMatches
        # with null theoretical group matches
        move_stmt = data_model_session.query(
            *[getattr(TempPeakGroupMatch, label) for label in peak_group_labels]).filter(
            ~TempPeakGroupMatch.matched).selectable

        def transform(row):
            out = dict(zip(peak_group_labels, row))
            peak_group_match_id = out.pop('id')
            out["theoretical_match_type"] = None
            out['matched'] = False
            out['hypothesis_sample_match_id'] = self.hypothesis_sample_match_id
            out['peak_group_id'] = peak_group_match_id
            assert out['weighted_monoisotopic_mass'] is not None
            return out

        conn = data_model_session.connection()
        batch = conn.execute(move_stmt)
        while True:
            items = batch.fetchmany(10000)
            if len(items) == 0:
                break
            conn.execute(TPeakGroupMatch.insert(), list(map(transform, items)))
        data_model_session.commit()

    def fit_regression(self):
        """Fit the L2 Logistic Regression Model against the temporary peak group table.
        Computes scores for each TempPeakGroupMatch and maps them to the referent PeakGroupMatch.

        .. warning::
            The regression operation is carried out **in memory**, however space used is proportional
            to the number of |Decon2LSPeakGroup| records in this hypothesis-sample match, not the total
            number of matches + unmatched groups, which is almost certainly going to be much larger.

        Returns
        -------
        sklearn.linear_model.LogisticRegression : The fitted model
        """
        features = self.features
        label = [T_TempPeakGroupMatch.c.matched]
        ids = [T_TempPeakGroupMatch.c.id]

        data_model_session = self.database_manager.session()
        conn = data_model_session.connection()
        feature_matrix = np.array(conn.execute(select(ids + features)).fetchall(), dtype=np.float64)
        # Drop all rows containing nan
        mask = ~np.isnan(feature_matrix).any(axis=1)
        feature_matrix = feature_matrix[mask]
        label_vector = np.array(conn.execute(select(label)).fetchall())[mask]
        classifier = ClassifierType()
        if self.model_parameters is None:
            classifier.fit(feature_matrix[:, 1:], label_vector.ravel())
        else:
            classifier.coef_ = np.asarray(self.model_parameters)
        scores = classifier.predict_proba(feature_matrix[:, 1:])[:, 1]
        for group_id, score in itertools.izip(feature_matrix[:, 0], scores):
                conn.execute(
                    TPeakGroupMatch.update().where(
                        TPeakGroupMatch.c.peak_group_id == group_id).values(
                        ms1_score=float(score)))
        data_model_session.commit()
        return classifier

    def clear_peak_groups(self):
        """Delete all TempPeakGroupMatch rows.
        """
        data_model_session = self.database_manager.session()
        data_model_session.query(TempPeakGroupMatch).delete()
        data_model_session.commit()
        data_model_session.close()

    def run(self):
        self.transfer_peak_groups()
        self.classifier = self.fit_regression()
        logger.info("Classes: %r", self.classifier.classes_)
        logger.info("Coefficients: %r", ([c.name for c in self.features], self.classifier.coef_))
        self.clear_peak_groups()
