import csv
import os
import logging
import itertools
import functools
import multiprocessing
import operator
try:
    logger = logging.getLogger("peak_grouping")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

try:
    from cPickle import Pickler, Unpickler
except:
    from pickle import Pickler, Unpickler

from sqlalchemy.ext.baked import bakery
from sqlalchemy import func, bindparam, select

import numpy as np

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, HypothesisSampleMatch, Decon2LSPeakGroup, PipelineException,
    PeakGroupDatabase, PeakGroupMatch, TempPeakGroupMatch, JointPeakGroupMatch,
    PeakGroupMatchToJointPeakGroupMatch, PeakGroupScoringModel)

from glycresoft_sqlalchemy.scoring import logistic_scoring

from glycresoft_sqlalchemy.utils import get_scale
from glycresoft_sqlalchemy.utils.collectiontools import flatten


from .common import (
    ppm_error, centroid_scan_error_regression, expanding_window,
    expected_a_peak_regression)

T_TempPeakGroupMatch = TempPeakGroupMatch.__table__
TPeakGroupMatch = PeakGroupMatch.__table__
T_JointPeakGroupMatch = JointPeakGroupMatch.__table__


query_oven = bakery()
SCALE = get_scale()

ClassifierType = logistic_scoring.LogisticModelScorer


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
    try:
        maximum_volume = max(g.total_volume for g in group_matches)
    except ValueError:
        maximum_volume = 1.

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

    try:
        ppm_error = max(g.ppm_error for g in group_matches)
    except ValueError:
        ppm_error = None

    instance_dict = {
        "first_scan_id": min_scan,
        "last_scan_id": max_scan,
        "scan_density": scan_density,
        "ppm_error": ppm_error,
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


def _batch_merge_groups(id_bunches, database_manager, minimum_abundance_ratio):
    session = database_manager()
    results = []

    id_bunches = [[y for x in bunch for y in x] for bunch in id_bunches]
    bunches = [session.query(PeakGroupMatch).filter(
        PeakGroupMatch.id.in_(bunch)).all() for bunch in id_bunches]

    try:
        for group_matches in bunches:

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
            try:
                maximum_volume = max(g.total_volume for g in group_matches)
            except ValueError:
                maximum_volume = 1.

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

            try:
                ppm_error = max(g.ppm_error for g in group_matches)
            except ValueError:
                ppm_error = None

            instance_dict = {
                "first_scan_id": min_scan,
                "last_scan_id": max_scan,
                "scan_density": scan_density,
                "ppm_error": ppm_error,
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
                "theoretical_match_id": group_matches[0].theoretical_match_id,
                "theoretical_match_type": group_matches[0].theoretical_match_type,
                "matched": group_matches[0].matched
            }

            results.append((
                instance_dict, [p.id for p in group_matches]))
        conn = session.connection()
        relations = []
        for instance_dict, member_ids in results:
            joint_id = conn.execute(T_JointPeakGroupMatch.insert(), instance_dict).lastrowid
            relations.extend({"peak_group_id": i, "joint_group_id": joint_id} for i in member_ids)

        conn.execute(PeakGroupMatchToJointPeakGroupMatch.insert(), relations)
        session.commit()
        return len(results)
    except Exception, e:
        logging.exception("An exception occurred in _batch_merge_groups", exc_info=e)


def _exactly_one_group_joiner(ids, database_manager, minimum_abundance_ratio=None):
    session = database_manager()
    results = []
    try:
        ids = [bunch[0] for bunch in ids]
        items = [session.query(PeakGroupMatch).get(i) for i in ids]
        hsm_id = items[0].hypothesis_sample_match_id

        for peak_group in items:
            fingerprint = peak_group.fingerprint = ':'.join(map(str, peak_group.peak_data['scan_times']))
            instance_dict = {
                "first_scan_id": peak_group.first_scan_id,
                "last_scan_id": peak_group.last_scan_id,
                "scan_density": peak_group.scan_density,
                "ppm_error": peak_group.ppm_error,
                "centroid_scan_estimate": peak_group.centroid_scan_estimate,
                "average_a_to_a_plus_2_ratio": peak_group.average_a_to_a_plus_2_ratio,
                "average_signal_to_noise": peak_group.average_signal_to_noise,
                "charge_state_count": peak_group.charge_state_count,
                "modification_state_count": 1,
                "total_volume": peak_group.total_volume,
                "scan_count": peak_group.scan_count,
                "peak_data": peak_group.peak_data,
                "fingerprint": fingerprint,
                "weighted_monoisotopic_mass": peak_group.weighted_monoisotopic_mass,
                "hypothesis_sample_match_id": hsm_id,
                "theoretical_match_id": peak_group.theoretical_match_id,
                "theoretical_match_type": peak_group.theoretical_match_type,
                "matched": peak_group.matched
            }
            results.append(instance_dict)
        conn = session.connection()

        conn.execute(T_JointPeakGroupMatch.insert(), results)
        session.commit()
        products = []
        for peak_group in items:
            joint_id = session.query(JointPeakGroupMatch.id).filter(
                JointPeakGroupMatch.hypothesis_sample_match_id == hsm_id,
                JointPeakGroupMatch.theoretical_match_id == peak_group.theoretical_match_id,
                (JointPeakGroupMatch.fingerprint == peak_group.fingerprint
                 if peak_group.theoretical_match_id is None else True)
                ).first()
            products.append({"peak_group_id": peak_group.id, "joint_group_id": joint_id[0]})
        session.execute(PeakGroupMatchToJointPeakGroupMatch.insert(), products)
        session.commit()
    except Exception, e:
        logger.exception("An error occurred in _exactly_one_group_joiner", exc_info=e)
        raise e
    return len(ids)


class MatchJoiner(PipelineModule):
    def __init__(self, database_path, hypothesis_sample_match_id, minimum_abundance_ratio=0.01,
                 grouping_error_tolerance=2e-5, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.minimum_abundance_ratio = minimum_abundance_ratio
        self.grouping_error_tolerance = grouping_error_tolerance
        self.n_processes = n_processes
        self.mass_shift_map = self.manager().query(
            HypothesisSampleMatch).get(
            hypothesis_sample_match_id).parameters['mass_shift_map']

    def stream_matched_ids(self, chunk_size=2500):
        chunk_size *= SCALE

        session = self.manager()

        gen = itertools.groupby(session.query(PeakGroupMatch.id, PeakGroupMatch.theoretical_match_id).filter(
            PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            PeakGroupMatch.matched).order_by(PeakGroupMatch.theoretical_match_id).all(), operator.itemgetter(1))

        getter = operator.itemgetter(0)

        batch = []
        i = 0
        for key, group in gen:
            batch.append([(getter(o),) for o in group])
            i += 1

            if i > chunk_size:
                yield batch
                batch = []
                i = 0
        yield batch

        session.close()

    def stream_unmatched_ids(self, chunk_size=2500):
        chunk_size *= SCALE

        session = self.manager()
        unmatched = session.query(PeakGroupMatch).filter(
            PeakGroupMatch.theoretical_match_id == None,
            PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).order_by(
            PeakGroupMatch.weighted_monoisotopic_mass.desc()).all()
        if len(unmatched) == 0:
            raise StopIteration()

        batch = []
        i = 0
        # mass_shift_map = unmatched[0].hypothesis_sample_match.parameters['mass_shift_map']
        for bunch in _group_unmatched_peak_groups_by_shifts(
                unmatched, self.mass_shift_map, self.grouping_error_tolerance):
            if len(bunch) == 0:
                continue
            batch.append([(p.id,) for p in bunch])
            i += 1
            if i > chunk_size:
                yield batch
                batch = []
                i = 0
        yield batch
        session.close()

    def prepare_task_fn(self):
        if len(self.mass_shift_map) == 1:
            return functools.partial(
                _exactly_one_group_joiner,
                database_manager=self.manager)
        else:
            return functools.partial(
                _batch_merge_groups,
                database_manager=self.manager,
                minimum_abundance_ratio=self.minimum_abundance_ratio)

    def run(self):
        cntr = 0
        last = 0
        task_fn = self.prepare_task_fn()
        if self.n_processes > 1:
            self.inform("Merging Matched (Concurrent)")
            pool = multiprocessing.Pool(self.n_processes)
            for increment in pool.imap_unordered(task_fn, self.stream_matched_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d groups merged", cntr)
                    last = cntr

            self.inform("Merging Unmatched (Concurrent)")
            for increment in pool.imap_unordered(task_fn, self.stream_unmatched_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d groups merged", cntr)
                    last = cntr

            pool.close()
            pool.terminate()

        else:
            self.inform("Merging Matched (Sequential)")
            for increment in itertools.imap(task_fn, self.stream_matched_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d groups merged", cntr)
                    last = cntr

            self.inform("Merging Unmatched (Sequential)")
            for increment in itertools.imap(task_fn, self.stream_unmatched_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d groups merged", cntr)
                    last = cntr


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
        task = MatchJoiner(
            self.manager.path, self.hypothesis_sample_match_id,
            self.minimum_abundance_ratio, self.match_tolerance,
            self.n_processes)
        task.start()

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

    def construct_model_matrix(self, session=None):
        features = self.features
        ids = self.ids
        label = self.label

        if session is None:
            data_model_session = self.manager.session()
        else:
            data_model_session = session
        conn = data_model_session.connection()

        id_vec = flatten(conn.execute(select(ids).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall())
        feature_matrix = np.array(conn.execute(select(features).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall(),
            dtype=np.float64)
        label_vector = np.array(conn.execute(select(label).where(
            T_JointPeakGroupMatch.c.hypothesis_sample_match_id == self.hypothesis_sample_match_id)).fetchall())

        return [label_vector, id_vec, feature_matrix]

    def fit_and_score(self):
        session = self.manager.session()
        label_vector, id_vec, feature_matrix = self.construct_model_matrix(session)
        if self.use_legacy_coefficients:
            model = session.query(PeakGroupScoringModel).filter_by(
                name=PeakGroupScoringModel.GENERIC_MODEL_NAME).first()
            if model is None:
                raise PipelineException("Generic PeakGroupScoringModel Not Found")
            classifier = logistic_scoring.from_peak_group_scoring_model(model)
        else:
            classifier = ClassifierType()
            classifier.fit(feature_matrix, label_vector.ravel())
        self.classifier = classifier
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

        if not len(data_model_session.connection().execute(id_stmt).fetchmany(2)) == 2:
            raise PipelineException("Hypothesis-Sample Match ID matches maps no PeakGroupMatches")

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
