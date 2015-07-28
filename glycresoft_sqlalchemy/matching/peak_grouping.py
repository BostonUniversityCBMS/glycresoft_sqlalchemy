import logging
import multiprocessing
import functools
import itertools
try:
    logger = logging.getLogger("peak_grouping")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

from sqlalchemy import func, bindparam, select
from sqlalchemy.orm import make_transient

import numpy as np

from sklearn.linear_model import LogisticRegression

from ..data_model import (Decon2LSPeak, Decon2LSPeakGroup, Decon2LSPeakToPeakGroupMap, PipelineModule,
                          SampleRun, Hypothesis,
                          TheoreticalCompositionMap, MassShift, HypothesisSampleMatch,
                          PeakGroupDatabase, PeakGroupMatch, TempPeakGroupMatch)

from ..utils.database_utils import get_or_create
from ..utils import pickle
from ..report import chromatogram

TDecon2LSPeakGroup = Decon2LSPeakGroup.__table__
T_TempPeakGroupMatch = TempPeakGroupMatch.__table__
TPeakGroupMatch = PeakGroupMatch.__table__


class Decon2LSPeakGrouper(PipelineModule):
    '''
    Pipeline Step to post-process Decon2LSPeaks, clustering them by mass
    and calculating trends across groups.
    '''
    def __init__(self, database_path, sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, max_charge_state=8,
                 minimum_abundance_ratio=0.01, minimum_mass=1200., maximum_mass=15000.,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        self.grouping_error_tolerance = grouping_error_tolerance
        self.minimum_scan_count = minimum_scan_count
        self.n_processes = n_processes
        self.sample_run_id = sample_run_id
        self.max_charge_state = max_charge_state
        self.minimum_mass = minimum_mass
        self.maximum_mass = maximum_mass
        self.minimum_abundance_ratio = minimum_abundance_ratio

    def group_peaks(self):
        '''
        Perform incremental clustering of similar peaks along the chromatographic dimension.
        Use these peak groups for feature fitting in downstream operations
        '''
        logger.info("Grouping Peaks")
        session = self.manager.session()

        map_insert = Decon2LSPeakToPeakGroupMap.insert()
        group_insert = TDecon2LSPeakGroup.insert()

        group_count = 0
        sample_run_id = self.sample_run_id
        id_query = session.query(Decon2LSPeakGroup.id).filter(
            Decon2LSPeakGroup.sample_run_id == sample_run_id).selectable

        clear_relation_query = Decon2LSPeakToPeakGroupMap.delete().where(
            Decon2LSPeakToPeakGroupMap.c.group_id.in_(id_query))

        session.execute(clear_relation_query)
        session.commit()
        session.query(Decon2LSPeakGroup).filter(
            Decon2LSPeakGroup.sample_run_id == sample_run_id).delete(
            synchronize_session=False)
        session.commit()
        logger.info("Cleared? %r", session.query(Decon2LSPeakGroup).count())
        conn = session.connection()
        for count, row in enumerate(Decon2LSPeak.from_sample_run(session.query(
                Decon2LSPeak.id, Decon2LSPeak.monoisotopic_mass).filter(
                Decon2LSPeak.charge <= self.max_charge_state,
                Decon2LSPeak.monoisotopic_mass.between(self.minimum_mass, self.maximum_mass)
                ).order_by(
                Decon2LSPeak.intensity.desc()), self.sample_run_id)):
            peak_id, monoisotopic_mass = row

            # Calculate the window around which to group peaks
            window_radius = monoisotopic_mass * self.grouping_error_tolerance
            window_min = monoisotopic_mass - window_radius
            window_max = monoisotopic_mass + window_radius

            group = session.query(Decon2LSPeakGroup.id).filter(
                Decon2LSPeakGroup.weighted_monoisotopic_mass.between(
                    window_min, window_max), Decon2LSPeakGroup.sample_run_id == sample_run_id).first()
            # If no group is within the window of acceptable masses, create a new group
            # with this peak as the seed
            if group is None:
                result = conn.execute(group_insert, {
                    "weighted_monoisotopic_mass": monoisotopic_mass,
                    "sample_run_id": sample_run_id
                    })
                group_count += 1
                group_id = result.lastrowid
            # Otherwise take the first group that matched.
            else:
                group_id = group[0]

            # Add this peak to the many-to-many mapping between peaks and groups with respect to this group.
            # Many-to-many relationship allows a
            conn.execute(map_insert, {"peak_id": peak_id, "group_id": group_id})
            if count % 1000 == 0:
                logger.info("%d peaks clustered, %d groups", count, group_count)
        session.commit()
        session.close()

    def stream_group_ids(self):
        '''
        Session-bounded generator of Decon2LSPeakGroup ids
        '''
        session = self.manager.session()
        try:
            for gid in session.query(Decon2LSPeakGroup.id):
                yield gid[0]
        except Exception, e:
            logger.exception("An error occurred while streaming ids, %r", locals(), exc_info=e)
            raise e
        finally:
            session.close()

    def update_groups(self):
        '''
        Once peaks have been completely assigned to clusters, calculate peak group
        features for each group.
        '''
        logger.info("Updating Groups")
        session = self.manager.session()
        task_fn = functools.partial(
            fill_out_group, database_manager=self.manager, minimum_scan_count=self.minimum_scan_count)
        accumulator = []
        if self.n_processes > 1:
            logger.info("Running concurrently")
            pool = multiprocessing.Pool(self.n_processes)
            count = 0
            for group in pool.imap_unordered(task_fn, self.stream_group_ids(), chunksize=25):
                count += 1
                if group is not None:
                    accumulator.append(group)
                if count % 10000 == 0:
                    session.bulk_update_mappings(Decon2LSPeakGroup, accumulator)
                    session.commit()
                    accumulator = []
                    logger.info("%d groups completed", count)
            pool.terminate()
        else:
            for count, group_id in enumerate(session.query(Decon2LSPeakGroup.id)):
                group = task_fn(group_id)
                if group is not None:
                    accumulator.append(group)
                if count % 10000 == 0:
                    session.bulk_update_mappings(Decon2LSPeakGroup, accumulator)
                    session.commit()
                    accumulator = []
                    logger.info("%d groups completed", count)

        session.bulk_update_mappings(Decon2LSPeakGroup, accumulator)
        session.commit()
        session.close()

    def estimate_trends(self):
        '''
        After assigning peak group features, impute the global
        trend for peak and scan shapes
        '''
        logger.info("Estimating peak trends")
        session = self.manager.session()
        engine = self.manager.connect()

        conn = engine.connect()
        cen_alpha, cen_beta = centroid_scan_error_regression(session)
        expected_a_alpha, expected_a_beta = expected_a_peak_regression(session)

        update_expr = TDecon2LSPeakGroup.update().values(
            centroid_scan_error=func.abs(
                TDecon2LSPeakGroup.c.centroid_scan_estimate - (
                    cen_alpha + cen_beta * TDecon2LSPeakGroup.c.weighted_monoisotopic_mass)),
            a_peak_intensity_error=func.abs(
                TDecon2LSPeakGroup.c.average_a_to_a_plus_2_ratio - (
                    expected_a_alpha + expected_a_beta * TDecon2LSPeakGroup.c.weighted_monoisotopic_mass))
                )
        max_weight = conn.execute(select([func.max(TDecon2LSPeakGroup.c.weighted_monoisotopic_mass)])).scalar()
        slices = [max_weight * float(i)/10. for i in range(1, 11)]
        for i in range(1, len(slices)):
            transaction = conn.begin()
            lower = slices[i - 1]
            upper = slices[i]
            logger.info("Updating slice %f-%f", lower, upper)
            step = update_expr.where(
                TDecon2LSPeakGroup.c.weighted_monoisotopic_mass.between(
                    lower, upper))
            conn.execute(step)
            transaction.commit()
        conn.close()

    def run(self):
        self.group_peaks()
        self.update_groups()
        self.estimate_trends()


def update_fit(group_id, database_manager, cen_alpha, cen_beta, expected_a_alpha, expected_a_beta):
    '''
    Given linear predicters, compute fitted values for centroid scan and
    expected monoisotopic intensity, and report the magnitude of the difference
    between the observed and expected.

    Parameters
    ----------
    group_id: int
        The primary key of a :class:`Decon2LSPeakGroup` record to predict on
    database_manager: DatabaseManager
        Provide connection to the database
    cen_alpha, cen_beta: float
        Linear regression parameters for centroid scan fit by weighted_monoisotopic_mass
    expected_a_alpha, expected_a_beta: float
        Linear regression parameters for expected monoisotopic intensity
        by weighted_monoisotopic_mass

    Returns
    -------
    dict: Update parameters for the :class:`Decon2LSPeakGroup`.
    '''
    session = database_manager.session()
    try:
        group = session.query(
            Decon2LSPeakGroup.id, Decon2LSPeakGroup.centroid_scan_estimate,
            Decon2LSPeakGroup.average_a_to_a_plus_2_ratio,
            Decon2LSPeakGroup.weighted_monoisotopic_mass).filter(
            Decon2LSPeakGroup.id == group_id).first()
        update = dict(
            centroid_scan_error=abs(group[1] - (cen_alpha + cen_beta * group[3])),
            a_peak_intensity_error=abs(group[2] - (expected_a_alpha + expected_a_beta * group[3])),
            gid=group_id)
        session.close()
        return update
    except Exception, e:
        logger.exception("An error occured. %r", locals(), exc_info=e)
        session.close()


def fill_out_group(group_id, database_manager, minimum_scan_count=1,
                   minimum_abundance_ratio=0.01):
    """Calculate peak group statistics for a given uninitialized :class:`Decon2LSPeakGroup`.

    Parameters
    ----------
    group_id : int
        The primary key of :class:`Decon2LSPeakGroup` to update
    database_manager : DatabaseManager
        Provides database access
    minimum_scan_count : int, optional
        Minimum number of scans to observe in order to accept a group

    Returns
    -------
    int : 0 if the group did not pass criteria, 1 otherwise
    """
    session = database_manager.session()
    try:
        group = session.query(Decon2LSPeakGroup).get(group_id)

        peaks = group.peaks.order_by(Decon2LSPeak.scan_id.asc()).all()
        for peak in peaks:
            peak.scan_time = peak.scan.time
        make_transient(group)

        charge_states = set()
        scan_times = set()
        peak_ids = []
        intensities = []
        full_width_half_maxes = []
        monoisotopic_masses = []
        signal_noises = []
        a_to_a_plus_2_ratios = []

        max_intensity = float(max(p.intensity for p in peaks))

        remove_peaks = set()

        for peak in peaks:
            if minimum_abundance_ratio > peak.intensity / max_intensity:
                remove_peaks.add(peak.id)
                continue
            charge_states.add(peak.charge)
            scan_times.add(peak.scan_time)
            peak_ids.append(peak.id)
            intensities.append(peak.intensity)
            full_width_half_maxes.append(peak.full_width_half_max)
            monoisotopic_masses.append(peak.monoisotopic_mass)
            signal_noises.append(peak.signal_to_noise)
            a_to_a_plus_2_ratios.append(
                peak.monoisotopic_intensity / float(peak.monoisotopic_plus_2_intensity)
                if peak.monoisotopic_plus_2_intensity > 0
                else 0
            )
        scan_count = len(scan_times)
        if scan_count < minimum_scan_count:
            logger.info("Deleting %r, with %d scans", group, scan_count)
            session.delete(group)
            session.commit()
            return None

        peaks = [p for p in peaks if p.id not in remove_peaks]
        scan_times = sorted(scan_times)
        count = len(monoisotopic_masses)
        fcount = float(count)

        min_scan = min(scan_times)
        max_scan = max(scan_times)

        average_signal_to_noise = sum(signal_noises) / fcount
        average_a_to_a_plus_2_ratio = sum(a_to_a_plus_2_ratios) / fcount

        total_volume = sum(i * f for i, f in zip(intensities, full_width_half_maxes))

        weighted_monoisotopic_mass = sum(m * i for m, i in zip(
            monoisotopic_masses, intensities)) / float(sum(intensities))

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
            scan_density = sum(window_densities) / float(len(window_densities))
        else:
            scan_density = 0.

        return {
            "id": group.id,
            "scan_density": scan_density,
            "average_signal_to_noise": average_signal_to_noise,
            "average_a_to_a_plus_2_ratio": average_a_to_a_plus_2_ratio,
            "centroid_scan_estimate": sum(scan_times) / fcount,
            "total_volume": total_volume,
            "weighted_monoisotopic_mass": weighted_monoisotopic_mass,
            "scan_count": scan_count,
            "first_scan_id": min_scan,
            "last_scan_id": max_scan,
            "charge_state_count": len(charge_states),
            "peak_data": {
                "scan_times": [p.scan_time for p in peaks],
                "intensities": [p.scan_id for p in peaks],
                "peak_ids": [p.id for p in peaks],
            }
        }
        # session.close()
        # return 1
    except Exception, e:
        logger.exception("An error occured. %r", 1, exc_info=e)
        session.close()


def expanding_window(series, threshold_gap_size=100):
    windows = []
    current_window = []
    last_item = 0
    for item in series:
        if item - last_item > threshold_gap_size:
            if len(current_window) > 0:
                windows.append(current_window)
            current_window = []
        current_window.append(item)
        last_item = item
    if len(current_window) > 0:
        windows.append(current_window)
    return windows


def centroid_scan_error_regression(session, minimum_abundance=250):
    '''
    alpha, beta = centroid_scan_error_regression(session)
    r = session.execute(select([func.abs(
            Decon2LSPeakGroup.centroid_scan_estimate - (
                alpha + Decon2LSPeakGroup.weighted_monoisotopic_mass * beta))]))
    r.fetchmany(200)
    '''
    mean_centroid_scan_estimate = session.query(func.avg(Decon2LSPeakGroup.centroid_scan_estimate)).filter(
            Decon2LSPeakGroup.total_volume > minimum_abundance).first()[0]
    mean_weighted_mass = session.query(func.avg(Decon2LSPeakGroup.weighted_monoisotopic_mass)).filter(
            Decon2LSPeakGroup.total_volume > minimum_abundance).first()[0]
    beta = float(
        session.query((func.sum(
                Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass) *
                      Decon2LSPeakGroup.centroid_scan_estimate)
                  / func.sum(
                ((Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass) *
                 (Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass))
            )).filter(Decon2LSPeakGroup.total_volume > minimum_abundance).first()[0])
    alpha = mean_centroid_scan_estimate - beta * mean_weighted_mass
    return alpha, beta


def expected_a_peak_regression(session):
    '''
    alpha, beta = expected_a_peak_regression(session)
    r = session.execute(select([func.abs(
            Decon2LSPeakGroup.average_a_to_a_plus_2_ratio - (
                alpha + Decon2LSPeakGroup.weighted_monoisotopic_mass * beta))]))
    r.fetchmany(200)
    '''
    mean_weighted_mass = session.query(func.avg(Decon2LSPeakGroup.weighted_monoisotopic_mass)).first()[0]
    mean_a_to_a_plus_2_ratio = session.query(func.avg(Decon2LSPeakGroup.average_a_to_a_plus_2_ratio)).first()[0]
    beta = float(
        session.query((func.sum(
                    (Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass) *
                    Decon2LSPeakGroup.average_a_to_a_plus_2_ratio)
                ) / func.sum(
                (Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass) *
                (Decon2LSPeakGroup.weighted_monoisotopic_mass - mean_weighted_mass)
        )).first()[0]
    )
    alpha = mean_a_to_a_plus_2_ratio - beta * mean_weighted_mass
    return alpha, beta


default_match_tolerance = 2e-5


def ppm_error(x, y):
    return (x - y) / y


def match_peak_group(search_id, search_type, database_manager, observed_ions_manager,
                     matching_tolerance, mass_shift_map, sample_run_id, hypothesis_sample_match_id):
    session = database_manager.session()
    try:
        search_target = session.query(search_type).get(search_id)

        base_mass = search_target.calculated_mass
        matches = []
        for mass_shift, count_range in mass_shift_map.items():
            shift = mass_shift.mass
            for shift_count in range(count_range):
                total_mass = base_mass + shift * shift_count
                for mass_match in observed_ions_manager.ppm_match_tolerance_search(
                        total_mass, matching_tolerance, sample_run_id):
                    mass_error = ppm_error(mass_match.weighted_monoisotopic_mass, total_mass)
                    matches.append((mass_match, mass_error, mass_shift, shift_count))
        params = []
        for mass_match, mass_error, mass_shift, shift_count in matches:
            # logger.debug("Peak data: %r", mass_match.peak_data)
            case = {
                "hypothesis_sample_match_id": hypothesis_sample_match_id,
                "theoretical_match_type": search_type.__name__,
                "theoretical_match_id": search_id,
                "charge_state_count": mass_match.charge_state_count,
                "scan_count": mass_match.scan_count,
                "first_scan_id": mass_match.first_scan_id,
                "last_scan_id": mass_match.last_scan_id,
                "ppm_error": mass_error,
                "scan_density": mass_match.scan_density,
                "weighted_monoisotopic_mass": mass_match.weighted_monoisotopic_mass + mass_shift.mass * shift_count,
                "total_volume": mass_match.total_volume,
                "average_a_to_a_plus_2_ratio": mass_match.average_a_to_a_plus_2_ratio,
                "a_peak_intensity_error": mass_match.a_peak_intensity_error,
                "centroid_scan_estimate": mass_match.centroid_scan_estimate,
                "centroid_scan_error": mass_match.centroid_scan_error,
                "average_signal_to_noise": mass_match.average_signal_to_noise,
                "peak_data": mass_match.peak_data,
                "matched": True,
                "mass_shift_type": mass_shift.name,
                "mass_shift_count": shift_count,
                "peak_group_id": mass_match.id,
                "peak_ids": mass_match.peak_ids
            }
            params.append(case)
        # session.bulk_insert_mappings(PeakGroupMatch, params)
        # session.commit()
        session.close()
        return params
    except Exception, e:
        logger.exception("An exception occurred in match_peak_group, %r", locals(), exc_info=e)
        raise e
    finally:
        session.close()


class PeakGroupMatching(PipelineModule):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=None, hypothesis_sample_match_id=None,
                 search_type="TheoreticalGlycanComposition",
                 match_tolerance=2e-5, mass_shift_map=None,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        session = self.manager.session()
        no_shift, created = get_or_create(session, MassShift, mass=0.0, name=u"NoShift")
        session.expunge(no_shift)
        lcms_database = PeakGroupDatabase(observed_ions_path)
        lcms_session = lcms_database.session()
        sample_run = lcms_session.query(SampleRun).get(sample_run_id or 1)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id

        hypothesis_sample_match = session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
        hypothesis_sample_match.sample_run_name = sample_run.name
        hypothesis_sample_match.parameters = {}
        hypothesis_sample_match.parameters['mass_shift_map'] = mass_shift_map
        session.add(hypothesis_sample_match)
        session.commit()

        self.hypothesis_sample_match_id = hypothesis_sample_match.id
        if mass_shift_map is None:
            mass_shift_map = {no_shift: 1}
        else:
            mass_shift_map[no_shift] = 1
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.mass_shift_map = mass_shift_map
        self.match_tolerance = match_tolerance
        self.n_processes = n_processes
        self.lcms_database = lcms_database
        self.search_type = TheoreticalCompositionMap[search_type]
        session.close()

    def stream_ids(self):
        session = self.manager.session()
        try:
            for gid in session.query(self.search_type.id):
                yield gid[0]
        except Exception, e:
            logger.info("An exception occurred while streaming ids", exc_info=e)
            raise e
        finally:
            session.close()

    def prepare_task_fn(self):
        fn = functools.partial(
            match_peak_group,
            search_type=self.search_type,
            database_manager=self.manager,
            observed_ions_manager=self.lcms_database,
            matching_tolerance=self.match_tolerance,
            mass_shift_map=self.mass_shift_map,
            sample_run_id=self.sample_run_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id)
        return fn

    def run(self):
        session = self.manager.session()

        task_fn = self.prepare_task_fn()
        counter = 0
        accumulator = []
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_ids(), chunksize=500):
                counter += 1
                if res is not None:
                    accumulator.extend(res)
                if counter % 1000 == 0:
                    logger.info("%d masses searched", counter)
                if len(accumulator) > 1000:
                    session.bulk_insert_mappings(PeakGroupMatch, accumulator)
                    session.commit()
                    accumulator = []
        else:
            for res in itertools.imap(task_fn, self.stream_ids()):
                counter += 1
                if res is not None:
                    accumulator.extend(res)
                if counter % 1000 == 0:
                    logger.info("%d masses searched", counter)
                if len(accumulator) > 10000:
                    session.bulk_insert_mappings(PeakGroupMatch, accumulator)
                    session.commit()
                    accumulator = []

        session.bulk_insert_mappings(PeakGroupMatch, accumulator)
        session.commit()
        session.close()


ClassifierType = LogisticRegression


class PeakGroupClassification(PipelineModule):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=None, hypothesis_sample_match_id=None,
                 model_parameters=None, minimum_mass=1500, maximum_mass=15000):
        self.database_manager = self.manager_type(database_path)
        self.lcms_database = PeakGroupDatabase(observed_ions_path)
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.model_parameters = model_parameters
        self.classifier = None
        self.minimum_mass = minimum_mass
        self.maximum_mass = maximum_mass

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
        labels = [
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
            'peak_ids',
            'peak_data',
            'ms1_score',
            'matched'
        ]
        stmt = lcms_database_session.query(
                Decon2LSPeakGroup).filter(
                Decon2LSPeakGroup.sample_run_id == self.sample_run_id,
                Decon2LSPeakGroup.weighted_monoisotopic_mass.between(self.minimum_mass, self.maximum_mass))
        batch = lcms_database_session.connection().execute(stmt.selectable)

        conn = data_model_session.connection()

        # Move all Decon2LSPeakGroups, regardless of whether or not they matched to
        # the temporary table.
        while True:
            items = batch.fetchmany(10000)
            if len(items) == 0:
                break
            self.buffer = items[0]
            conn.execute(T_TempPeakGroupMatch.insert(), [dict(zip(labels, row)) for row in items])
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
            TempPeakGroupMatch).filter(
            ~TempPeakGroupMatch.matched).selectable

        def transform(row):
            out = dict(zip(labels, row))
            peak_group_match_id = out.pop('id')
            out["theoretical_match_type"] = None
            out['matched'] = False
            out['hypothesis_sample_match_id'] = self.hypothesis_sample_match_id
            out['peak_group_id'] = peak_group_match_id
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
        features = [
            T_TempPeakGroupMatch.c.charge_state_count,
            T_TempPeakGroupMatch.c.scan_density,
            T_TempPeakGroupMatch.c.scan_count,
            T_TempPeakGroupMatch.c.total_volume,
            T_TempPeakGroupMatch.c.a_peak_intensity_error,
            T_TempPeakGroupMatch.c.centroid_scan_error,
            T_TempPeakGroupMatch.c.average_signal_to_noise
        ]
        label = [T_TempPeakGroupMatch.c.matched]
        ids = [T_TempPeakGroupMatch.c.id]

        data_model_session = self.database_manager.session()
        conn = data_model_session.connection()
        feature_matrix = np.array(conn.execute(select(features)).fetchall(), dtype=np.float64)
        label_vector = np.array(conn.execute(select(label)).fetchall())
        classifier = ClassifierType()
        if self.model_parameters is None:
            classifier.fit(feature_matrix, label_vector.ravel())
        else:
            classifier.coef_ = np.asarray(self.model_parameters)
        scores = classifier.predict_proba(feature_matrix)[:, 1]
        i = 0
        for group_id, in conn.execute(select(ids)):
            conn.execute(
                TPeakGroupMatch.update().where(
                    TPeakGroupMatch.c.peak_group_id == group_id).values(
                    ms1_score=scores[i]))
            i += 1
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
        logger.info("Coefficients: %r", self.classifier.coef_)
        self.clear_peak_groups()


class LCMSPeakClusterSearch(PipelineModule):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, hypothesis_sample_match_id=None,
                 search_type="TheoreticalGlycanComposition",
                 match_tolerance=2e-5, mass_shift_map=None, regression_parameters=None,
                 n_processes=4, **kwargs):
        self.manager = self.manager_type(database_path)
        self.observed_ions_path = observed_ions_path
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.grouping_error_tolerance = grouping_error_tolerance
        self.minimum_scan_count = minimum_scan_count
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.search_type = search_type
        self.match_tolerance = match_tolerance
        self.mass_shift_map = mass_shift_map
        self.n_processes = n_processes
        self.regression_parameters = regression_parameters
        self.options = kwargs

    def run(self):
        grouper = Decon2LSPeakGrouper(
            self.observed_ions_path, self.sample_run_id, self.grouping_error_tolerance,
            self.minimum_scan_count, self.n_processes)
        if not self.options.get("skip_grouping", False):
            grouper.start()

        database_manager = self.manager
        session = database_manager.session()
        sample_name = grouper.manager.session().query(SampleRun).get(self.sample_run_id).name
        hypothesis_sample_match, created = get_or_create(
            session, HypothesisSampleMatch,
            id=self.hypothesis_sample_match_id,
            target_hypothesis_id=self.hypothesis_id)
        if self.hypothesis_sample_match_id is not None:
            if not self.options.get("skip_matching", False):
                hypothesis_sample_match.peak_group_matches.delete(synchronize_session=False)
                session.commit()
                logger.info("Cleared? %r", hypothesis_sample_match.peak_group_matches.count())
            else:
                hypothesis_sample_match.peak_group_matches.filter(
                    ~PeakGroupMatch.matched).delete(synchronize_session=False)
                session.commit()
                logger.info("Cleared? %r", hypothesis_sample_match.peak_group_matches.filter(
                    ~PeakGroupMatch.matched).count())
        hypothesis_sample_match.name = "{}_on_{}_ms1".format(
            session.query(Hypothesis).get(self.hypothesis_id).name,
            sample_name)
        session.add(hypothesis_sample_match)
        session.commit()

        self.hypothesis_sample_match_id = hypothesis_sample_match.id

        matcher = PeakGroupMatching(
            self.database_path, self.observed_ions_path, self.hypothesis_id,
            self.sample_run_id, self.hypothesis_sample_match_id,
            self.search_type, self.match_tolerance, self.mass_shift_map,
            self.n_processes)

        if not self.options.get("skip_matching", False):
            matcher.start()

        classifier = PeakGroupClassification(
            self.database_path, self.observed_ions_path, self.hypothesis_id,
            self.sample_run_id, self.hypothesis_sample_match_id, self.regression_parameters)

        classifier.start()
        hypothesis_sample_match.parameters['classifier'] = pickle.dumps(classifier.classifier)
