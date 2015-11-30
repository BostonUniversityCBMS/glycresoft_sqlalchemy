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

from sqlalchemy.ext.baked import bakery
from sqlalchemy import func, bindparam, select
from sqlalchemy.orm import make_transient

import numpy as np

from sklearn.linear_model import LogisticRegression

from glycresoft_sqlalchemy.data_model import (Decon2LSPeak, Decon2LSPeakGroup, Decon2LSPeakToPeakGroupMap, PipelineModule,
                          SampleRun, Hypothesis, MS1GlycanHypothesisSampleMatch, hypothesis_sample_match_type,
                          TheoreticalCompositionMap, MassShift, HypothesisSampleMatch,
                          PeakGroupDatabase, PeakGroupMatch, TempPeakGroupMatch, JointPeakGroupMatch,
                          PeakGroupMatchToJointPeakGroupMatch)

from glycresoft_sqlalchemy.utils.database_utils import get_or_create, toggle_indices
from glycresoft_sqlalchemy.utils import pickle
from glycresoft_sqlalchemy.utils.collectiontools import flatten

from .common import (
    expanding_window, expected_a_peak_regression, centroid_scan_error_regression)

TDecon2LSPeakGroup = Decon2LSPeakGroup.__table__
T_TempPeakGroupMatch = TempPeakGroupMatch.__table__
TPeakGroupMatch = PeakGroupMatch.__table__
T_JointPeakGroupMatch = JointPeakGroupMatch.__table__


query_oven = bakery()
get_group_id_by_mass_window = query_oven(lambda session: session.query(Decon2LSPeakGroup.id))
get_group_id_by_mass_window += lambda q: q.filter(Decon2LSPeakGroup.weighted_monoisotopic_mass.between(
            bindparam("lower"), bindparam("upper")))
get_group_id_by_mass_window += lambda q: q.filter(Decon2LSPeakGroup.sample_run_id == bindparam("sample_run_id"))


class Decon2LSPeakGrouper(PipelineModule):
    '''
    Pipeline Step to post-process Decon2LSPeaks, clustering them by mass
    and calculating trends across groups.
    '''
    def __init__(self, database_path, sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, max_charge_state=8,
                 minimum_abundance_ratio=0.01, minimum_mass=1200.,
                 maximum_mass=15000., minimum_signal_to_noise=1.,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        self.grouping_error_tolerance = grouping_error_tolerance
        self.minimum_scan_count = minimum_scan_count
        self.n_processes = n_processes
        self.sample_run_id = sample_run_id
        self.max_charge_state = max_charge_state
        self.minimum_signal_to_noise = minimum_signal_to_noise
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
        map_items = []
        for count, row in enumerate(Decon2LSPeak.from_sample_run(session.query(
                Decon2LSPeak.id, Decon2LSPeak.monoisotopic_mass).filter(
                Decon2LSPeak.charge <= self.max_charge_state,
                Decon2LSPeak.signal_to_noise >= self.minimum_signal_to_noise,
                Decon2LSPeak.monoisotopic_mass.between(self.minimum_mass, self.maximum_mass)
                ).order_by(
                Decon2LSPeak.intensity.desc()), self.sample_run_id)):
            peak_id, monoisotopic_mass = row

            # Calculate the window around which to group peaks
            window_radius = monoisotopic_mass * self.grouping_error_tolerance
            window_min = monoisotopic_mass - window_radius
            window_max = monoisotopic_mass + window_radius

            # group = session.query(Decon2LSPeakGroup.id).filter(
            #     Decon2LSPeakGroup.weighted_monoisotopic_mass.between(
            #         window_min, window_max), Decon2LSPeakGroup.sample_run_id == sample_run_id).first()
            group = get_group_id_by_mass_window(session).params(
                lower=window_min, upper=window_max, sample_run_id=sample_run_id).first()
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

            map_items.append({"peak_id": peak_id, "group_id": group_id})
            # Add this peak to the many-to-many mapping between peaks and groups with respect to this group.
            # Many-to-many relationship allows a peak to be assigned to a group without changing the peak.
            if count % 1000 == 0:
                logger.info("%d peaks clustered, %d groups", count, group_count)
                conn.execute(map_insert, map_items)
                map_items = []

        conn.execute(map_insert, map_items)
        map_items = []
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
        slices = [0] + [max_weight * float(i)/10. for i in range(1, 11)]
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

        transaction = conn.begin()
        lower = slices[len(slices) - 1]
        step = update_expr.where(
            TDecon2LSPeakGroup.c.weighted_monoisotopic_mass >= lower)
        conn.execute(step)
        transaction.commit()

        conn.close()

    def run(self):
        self.group_peaks()
        self.update_groups()
        self.estimate_trends()


def update_fit(group_id, database_manager, cen_alpha, cen_beta, expected_a_alpha, expected_a_beta,
               source_model=Decon2LSPeakGroup):
    '''
    Given linear predicters, compute fitted values for centroid scan and
    expected monoisotopic intensity, and report the magnitude of the difference
    between the observed and expected.

    Parameters
    ----------
    group_id: int
        The primary key of a :class:`source_model` record to predict on
    database_manager: DatabaseManager
        Provide connection to the database
    cen_alpha, cen_beta: float
        Linear regression parameters for centroid scan fit by weighted_monoisotopic_mass
    expected_a_alpha, expected_a_beta: float
        Linear regression parameters for expected monoisotopic intensity
        by weighted_monoisotopic_mass

    Returns
    -------
    dict: Update parameters for the :class:`source_model`.
    '''
    session = database_manager.session()
    try:
        group = session.query(
            source_model.id, source_model.centroid_scan_estimate,
            source_model.average_a_to_a_plus_2_ratio,
            source_model.weighted_monoisotopic_mass).filter(
            source_model.id == group_id).first()
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
                "charge_states": charge_states,
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
