import logging
import multiprocessing
import functools
import itertools
try:
    logger = logging.getLogger("peak_grouping")
except:
    pass

from sqlalchemy import func, bindparam
from ..data_model import (Decon2LSPeak, Decon2LSPeakGroup, Decon2LSPeakToPeakGroupMap, PipelineModule,
                          TheoreticalGlycanComposition, TheoreticalGlycopeptideComposition, SampleRun,
                          TheoreticalCompositionMap, MassShift, HypothesisSampleMatch, DatabaseManager,
                          PeakGroupDatabase, PeakGroupMatch)

from ..utils.database_utils import get_or_create

TDecon2LSPeakGroup = Decon2LSPeakGroup.__table__


class Decon2LSPeakGrouper(PipelineModule):
    '''
    Pipeline Step to post-process Decon2LSPeaks, clustering them by mass
    and calculating trends across groups.
    '''
    def __init__(self, database_path, sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.grouping_error_tolerance = grouping_error_tolerance
        self.minimum_scan_count = minimum_scan_count
        self.n_processes = n_processes
        self.sample_run_id = sample_run_id

    def group_peaks(self):
        '''
        Perform incremental clustering of similar peaks along the chromatographic dimension.
        Use these peak groups for feature fitting in downstream operations
        '''
        logger.info("Grouping Peaks")
        session = self.manager.session()
        conn = session.connection()
        map_insert = Decon2LSPeakToPeakGroupMap.insert()
        group_insert = TDecon2LSPeakGroup.insert()
        for count, row in enumerate(session.query(Decon2LSPeak.id, Decon2LSPeak.monoisotopic_mass).order_by(
                Decon2LSPeak.intensity.desc()).filter(Decon2LSPeak.from_sample_run(self.sample_run_id))):
            peak_id, monoisotopic_mass = row

            # Calculate the window around which to group peaks
            window_radius = monoisotopic_mass * self.grouping_error_tolerance
            window_min = monoisotopic_mass - window_radius
            window_max = monoisotopic_mass + window_radius

            group = session.query(Decon2LSPeakGroup.id).filter(
                Decon2LSPeakGroup.weighted_monoisotopic_mass.between(
                    window_min, window_max)).first()
            # If no group is within the window of acceptable masses, create a new group
            # with this peak as the seed
            if group is None:
                result = conn.execute(group_insert, {"weighted_monoisotopic_mass": monoisotopic_mass})
                group_id = result.lastrowid
            # Otherwise take the first group that matched.
            else:
                group_id = group[0]

            # Add this peak to the many-to-many mapping between peaks and groups with respect to this group.
            # Many-to-many relationship allows a
            conn.execute(map_insert, {"peak_id": peak_id, "group_id": group_id})
            if count % 1000 == 0:
                logger.info("%d peaks clustered", count)
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
        if self.n_processes > 1:
            logger.info("Running concurrently")
            pool = multiprocessing.Pool(self.n_processes)
            count = 0
            for group in pool.imap_unordered(task_fn, self.stream_group_ids(), chunksize=25):
                # if group is not None:
                #     session.add(group)
                count += 1
                if count % 1000 == 0:
                    logger.info("%d groups completed", count)

        else:
            for count, group_id in enumerate(session.query(Decon2LSPeakGroup.id)):
                task_fn(group_id)
                # if group is not None:
                #     session.add(group)
                if count % 1000 == 0:
                    logger.info("%d groups completed", count)

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

        task_fn = functools.partial(
            update_fit, database_manager=self.manager, cen_alpha=cen_alpha,
            cen_beta=cen_beta, expected_a_alpha=expected_a_alpha,
            expected_a_beta=expected_a_beta)
        count = 0
        accumulator = []
        stmt = TDecon2LSPeakGroup.update().where(TDecon2LSPeakGroup.c.id == bindparam("gid")).\
            values(centroid_scan_error=bindparam("centroid_scan_error"),
                   a_peak_intensity_error=bindparam("a_peak_intensity_error"))
        transaction = conn.begin()
        if self.n_processes > 1:
            logger.info("Updating peaks concurrently")
            pool = multiprocessing.Pool(self.n_processes)
            for update in pool.imap_unordered(task_fn, self.stream_group_ids(), chunksize=50):
                accumulator.append(update)
                count += 1
                if count % 100000 == 0:
                    conn.execute(stmt, accumulator)
                    accumulator = []
                    transaction.commit()
                    transaction = conn.begin()
            pool.close()
        else:
            logger.info("Updating peaks sequentially")
            for gid in self.stream_group_ids():
                update = task_fn(gid)
                conn.execute(stmt, update)
                count += 1
                if count % 1000 == 0:
                    pass
                    # conn = session.connection()
        conn.execute(stmt, accumulator)
        accumulator = []
        transaction.commit()
        transaction.close()
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
        return update
    except Exception, e:
        logger.exception("An error occured. %r", locals(), exc_info=e)
    finally:
        session.close()


def fill_out_group(group_id, database_manager, minimum_scan_count=1):
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
    TYPE : Description
    """
    session = database_manager.session()
    try:
        group = session.query(Decon2LSPeakGroup).get(group_id)

        charge_states = set()
        scan_ids = set()
        peak_ids = []
        intensities = []
        full_width_half_maxes = []
        monoisotopic_masses = []
        signal_noises = []
        a_to_a_plus_2_ratios = []
        for peak in group.peaks:
                charge_states.add(peak.charge)
                scan_ids.add(peak.scan_id)
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
        scan_count = len(scan_ids)
        if scan_count < minimum_scan_count:
            logger.info("Deleting %r, with %d scans", group, scan_count)
            session.delete(group)
            return 0

        scan_ids = list(scan_ids)
        count = len(monoisotopic_masses)
        fcount = float(count)

        min_scan = min(scan_ids)
        max_scan = max(scan_ids)

        average_signal_to_noise = sum(signal_noises) / fcount
        average_a_to_a_plus_2_ratio = sum(a_to_a_plus_2_ratios) / fcount

        total_volume = sum(i * f for i, f in zip(intensities, full_width_half_maxes))

        weighted_monoisotopic_mass = sum(m * i for m, i in zip(
            monoisotopic_masses, intensities)) / float(sum(intensities))

        scan_density = scan_count / float(max_scan - min_scan) if scan_count > 1 else 0

        group.scan_density = scan_density
        group.most_abundant = True
        group.average_signal_to_noise = average_signal_to_noise
        group.average_a_to_a_plus_2_ratio = average_a_to_a_plus_2_ratio
        group.centroid_scan_estimate = sum(scan_ids) / fcount
        group.total_volume = total_volume
        group.weighted_monoisotopic_mass = weighted_monoisotopic_mass
        group.scan_count = scan_count
        group.first_scan_id = min_scan
        group.last_scan_id = max_scan
        group.charge_state_count = len(charge_states)
        group.peak_data = {
            "scan_ids": scan_ids,
            "signal_to_noises": signal_noises,
            "a_to_a_plus_2_ratios": a_to_a_plus_2_ratios,
            "monoisotopic_masses": monoisotopic_masses,
            "intensities": intensities,
            "peak_ids": peak_ids
        }
        session.add(group)
        session.commit()
        return 1
    except Exception, e:
        logger.exception("An error occured. %r", locals(), exc_info=e)
    finally:
        session.close()


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

    # centroid_scan_estimate_sum_of_squares = session.query(
    #     (func.sum(
    #         (Decon2LSPeakGroup.centroid_scan_estimate - mean_centroid_scan_estimate) *
    #         (Decon2LSPeakGroup.centroid_scan_estimate - mean_centroid_scan_estimate)))).first()[0]
    # residuals = session.execute(select([func.abs(
    #         Decon2LSPeakGroup.centroid_scan_estimate - (
    #             alpha + Decon2LSPeakGroup.weighted_monoisotopic_mass * beta))]))
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
                     matching_tolerance, mass_shift_map, sample_run_id, hypothesis_match_id):
    session = database_manager.session()
    search_target = session.query(search_type).get(search_id)

    base_mass = search_target.mass
    matches = []
    for mass_shift, count_range in mass_shift_map.items():
        shift = mass_shift.mass
        for shift_count in range(count_range):
            total_mass = base_mass + shift * shift_count
            for mass_match in observed_ions_manager.ppm_match_tolerance_search(
                    total_mass, matching_tolerance, sample_run_id):
                mass_error = ppm_error(mass_match, total_mass)
                matches.append(mass_match, mass_error, mass_shift, shift_count)
    params = []
    for mass_match, mass_error, mass_shift, shift_count in matches:

        case = {
            "hypothesis_match_id": hypothesis_match_id,
            "theoretical_match_type": search_type,
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
            "peak_ids": mass_match.peak_ids,
            "peak_data": mass_match.peak_data,
            "matched": True,
            "mass_shift_type": mass_shift,
            "mass_shift_count": shift_count
        }
        params.append(case)
        mass_match.matched = True
        session.add(mass_match)
    session.bulk_insert_mappings(PeakGroupMatch, params)
    session.commit()
    return len(params)


class PeakGroupMatching(PipelineModule):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=None, hypothesis_sample_match_id=None,
                 search_type="TheoreticalGlycanComposition",
                 match_tolernace=2e-5, mass_shift_map=None,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        session = self.manager.session()
        no_shift = get_or_create(session, MassShift, mass=0.0, name=u"NoShift")
        lcms_database = PeakGroupDatabase(observed_ions_path)
        lcms_session = lcms_database.session()
        sample_run = lcms_session.query(SampleRun).get(sample_run_id or 1)

        hypothesis_sample_match = get_or_create(session, HypothesisSampleMatch,
                                                id=hypothesis_sample_match_id,
                                                target_hypothesis_id=hypothesis_id,
                                                sample_run_name=sample_run.name)
        session.add(hypothesis_sample_match)
        session.commit()

        self.hypothesis_sample_match_id = hypothesis_sample_match.id
        if mass_shift_map is None:
            mass_shift_map = {no_shift: 1}
        else:
            mass_shift_map[no_shift] = 1
        self.mass_shift_map = mass_shift_map
        self.n_processes = n_processes
        self.lcms_database = lcms_database
        self.search_type = TheoreticalCompositionMap[search_type]

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

    def run(self):
        session = self.manager.session()

        task_fn = callable

        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_ids(), chunksize=500):
                pass
        else:
            for res in itertools.imap(task_fn, self.stream_ids()):
                pass
        session.commit()
        session.close()
