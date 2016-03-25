import os
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

from glycresoft_sqlalchemy.data_model import (
    PipelineModule,
    SampleRun, Decon2LSPeakGroup,
    TheoreticalCompositionMap, MassShift, HypothesisSampleMatch,
    PeakGroupDatabase, PeakGroupMatch, JointPeakGroupMatch
)

from .common import ppm_error

from glycresoft_sqlalchemy.utils import get_scale
from glycresoft_sqlalchemy.utils.database_utils import get_or_create, toggle_indices

TPeakGroupMatch = PeakGroupMatch.__table__


SCALE = get_scale()

query_oven = bakery()


def yield_ids(session, theoretical_type, hypothesis_id, chunk_size=100, filter=lambda q: q):
    base_query = filter(session.query(theoretical_type.id).filter(
        theoretical_type.from_hypothesis(hypothesis_id))).all()
    last = 0
    final = len(base_query)
    while 1:
        next_chunk = base_query[last:(last + chunk_size)]
        if last <= final:
            yield next_chunk
            last += chunk_size
        else:
            break


def yield_peak_group_ids(session, group_type, sample_run_id, chunk_size=100, filter=lambda q: q):
    base_query = filter(session.query(
        group_type.id).filter(group_type.sample_run_id == sample_run_id)).all()
    last = 0
    final = len(base_query)
    while 1:
        next_chunk = base_query[last:(last + chunk_size)]
        if last <= final:
            yield next_chunk
            last += chunk_size
        else:
            break


def batch_match_theoretical_composition(
        peak_group_ids, search_type, database_manager, observed_ions_manager,
        matching_tolerance, mass_shift_map, hypothesis_id,
        hypothesis_sample_match_id):
    session = database_manager()
    ions_session = observed_ions_manager()
    params = []
    try:
        for peak_group_id in peak_group_ids:
            peak_group = ions_session.query(Decon2LSPeakGroup).get(peak_group_id[0])
            base_mass = peak_group.weighted_monoisotopic_mass
            matches = []
            for mass_shift, count_range in mass_shift_map.items():
                shift = mass_shift.mass
                for shift_count in range(1, count_range + 1):
                    total_mass = base_mass + (shift * shift_count)
                    for mass_match in search_type.ppm_error_tolerance_search(
                            session=session,
                            mass=total_mass, tolerance=matching_tolerance,
                            hypothesis_id=hypothesis_id):
                        mass_error = ppm_error(mass_match.calculated_mass, total_mass)
                        matches.append((mass_match.id, mass_error, mass_shift, shift_count))

            for mass_match_id, mass_error, mass_shift, shift_count in matches:
                case = {
                    "hypothesis_sample_match_id": hypothesis_sample_match_id,
                    "theoretical_match_type": search_type.__name__,
                    "theoretical_match_id": mass_match_id,
                    "charge_state_count": peak_group.charge_state_count,
                    "scan_count": peak_group.scan_count,
                    "first_scan_id": peak_group.first_scan_id,
                    "last_scan_id": peak_group.last_scan_id,
                    "ppm_error": mass_error,
                    "scan_density": peak_group.scan_density,
                    "weighted_monoisotopic_mass": peak_group.weighted_monoisotopic_mass,
                    "total_volume": peak_group.total_volume,
                    "average_a_to_a_plus_2_ratio": peak_group.average_a_to_a_plus_2_ratio,
                    "a_peak_intensity_error": peak_group.a_peak_intensity_error,
                    "centroid_scan_estimate": peak_group.centroid_scan_estimate,
                    "centroid_scan_error": peak_group.centroid_scan_error,
                    "average_signal_to_noise": peak_group.average_signal_to_noise,
                    "peak_data": peak_group.peak_data,
                    "matched": True,
                    "mass_shift_type": mass_shift.id,
                    "mass_shift_count": shift_count,
                    "peak_group_id": peak_group.id,
                }
                params.append(case)
        session.bulk_insert_mappings(PeakGroupMatch, params)
        session.commit()

    except Exception, e:
        logger.exception("An exception occurred in match_peak_group, %r", locals(), exc_info=e)
        raise e
    finally:
        session.close()
        return len(peak_group_ids)


def batch_match_peak_group(search_ids, search_type, database_manager, observed_ions_manager,
                           matching_tolerance, mass_shift_map, sample_run_id,
                           hypothesis_sample_match_id):
    session = database_manager.session()
    params = []
    try:
        for search_id in search_ids:
            search_target = session.query(search_type).get(search_id)
            base_mass = search_target.calculated_mass
            matches = []
            for mass_shift, count_range in mass_shift_map.items():
                shift = mass_shift.mass
                for shift_count in range(1, count_range + 1):
                    total_mass = base_mass + (shift * shift_count)

                    for mass_match in observed_ions_manager.ppm_match_tolerance_search(
                            total_mass, matching_tolerance, sample_run_id):
                        mass_error = ppm_error(mass_match.weighted_monoisotopic_mass, total_mass)
                        matches.append((mass_match, mass_error, mass_shift, shift_count))
            for mass_match, mass_error, mass_shift, shift_count in matches:
                case = {
                    "hypothesis_sample_match_id": hypothesis_sample_match_id,
                    "theoretical_match_type": search_type.__name__,
                    "theoretical_match_id": search_id[0],
                    "charge_state_count": mass_match.charge_state_count,
                    "scan_count": mass_match.scan_count,
                    "first_scan_id": mass_match.first_scan_id,
                    "last_scan_id": mass_match.last_scan_id,
                    "ppm_error": mass_error,
                    "scan_density": mass_match.scan_density,
                    "weighted_monoisotopic_mass": mass_match.weighted_monoisotopic_mass,
                    "total_volume": mass_match.total_volume,
                    "average_a_to_a_plus_2_ratio": mass_match.average_a_to_a_plus_2_ratio,
                    "a_peak_intensity_error": mass_match.a_peak_intensity_error,
                    "centroid_scan_estimate": mass_match.centroid_scan_estimate,
                    "centroid_scan_error": mass_match.centroid_scan_error,
                    "average_signal_to_noise": mass_match.average_signal_to_noise,
                    "peak_data": mass_match.peak_data,
                    "matched": True,
                    "mass_shift_type": mass_shift.id,
                    "mass_shift_count": shift_count,
                    "peak_group_id": mass_match.id,
                }
                params.append(case)
        session.bulk_insert_mappings(PeakGroupMatch, params)
        session.commit()
    except Exception, e:
        logger.exception("An exception occurred in match_peak_group, %r", locals(), exc_info=e)
        raise e
    finally:
        session.close()
        return len(search_ids)


def match_peak_group(search_id, search_type, database_manager, observed_ions_manager,
                     matching_tolerance, mass_shift_map, sample_run_id, hypothesis_sample_match_id):
    session = database_manager.session()
    try:
        search_target = session.query(search_type).get(search_id)

        base_mass = search_target.calculated_mass
        matches = []
        for mass_shift, count_range in mass_shift_map.items():
            shift = mass_shift.mass
            for shift_count in range(1, count_range + 1):
                total_mass = base_mass + (shift * shift_count)

                for mass_match in observed_ions_manager.ppm_match_tolerance_search(
                        total_mass, matching_tolerance, sample_run_id):
                    mass_error = ppm_error(mass_match.weighted_monoisotopic_mass, total_mass)
                    matches.append((mass_match, mass_error, mass_shift, shift_count))
        params = []
        for mass_match, mass_error, mass_shift, shift_count in matches:
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
                "weighted_monoisotopic_mass": mass_match.weighted_monoisotopic_mass,
                "total_volume": mass_match.total_volume,
                "average_a_to_a_plus_2_ratio": mass_match.average_a_to_a_plus_2_ratio,
                "a_peak_intensity_error": mass_match.a_peak_intensity_error,
                "centroid_scan_estimate": mass_match.centroid_scan_estimate,
                "centroid_scan_error": mass_match.centroid_scan_error,
                "average_signal_to_noise": mass_match.average_signal_to_noise,
                "peak_data": mass_match.peak_data,
                "matched": True,
                "mass_shift_type": mass_shift.id,
                "mass_shift_count": shift_count,
                "peak_group_id": mass_match.id,
            }
            params.append(case)

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
        hypothesis_sample_match.parameters = hypothesis_sample_match.parameters or {}
        hypothesis_sample_match.parameters['mass_shift_map'] = mass_shift_map
        session.add(hypothesis_sample_match)
        session.commit()

        self.hypothesis_sample_match_id = hypothesis_sample_match.id
        if mass_shift_map is None:
            mass_shift_map = {no_shift: 1}
        else:
            mass_shift_map[no_shift] = 1
        logger.info("Mass Shift Map: %r", mass_shift_map)
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.mass_shift_map = mass_shift_map
        self.match_tolerance = match_tolerance
        self.n_processes = n_processes
        self.lcms_database = lcms_database
        self.search_type = TheoreticalCompositionMap.get(search_type, search_type)
        if not isinstance(self.search_type, type):
            raise TypeError("{} is not a type".format(self.search_type))
        session.close()

    def stream_ids(self):
        session = self.manager.session()
        try:
            for gid in session.query(self.search_type.id).filter(
                    self.search_type.from_hypothesis(self.hypothesis_id)):
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
        toggler = toggle_indices(session, PeakGroupMatch)
        toggler.drop()
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
            pool.terminate()
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
        logger.info("Search Complete.")
        toggler.create()
        session.close()


class BatchPeakGroupMatching(PeakGroupMatching):

    def __init__(self, *args, **kwargs):
        super(BatchPeakGroupMatching, self).__init__(*args, **kwargs)

    def stream_ids(self, chunk_size=400):
        chunk_size *= SCALE
        session = self.manager.session()
        try:
            for gids in yield_ids(session, self.search_type, self.hypothesis_id, chunk_size=chunk_size):
                yield gids
        except Exception, e:
            logger.info("An exception occurred while streaming ids", exc_info=e)
            raise e
        finally:
            session.close()

    def prepare_task_fn(self):
        fn = functools.partial(
            batch_match_peak_group,
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

        counter = 0
        last = 0
        step = 1000

        task_fn = self.prepare_task_fn()
        toggler = toggle_indices(session, PeakGroupMatch)
        toggler.drop()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_ids()):
                counter += res
                if counter > (last + step):
                    last += step
                    logger.info("%d masses searched", counter)

        else:
            for res in itertools.imap(task_fn, self.stream_ids()):
                counter += res
                if counter > (last + step):
                    last += step
                    logger.info("%d masses searched", counter)
        logger.info("Search Complete.")
        toggler.create()
        session.close()


class BatchPeakGroupMatchingSearchGroups(PeakGroupMatching):
    def __init__(self, *args, **kwargs):
        super(BatchPeakGroupMatchingSearchGroups, self).__init__(*args, **kwargs)

    def stream_ids(self, chunk_size=200):
        chunk_size *= SCALE
        session = self.lcms_database()
        try:
            for gids in yield_peak_group_ids(session, Decon2LSPeakGroup, self.sample_run_id, chunk_size=chunk_size):
                yield gids
        except Exception, e:
            logger.info("An exception occurred while streaming ids", exc_info=e)
            raise e
        finally:
            session.close()

    def prepare_task_fn(self):
        fn = functools.partial(
            batch_match_theoretical_composition,
            search_type=self.search_type,
            database_manager=self.manager,
            observed_ions_manager=self.lcms_database,
            matching_tolerance=self.match_tolerance,
            mass_shift_map=self.mass_shift_map,
            hypothesis_id=self.hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id)
        return fn

    def run(self):
        session = self.manager.session()

        counter = 0
        last = 0
        step = 1000

        task_fn = self.prepare_task_fn()
        toggler = toggle_indices(session, PeakGroupMatch)
        toggler.drop()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_ids()):
                counter += res
                if counter > (last + step):
                    last += step
                    logger.info("%d masses searched", counter)
            pool.terminate()
        else:
            for res in itertools.imap(task_fn, self.stream_ids()):
                counter += res
                if counter > (last + step):
                    last += step
                    logger.info("%d masses searched", counter)
        logger.info("Search Complete.")
        toggler.create()
        session.close()


BatchPeakGroupMatching = BatchPeakGroupMatchingSearchGroups
