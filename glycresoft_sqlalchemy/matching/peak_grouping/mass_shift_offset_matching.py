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
    SampleRun,
    TheoreticalCompositionMap, MassShift, HypothesisSampleMatch,
    PeakGroupDatabase, PeakGroupMatch,
)

from .common import ppm_error

from glycresoft_sqlalchemy.utils.database_utils import get_or_create, toggle_indices

TPeakGroupMatch = PeakGroupMatch.__table__


query_oven = bakery()


def yield_ids(session, theoretical_type, hypothesis_id, chunk_size=1000, filter=lambda q: q):
    base_query = filter(session.query(theoretical_type.id).filter(theoretical_type.from_hypothesis(hypothesis_id)))
    chunk = []

    for item in base_query:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
    yield chunk


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
                    # This query is recompiled every time, using up the majority of the time spent
                    # per call. Look into cached compilation with bind parameters and the BakedQuery
                    # pattern.
                    for mass_match in observed_ions_manager.ppm_match_tolerance_search(
                            total_mass, matching_tolerance, sample_run_id):
                        mass_error = ppm_error(mass_match.weighted_monoisotopic_mass, total_mass)
                        matches.append((mass_match, mass_error, mass_shift, shift_count))
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
        return len(params)


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
                # This query is recompiled every time, using up the majority of the time spent
                # per call. Look into cached compilation with bind parameters and the BakedQuery
                # pattern.
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
        toggler.create()
        session.close()


class BatchPeakGroupMatching(PeakGroupMatching):

    def __init__(self, *args, **kwargs):
        super(BatchPeakGroupMatching, self).__init__(*args, **kwargs)

    def stream_ids(self):
        session = self.manager.session()
        try:
            for gids in yield_ids(session, self.theoretical_type, self.hypothesis_id):
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

        task_fn = self.prepare_task_fn()
        counter = 0
        toggler = toggle_indices(session, PeakGroupMatch)
        toggler.drop()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_ids(), chunksize=500):
                counter += res
                if counter % 1000 == 0:
                    logger.info("%d masses searched", counter)

        else:
            for res in itertools.imap(task_fn, self.stream_ids()):
                counter += res
                if counter % 1000 == 0:
                    logger.info("%d masses searched", counter)
        toggler.create()
        session.close()