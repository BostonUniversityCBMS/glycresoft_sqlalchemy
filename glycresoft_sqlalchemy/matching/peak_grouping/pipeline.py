import logging
try:
    logger = logging.getLogger("peak_grouping")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

from sqlalchemy.ext.baked import bakery


from glycresoft_sqlalchemy.data_model import (
    PipelineModule,
    SampleRun, Hypothesis,
    TheoreticalCompositionMap, HypothesisSampleMatch,
    PeakGroupMatch)

from glycresoft_sqlalchemy.utils.database_utils import get_or_create
from glycresoft_sqlalchemy.utils import pickle

from .grouper import Decon2LSPeakGrouper
from .mass_shift_offset_matching import PeakGroupMatching, BatchPeakGroupMatching
from .classification import PeakGroupMassShiftJoiningClassifier, PeakGroupClassification


query_oven = bakery()
PeakGroupType = PeakGroupMatch


class LCMSPeakClusterSearch(PipelineModule):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, hypothesis_sample_match_id=None,
                 search_type="TheoreticalGlycanComposition",
                 match_tolerance=2e-5, mass_shift_map=None,
                 regression_parameters=None, minimum_mass=None,
                 maximum_mass=None,
                 n_processes=4, **kwargs):
        self.manager = self.manager_type(database_path)
        self.observed_ions_path = observed_ions_path
        self.hypothesis_id = hypothesis_id
        self.sample_run_id = sample_run_id
        self.grouping_error_tolerance = grouping_error_tolerance
        self.minimum_scan_count = minimum_scan_count
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.search_type = TheoreticalCompositionMap[search_type]
        self.match_tolerance = match_tolerance
        self.mass_shift_map = mass_shift_map
        self.n_processes = n_processes
        self.regression_parameters = regression_parameters
        self.minimum_mass = minimum_mass
        self.maximum_mass = maximum_mass
        self.options = kwargs
        session = self.manager.session()
        hypothesis = session.query(Hypothesis).get(self.hypothesis_id)
        self.hypothesis_sample_match_type = HypothesisSampleMatch.hierarchy_root[hypothesis.__class__]
        session.close()

    def do_matching(self):
        matcher = BatchPeakGroupMatching(
            self.database_path, self.observed_ions_path, self.hypothesis_id,
            self.sample_run_id, self.hypothesis_sample_match_id,
            self.search_type, self.match_tolerance, self.mass_shift_map,
            self.n_processes)

        if not self.options.get("skip_matching", False):
            matcher.start()

    def do_grouping(self):
        grouper = Decon2LSPeakGrouper(
            self.observed_ions_path, self.sample_run_id, self.grouping_error_tolerance,
            minimum_mass=self.minimum_mass, maximum_mass=self.maximum_mass,
            minimum_scan_count=self.minimum_scan_count, n_processes=self.n_processes)
        if not self.options.get("skip_grouping", False):
            grouper.start()
        return grouper

    def prepare_hypothesis_sample_match(self, session, sample_name):
        hypothesis_sample_match, created = get_or_create(
            session, self.hypothesis_sample_match_type,
            id=self.hypothesis_sample_match_id,
            target_hypothesis_id=self.hypothesis_id)
        if self.hypothesis_sample_match_id is not None:
            if not self.options.get("skip_matching", False):
                hypothesis_sample_match.peak_group_matches.delete(synchronize_session=False)
                session.commit()
                logger.info("Cleared? %r", hypothesis_sample_match.peak_group_matches.count())
            else:
                hypothesis_sample_match.peak_group_matches.filter(
                    ~PeakGroupType.matched).delete(synchronize_session=False)
                session.commit()
                logger.info("Cleared? %r", hypothesis_sample_match.peak_group_matches.filter(
                    ~PeakGroupType.matched).count())
        self.options.setdefault("hypothesis_sample_match_name", "{} @ {}_ms1".format(
            session.query(Hypothesis).get(self.hypothesis_id).name,
            sample_name))

        hypothesis_sample_match.name = HypothesisSampleMatch.make_unique_name(
            session, self.options["hypothesis_sample_match_name"])

        hypothesis_sample_match.parameters = dict(
            observed_ions_path=self.observed_ions_path,
            hypothesis_id=self.hypothesis_id,
            sample_run_id=self.sample_run_id,
            grouping_error_tolerance=self.grouping_error_tolerance,
            minimum_scan_count=self.minimum_scan_count,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            search_type=self.search_type,
            match_tolerance=self.match_tolerance,
            mass_shift_map=self.mass_shift_map,
            n_processes=self.n_processes,
            minimum_mass=self.minimum_mass,
            maximum_mass=self.maximum_mass
        )
        session.add(hypothesis_sample_match)
        session.commit()

        self.hypothesis_sample_match_id = hypothesis_sample_match.id
        return hypothesis_sample_match

    def do_classification(self):
        classifier = PeakGroupMassShiftJoiningClassifier(
            self.database_path, self.observed_ions_path,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            sample_run_id=self.sample_run_id,
            n_processes=self.n_processes)

        classifier.start()
        return classifier

    def run(self):
        grouper = self.do_grouping()

        database_manager = self.manager
        session = database_manager.session()
        sample_name = grouper.manager.session().query(SampleRun).get(self.sample_run_id).name

        hypothesis_sample_match = self.prepare_hypothesis_sample_match(session, sample_name)

        self.do_matching()

        classifier = self.do_classification()
        hypothesis_sample_match.parameters['classifier'] = pickle.dumps(classifier.classifier)
        session.add(hypothesis_sample_match)
        session.commit()
