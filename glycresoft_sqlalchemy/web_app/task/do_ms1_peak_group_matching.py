from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptideComposition,
    MS1GlycopeptideHypothesisSampleMatch, SampleRun, make_transient, MassShift)

from glycresoft_sqlalchemy.matching.peak_grouping import (
    LCMSPeakClusterSearch)

import os
import pickle

from .task_process import NullPipe, Message, Task


class CommunicativeLCMSPeakClusterSearch(LCMSPeakClusterSearch):
    def __init__(self, database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1, hypothesis_sample_match_id=None,
                 search_type="TheoreticalGlycanComposition",
                 match_tolerance=2e-5, mass_shift_map=None, regression_parameters=None,
                 comm=NullPipe(), **kwargs):
        kwargs.setdefault("n_processes", 4)
        dbm = self.manager_type(database_path)
        session = dbm.session()
        new_map = {}
        for k, v in mass_shift_map.items():
            s = session.query(MassShift).get(k)
            new_map[s] = v
        mass_shift_map = new_map
        session.close()

        super(CommunicativeLCMSPeakClusterSearch, self).__init__(
            database_path, observed_ions_path, hypothesis_id, sample_run_id,
            grouping_error_tolerance, minimum_scan_count, hypothesis_sample_match_id,
            search_type, match_tolerance, mass_shift_map, regression_parameters, **kwargs)
        self.comm = comm

    def run(self):
        self.comm.send(Message("Begin grouping peaks", "update"))
        grouper = self.do_grouping()

        database_manager = self.manager
        session = database_manager.session()
        sample_name = grouper.manager.session().query(SampleRun).get(self.sample_run_id).name

        hypothesis_sample_match = self.prepare_hypothesis_sample_match(session, sample_name)

        self.comm.send(Message("Begin matching masses", "update"))
        self.do_matching()

        self.comm.send(Message("Begin peak group scoring", "update"))
        classifier = self.do_classification()
        hypothesis_sample_match.parameters['classifier'] = pickle.dumps(classifier.classifier)

        self.comm.send(Message("Peak Cluster Search Complete", "update"))

        self.comm.send(Message(hypothesis_sample_match.to_json(), "new-hypothesis-sample-match"))


def taskmain(*args, **kwargs):
    return CommunicativeLCMSPeakClusterSearch(*args, **kwargs).start()


class LCMSSearchTask(Task):
    def __init__(self,
                 database_path, observed_ions_path, hypothesis_id,
                 sample_run_id=1, grouping_error_tolerance=8e-5,
                 minimum_scan_count=1,
                 search_type="TheoreticalGlycanComposition",
                 match_tolerance=2e-5, mass_shift_map=None, regression_parameters=None,
                 callback=lambda: 0, **kwargs):
        args = (
            database_path, observed_ions_path, hypothesis_id,
            sample_run_id, grouping_error_tolerance,
            minimum_scan_count, None,
            search_type, match_tolerance, mass_shift_map, regression_parameters)
        name = "LC-MS Search {} @ {}".format(hypothesis_id, os.path.basename(observed_ions_path))
        kwargs.setdefault('name', name)
        super(LCMSSearchTask, self).__init__(taskmain, args, callback, **kwargs)
