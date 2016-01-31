from glycresoft_sqlalchemy.data_model import (
    PipelineModule, MS2GlycopeptideHypothesisSampleMatch,
    SampleRun, Hypothesis, HypothesisSampleMatch)
from glycresoft_sqlalchemy.matching.glycopeptide.fragment_matching import IonMatching
from glycresoft_sqlalchemy.scoring import score_spectrum_matches, target_decoy


class GlycopeptideFragmentMatchingPipeline(PipelineModule):

    def __init__(self, database_path, observed_ions_path,
                 target_hypothesis_id, decoy_hypothesis_id, hypothesis_sample_match_id=None,
                 sample_run_id=None, scorer=None, ms1_tolerance=1e-5,
                 ms2_tolerance=2e-5, intensity_threshold=150., n_processes=4,
                 **kwargs):
        if sample_run_id is None:
            sample_run_id = 1
        self.manager = self.manager_type(database_path)
        self.observed_ions_path = observed_ions_path
        self.target_hypothesis_id = target_hypothesis_id
        self.decoy_hypothesis_id = decoy_hypothesis_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.sample_run_id = sample_run_id
        self.ms1_tolerance = ms1_tolerance
        self.ms2_tolerance = ms2_tolerance
        self.intensity_threshold = intensity_threshold
        self.n_processes = n_processes
        self.options = kwargs
        self.scorer = scorer

    def prepare_hypothesis_sample_match(self):
        session = self.manager()
        if self.hypothesis_sample_match_id is None:
            ion_session = self.manager_type(self.observed_ions_path)
            sample_run = ion_session.query(SampleRun).get(self.sample_run_id)
            target_hypothesis = session.query(Hypothesis).get(self.target_hypothesis_id)

            count = session.query(MS2GlycopeptideHypothesisSampleMatch).filter_by(
                target_hypothesis_id=target_hypothesis.id, sample_run_name=sample_run.name).count()
            if count > 1:
                count_part = " (%d)" % (count)
            else:
                count_part = ""
            hsm = MS2GlycopeptideHypothesisSampleMatch(
                target_hypothesis_id=self.target_hypothesis_id,
                decoy_hypothesis_id=self.decoy_hypothesis_id,
                sample_run_name=sample_run.name,
                name=self.options.get("hypothesis_sample_match_name",
                    "%s @ %s%s" % (
                        target_hypothesis.name, sample_run.name, count_part)))

            session.add(hsm)
            session.commit()
            self.hypothesis_sample_match_id = hsm.id
            hsm.parameters.update({
                "observed_ions_path": self.observed_ions_path,
                "target_hypothesis_id": self.target_hypothesis_id,
                "decoy_hypothesis_id": self.decoy_hypothesis_id,
                "hypothesis_sample_match_id": self.hypothesis_sample_match_id,
                "sample_run_id": self.sample_run_id,
                "ms1_tolerance": self.ms1_tolerance,
                "ms2_tolerance": self.ms2_tolerance,
                "intensity_threshold": self.intensity_threshold,
                "n_processes": self.n_processes,
                "options": self.options,
            })
        else:
            hsm = session.query(
                MS2GlycopeptideHypothesisSampleMatch).get(self.hypothesis_sample_match_id)
            # Do clean up of previous matches?
            hsm.parameters.update({
                "observed_ions_path": self.observed_ions_path,
                "target_hypothesis_id": self.target_hypothesis_id,
                "decoy_hypothesis_id": self.decoy_hypothesis_id,
                "hypothesis_sample_match_id": self.hypothesis_sample_match_id,
                "sample_run_id": self.sample_run_id,
                "ms1_tolerance": self.ms1_tolerance,
                "ms2_tolerance": self.ms2_tolerance,
                "intensity_threshold": self.intensity_threshold,
                "n_processes": self.n_processes,
                "options": self.options,
            })
        session.add(hsm)
        session.commit()

    def do_target_matching(self):
        task = IonMatching(
            self.database_path,
            self.target_hypothesis_id,
            self.observed_ions_path,
            observed_ions_type='db',
            sample_run_id=self.sample_run_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            intensity_threshold=self.intensity_threshold,
            ms1_tolerance=self.ms1_tolerance,
            ms2_tolerance=self.ms2_tolerance,
            n_processes=self.n_processes
            )
        task.start()

    def do_decoy_matching(self):
        task = IonMatching(
            self.database_path,
            self.decoy_hypothesis_id,
            self.observed_ions_path,
            observed_ions_type='db',
            sample_run_id=self.sample_run_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            intensity_threshold=self.intensity_threshold,
            ms1_tolerance=self.ms1_tolerance,
            ms2_tolerance=self.ms2_tolerance,
            n_processes=self.n_processes
            )
        task.start()

    def do_matching(self):
        self.do_target_matching()
        self.do_decoy_matching()

    def do_target_spectrum_assignment(self):
        task = score_spectrum_matches.SimpleSpectrumAssignment(
            self.database_path,
            hypothesis_id=self.target_hypothesis_id,
            scorer=self.scorer,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            n_processes=self.n_processes)
        task.start()

    def do_decoy_spectrum_assignment(self):
        task = score_spectrum_matches.SimpleSpectrumAssignment(
            self.database_path,
            hypothesis_id=self.decoy_hypothesis_id,
            scorer=self.scorer,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            n_processes=self.n_processes)
        task.start()

    def do_spectrum_assignment(self):
        self.do_target_spectrum_assignment()
        self.do_decoy_spectrum_assignment()

    def do_target_decoy_fdr_estimation(self):
        tda = target_decoy.TargetDecoyAnalyzer(
            self.database_path,
            self.target_hypothesis_id,
            self.decoy_hypothesis_id,
            self.hypothesis_sample_match_id)
        tda.start()

    def do_store_observed_data(self):
        session = self.manager()
        ion_session = self.manager_type(self.observed_ions_path)
        sample_run = ion_session.query(SampleRun).get(self.sample_run_id)
        hypothesis_sample_match = session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
        is_loaded = session.query(SampleRun).filter_by(uuid=sample_run.uuid).first() is None
        if not is_loaded:
            hypothesis_sample_match.copy_tandem_sample_run(sample_run)

    def run(self):
        self.prepare_hypothesis_sample_match()
        self.do_matching()
        self.do_spectrum_assignment()
        self.do_target_decoy_fdr_estimation()
