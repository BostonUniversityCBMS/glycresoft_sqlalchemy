from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch,
    MS2GlycopeptideHypothesisSampleMatch, SampleRun, HypothesisSampleMatch)

from glycresoft_sqlalchemy.matching.glycopeptide.fragment_matching import (
    IonMatching, ms1_tolerance_default, ms2_tolerance_default)
from glycresoft_sqlalchemy.matching.glycopeptide import pipeline
# from glycresoft_sqlalchemy.matching.matching import IonMatching, ms1_tolerance_default, ms2_tolerance_default
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from glycresoft_sqlalchemy.scoring import target_decoy, score_spectrum_matches
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder import (
    search_space_builder as search_space_builder_lib, make_decoys)
import os

from .task_process import NullPipe, Message, Task


class CommunicativeGlycopeptideFragmentMatchingPipeline(pipeline.GlycopeptideFragmentMatchingPipeline):
    def __init__(self, database_path, observed_ions_path,
                 target_hypothesis_id, decoy_hypothesis_id, hypothesis_sample_match_id=None,
                 sample_run_id=None, scorer=None, ms1_tolerance=1e-5,
                 ms2_tolerance=2e-5, intensity_threshold=150., n_processes=4,
                 comm=NullPipe(), **kwargs):
        super(CommunicativeGlycopeptideFragmentMatchingPipeline, self).__init__(
            database_path=database_path, observed_ions_path=observed_ions_path,
            target_hypothesis_id=target_hypothesis_id, decoy_hypothesis_id=decoy_hypothesis_id,
            hypothesis_sample_match_id=hypothesis_sample_match_id, sample_run_id=sample_run_id,
            scorer=scorer, ms1_tolerance=ms1_tolerance, ms2_tolerance=ms2_tolerance,
            intensity_threshold=intensity_threshold, n_processes=n_processes, **kwargs)
        self.comm = comm

    def do_target_matching(self):
        self.comm.send(
            Message("Begin Matching for target %d" % self.target_hypothesis_id, "update"))
        super(CommunicativeGlycopeptideFragmentMatchingPipeline, self).do_target_matching()

    def do_decoy_matching(self):
        self.comm.send(
            Message("Begin Matching for decoy %d" % self.decoy_hypothesis_id, "update"))
        super(CommunicativeGlycopeptideFragmentMatchingPipeline, self).do_decoy_matching()

    def do_target_decoy_fdr_estimation(self):
        self.comm.send(Message("Begin TDA", "update"))
        super(CommunicativeGlycopeptideFragmentMatchingPipeline, self).do_target_decoy_fdr_estimation()


def taskmain(
        database_path, observed_ions_path, target_hypothesis_id=None,
        decoy_hypothesis_id=None, source_hypothesis_sample_match_id=None,
        observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default,
        comm=NullPipe(),
        **kwargs):

    manager = DatabaseManager(database_path)
    manager.initialize()
    if target_hypothesis_id is None:
        if source_hypothesis_sample_match_id is None:
            raise Exception("No target hypotesis or hypotesis sample match given")

        session = manager.session()
        source_hsm = session.query(HypothesisSampleMatch).get(source_hypothesis_sample_match_id)
        search_space_builder = search_space_builder_lib.constructs[source_hsm.__class__]
        builder = search_space_builder.from_hypothesis_sample_match(
            database_path, source_hsm.id, n_processes=kwargs.get("n_processes", 4))
        target_hypothesis_id = builder.start()
        builder = make_decoys.BatchingDecoySearchSpaceBuilder(
            database_path, hypothesis_ids=[target_hypothesis_id], n_processes=kwargs.get("n_processes", 4))
        decoy_hypothesis_id = builder.start()
        decoy_hypothesis_id = decoy_hypothesis_id[0]
    else:
        target_hypothesis = manager.query(Hypothesis).get(target_hypothesis_id)

        decoy_hypothesis_id = target_hypothesis.decoy_hypothesis_id()
        if decoy_hypothesis_id is None:
            builder = make_decoys.BatchingDecoySearchSpaceBuilder(
                database_path, hypothesis_ids=[target_hypothesis_id], n_processes=kwargs.get("n_processes", 4))
            decoy_hypothesis_id = builder.start()
            decoy_hypothesis_id = decoy_hypothesis_id[0]

    if observed_ions_type == 'bupid_yaml' and observed_ions_path[-3:] != '.db':
        comm.send(Message("Converting %s to db" % observed_ions_path))
        parser = BUPIDMSMSYamlParser(observed_ions_path, manager.bridge_address())
        observed_ions_path = parser.manager.path
        observed_ions_type = 'db'
        sample_name = parser.sample_run_name
    else:
        sample_name = DatabaseManager(
            observed_ions_path).session().query(SampleRun).get(sample_run_id).name
    session = manager.session()
    target_name = manager.session().query(
                Hypothesis).get(target_hypothesis_id).name
    # if decoy_hypothesis_id is not None:
    #     hsm = MS2GlycopeptideHypothesisSampleMatch(
    #         name="{} @ {}".format(target_name, sample_name),
    #         target_hypothesis_id=target_hypothesis_id,
    #         decoy_hypothesis_id=decoy_hypothesis_id,
    #         sample_run_name=sample_name
    #     )
    #     session.add(hsm)
    #     session.commit()
    #     comm.send(Message("Created %r" % hsm))
    #     hsm_id = hsm.id
    # else:
    #     comm.send(Message("No MS2GlycopeptideHypothesisSampleMatch used"))
    #     hsm_id = None
    job = CommunicativeGlycopeptideFragmentMatchingPipeline(
        database_path, observed_ions_path, target_hypothesis_id=target_hypothesis_id,
        decoy_hypothesis_id=decoy_hypothesis_id, hypothesis_sample_match_id=None,
        sample_run_id=sample_run_id, scorer=None, ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance, intensity_threshold=float(kwargs.get("intensity_threshold", 150.)),
        n_processes=kwargs.get('n_processes', 4), comm=comm)
    job.start()
    hsm_id = job.hypothesis_sample_match_id
    hsm = session.query(HypothesisSampleMatch).get(hsm_id)
    # comm.send(Message("Begin Matching for target %d" % target_hypothesis_id, "update"))
    # job = IonMatching(
    #     database_path,
    #     hypothesis_id=target_hypothesis_id,
    #     observed_ions_path=observed_ions_path,
    #     observed_ions_type=observed_ions_type,
    #     hypothesis_sample_match_id=hsm_id,
    #     sample_run_id=sample_run_id,
    #     ms1_tolerance=ms1_tolerance,
    #     ms2_tolerance=ms2_tolerance,
    #     n_processes=kwargs.get("n_processes", 4))
    # job.start()
    # comm.send(Message("End Matching for target %d. %s, %r" % (target_hypothesis_id, job.status, job.error)))
    # if job.error is not None:
    #     comm.send(Message(job.error, 'error'))
    # if decoy_hypothesis_id is None:
    #     return
    # comm.send(Message("Begin Matching for decoy %d" % decoy_hypothesis_id, "update"))
    # job = IonMatching(
    #     database_path,
    #     hypothesis_id=decoy_hypothesis_id,
    #     observed_ions_path=observed_ions_path,
    #     observed_ions_type=observed_ions_type,
    #     hypothesis_sample_match_id=hsm_id,
    #     sample_run_id=sample_run_id,
    #     ms1_tolerance=ms1_tolerance,
    #     ms2_tolerance=ms2_tolerance,
    #     n_processes=kwargs.get("n_processes", 4))
    # job.start()

    # job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, target_hypothesis_id, hsm_id)
    # job.start()

    # job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, decoy_hypothesis_id, hsm_id)
    # job.start()

    # comm.send(Message("Begin TDA", "update"))
    # job = target_decoy.TargetDecoyAnalyzer(database_path, target_hypothesis_id, decoy_hypothesis_id, hsm_id)
    # job.start()
    if hsm_id is not None:
        comm.send(Message(hsm.to_json(), "new-hypothesis-sample-match"))
    return


class TandemMSGlycoproteomicsSearchTask(Task):
    def __init__(
            self, database_path, observed_ions_path, target_hypothesis_id,
            decoy_hypothesis_id, source_hypothesis_sample_match_id,
            observed_ions_type, sample_run_id,
            ms1_tolerance, ms2_tolerance, callback=lambda: 0, **kwargs):
        args = (
            database_path, observed_ions_path, target_hypothesis_id,
            decoy_hypothesis_id, source_hypothesis_sample_match_id,
            observed_ions_type, sample_run_id,
            ms1_tolerance, ms2_tolerance)
        name = "Tandem MS Glycoproteomics Search {} @ {}".format(target_hypothesis_id, os.path.basename(observed_ions_path))
        kwargs.setdefault('name', name)
        super(TandemMSGlycoproteomicsSearchTask, self).__init__(taskmain, args, callback, **kwargs)
