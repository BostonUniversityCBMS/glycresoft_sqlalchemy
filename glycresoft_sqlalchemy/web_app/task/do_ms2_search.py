from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch,
    MS2GlycopeptideHypothesisSampleMatch, SampleRun)

from glycresoft_sqlalchemy.matching.matching import IonMatching, ms1_tolerance_default, ms2_tolerance_default
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from glycresoft_sqlalchemy.scoring import target_decoy, score_spectrum_matches

import os

from .task_process import NullPipe, Message, Task


def taskmain(
        database_path, observed_ions_path, target_hypothesis_id=None,
        decoy_hypothesis_id=None, observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default,
        comm=NullPipe(),
        **kwargs):
    if target_hypothesis_id is None:
        target_hypothesis_id = 1

    manager = DatabaseManager(database_path)
    manager.initialize()

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
    if decoy_hypothesis_id is not None:
        hsm = MS2GlycopeptideHypothesisSampleMatch(
            name="{} @ {}".format(target_name, sample_name),
            target_hypothesis_id=target_hypothesis_id,
            decoy_hypothesis_id=decoy_hypothesis_id,
            sample_run_name=sample_name
        )
        session.add(hsm)
        session.commit()
        comm.send(Message("Created %r" % hsm))
        hsm_id = hsm.id
    else:
        comm.send(Message("No MS2GlycopeptideHypothesisSampleMatch used"))
        hsm_id = None

    comm.send(Message("Begin Matching for target %d" % target_hypothesis_id, "update"))
    job = IonMatching(
        database_path,
        hypothesis_id=target_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_sample_match_id=hsm_id,
        sample_run_id=sample_run_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()
    comm.send(Message("End Matching for target %d. %s, %r" % (target_hypothesis_id, job.status, job.error)))
    if job.error is not None:
        comm.send(Message(job.error, 'error'))
    if decoy_hypothesis_id is None:
        return
    comm.send(Message("Begin Matching for decoy %d" % decoy_hypothesis_id, "update"))
    job = IonMatching(
        database_path,
        hypothesis_id=decoy_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_sample_match_id=hsm_id,
        sample_run_id=sample_run_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, target_hypothesis_id, hsm_id)
    job.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, decoy_hypothesis_id, hsm_id)
    job.start()

    comm.send(Message("Begin TDA", "update"))
    job = target_decoy.TargetDecoyAnalyzer(database_path, target_hypothesis_id, decoy_hypothesis_id, hsm_id)
    job.start()
    if hsm_id is not None:
        comm.send(Message(hsm.to_json(), "new-hypothesis-sample-match"))
    return


class TandemMSGlycoproteomicsSearchTask(Task):
    def __init__(
            self, database_path, observed_ions_path, target_hypothesis_id,
            decoy_hypothesis_id, observed_ions_type, sample_run_id,
            ms1_tolerance, ms2_tolerance, callback=lambda: 0, **kwargs):
        args = (
            database_path, observed_ions_path, target_hypothesis_id,
            decoy_hypothesis_id, observed_ions_type, sample_run_id,
            ms1_tolerance, ms2_tolerance)
        name = "Tandem MS Glycoproteomics Search {} @ {}".format(target_hypothesis_id, os.path.basename(observed_ions_path))
        kwargs.setdefault('name', name)
        super(TandemMSGlycoproteomicsSearchTask, self).__init__(taskmain, args, callback, **kwargs)
