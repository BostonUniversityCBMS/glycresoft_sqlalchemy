from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch,
    HypothesisSampleMatch, SampleRun)

from glycresoft_sqlalchemy.matching.matching import IonMatching, ms1_tolerance_default, ms2_tolerance_default
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from glycresoft_sqlalchemy.scoring import target_decoy

from .task_process import NullPipe, Message


def taskmain(
        database_path, observed_ions_path, target_hypothesis_id=None,
        decoy_hypothesis_id=None, observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default,
        comm=None,
        **kwargs):
    if comm is None:
        comm = NullPipe()
    if target_hypothesis_id is None:
        target_hypothesis_id = 1

    comm.send(Message("Starting do_search:taskmain"))

    manager = DatabaseManager(database_path)
    manager.initialize()

    if observed_ions_type == 'bupid_yaml' and observed_ions_path[-3:] != '.db':
        comm.send(Message("Converting %s to db" % observed_ions_path))
        parser = BUPIDMSMSYamlParser(observed_ions_path)
        observed_ions_path = parser.manager.path
        observed_ions_type = 'db'
        sample_name = parser.sample_run_name
    else:
        sample_name = ','.join(x[0] for x in DatabaseManager(observed_ions_path).session().query(SampleRun.name).all())
    session = manager.session()
    if decoy_hypothesis_id is not None:
        hsm = HypothesisSampleMatch(
            target_hypothesis_id=target_hypothesis_id,
            decoy_hypothesis_id=decoy_hypothesis_id,
            sample_run_name=sample_name
        )
        session.add(hsm)
        session.commit()
        comm.send(Message("Created %r" % hsm))
        hsm_id = hsm.id
    else:
        comm.send(Message("No HypothesisSampleMatch used"))
        hsm_id = None

    comm.send(Message("Begin Matching for target %d" % target_hypothesis_id))
    job = IonMatching(
        database_path,
        hypothesis_id=target_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_sample_match_id=hsm_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()
    comm.send(Message("End Matching for target %d. %s, %r" % (target_hypothesis_id, job.status, job.error)))
    if job.error is not None:
        comm.send(Message(job.error, 'error'))
    if decoy_hypothesis_id is None:
        return
    comm.send(Message("Begin Matching for decoy %d" % decoy_hypothesis_id))
    job = IonMatching(
        database_path,
        hypothesis_id=decoy_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_sample_match_id=hsm_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()
    comm.send(Message("Begin Matching for decoy %d" % decoy_hypothesis_id))

    job = target_decoy.TargetDecoyAnalyzer(database_path, target_hypothesis_id, decoy_hypothesis_id)
    job.start()
    return
