from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch

from glycresoft_sqlalchemy.matching.matching import IonMatching, ms1_tolerance_default, ms2_tolerance_default
from glycresoft_sqlalchemy.scoring import target_decoy


def taskmain(database_path, observed_ions_path, observed_ions_type='bupid_yaml',
             hypothesis_ids=None,
             ms1_tolerance=ms1_tolerance_default,
             ms2_tolerance=ms2_tolerance_default,
             hypothesis_match_id=None,
             sample_run_id=None,
             n_processes=4):
    for hypothesis_id in hypothesis_ids:
        job = IonMatching(
            database_path=database_path,
            observed_ions_path=observed_ions_path,
            observed_ions_type=observed_ions_type,
            hypothesis_id=hypothesis_id,
            ms1_tolerance=ms1_tolerance,
            ms2_tolerance=ms2_tolerance,
            )
        job.start()

    manager = DatabaseManager(database_path)
    session = manager.session()
    if hypothesis_ids is None:
        for hypothesis in session.query(Hypothesis).filter(~Hypothesis.is_decoy):
            for decoy_hypothesis in hypothesis.parameters.get("decoy", []):
                decoy_id = decoy_hypothesis["hypothesis_id"]
                job = target_decoy.TargetDecoyAnalyzer(database_path, hypothesis.id, decoy_id)
                job.start()
    return
