import argparse
import logging
try:
    logger = logging.getLogger("run_search")
except:
    pass
from glycresoft_sqlalchemy.matching import matching
from glycresoft_sqlalchemy.scoring import target_decoy
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from glycresoft_sqlalchemy.data_model import DatabaseManager, HypothesisSampleMatch, SampleRun, Hypothesis

import summarize

ms1_tolerance_default = matching.ms1_tolerance_default
ms2_tolerance_default = matching.ms2_tolerance_default


def run_search(
        database_path, observed_ions_path, target_hypothesis_id=None,
        decoy_hypothesis_id=None, observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default, **kwargs):
    if target_hypothesis_id is None:
        target_hypothesis_id = 1
    manager = DatabaseManager(database_path)
    manager.initialize()
    if observed_ions_type == 'bupid_yaml' and observed_ions_path[-3:] != '.db':
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
        hsm_id = hsm.id
    else:
        hsm_id = None

    job = matching.IonMatching(
        database_path,
        hypothesis_id=target_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_match_id=hsm_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()

    if decoy_hypothesis_id is None:
        return

    job = matching.IonMatching(
        database_path,
        hypothesis_id=decoy_hypothesis_id,
        observed_ions_path=observed_ions_path,
        observed_ions_type=observed_ions_type,
        hypothesis_match_id=hsm_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()

    job = target_decoy.TargetDecoyAnalyzer(database_path, target_hypothesis_id, decoy_hypothesis_id)
    job.start()

    summarize.main(database_path)


app = argparse.ArgumentParser('database-search')
app.add_argument("database_path")
app.add_argument("target_hypothesis_id")
app.add_argument("-n", "--n-processes", default=4, required=False, type=int)
app.add_argument("-i", "--observed-ions-path")
app.add_argument("-p", "--observed-ions-type", default='bupid_yaml', choices=["bupid_yaml", "db"])
app.add_argument("-d", "--decoy-hypothesis-id", type=int, default=None, required=False)
app.add_argument("-t1", "--ms1-tolerance", default=ms1_tolerance_default, required=False, type=float)
app.add_argument("-t2", "--ms2-tolerance", default=ms2_tolerance_default, required=False, type=float)


def main():
    args = app.parse_args()
    logger.debug("Arguments %r", args)
    run_search(**args.__dict__)

if __name__ == '__main__':
    main()
