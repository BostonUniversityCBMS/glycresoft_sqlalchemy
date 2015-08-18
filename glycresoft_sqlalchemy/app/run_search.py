import argparse
import logging
try:
    logger = logging.getLogger("run_search")
except:
    pass
from glycresoft_sqlalchemy.matching import matching, peak_grouping
from glycresoft_sqlalchemy.scoring import target_decoy, score_spectrum_matches
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from glycresoft_sqlalchemy.spectra.decon2ls_sa import Decon2LSIsosParser
from glycresoft_sqlalchemy.data_model import DatabaseManager, MS2GlycopeptideHypothesisSampleMatch, SampleRun, Hypothesis

import summarize

ms1_tolerance_default = matching.ms1_tolerance_default
ms2_tolerance_default = matching.ms2_tolerance_default


def run_ms2_glycoproteomics_search(
        database_path, observed_ions_path, target_hypothesis_id=None,
        decoy_hypothesis_id=None, observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default, **kwargs):
    if target_hypothesis_id is None:
        target_hypothesis_id = 1
    manager = DatabaseManager(database_path)
    manager.initialize()
    if observed_ions_type == 'bupid_yaml' and observed_ions_path[-3:] != '.db':
        parser = BUPIDMSMSYamlParser(observed_ions_path, manager.bridge_address())
        observed_ions_path = parser.manager.path
        observed_ions_type = 'db'
        sample_name = parser.sample_run_name
    else:
        sample_name = ','.join(x[0] for x in DatabaseManager(observed_ions_path).session().query(SampleRun.name).all())
    session = manager.session()
    if decoy_hypothesis_id is not None:
        hsm = MS2GlycopeptideHypothesisSampleMatch(
            target_hypothesis_id=target_hypothesis_id,
            decoy_hypothesis_id=decoy_hypothesis_id,
            sample_run_name=sample_name,
            name="{hypothesis.name}_on_{sample_name}".format(hypothesis=session.query(
                Hypothesis).get(target_hypothesis_id), sample_name=sample_name)
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
        hypothesis_sample_match_id=hsm_id,
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
        hypothesis_sample_match_id=hsm_id,
        ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance,
        n_processes=kwargs.get("n_processes", 4))
    job.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, target_hypothesis_id, hsm_id)
    job.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(database_path, decoy_hypothesis_id, hsm_id)
    job.start()

    job = target_decoy.TargetDecoyAnalyzer(database_path, target_hypothesis_id, decoy_hypothesis_id, hsm_id)
    job.start()

    summarize.main(database_path)


def run_ms1_search(
        database_path, observed_ions_path, hypothesis_id=None,
        observed_ions_type='isos', sample_run_id=None,
        grouping_tolerance=8e-5, search_type=None,
        match_tolerance=1e-5, mass_shift=None, n_processes=4, **kwargs):
    search_type = {
        "glycopeptide": "TheoreticalGlycopeptideComposition",
        "glycan": "TheoreticalGlycanComposition"
    }[search_type]

    manager = DatabaseManager(database_path)

    if observed_ions_type == 'isos':
        parser = Decon2LSIsosParser(observed_ions_path, manager.bridge_address())
        observed_ions_path = parser.manager.path

    # mass_shift converted to dictionary for mass_shift_map

    pipeline = peak_grouping.LCMSPeakClusterSearch(
        database_path, observed_ions_path, hypothesis_id, sample_run_id=1,
        grouping_error_tolerance=grouping_tolerance,
        search_type=search_type, match_tolerance=match_tolerance, mass_shift_map=None,
        n_processes=n_processes, **kwargs)
    pipeline.start()


app = argparse.ArgumentParser('database-search')

subparsers = app.add_subparsers()

ms2_glycoproteomics_app = subparsers.add_parser("ms2-glycoproteomics")

ms2_glycoproteomics_app.add_argument("database_path")
ms2_glycoproteomics_app.add_argument("target_hypothesis_id")
ms2_glycoproteomics_app.add_argument("-n", "--n-processes", default=4, required=False, type=int)
ms2_glycoproteomics_app.add_argument("-i", "--observed-ions-path")
ms2_glycoproteomics_app.add_argument("-p", "--observed-ions-type", default='bupid_yaml', choices=["bupid_yaml", "db"])
ms2_glycoproteomics_app.add_argument("-d", "--decoy-hypothesis-id", type=int, default=None, required=False)
ms2_glycoproteomics_app.add_argument(
    "-t1", "--ms1-tolerance", default=ms1_tolerance_default, required=False, type=float)
ms2_glycoproteomics_app.add_argument(
    "-t2", "--ms2-tolerance", default=ms2_tolerance_default, required=False, type=float)
ms2_glycoproteomics_app.set_defaults(task=run_ms2_glycoproteomics_search)

ms1_app = subparsers.add_parser("ms1")
ms1_app.add_argument("database_path")
ms1_app.add_argument("hypothesis_id")
ms1_app.add_argument("-n", "--n-processes", default=4, required=False, type=int)
ms1_app.add_argument("-i", "--observed-ions-path")
ms1_app.add_argument("-p", "--observed-ions-type", default='isos', choices=["isos", "db"])
ms1_app.add_argument("-t", "--match-tolerance", default=1e-5, required=False, type=float)
ms1_app.add_argument("-g", "--grouping-tolerance", default=8e-5, required=False, type=float)
ms1_app.add_argument("-s", "--search-type", default='glycopeptide', choices=['glycan', 'glycopeptide'])
# ms1_app.add_argument("-m", "--mass-shift", action='append', nargs=4)
ms1_app.add_argument('--skip-grouping', action='store_true', required=False)
ms1_app.add_argument('--skip-matching', action='store_true', required=False)
ms1_app.add_argument('--hypothesis-sample-match-id', action='store', default=None, required=False)
ms1_app.set_defaults(task=run_ms1_search)


def main():
    args = app.parse_args()
    logger.debug("Arguments %r", args)
    task = args.task
    del args.task
    task(**args.__dict__)

if __name__ == '__main__':
    main()
