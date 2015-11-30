import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy import data_model
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder, pooling_make_decoys
from glycresoft_sqlalchemy.matching import matching
from glycresoft_sqlalchemy.scoring import target_decoy, score_spectrum_matches
from glycresoft_sqlalchemy.search_space_builder import integrated_omics

# import multiprocessing
# multiprocessing.log_to_stderr()


def test_main():
    db_file_name = "./datafiles/integrated_omics_simple.db"
    # try:
    #     os.remove(db_file_name)
    #     pass
    # except:
    #     pass

    # # db_file_name = "postgresql:///db"

    # job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
    #     db_file_name,
    #     protein_ids=["P02763|A1AG1_HUMAN", "P19652|A1AG2_HUMAN"],  # "P19652|A1AG2_HUMAN",
    #     mzid_path="datafiles/AGP_Proteomics2.mzid",
    #     glycomics_path="datafiles/human_n_glycans.txt",
    #     n_processes=4)
    # job.start()

    # ec = os.system(
    #     "glycresoft-database-search ms1 -n 5 -i datafiles/20140918_01_isos.db -p db -g 2e-5 --skip-grouping {db_file_name} 1".format(
    #         db_file_name=db_file_name))
    # assert ec == 0
    # job = exact_search_space_builder.ExactSearchSpaceBuilder.from_hypothesis_sample_match(
    #     db_file_name, 1, 4)
    # hypothesis_id = job.start()
    hypothesis_id = 3

    job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(db_file_name, hypothesis_ids=[hypothesis_id])
    decoy_hypothesis_id = job.start()[0]
    # decoy_hypothesis_id = 3
    manager = data_model.DatabaseManager(db_file_name)
    session = manager.session()
    hsm = data_model.MS2GlycopeptideHypothesisSampleMatch(
        target_hypothesis_id=hypothesis_id,
        decoy_hypothesis_id=decoy_hypothesis_id,
        sample_run_name="20140918_01.yaml",
        name="End-to-End AGP @ 20140918_01")
    session.add(hsm)
    session.commit()

    hsm_id = hsm.id

    matcher = matching.IonMatching(
        db_file_name, hypothesis_id, r"datafiles\20140918_01.db",
        "db", ms1_tolerance=1e-5, ms2_tolerance=2e-5,
        hypothesis_sample_match_id=hsm_id, sample_run_id=1, n_processes=6)
    matcher.start()

    matcher = matching.IonMatching(
        db_file_name, decoy_hypothesis_id, r"datafiles\20140918_01.db",
        "db", ms1_tolerance=1e-5, ms2_tolerance=2e-5,
        hypothesis_sample_match_id=hsm_id, sample_run_id=1, n_processes=6)

    matcher.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(db_file_name, hypothesis_id, hsm_id, n_processes=6)
    job.start()

    job = score_spectrum_matches.SimpleSpectrumAssignment(db_file_name, decoy_hypothesis_id, hsm_id, n_processes=6)
    job.start()

    tda = target_decoy.TargetDecoyAnalyzer(db_file_name, hypothesis_id, decoy_hypothesis_id)
    tda.start()


if __name__ == '__main__':
    test_main()
