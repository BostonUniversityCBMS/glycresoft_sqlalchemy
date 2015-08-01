import os
import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")

from glycresoft_sqlalchemy import data_model
from glycresoft_sqlalchemy.search_space_builder import naive_glycopeptide_hypothesis
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder, pooling_make_decoys
from glycresoft_sqlalchemy.matching import matching, peak_grouping
from glycresoft_sqlalchemy.scoring import target_decoy
from glycresoft_sqlalchemy.search_space_builder import integrated_omics


def test_main():
    db_file_name = "./datafiles/integrated_omics_simple.db"
    try:
        os.remove(db_file_name)
    except:
        pass

    i = integrated_omics.load_proteomics("datafiles/integrated_omics_simple.db", "datafiles/AGP_Proteomics2.mzid")
    integrated_omics.load_glycomics_naive("datafiles/integrated_omics_simple.db", "datafiles/human_n_glycan.csv", i)
    job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
        "./datafiles/integrated_omics_simple.db", i, n_processes=4)
    job.start()

    ec = os.system(
        r"glycresoft-database-search ms1 -n 6 -i datafiles\20140918_01_isos.db {db_file_name}\
         -p db -g 2e-5 --skip-grouping".format(db_file_name=db_file_name))
    assert ec == 0
    job = exact_search_space_builder.ExactSearchSpaceBuilder.from_hypothesis(
        "datafiles/naive_glycopeptide.db", 1, 6)
    hypothesis_id = job.start()
    job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(db_file_name, hypothesis_ids=[hypothesis_id])
    decoy_hypothesis_id = job.start()[0]
    manager = data_model.DatabaseManager(db_file_name)
    session = manager.session()
    hsm = data_model.HypothesisSampleMatch(
        target_hypothesis_id=hypothesis_id,
        decoy_hypothesis_id=decoy_hypothesis_id,
        sample_run_name="20140918_01.yaml",
        name="End-to-End AGP @ 20140918_01")
    session.add(hsm)
    session.commit()

    hsm_id = hsm.id

    matcher = matching.IonMatching(db_file_name, hypothesis_id, r"datafiles\20140918_01.db",
                               "db", ms1_tolerance=1e-5, ms2_tolerance=2e-5,
                               hypothesis_sample_match_id=hsm_id, sample_run_id=1, n_processes=8)
    matcher.start()
    matcher = matching.IonMatching(db_file_name, decoy_hypothesis_id, r"datafiles\20140918_01.db",
                               "db", ms1_tolerance=1e-5, ms2_tolerance=2e-5,
                               hypothesis_sample_match_id=hsm_id, sample_run_id=1, n_processes=8)
    matcher.start()
    tda = target_decoy.TargetDecoyAnalyzer(db_file_name, hypothesis_id, decoy_hypothesis_id)
    tda.start()


if __name__ == '__main__':
    test_main()
