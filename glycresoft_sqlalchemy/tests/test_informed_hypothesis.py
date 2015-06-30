import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.search_space_builder import integrated_omics
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder
from glycresoft_sqlalchemy.search_space_builder import make_decoys
from glycresoft_sqlalchemy.matching import matching
from glycresoft_sqlalchemy.scoring import target_decoy
import summarize

db_file_name = "./datafiles/integrated_omics_simple.db"


def test_main():
    i = 1
    try:
        os.remove(db_file_name)
    except:
        pass
    i = integrated_omics.load_proteomics(db_file_name, "datafiles/AGP_Proteomics2.mzid")
    integrated_omics.load_glycomics_naive(db_file_name, "datafiles/human_n_glycan.csv", i)
    job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
        db_file_name, i, n_processes=4)
    job.start()
    job = integrated_omics.IntegratedOmicsMS1LegacyCSV(db_file_name, i)
    job.start()

    enzyme = []
    assert os.path.exists(db_file_name)
    job = exact_search_space_builder.ExactSearchSpaceBuilder("./datafiles/ResultOf20140918_01_isos_Informed.csv",
                                                             db_file_name, i, enzyme=enzyme)
    job.start()

    job = make_decoys.DecoySearchSpaceBuilder(db_file_name, hypothesis_ids=[i], n_processes=4)
    job.start()
    job = matching.IonMatching(db_file_name, "./datafiles/20140918_01.yaml.db",
                               hypothesis_ids=None, n_processes=6)
    job.start()
    tda = target_decoy.TargetDecoyAnalyzer(db_file_name, i, i + 1)
    tda.start()
    summarize.main(db_file_name)

if __name__ == '__main__':
    test_main()
