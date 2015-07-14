import os
import logging
# try:
#     logging.basicConfig(level=logging.DEBUG, filemode='w',
#                         format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
#                         datefmt="%H:%M:%S")
# except:
#     pass

from glycresoft_sqlalchemy.search_space_builder import integrated_omics
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder
from glycresoft_sqlalchemy.search_space_builder import make_decoys
from glycresoft_sqlalchemy.app.run_search import run_ms2_glycoproteomics_search

db_file_name = "./datafiles/integrated_omics_simple.db"


def test_main():
    i = 1
    try:
        os.remove(db_file_name)
    except:
        print "Did not clear %s" % db_file_name
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
    run_ms2_glycoproteomics_search(
        db_file_name, "./datafiles/20140918_01.db", target_hypothesis_id=i,
        decoy_hypothesis_id=i + 1, observed_ions_type='bupid_yaml', n_processes=6)

if __name__ == '__main__':
    test_main()
