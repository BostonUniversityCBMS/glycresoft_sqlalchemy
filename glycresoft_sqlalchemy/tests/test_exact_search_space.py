import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder


def test_main():
    enzyme = ["trypsin"]
    assert os.path.exists("./datafiles/integrated_omics_simple.db")
    job = exact_search_space_builder.ExactSearchSpaceBuilder("./datafiles/ResultOf20140918_01_isos_Informed.csv",
                                                             "./datafiles/integrated_omics_simple.db",
                                                             site_list="./datafiles/sitelist.txt", enzyme=enzyme,
                                                             n_processes=1)
    job.start()

if __name__ == '__main__':
    test_main()
