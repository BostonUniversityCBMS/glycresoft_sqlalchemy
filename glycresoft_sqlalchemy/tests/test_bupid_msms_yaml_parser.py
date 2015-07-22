import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.spectra import bupid_topdown_deconvoluter_sa

import os


def test_main():
    try:
        os.remove("datafiles/20140918_01.db")
    except:
        pass
    f = "datafiles/20140918_01.yaml"
    bupid_topdown_deconvoluter_sa.process_data_file(f)

if __name__ == '__main__':
    test_main()
