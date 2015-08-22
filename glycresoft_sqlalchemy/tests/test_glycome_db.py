import os
import logging
try:
    logging.basicConfig(level=logging.INFO, filemode='w', filename="testlog",
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy import data_model
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import glycomedb_utils

def test_main():
    try:
        os.remove("datafiles/Glycome-DB.db")
    except:
        pass
    job = glycomedb_utils.GlycomeDBDownloader("datafiles/Glycome-DB.db")
    job.start()

if __name__ == '__main__':
    test_main()
