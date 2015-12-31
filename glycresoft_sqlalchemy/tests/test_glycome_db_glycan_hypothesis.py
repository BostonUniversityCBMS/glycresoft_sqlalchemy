import os
import logging
try:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.search_space_builder.glycan_builder import glycomedb_utils


db_file = "datafiles/human_glycan_structure_hypothesis.db"


def test_main():
    try:
        os.remove(db_file)
    except:
        pass
    job = glycomedb_utils.GlycomeDBHypothesis(
        db_file, taxa_ids=[9606], glycomedb_path="datafiles/Glycome-DB.db", motif_family='N-Linked Glycans')
    job.start()

if __name__ == '__main__':
    test_main()
