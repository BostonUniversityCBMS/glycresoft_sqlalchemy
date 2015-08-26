import os
import logging
try:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.search_space_builder.glycan_builder import structure_fragmentation


db_file = "datafiles/human_glycan_structure_hypothesis.db"


def test_main():
    job = structure_fragmentation.GlycanStructureFragmenter(db_file, 1, kind="ABCXYZ", n_processes=4)
    job.start()

if __name__ == '__main__':
    test_main()
