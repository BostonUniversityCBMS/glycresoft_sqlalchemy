import os

from glycresoft_sqlalchemy.data_model import Hypothesis, DatabaseManager, get_or_create
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms1 import include_glycomics
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder import glycan_utilities


def testmain():
    db_file_name = "./datafiles/glycan_combination_test.db"
    glycan_file_path = "./datafiles/human_n_glycans.txt"
    try:
        os.remove(db_file_name)
    except:
        pass

    manager = DatabaseManager(db_file_name)
    manager.initialize()
    session = manager()
    hypothesis, _ = get_or_create(session, Hypothesis, name="Test Hypothesis")
    session.commit()

    include_glycomics.MS1GlycanImporter(db_file_name, glycan_file_path, hypothesis_id=hypothesis.id).start()
    glycan_utilities.create_combinations(session, 2, hypothesis.id)


if __name__ == '__main__':
    testmain()
