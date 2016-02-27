import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics


def testmain():
    db_file_name = "./datafiles/combinatorial_glycome.db"
    try:
        os.remove(db_file_name)
    except:
        pass
    # rules_table = {
    #     "Hex": (3, 10),
    #     "HexNAc": (2, 8),
    #     "Fuc": (0, 5),
    #     "NeuAc": (0, 4)
    # }
    rules_table, constraints = constrained_combinatorics.parse_rules_from_file("./datafiles/combinatorial_rules.txt")
    job = constrained_combinatorics.ConstrainedCombinatoricsGlycanHypothesisBuilder(
        db_file_name, rules_table=rules_table, constraints_list=constraints)
    job.start()


if __name__ == '__main__':
    testmain()
