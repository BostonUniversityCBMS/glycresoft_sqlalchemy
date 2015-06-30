import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s:%(processName)s- %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy.search_space_builder import naive_glycopeptide_hypothesis


def test_main():
    try:
        os.remove("./datafiles/naive_glycopeptide.db")
    except:
        pass
    constant_mods, variable_mods = (["Carbamidomethyl (C)"], ["Deamidated (Q)", "Deamidated (N)"])
    enzyme = 'trypsin'
    job = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesisBuilder(
        "./datafiles/naive_glycopeptide.db", "test", "./datafiles/proteins.fasta",
        None, "./datafiles/human_n_glycan.csv", 'csv', constant_mods,
        variable_mods, enzyme, n_processes=6)
    job.start()


if __name__ == '__main__':
    test_main()
