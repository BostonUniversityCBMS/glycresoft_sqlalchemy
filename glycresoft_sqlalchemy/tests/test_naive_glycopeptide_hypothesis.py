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
    constant_mods, variable_mods = (["Carbamidomethyl (C)"], ["Deamidated"])
    enzyme = 'trypsin'
    job = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesisBuilder(
        "./datafiles/naive_glycopeptide.db", "test", "./datafiles/proteins_agp_only.fasta",
        None, "./datafiles/human_n_glycans.txt", 'txt', constant_mods,
        variable_mods, enzyme, max_missed_cleavages=1, maximum_glycosylation_sites=1, n_processes=5)
    job.start()


if __name__ == '__main__':
    test_main()
