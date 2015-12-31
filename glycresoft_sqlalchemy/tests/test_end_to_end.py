import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.search_space_builder import naive_glycopeptide_hypothesis
from glycresoft_sqlalchemy.search_space_builder import search_space_builder, pooling_make_decoys
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics
from glycresoft_sqlalchemy.matching.glycopeptide.pipeline import GlycopeptideFragmentMatchingPipeline


def test_main():
    db_file_name = "./datafiles/naive_glycopeptide.db"
    try:
        os.remove(db_file_name)
    except:
        pass

    rules_table = {
        "Hex": (3, 8),
        "HexNAc": (2, 8),
        "Fuc": (0, 5),
        "NeuAc": (0, 4)
    }

    job = constrained_combinatorics.ConstrainedCombinatoricsGlycanHypothesisBuilder(
        db_file_name, rules_table=rules_table, constraints_list=[])
    combn_glycan_hypothesis_id = job.start()
    print "combn_glycan_hypothesis_id", combn_glycan_hypothesis_id

    constant_mods, variable_mods = (["Carbamidomethyl (C)"], ["Deamidated (N)"])
    enzyme = 'trypsin'
    job = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesisBuilder(
        database_path=db_file_name,
        hypothesis_name="test",
        protein_file="./datafiles/proteins_agp_only.fasta",
        site_list_file=None,
        constant_modifications=constant_mods,
        variable_modifications=variable_mods,
        enzyme=enzyme,
        glycomics_path=db_file_name,
        glycomics_format='hypothesis',
        source_hypothesis_id=combn_glycan_hypothesis_id,
        # glycomics_path='./datafiles/human_n_glycans.txt',
        # glycomics_format='txt',
        maximum_glycosylation_sites=1,
        max_missed_cleavages=1,
        n_processes=6)
    glycopeptide_hypothesis_id = job.start()

    return 1

    ec = os.system(
        ("glycresoft-database-search ms1 -n 6 -i datafiles/20140918_01_isos.db datafiles/naive_glycopeptide.db %d "
         "-p db -g 2e-5 --skip-grouping") % glycopeptide_hypothesis_id)
    assert ec == 0
    job = search_space_builder.BatchingTheoreticalSearchSpaceBuilder.from_hypothesis_sample_match(
        "datafiles/naive_glycopeptide.db", 1, 6)
    hypothesis_id = job.start()
    job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(db_file_name, hypothesis_ids=[hypothesis_id])
    decoy_hypothesis_id = job.start()[0]

    job = GlycopeptideFragmentMatchingPipeline(
        db_file_name, r"datafiles\20140918_01.db",
        target_hypothesis_id=hypothesis_id,
        decoy_hypothesis_id=decoy_hypothesis_id,
        sample_run_name="20140918_01.yaml",
        hypothesis_sample_match_name="End-to-End AGP @ 20140918_01")
    job.start()

if __name__ == '__main__':
    test_main()
