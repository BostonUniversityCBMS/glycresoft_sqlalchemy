import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy.search_space_builder import integrated_omics
from glycresoft_sqlalchemy.data_model import TheoreticalGlycopeptideComposition, func
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics

def test_main():
    db_file = "./datafiles/build_informed_hypothesis_test.db"
    os.remove(db_file)

    rules_table = {
        "Hex": (3, 8),
        "HexNAc": (2, 8),
        "Fuc": (0, 5),
        "NeuAc": (0, 4)
    }

    job = constrained_combinatorics.ConstrainedCombinatoricsGlycanHypothesisBuilder(
        db_file, rules_table=rules_table, constraints_list=[])
    job.start()

    n_processes = 6

    simple_proteome = [
        "P02763|A1AG1_HUMAN", "P19652|A1AG2_HUMAN"
    ]

    job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
        db_file,
        protein_ids=simple_proteome,
        mzid_path="datafiles/AGP_Proteomics2.mzid",
        glycomics_path='./datafiles/human_n_glycans.txt',
        glycomics_format='txt',
        # glycomics_path=db_file,
        # glycomics_format='hypothesis',
        # source_hypothesis_id=1,
        maximum_glycosylation_sites=1,
        include_all_baseline=True,
        n_processes=n_processes)

    job.start()

    db = job.manager.session()
    dups = db.query(
        TheoreticalGlycopeptideComposition.glycopeptide_sequence,
        func.count(TheoreticalGlycopeptideComposition.glycopeptide_sequence)).group_by(
        TheoreticalGlycopeptideComposition.glycopeptide_sequence,
        TheoreticalGlycopeptideComposition.protein_id).having(
        func.count(TheoreticalGlycopeptideComposition.glycopeptide_sequence) > 1).all()
    assert len(dups) == 0, "Duplicate Sequences Found"
if __name__ == '__main__':
    test_main()
