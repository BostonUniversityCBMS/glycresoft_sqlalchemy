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


def test_main():
    i = 1
    try:
        os.remove("./datafiles/integrated_omics.db")
    except:
        pass
    i = integrated_omics.load_proteomics("datafiles/integrated_omics.db", "datafiles/AGP_Transferin_Proteomics.mzid")
    integrated_omics.load_glycomics_naive("datafiles/integrated_omics.db", "datafiles/human_n_glycan.csv", i)
    job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
        "./datafiles/integrated_omics.db", i, protein_ids=[2, 4, 5, 6, 8, 9, 15, 20], n_processes=4)
    job.start()
    job = integrated_omics.IntegratedOmicsMS1LegacyCSV("datafiles/integrated_omics.db", i)
    job.start()

    db = job.manager.session()
    dups = db.query(
        TheoreticalGlycopeptideComposition.glycopeptide_sequence,
        func.count(TheoreticalGlycopeptideComposition.glycopeptide_sequence)).group_by(
        TheoreticalGlycopeptideComposition.glycopeptide_sequence).having(
        func.count(TheoreticalGlycopeptideComposition.glycopeptide_sequence) > 1).all()
    assert len(dups) == 0, "Duplicate Sequences Found"

if __name__ == '__main__':
    test_main()
