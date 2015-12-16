from glycresoft_sqlalchemy.data_model import DatabaseManager, TheoreticalGlycopeptide, Protein
from glycresoft_sqlalchemy.search_space_builder import make_decoys

import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")


def test_main():
    db_file_name = "datafiles/integrated_omics_simple.db"
    dbm = DatabaseManager(db_file_name)
    s = dbm()
    print s.query(TheoreticalGlycopeptide).join(Protein).filter(Protein.hypothesis_id == 3).count()
    job = make_decoys.BatchingDecoySearchSpaceBuilder(db_file_name, hypothesis_ids=[3], n_processes=4)
    # job.session.query(make_decoys.Hypothesis.id).filter(make_decoys.Hypothesis.name.like("decoy")).delete('fetch')
    # job.session.commit()
    job.start()

if __name__ == '__main__':
    test_main()
