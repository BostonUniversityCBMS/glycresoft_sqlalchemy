from glycresoft_sqlalchemy.data_model import DatabaseManager, GlycopeptideMatch, Experiment, TheoreticalGlycopeptide
from sqlalchemy.sql.functions import max as sql_max


def main(database_path):
    session = DatabaseManager(database_path).session()
    for experiment in session.query(Experiment):
        print experiment
        for protein in experiment.proteins.values():
            print protein
            print "Top Score: ", session.query(
                sql_max(GlycopeptideMatch.ms2_score)).filter(GlycopeptideMatch.protein_id == protein.id).first()
            print "Number of matches: ", session.query(GlycopeptideMatch).filter(
                GlycopeptideMatch.protein_id == protein.id).count()


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
