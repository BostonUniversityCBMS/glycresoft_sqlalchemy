from ..data_model import DatabaseManager, TheoreticalGlycopeptide, GlycopeptideMatch
from .score_matches import evaluate


def rescore(matched, theoretical):
    evaluate(matched, theoretical)


class MatchScorer(object):

    def __init__(self, database_path):
        self.manager = DatabaseManager(database_path)

    def run(self):
        session = self.manager.session()
        for i in session.query(GlycopeptideMatch.id):
            i = i[0]
            q = session.query(GlycopeptideMatch, TheoreticalGlycopeptide).filter(
                GlycopeptideMatch.id == i, GlycopeptideMatch.theoretical_glycopeptide_id == TheoreticalGlycopeptide.id)
            for match, theoretical in q:
                rescore(match, theoretical)
                session.add(match)

        session.commit()


def run(session):
        for i in session.query(GlycopeptideMatch.id):
            i = i[0]
            q = session.query(GlycopeptideMatch, TheoreticalGlycopeptide).filter(
                GlycopeptideMatch.id == i, GlycopeptideMatch.theoretical_glycopeptide_id == TheoreticalGlycopeptide.id)
            for match, theoretical in q:
                rescore(match, theoretical)
                session.add(match)
        session.commit()
