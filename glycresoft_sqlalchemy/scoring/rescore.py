from ..data_model import DatabaseManager, TheoreticalGlycopeptide, GlycopeptideMatch
from .score_matches import apply


def rescore(matched, theoretical):
    apply(matched, theoretical)


class MatchScorer(object):

    def __init__(self, database_path):
        self.manager = DatabaseManager(database_path)

    def run(self):
        session = self.manager.session()
        for i in session.query(GlycopeptideMatch.id):
            i = i[0]
            q = session.query(GlycopeptideMatch, TheoreticalGlycopeptide).filter(
                GlycopeptideMatch.id == i, TheoreticalGlycopeptide.id == i)
            for match, theoretical in q:
                rescore(match, theoretical)
                session.add(match)

        session.commit()
