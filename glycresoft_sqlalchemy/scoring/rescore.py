from ..data_model import DatabaseManager, TheoreticalGlycopeptide, GlycopeptideMatch, PipelineModule
from .score_matches import evaluate


def rescore(matched, theoretical):
    evaluate(matched, theoretical)


class MatchScorer(PipelineModule):
    def __init__(self, database_path):
        self.manager = self.manager_type(database_path)

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
