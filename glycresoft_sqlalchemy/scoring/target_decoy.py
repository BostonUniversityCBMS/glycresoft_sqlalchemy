from collections import namedtuple

from sqlalchemy import func, distinct
from ..data_model import DatabaseManager, Experiment, GlycopeptideMatch, Protein


Threshold = namedtuple("Threshold", ("score", "targets", "decoys", "fdr"))


class TargetDecoyAnalyzer(object):
    manager_type = DatabaseManager

    def __init__(self, database_path, target_experiment_id=None, decoy_experiment_id=None):
        self.manager = self.manager_type(database_path)
        self.target_id = target_experiment_id
        self.decoy_id = decoy_experiment_id

    def thresholds(self):
        session = self.manager.session()
        thresholds = session.query(distinct(func.round(GlycopeptideMatch.ms2_score, 2)))
        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id, Protein.id == self.target_id)
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id, Protein.id == self.decoy_id)

        results = {}

        for score in thresholds:
            decoys_at, = list(dq.filter(GlycopeptideMatch.ms2_score >= score).count())
            targets_at, = list(tq.filter(GlycopeptideMatch.ms2_score >= score).count())
            print targets_at
            ratio = 0 # decoys_at / float(targets_at)

            result = Threshold(score, targets_at, decoys_at, ratio)
            results[score] = result

        return results

    def run(self):
        thresholds = self.thresholds()
        return thresholds
