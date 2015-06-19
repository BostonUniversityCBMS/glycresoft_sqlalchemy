from itertools import chain
from collections import namedtuple

from sqlalchemy import func, distinct
from ..data_model import DatabaseManager, GlycopeptideMatch, Protein


Threshold = namedtuple("Threshold", ("score", "targets", "decoys", "fdr"))

ms2_score = GlycopeptideMatch.ms2_score
p_value = GlycopeptideMatch.p_value


class TargetDecoyAnalyzer(object):
    manager_type = DatabaseManager

    def __init__(self, database_path, target_experiment_id=None, decoy_experiment_id=None):
        self.manager = self.manager_type(database_path)
        self.target_id = target_experiment_id
        self.decoy_id = decoy_experiment_id

    def target_decoy_ratio(self, cutoff, score=ms2_score):
        session = self.manager.session()
        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id)
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.decoy_id)

        decoys_at = dq.filter(GlycopeptideMatch.ms2_score >= cutoff).count()
        targets_at = tq.filter(GlycopeptideMatch.ms2_score >= cutoff).count()
        try:
            ratio = decoys_at / float(targets_at)
        except ZeroDivisionError:
            ratio = 1.

        return ratio, targets_at, decoys_at

    def global_thresholds(self):
        session = self.manager.session()

        thresholds = session.query(distinct(func.round(GlycopeptideMatch.ms2_score, 2)))

        results = {}

        for score in thresholds:
            score = score[0]
            ratio, targets_at, decoys_at = self.target_decoy_ratio(score)
            if ratio >= 0.5:
                continue
            result = Threshold(score, targets_at, decoys_at, ratio)
            results[score] = result

        session.close()
        return results

    def glycoform_thresholds(self):
        session = self.manager.session()

        thresholds = session.query(distinct(func.round(GlycopeptideMatch.ms2_score, 2)))
        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id).group_by(
            GlycopeptideMatch.base_peptide_sequence, GlycopeptideMatch.glycan_composition_str).order_by(
            GlycopeptideMatch.ms2_score.desc())
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id, Protein.experiment_id == self.decoy_id).group_by(
            GlycopeptideMatch.base_peptide_sequence, GlycopeptideMatch.glycan_composition_str).order_by(
            GlycopeptideMatch.ms2_score.desc())

        results = {}

        for score in thresholds:
            score = score[0]
            decoys_at = dq.filter(GlycopeptideMatch.ms2_score >= score).count()
            targets_at = tq.filter(GlycopeptideMatch.ms2_score >= score).count()
            try:
                ratio = decoys_at / float(targets_at)
            except ZeroDivisionError:
                ratio = 0
            if ratio >= 0.5:
                continue
            result = Threshold(score, targets_at, decoys_at, ratio)
            results[score] = result

        session.close()
        return results

    def p_values(self):
        session = self.manager.session()

        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id).order_by(
            GlycopeptideMatch.ms2_score.desc())
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.decoy_id)

        total_decoys = float(dq.count())

        last_score = 0
        last_p_value = 0
        for target in tq:
            if target.ms2_score == last_score:
                target.p_value = last_p_value
                session.add(target)
            else:
                session.commit()
                decoys_at = dq.filter(GlycopeptideMatch.ms2_score >= target.ms2_score).count()
                last_score = target.ms2_score
                last_p_value = decoys_at / total_decoys
                target.p_value = last_p_value
                session.add(target)
        session.commit()
        session.close()

    def estimate_percent_incorrect_targets(self, cutoff, score=ms2_score):
        session = self.manager.session()

        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id)

        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.decoy_id)

        target_cut = tq.filter(score >= 0, score < cutoff).count()
        decoy_cut = dq.filter(score >= 0, score < cutoff).count()
        percent_incorrect_targets = target_cut / float(decoy_cut)
        return percent_incorrect_targets

    def fdr_with_percent_incorrect_targets(self, cutoff):
        percent_incorrect_targets = self.estimate_percent_incorrect_targets(cutoff)
        return percent_incorrect_targets * self.target_decoy_ratio(cutoff)[0]

    def _calculate_q_values(self):
        session = self.manager.session()
        thresholds = chain.from_iterable(session.query(distinct(GlycopeptideMatch.ms2_score)).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id).order_by(
            GlycopeptideMatch.ms2_score.asc()))

        mapping = {}
        last_score = 1
        last_q_value = 0
        for threshold in thresholds:
            try:
                q_value = self.fdr_with_percent_incorrect_targets(threshold)
                # If a worse score has a higher q-value than a better score, use that q-value
                # instead.
                if last_q_value < q_value and last_score < threshold:
                    q_value = last_q_value
                last_q_value = q_value
                last_score = threshold
                mapping[threshold] = q_value
            except ZeroDivisionError:
                mapping[threshold] = 1.
        return mapping

    def q_values(self):
        session = self.manager.session()

        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.experiment_id == self.target_id)

        q_map = self._calculate_q_values()

        for target in tq:
            target.q_value = q_map[target.ms2_score]
            session.add(target)
        session.commit()
        session.close()

    def run(self):
        thresholds = self.global_thresholds()
        self.p_values()
        self.q_values()
        return thresholds
