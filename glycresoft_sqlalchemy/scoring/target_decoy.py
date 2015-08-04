import logging
try:
    logger = logging.getLogger("target_decoy")
except:
    pass
from itertools import chain
from collections import namedtuple, Counter, defaultdict

from sqlalchemy import func, distinct
from ..data_model import DatabaseManager, GlycopeptideMatch, Protein
from ..data_model import PipelineModule


Threshold = namedtuple("Threshold", ("score", "targets", "decoys", "fdr"))

ms2_score = GlycopeptideMatch.ms2_score
p_value = GlycopeptideMatch.p_value


class RangeCounter(Counter):
    def add_below(self, key, value):
        for pkey in list(self.keys()):
            if pkey <= key:
                self[pkey] += value

    def add_above(self, key, value):
        for pkey in list(self.keys()):
            if pkey > key:
                self[pkey] += value


class TargetDecoyAnalyzer(PipelineModule):
    manager_type = DatabaseManager

    def __init__(self, database_path, target_hypothesis_id=None, decoy_hypothesis_id=None,
                 hypothesis_sample_match_id=None):
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.target_id = target_hypothesis_id
        self.decoy_id = decoy_hypothesis_id

        session = self.manager.session()

        self.target_count = session.query(
            GlycopeptideMatch.ms2_score).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id).count()

        self.decoy_count = session.query(
            GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id).count()

        session.close()
        self.n_targets_at = {}
        self.n_decoys_at = {}

    def calculate_thresholds(self):
        session = self.manager.session()
        targets = session.query(GlycopeptideMatch.ms2_score, func.count(GlycopeptideMatch.ms2_score)).group_by(
            GlycopeptideMatch.ms2_score).order_by(
            GlycopeptideMatch.ms2_score.asc()).filter(
            (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id) if
            self.hypothesis_sample_match_id is not None else True,
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id).all()

        running_total = 0
        targets_above = {}
        for score, at in targets:
            targets_above[score] = running_total
            running_total += at

        for score, at in list(targets_above.items()):
            targets_above[score] = running_total - targets_above[score]

        decoys = session.query(GlycopeptideMatch.ms2_score, func.count(GlycopeptideMatch.ms2_score)).group_by(
            GlycopeptideMatch.ms2_score).order_by(
            GlycopeptideMatch.ms2_score.asc()).filter(
            (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id) if
            self.hypothesis_sample_match_id is not None else True,
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id).all()

        running_total = 0
        decoys_above = {}
        for score, at in decoys:
            decoys_above[score] = running_total
            running_total += at

        for score, at in list(decoys_above.items()):
            decoys_above[score] = running_total - decoys_above[score]

        self.n_targets_at = targets_above
        self.n_decoys_at = decoys_above

        return targets_above, decoys_above

    def calculate_n_decoys_at(self, threshold):
        if threshold in self.n_decoys_at:
            return self.n_decoys_at[threshold]
        else:
            session = self.manager.session()
            self.n_decoys_at[threshold] = session.query(
                GlycopeptideMatch).filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == self.decoy_id).filter(
                GlycopeptideMatch.ms2_score >= threshold,
                (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id) if
                self.hypothesis_sample_match_id is not None else True).count()
            return self.n_decoys_at[threshold]

    def calculate_n_targets_at(self, threshold):
        if threshold in self.n_targets_at:
            return self.n_targets_at[threshold]
        else:
            session = self.manager.session()
            self.n_targets_at[threshold] = session.query(
                GlycopeptideMatch).filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == self.target_id).filter(
                GlycopeptideMatch.ms2_score >= threshold,
                (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id) if
                self.hypothesis_sample_match_id is not None else True).count()
            return self.n_targets_at[threshold]

    def target_decoy_ratio(self, cutoff, score=ms2_score):

        decoys_at = self.calculate_n_decoys_at(cutoff)
        targets_at = self.calculate_n_targets_at(cutoff)
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

    def p_values(self):
        logger.info("Computing p-values")
        session = self.manager.session()

        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id).order_by(
            GlycopeptideMatch.ms2_score.desc())
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id)

        total_decoys = float(dq.count())
        if total_decoys == 0:
            raise ValueError("No decoy matches found")
        last_score = 0
        last_p_value = 0
        for target in tq:
            if target.ms2_score == last_score:
                target.p_value = last_p_value
                session.add(target)
            else:
                session.commit()
                decoys_at = self.calculate_n_decoys_at(target.ms2_score)
                last_score = target.ms2_score

                last_p_value = decoys_at / total_decoys
                target.p_value = last_p_value
                session.add(target)
        session.commit()
        session.close()

    def estimate_percent_incorrect_targets(self, cutoff, score=ms2_score):
        session = self.manager.session()

        target_cut = self.target_count - self.calculate_n_targets_at(cutoff)
        decoy_cut = self.decoy_count - self.calculate_n_decoys_at(cutoff)
        percent_incorrect_targets = target_cut / float(decoy_cut)
        session.close()
        return percent_incorrect_targets

    def fdr_with_percent_incorrect_targets(self, cutoff):
        percent_incorrect_targets = self.estimate_percent_incorrect_targets(cutoff)
        return percent_incorrect_targets * self.target_decoy_ratio(cutoff)[0]

    def _calculate_q_values(self):
        session = self.manager.session()
        thresholds = chain.from_iterable(session.query(distinct(GlycopeptideMatch.ms2_score)).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id).order_by(
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
        session.close()
        return mapping

    def q_values(self):
        logger.info("Computing q-values")
        session = self.manager.session()

        q_map = self._calculate_q_values()
        for k in q_map:
            logger.info("Updating entries with score %f -> %f", k, q_map[k])
            session.query(GlycopeptideMatch).filter(
                GlycopeptideMatch.ms2_score == k).update(
                {"q_value": q_map[k]}, synchronize_session=False)

        session.commit()
        session.close()

    def run(self):
        self.calculate_thresholds()
        self.q_values()
