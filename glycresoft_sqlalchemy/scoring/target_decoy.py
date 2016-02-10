import logging
try:
    logger = logging.getLogger("target_decoy")
except:
    pass
from itertools import chain
from collections import OrderedDict, defaultdict

from sqlalchemy import func, distinct
from ..data_model import (
    GlycopeptideMatch, Protein, PipelineModule, HypothesisSampleMatch,
    GlycopeptideSpectrumMatch, GlycopeptideSpectrumMatchScore)


ms2_score = GlycopeptideMatch.ms2_score


class TargetDecoyAnalyzer(PipelineModule):
    @classmethod
    def from_hypothesis_sample_match(cls, database_path, hsm, **kwargs):
        return cls(database_path, hsm.target_hypothesis_id, hsm.decoy_hypothesis_id, hsm.id, **kwargs)

    def __init__(self, database_path, target_hypothesis_id=None, decoy_hypothesis_id=None,
                 hypothesis_sample_match_id=None, with_pit=False, **kwargs):
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.target_id = target_hypothesis_id
        self.decoy_id = decoy_hypothesis_id
        self.with_pit = with_pit

        session = self.manager.session()

        self.target_count = session.query(
            GlycopeptideMatch.ms2_score).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id,
            GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).count()

        self.decoy_count = session.query(
            GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id,
            GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).count()

        hsm = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
        if hsm is not None:
            hsm.parameters["target_decoy"] = {
                "with_pit": with_pit,
                "target_hypothesis_id": target_hypothesis_id,
                "decoy_hypothesis_id": decoy_hypothesis_id,
            }
            session.add(hsm)
            session.commit()

        session.close()
        self.n_targets_at = {}
        self.n_decoys_at = {}

    def calculate_thresholds(self):
        session = self.manager.session()
        targets = session.query(GlycopeptideMatch.ms2_score, func.count(GlycopeptideMatch.ms2_score)).group_by(
            GlycopeptideMatch.ms2_score).order_by(
            GlycopeptideMatch.ms2_score.asc()).filter(
            (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id),
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id).all()

        running_total = 0
        targets_above = OrderedDict()
        for score, at in targets:
            targets_above[score] = running_total
            running_total += at

        for score, at in list(targets_above.items()):
            targets_above[score] = running_total - targets_above[score]

        decoys = session.query(GlycopeptideMatch.ms2_score, func.count(GlycopeptideMatch.ms2_score)).group_by(
            GlycopeptideMatch.ms2_score).order_by(
            GlycopeptideMatch.ms2_score.asc()).filter(
            (GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id),
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id).all()

        running_total = 0
        decoys_above = OrderedDict()
        for score, at in decoys:
            decoys_above[score] = running_total
            running_total += at

        for score, at in list(decoys_above.items()):
            decoys_above[score] = running_total - decoys_above[score]

        self.n_targets_at = targets_above
        self.n_decoys_at = decoys_above

        return targets_above, decoys_above

    def n_decoys_above_threshold(self, threshold):
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

    def n_targets_above_threshold(self, threshold):
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

    def target_decoy_ratio(self, cutoff):

        decoys_at = self.n_decoys_above_threshold(cutoff)
        targets_at = self.n_targets_above_threshold(cutoff)
        try:
            ratio = decoys_at / float(targets_at)
        except ZeroDivisionError:
            ratio = 1.
        return ratio, targets_at, decoys_at

    def p_values(self):
        logger.info("Computing p-values")
        session = self.manager.session()

        tq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id,
            GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).order_by(
            GlycopeptideMatch.ms2_score.desc())
        dq = session.query(GlycopeptideMatch).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.decoy_id,
            GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id)

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
                decoys_at = self.n_decoys_above_threshold(target.ms2_score)
                last_score = target.ms2_score

                last_p_value = decoys_at / total_decoys
                target.p_value = last_p_value
                session.add(target)
        session.commit()
        session.close()

    def estimate_percent_incorrect_targets(self, cutoff):
        target_cut = self.target_count - self.n_targets_above_threshold(cutoff)
        decoy_cut = self.decoy_count - self.n_decoys_above_threshold(cutoff)
        percent_incorrect_targets = target_cut / float(decoy_cut)
        return percent_incorrect_targets

    def fdr_with_percent_incorrect_targets(self, cutoff):
        if self.with_pit:
            percent_incorrect_targets = self.estimate_percent_incorrect_targets(cutoff)
        else:
            percent_incorrect_targets = 1.0
        return percent_incorrect_targets * self.target_decoy_ratio(cutoff)[0]

    def _calculate_q_values(self):
        session = self.manager.session()
        thresholds = chain.from_iterable(session.query(distinct(GlycopeptideMatch.ms2_score)).filter(
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == self.target_id,
            GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).order_by(
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
        for k in sorted(q_map):
            logger.info("Updating entries with score %r -> %f", k, q_map[k])
            session.query(GlycopeptideMatch).filter(
                GlycopeptideMatch.ms2_score == k,
                GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).update(
                {"q_value": q_map[k]}, synchronize_session=False)

        session.commit()
        session.close()

    def run(self):
        self.calculate_thresholds()
        self.q_values()


class TargetDecoySpectrumMatchAnalyzer(TargetDecoyAnalyzer):
    def __init__(self, database_path, target_hypothesis_id=None, decoy_hypothesis_id=None,
                 hypothesis_sample_match_id=None, with_pit=False, score="simple_ms2_score", **kwargs):
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.target_id = target_hypothesis_id
        self.decoy_id = decoy_hypothesis_id
        self.with_pit = with_pit

        self.score = score

        session = self.manager.session()

        self.target_count = session.query(
            GlycopeptideSpectrumMatch.id).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.target_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).count()

        self.decoy_count = session.query(
            GlycopeptideSpectrumMatch.id).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.decoy_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).count()

        hsm = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
        if hsm is not None:
            hsm.parameters["target_decoy"] = {
                "with_pit": with_pit,
                "target_hypothesis_id": target_hypothesis_id,
                "decoy_hypothesis_id": decoy_hypothesis_id,
                "score": score
            }
            session.add(hsm)
            session.commit()

        session.close()
        self.n_targets_at = {}
        self.n_decoys_at = {}

    def calculate_thresholds(self):
        session = self.manager.session()
        targets = session.query(GlycopeptideSpectrumMatchScore.value, func.count(
            GlycopeptideSpectrumMatchScore.value)).join(
            GlycopeptideSpectrumMatch).group_by(
            GlycopeptideSpectrumMatchScore.value).order_by(
            GlycopeptideSpectrumMatchScore.value.asc()).filter(
            GlycopeptideSpectrumMatch.best_match,
            GlycopeptideSpectrumMatchScore.name == self.score,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.hypothesis_id == self.target_id).all()

        running_total = 0
        targets_above = OrderedDict()
        for score, at in targets:
            targets_above[score] = running_total
            running_total += at

        for score, at in list(targets_above.items()):
            targets_above[score] = running_total - targets_above[score]

        decoys = session.query(GlycopeptideSpectrumMatchScore.value, func.count(
            GlycopeptideSpectrumMatchScore.value)).join(
            GlycopeptideSpectrumMatch).group_by(
            GlycopeptideSpectrumMatchScore.value).order_by(
            GlycopeptideSpectrumMatchScore.value.asc()).filter(
            GlycopeptideSpectrumMatch.best_match,
            GlycopeptideSpectrumMatchScore.name == self.score,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.hypothesis_id == self.decoy_id).all()

        running_total = 0
        decoys_above = OrderedDict()
        for score, at in decoys:
            decoys_above[score] = running_total
            running_total += at

        for score, at in list(decoys_above.items()):
            decoys_above[score] = running_total - decoys_above[score]

        self.n_targets_at = targets_above
        self.n_decoys_at = decoys_above

        return targets_above, decoys_above

    def n_decoys_above_threshold(self, threshold):
        if threshold in self.n_decoys_at:
            return self.n_decoys_at[threshold]
        else:
            session = self.manager.session()
            self.n_decoys_at[threshold] = session.query(
                GlycopeptideSpectrumMatch.id).join(GlycopeptideSpectrumMatchScore).filter(
                GlycopeptideSpectrumMatch.best_match,
                GlycopeptideSpectrumMatchScore.name == self.score,
                GlycopeptideSpectrumMatch.hypothesis_id == self.decoy_id).filter(
                GlycopeptideSpectrumMatchScore.value >= threshold,
                (GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id)
                ).count()
            return self.n_decoys_at[threshold]

    def n_targets_above_threshold(self, threshold):
        if threshold in self.n_targets_at:
            return self.n_targets_at[threshold]
        else:
            session = self.manager.session()
            self.n_targets_at[threshold] = session.query(
                GlycopeptideSpectrumMatch.id).join(GlycopeptideSpectrumMatchScore).filter(
                GlycopeptideSpectrumMatch.best_match,
                GlycopeptideSpectrumMatchScore.name == self.score,
                GlycopeptideSpectrumMatch.hypothesis_id == self.target_id).filter(
                GlycopeptideSpectrumMatchScore.value >= threshold,
                (GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id)
                ).count()
            return self.n_targets_at[threshold]

    def target_decoy_ratio(self, cutoff):

        decoys_at = self.n_decoys_above_threshold(cutoff)
        targets_at = self.n_targets_above_threshold(cutoff)
        try:
            ratio = decoys_at / float(targets_at)
        except ZeroDivisionError:
            ratio = 1.
        return ratio, targets_at, decoys_at

    def estimate_percent_incorrect_targets(self, cutoff):
        target_cut = self.target_count - self.n_targets_above_threshold(cutoff)
        decoy_cut = self.decoy_count - self.n_decoys_above_threshold(cutoff)
        percent_incorrect_targets = target_cut / float(decoy_cut)
        return percent_incorrect_targets

    def fdr_with_percent_incorrect_targets(self, cutoff):
        if self.with_pit:
            percent_incorrect_targets = self.estimate_percent_incorrect_targets(cutoff)
        else:
            percent_incorrect_targets = 1.0
        return percent_incorrect_targets * self.target_decoy_ratio(cutoff)[0]

    def _calculate_q_values(self):
        session = self.manager.session()
        thresholds = chain.from_iterable(session.query(distinct(
            GlycopeptideSpectrumMatchScore.value)).join(
            GlycopeptideSpectrumMatch).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.target_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            GlycopeptideSpectrumMatchScore.name == self.score).order_by(
            GlycopeptideSpectrumMatchScore.value.asc()))

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
        accumulator = []
        for k in sorted(q_map):
            ids = [x[0] for x in session.query(GlycopeptideSpectrumMatch.id).join(
                GlycopeptideSpectrumMatchScore).filter(
                GlycopeptideSpectrumMatchScore.value == k,
                GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                GlycopeptideSpectrumMatchScore.name == self.score)]
            accumulator.extend(
                GlycopeptideSpectrumMatchScore(
                    name='q_value', value=q_map[k], spectrum_match_id=i) for i in ids)
            if len(accumulator) > 200:
                logger.info("Updating entries with score %f -> %f", k, q_map[k])
                session.bulk_save_objects(accumulator)
                accumulator = []

        session.bulk_save_objects(accumulator)
        accumulator = []

        session.commit()
        session.close()

    def run(self):
        self.calculate_thresholds()
        self.q_values()


class InMemoryTargetDecoyAnalyzer(object):
    '''
    A work in progress
    '''
    def __init__(self, target_series, decoy_series, with_pit=True):
        self.targets = target_series
        self.decoys = decoy_series
        self.target_count = len(target_series)
        self.decoy_count = len(decoy_series)
        self.with_pit = with_pit
        self.calculate_thresholds()

    def calculate_thresholds(self):
        self.n_targets_at = {}
        self.n_decoys_at = {}

        target_series = self.targets
        decoy_series = self.decoys

        thresholds = sorted({case.ms2_score for case in target_series} | {case.ms2_score for case in decoy_series})
        g = iter(thresholds)

        target_counts = defaultdict(int)
        decoy_counts = defaultdict(int)

        count = 0
        current_threshold = g.next()
        for series, counter in [(target_series, target_counts), (decoy_series, decoy_counts)]:
            g = iter(thresholds)
            count = 0
            current_threshold = g.next()
            for case in series:
                if case.ms2_score <= current_threshold:
                    count += 1
                else:
                    counter[current_threshold] = count
                    current_threshold = g.next()
                    while case.ms2_score > current_threshold:
                        counter[current_threshold] = count
                        current_threshold = g.next()
                    count += 1
            for t in g:
                counter[current_threshold] = count

        self.n_targets_at = {}
        for t, c in target_counts.items():
            self.n_targets_at[t] = self.target_count - c

        self.n_decoys_at = {}
        for t, c in decoy_counts.items():
            self.n_decoys_at[t] = self.decoy_count - c

        self.thresholds = thresholds

    def n_decoys_above_threshold(self, threshold):
        return self.n_decoys_at.get(threshold, 0.)

    def n_targets_above_threshold(self, threshold):
        return self.n_targets_at.get(threshold, 0.)

    def target_decoy_ratio(self, cutoff):

        decoys_at = self.n_decoys_above_threshold(cutoff)
        targets_at = self.n_targets_above_threshold(cutoff)
        try:
            ratio = decoys_at / float(targets_at)
        except ZeroDivisionError:
            ratio = 1.
        return ratio, targets_at, decoys_at

    def estimate_percent_incorrect_targets(self, cutoff):
        target_cut = self.target_count - self.n_targets_above_threshold(cutoff)
        decoy_cut = self.decoy_count - self.n_decoys_above_threshold(cutoff)
        percent_incorrect_targets = target_cut / float(decoy_cut)

        return percent_incorrect_targets

    def fdr_with_percent_incorrect_targets(self, cutoff):
        if self.with_pit:
            percent_incorrect_targets = self.estimate_percent_incorrect_targets(cutoff)
        else:
            percent_incorrect_targets = 1.0
        return percent_incorrect_targets * self.target_decoy_ratio(cutoff)[0]

    def _calculate_q_values(self):
        thresholds = self.thresholds
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
        q_map = self._calculate_q_values()
        for target in self.targets:
            target.q_value2 = q_map[target.ms2_score]
        for decoy in self.decoys:
            decoy.q_value2 = q_map[decoy.ms2_score]
