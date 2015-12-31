from collections import OrderedDict

from glycresoft_sqlalchemy.data_model import (
    GlycopeptideMatch, HypothesisSampleMatch, PipelineModule,
    Protein, object_session, GlycopeptideSpectrumMatchScore,
    GlycopeptideSpectrumMatch)

from sqlalchemy import func
from matplotlib import pyplot as plt
import numpy as np

from scipy.interpolate import interp1d
from sklearn.metrics import auc


def count_at_score_glycopeptide_match(session, hypothesis_id, hypothesis_sample_match_id):
    targets = session.query(GlycopeptideMatch.q_value, func.count(GlycopeptideMatch.q_value)).group_by(
            GlycopeptideMatch.q_value).order_by(
            GlycopeptideMatch.q_value.desc()).filter(
            (GlycopeptideMatch.hypothesis_sample_match_id == hypothesis_sample_match_id),
            GlycopeptideMatch.protein_id == Protein.id,
            Protein.hypothesis_id == hypothesis_id).all()

    running_total = 0
    targets_above = OrderedDict()
    for score, at in targets:
        targets_above[score] = running_total
        running_total += at

    for score, at in list(targets_above.items()):
        targets_above[score] = running_total - targets_above[score]
    return targets_above


def count_at_score_glycopeptide_spectrum_match(session, hypothesis_id, hypothesis_sample_match_id):
    targets = session.query(GlycopeptideSpectrumMatchScore.value, func.count(
            GlycopeptideSpectrumMatchScore.value)).join(
            GlycopeptideSpectrumMatch).group_by(
            GlycopeptideSpectrumMatchScore.value).order_by(
            GlycopeptideSpectrumMatchScore.value.desc()).filter(
            GlycopeptideSpectrumMatchScore.name == "q_value",
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.hypothesis_id == hypothesis_id).all()

    running_total = 0
    targets_above = OrderedDict()
    for score, at in targets:
        targets_above[score] = running_total
        running_total += at

    for score, at in list(targets_above.items()):
        targets_above[score] = running_total - targets_above[score]
    return targets_above


def detect_rapid_drop(counts, step_size=1, threshold=0.5):
    i = 0
    size = len(counts)
    while i + step_size < size:
        ratio = counts[i + step_size] / float(counts[i])
        if ratio < threshold:
            return i + step_size
        i += 1


def curve_for(score_count_pairs, ax=None, threshold=True, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)

    if threshold is True:
        threshold = detect_rapid_drop(score_count_pairs.values())
    pairs = list(score_count_pairs.items())[threshold:]
    ax.plot(*reversed(zip(*pairs)), **kwargs)
    return ax


class HypothesisSampleMatchComparer(PipelineModule):
    def __init__(self, database_path, *ids):
        self.manager = self.manager_type(database_path)
        self.ids = ids
        self.index = OrderedDict()

    def count_at_score(self, session, hypothesis_id, hypothesis_sample_match_id):
        counts = count_at_score_glycopeptide_match(session, hypothesis_id, hypothesis_sample_match_id)
        return counts

    def add_hypothesis_sample_match(self, hypothesis_sample_match):
        session = object_session(hypothesis_sample_match)
        self.index[hypothesis_sample_match.name] = self.count_at_score(
            session, hypothesis_sample_match.target_hypothesis_id, hypothesis_sample_match.id)

    def add_hypothesis_sample_match_by_id(self, session, id):
        hsm = session.query(HypothesisSampleMatch).get(id)
        self.add_hypothesis_sample_match(hsm)

    def plot(self):
        ax = None
        for series, counts in self.index.items():
            ax = curve_for(counts, ax, label=series, alpha=0.5)
        ax.set_ylim(0, 0.15)
        ticks = np.arange(0, 0.15, 0.01)
        ax.yaxis.set_ticks(ticks)

        maximum = ax.get_xlim()[1]
        step = maximum * 0.025
        ticks = list(map(int, np.arange(0, maximum, step)))

        ax.xaxis.set_ticks(ticks)
        ax.set_xticklabels(ticks, rotation=90, ha='center')

        ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("GSMs Retained")
        ax.set_ylabel("FDR")
        ax.set_title("Analysis Comparison")

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        return ax

    def run(self):
        session = self.manager()
        for id in self.ids:
            self.add_hypothesis_sample_match_by_id(session, id)
        return self.plot()


def _name_to_id(session, *names):
    return [x[0] for x in session.query(HypothesisSampleMatch.id).filter(
        HypothesisSampleMatch.name.in_(list(names))).all()]


class ROCCurvePlot(HypothesisSampleMatchComparer):
    x_axis_max = None
    y_axis_max = None

    def add_hypothesis_sample_match(self, hypothesis_sample_match):
        session = object_session(hypothesis_sample_match)
        self.index[hypothesis_sample_match.name] = ROCSeries(
            self.count_at_score(session, hypothesis_sample_match.decoy_hypothesis_id, hypothesis_sample_match.id),
            self.count_at_score(session, hypothesis_sample_match.target_hypothesis_id, hypothesis_sample_match.id)
            )

    def plot(self):
        ax = None
        # false_max = max(self.index.values(), key=lambda x: x.xmax).xmax
        # true_max = max(self.index.values(), key=lambda x: x.ymax).ymax

        false_max = None
        true_max = None

        fig, ax = plt.subplots(1)

        for title, series in self.index.items():
            curve_path, auc_score = series(xmax=false_max, ymax=true_max)
            ax.plot(*zip(*curve_path), label="%s (%0.2f%%)" % (title, 100 * auc_score), alpha=0.5, lw=2)

        cap = 1.01

        ax.set_xlim(0, cap)
        ax.set_ylim(0, cap)
        ax.plot((0, cap), (0, cap), ls='--', c='black')

        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.set_title("ROC Curve")
        ax.legend(frameon=False, fontsize=8, loc='lower right')
        return ax


class ROCSeries(object):
    def __init__(self, false_positives, true_positives):
        self.x = false_positives
        self.y = true_positives
        self.offset = max(detect_rapid_drop(true_positives.values()), detect_rapid_drop(false_positives.values()))
        self.xmax = float(false_positives.values()[self.offset + 1])
        self.ymax = float(true_positives.values()[self.offset + 1])

    def __repr__(self):
        return "ROCSeries({}, {}, {})".format(self.xmax, self.ymax, self.offset)

    def path(self, xmax=None, ymax=None):
        xf = interp1d(*zip(*self.x.items()))
        yf = interp1d(*zip(*self.y.items()))

        path = []
        # offset = self.offset

        if xmax is None:
            xmax = self.xmax

        if ymax is None:
            ymax = self.ymax

        last_xi = 0
        last_yi = 0

        # ykeys = set(self.y.keys()[offset:])
        # xkeys = set(self.x.keys()[offset:])
        # all_keys = ykeys | xkeys
        # all_keys.remove(None)
        all_keys = np.arange(0, 1., 0.01)

        for i in (all_keys):
            try:
                xi = int(xf(i))
                last_xi = xi
            except TypeError:
                xi = last_xi
            try:
                yi = int(yf(i))
                last_yi = yi
            except TypeError:
                yi = last_yi
            path.append((min(xi/xmax, 1), min(yi/ymax, 1)))
        auc_score = auc(*zip(*path))
        return path, auc_score

    __call__ = path
