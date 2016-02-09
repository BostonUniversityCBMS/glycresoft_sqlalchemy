'''
Derived from the code sklearn.linear_model.

Used to avoid needing to reconstitute the complete `LogisticRegression` instance
to score data based on previously fitted results, as the sklearn implementation
requires some state configuration done in the fitting process in order to carry out
the classification on any data.
'''
import numpy as np
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.linear_model import LogisticRegression

from glycresoft_sqlalchemy.data_model import PeakGroupScoringModel


class LogisticModelScorer(object):
    def __init__(self, coef=None, intercept=None, labels=None):
        self.coef = coef
        self.intercept = intercept
        self.labels = labels

    # Much taken from sklearn's LogisticRegression

    def decision_function(self, X):
        if X.shape[1] != self.coef.shape[1]:
            raise ValueError(
                "X has %d features per sample; expecting %d" % (X.shape[1], self.coef.shape[1]))
        scores = safe_sparse_dot(X, self.coef.T, dense_output=True) + self.intercept
        return scores.ravel() if scores.shape[1] == 1 else scores

    def _predict_proba_lr(self, X):
        X = np.asanyarray(X)
        prob = self.decision_function(X)
        prob *= -1
        np.exp(prob, prob)
        prob += 1
        np.reciprocal(prob, prob)
        if len(prob.shape) == 1:
            return np.vstack([1 - prob, prob]).T
        else:
            # OvR normalization, like LibLinear's predict_probability
            prob /= prob.sum(axis=1).reshape((prob.shape[0], -1))
            return prob

    def predict_proba(self, X):
        return self._predict_proba_lr(X)

    def fit(self, X, y, **kwargs):
        lr = LogisticRegression(**kwargs)
        lr.fit(X, y)
        self.intercept = lr.intercept_
        self.coef = lr.coef_

    __call__ = predict_proba

    evaluate = predict_proba


def from_peak_group_scoring_model(model):
    return LogisticModelScorer(model.to_parameter_vector(), model.intercept)


def to_peak_group_scoring_model(scorer):
    inst = PeakGroupScoringModel.from_parameter_vector(scorer.coef)
    inst.intercept = scorer.intercept
    return inst
