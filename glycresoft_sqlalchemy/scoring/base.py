import abc

from glycresoft_sqlalchemy.utils import simple_repr
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatchScore


class ScorerBase(object):
    __meta__ = abc.ABCMeta

    def __init__(self, score_name, attr_name=None):
        self.score_name = score_name
        if attr_name is None:
            attr_name = score_name
        self.attr_name = attr_name

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    __repr__ = simple_repr

    @abc.abstractmethod
    def evaluate(self, *args, **kwargs):
        raise NotImplementedError()


class GlycopeptideSpectrumMatchScorer(ScorerBase):
    def __call__(self, *args, **kwargs):
        value = self.evaluate(*args, **kwargs)
        if value is None:
            value = getattr(args[0], self.attr_name)
        return GlycopeptideSpectrumMatchScore(name=self.score_name, value=value, spectrum_match_id=args[0].id)
