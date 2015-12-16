from glycresoft_sqlalchemy.utils import simple_repr
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatchScore


class ScorerBase(object):
    def __init__(self, scorer, score_name, attr_name=None):
        self.scorer = scorer
        self.score_name = score_name
        if attr_name is None:
            attr_name = score_name
        self.attr_name = attr_name

    def __call__(self, *args, **kwargs):
        return self.scorer(*args, **kwargs)

    __repr__ = simple_repr


class GlycopeptideSpectrumMatchScorer(ScorerBase):
    def __call__(self, *args, **kwargs):
        value = self.scorer(*args, **kwargs)
        if value is None:
            value = getattr(args[0], self.attr_name)
        return GlycopeptideSpectrumMatchScore(name=self.score_name, value=value, spectrum_match_id=args[0].id)
