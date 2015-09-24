import operator
import itertools
from glycresoft_sqlalchemy.utils import common_math
from glycresoft_sqlalchemy.structure import sequence

Sequence = sequence.Sequence
chain_iterable = itertools.chain.from_iterable
get_intensity = operator.attrgetter("intensity")

DPeak = common_math.DPeak
intensity_ratio_function = common_math.intensity_ratio_function
intensity_rank = common_math.intensity_rank
ppm_error = common_math.ppm_error
MassOffsetFeature = common_math.MassOffsetFeature
search_spectrum = common_math.search_spectrum
MatchedSpectrum = common_math.MatchedSpectrum
