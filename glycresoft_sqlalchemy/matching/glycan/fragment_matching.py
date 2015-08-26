import math
import multiprocessing
import functools
import itertools
import logging
import operator
try:
    logger = logging.getLogger(__name__)
except:
    pass
from os.path import splitext
from collections import defaultdict

from ..spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser

from ..data_model import (
    MS2GlycanHypothesis, MS2GlycanHypothesisSampleMatch, PipelineModule, MSMSSqlDB, )


neutral_mass_getter = operator.attrgetter("neutral_mass")
key_getter = operator.itemgetter('key')

decon_format_lookup = {
    "bupid_yaml": BUPIDMSMSYamlParser,
    "db": MSMSSqlDB,
}

PROTON = 1.007276035

ms1_tolerance_default = 1e-5
ms2_tolerance_default = 2e-5


def ppm_error(x, y):
    return (x - y) / y
