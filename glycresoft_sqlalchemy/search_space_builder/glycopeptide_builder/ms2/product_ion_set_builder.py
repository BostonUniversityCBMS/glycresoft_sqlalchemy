import re
import itertools
import functools
import multiprocessing
import operator
import logging

from glycresoft_sqlalchemy.structure import sequence, constants
from glycresoft_sqlalchemy.utils import collectiontools
from glycresoft_sqlalchemy.data_model import (
    TheoreticalPeptideProductIon, TheoreticalGlycopeptideStubIon,
    PipelineModule, TheoreticalGlycopeptide, ProductIonSet)

Sequence = sequence.Sequence
kind_getter = operator.attrgetter("kind")

logger = logging.getLogger("fragment_generation")

