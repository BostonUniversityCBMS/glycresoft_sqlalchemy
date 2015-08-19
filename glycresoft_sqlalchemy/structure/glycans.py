'''
Deprecated Simplisitic Glycan Implementation
'''

import re

from collections import OrderedDict
from .composition import Composition
from .modification import AnonymousModificationRule
import glypy

oxonium_ions = {
    "HexNAc": 204.0864,
    "HexNAc-H2O": 186.0754,
    "HexNAc-2H2O": 168.0650,
    "FragmentOfHexNAc": 138.0542,
    "dHex": 146.05791,
    "Hex": 163.0601,
    "HexHexNAc": 366.1394,
    "NeuAc": 292.1026,
    "Neu5Ac-H2O": 274.0920,
}


class SimpleGlycan(glypy.GlycanComposition):

    def __init__(self, *args, **kwargs):
        super(SimpleGlycan).__init__(self, *args, **kwargs)

    @property
    def molecular_weight(self):
        return self.mass

    @property
    def mass(self):
        return super(SimpleGlycan, self).mass()

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(self.serialize(), self.mass)))
