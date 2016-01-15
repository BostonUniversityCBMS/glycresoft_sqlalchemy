'''
Deprecated Simplisitic Glycan Implementation
'''

from .modification import AnonymousModificationRule
import glypy
from glypy import Composition

oxonium_ions = {
    "HexNAc": 204.0864,
    "HexNAc-H2O": 186.0754,
    "HexNAc-2H2O": 168.0650,
    "FragmentOfHexNAc": 138.0542,
    "dHex": 146.05791,
    "Hex": 163.0601,
    "HexHexNAc": 366.1394,
    "NeuAc": 292.1026,
    "NeuAc-H2O": 274.0920,
}


class SimpleGlycan(glypy.GlycanComposition):

    def __init__(self, *args, **kwargs):
        super(SimpleGlycan).__init__(self, *args, **kwargs)

    @property
    def mass(self):
        return super(SimpleGlycan, self).mass()

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(self.serialize(), self.mass)))


class GlycosylationSite(object):

    def __init__(self, parent, position, glycosylation=None):
        self.parent = parent
        self.position = position
        self.glycosylation = glycosylation
        self.mass = glycosylation.mass()

    def __repr__(self):
        rep = 'GlycosylationSite(%d, %s)' % (self.position, self.glycosylation)
        return rep

    def __eq__(self, other):
        try:
            return (self.position == other.position) and (self.glycosylation == other.glycosylation)
        except:
            return False

    def __ne__(self, other):
        try:
            return (self.position != other.position) and (self.glycosylation != other.glycosylation)
        except:
            return True

    def __hash__(self):
        return hash((self.position, self.glycosylation))


class MassableValueDict(dict):
    def total_composition(self):
        total = Composition()
        for key, value in self.items():
            total += value.total_composition()
        return total

    def mass(self, charge=0, average=False, mass_data=None):
        if charge == 0:
            total = 0.
            for key, value in self.items():
                total += value.mass(average=average, mass_data=mass_data)
            return total
        else:
            return self.total_composition().calc_mass(
                charge=charge, average=average, mass_data=mass_data)


class GlycosylationManager(MassableValueDict):
    def __init__(self, parent, *glycosites):
        self.parent = parent
        for glycosite in glycosites:
            self[glycosite.position] = glycosite
