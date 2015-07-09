'''
Deprecated Simplisitic Glycan Implementation
'''

import csv
import re

from collections import OrderedDict
from .composition import Composition
from .modification import AnonymousModificationRule


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


class SimpleGlycan(object):

    @staticmethod
    def glycan_composition_string(count_list):
        return "[{0}]".format(";".join(map(lambda x: str(int(x)), count_list)))

    def __init__(self, mass, composition, glycan_identities=None, **kwargs):
        if isinstance(composition, basestring):
            composition = map(int, composition[1:-1].split(";"))
        self.composition = composition
        self.mass = mass
        self.glycan_identities = glycan_identities
        self.data = kwargs

    @property
    def molecular_weight(self):
        return self.mass

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(map(str, self.composition))), self.mass)

    def __repr__(self):
        return "Glycan({composition} {mass} {data})".format(**self.__dict__)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        comp_eq = all(i == j for i, j in zip(self.composition, other.composition))
        return comp_eq

    def __ne__(self, other):
        comp_eq = any(i != j for i, j in zip(self.composition, other.composition))
        return comp_eq

    # @property
    # def __dict__(self):
    #     return {n: getattr(self, n) for n in self.__slots__}

    @property
    def composition_dict(self):
        composition_dict = OrderedDict()
        iter_identities = iter(self.glycan_identities)
        for count in self.composition:
            composition_dict[iter_identities.next()] = count
        return composition_dict

    def dehydrate(self):
        if "Water" in self.glycan_identities:
            try:
                i = self.glycan_identities.index("Water")
                if self.composition[i] > 0:
                    self.composition[i] -= 1
            except Exception, e:
                print e
            self.mass -= Composition("H2O").mass
        return self

    def deadduct(self):
        if "adduct" in self.data:
            self.mass -= self.data['adduct'].mass
        return self


class Adduct(object):

    @classmethod
    def from_csv(cls, row):
        symbol, mass = row["Adduct/Replacement"].split("/")
        mass = float(mass)
        amount = int(row["Adduct Amount"])
        return cls(symbol, mass, amount)

    def __init__(self, symbol, singleton_mass, amount):
        self.symbol = symbol
        self.amount = amount
        self.mass = singleton_mass * amount


class GlycanHypothesis(list):
    def __init__(self, csv_path=None):
        self.glycans = []
        if csv_path is not None:
            self.reader = csv.DictReader(open(csv_path))
            self.glycan_identities = self.get_glycan_identities(self.reader.fieldnames)
            self.extend(self.reader)
        self.extract_glycans()

    def extract_glycans(self):
        self.glycans = list({SimpleGlycan(float(r["Molecular Weight"]),
                                    r["Compositions"],
                                    glycan_identities=self.glycan_identities) for r in self})
        return self.glycans

    def get_glycan_identities(self, colnames):
        glycan_identity = []
        extract_state = False
        for col in colnames:
            if col == "Compositions":
                extract_state = True
                continue
            elif col == "Adduct/Replacement":
                extract_state = False
                break
            elif extract_state:
                glycan_identity.append(col.replace("G:", ""))
        return glycan_identity


class GlycanMS1Results(list):
    def __init__(self, csv_path=None):
        self.glycans = []
        if csv_path is not None:
            self.reader = csv.DictReader(open(csv_path))
            self.glycan_identities = self.get_glycan_identities(self.reader.fieldnames)
            self.extend(self.reader)
        self.extract_glycans()

    def get_glycan_identities(self, colnames):
        glycan_identity = []
        extract_state = False
        for col in colnames:
            if col == "Hypothesis MW":
                extract_state = True
                continue
            elif col == "Adduct/Replacement":
                extract_state = False
                break
            elif extract_state:
                glycan_identity.append(col.replace("G:", ""))
        return glycan_identity

    def extract_glycans(self):
        self.glycans = list({SimpleGlycan(float(r["Hypothesis MW"]),
                                    r["Compositions"],
                                    score=float(r["Score"]),
                                    volume=float(r["Total Volume"]),
                                    glycan_identities=self.glycan_identities,
                                    adduct=Adduct.from_csv(r)) for r in self
                            if r["Compositions"] != "0"})
        return self.glycans


def pack_glycan_string(glycan_map, glycan_order):
    pack_string = "[{0}]".format(";".join([glycan_map[g] for g in glycan_order]))
    return pack_string


def parse_glycan_pack_string(pack_string, glycan_order):
    glycan_map = {g: int(c) for g, c in zip(glycan_order, re.findall('\\b\\d+\\b', pack_string))}
    return glycan_map
