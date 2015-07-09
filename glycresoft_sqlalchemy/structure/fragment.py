import re

from .modification import Modification
from .composition import Composition
from ..utils.collectiontools import descending_combination_counter

fragment_pairing = {
    "a": "x",
    "A": "X",
    "b": "y",
    "B": "Y",
    "c": "z",
    "C": "Z",
    "x": "a",
    "X": "A",
    "y": "b",
    "Y": "B",
    "z": "c",
    "Z": "C",
}

fragment_shift ={
    'b': Composition('H+').mass,
    'B': Composition('H+').mass,
    'y': Composition('H2O').mass + Composition('H+').mass,
    'Y': Composition('H2O').mass + Composition('H+').mass
}

fragment_direction = {
    "a": 1,
    "A": 1,
    "b": 1,
    "B": 1,
    "c": 1,
    "C": 1,
    "x": -1,
    "X": -1,
    "y": -1,
    "Y": -1,
    "z": -1,
    "Z": -1,
}


class Fragment(object):
    """Glycopeptide Fragment"""

    parser = re.compile("(?P<kind>[abcxyzABCXYZ])(?P<position>[0-9]+)(?P<modificaiton>\+.*)?")

    @classmethod
    def parse(cls, frag_name):
        matches = cls.parser.search(frag_name)
        data = matches.groupdict()
        return data

    @classmethod
    def modification_name_block(cls, mod_dict):
        name_block = []
        for mod_name in cls.concerned_mods:
            if mod_name in mod_dict:
                if mod_dict[mod_name] > 1:
                    name_block.extend(['+', str(mod_dict[mod_name]), mod_name])
                elif mod_dict[mod_name] == 1:
                    name_block.extend(['+', mod_name])

        return "".join(name_block)

    concerned_mods = ['HexNAc']

    def __init__(self, frag_type, pos, mod_dict, mass, golden_pairs=None):
        if golden_pairs is None:
            golden_pairs = []
        self.type = frag_type[0].lower() + frag_type[1:]
        # The mass value is the bare backbone's mass
        self.bare_mass = mass

        self.mass = mass

        self.pos = pos
        self.mod_dict = mod_dict

        for key, value in self.mod_dict.items():
            self.mass += Modification(key).mass * value

        self.golden_pairs = golden_pairs

    def get(self):
        """Simply return string like B2, Y3 with no modificaiton information."""
        fragment_name = []
        fragment_name.append(self.type)
        fragment_name.append(str(self.pos))
        return ''.join(fragment_name)

    def get_mass(self):
        return self.mass

    def get_modification_number(self, mod_name):
        if mod_name in self.mod_dict:
            return self.mod_dict[mod_name]
        else:
            return 0

    def get_modifications(self):
        for name, number in self.mod_dict.items():
            yield [name, number]

    def get_fragment_name(self):
        """Connect the information into string."""
        fragment_name = []
        fragment_name.append(self.type)
        fragment_name.append(str(self.pos))

        # Only concerned modifications are reported.
        concerned_mod = ['HexNAc']
        for mod_name in concerned_mod:
            if mod_name in self.mod_dict:
                if self.mod_dict[mod_name] > 1:
                    fragment_name.extend(['+', str(self.mod_dict[mod_name]), mod_name])
                elif self.mod_dict[mod_name] == 1:
                    fragment_name.extend(['+', mod_name])
                else:
                    pass

        return ''.join(fragment_name)

    def partial_loss(self, modifications=None):
        if modifications is None:
            modifications = self.concerned_mods
        mods = dict(self.mod_dict)
        mods_of_interest = {k: v for k, v in mods.items() if k in modifications}

        other_mods = {k: v for k, v in mods.items() if k not in modifications}
        for mod in descending_combination_counter(mods_of_interest):
            other_mods.update({k: v for k, v in mod.items() if v != 0})
            yield Fragment(self.type, self.pos, dict(other_mods), self.bare_mass, golden_pairs=self.golden_pairs)

    def to_json(self):
        d = dict(self.__dict__)
        d['key'] = self.name

    @property
    def name(self):
        return self.get_fragment_name()

    def __repr__(self):
        return "Fragment(%(type)s @ %(pos)s %(mass)s [%(mod_dict)s])" % self.__dict__
