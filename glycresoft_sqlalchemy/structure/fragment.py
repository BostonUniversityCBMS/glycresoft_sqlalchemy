import re
from collections import defaultdict

from .modification import Modification
from .composition import Composition
from ..utils.collectiontools import descending_combination_counter
from ..utils import simple_repr

fragment_pairing = {
    "a": "x",
    "b": "y",
    "c": "z",
    "x": "a",
    "y": "b",
    "z": "c",
}

fragment_shift = {
    'b': 0.,
    'y': Composition('H2O').mass
}

fragment_direction = {
    "a": 1,
    "b": 1,
    "c": 1,
    "x": -1,
    "y": -1,
    "z": -1,
}

generic_neutral_losses_composition = {
    "-NH3": -Composition("NH3"),
    "-H2O": -Composition("H2O"),
    "-NH3-NH3": -Composition("(NH3)2"),
    "-H2O-H2O": -Composition("(H2O)2"),
    "-NH3-H2O": -Composition("NH3H2O")
}


class NeutralLoss(object):
    def __init__(self, name, composition=None):
        if composition is None:
            composition = generic_neutral_losses_composition[name]
        self.name = name
        self.composition = composition
        self.mass = composition.mass

    def clone(self):
        return self.__class__(self.name, self.composition.clone())

    def __str__(self):
        return self.name

    def __repr__(self):
        return "NeutralLoss(name=%r)" % self.name


class FragmentBase(object):
    _neutral_loss = None
    _glycosylation = None
    _name = None

    def get_neutral_loss(self):
        return self._neutral_loss

    def set_neutral_loss(self, neutral_loss):
        if self._neutral_loss is not None:
            self.mass -= self._neutral_loss.mass
        self._neutral_loss = neutral_loss
        if neutral_loss is not None:
            self.mass += neutral_loss.mass

    def get_fragment_name(self):
        parts = [self._name]
        neutral_loss = self.neutral_loss
        if neutral_loss is not None:
            parts.append(str(neutral_loss))
        return ''.join(parts)

    def set_name(self, name):
        self._name = name

    name = property(get_fragment_name, set_name)
    neutral_loss = property(get_neutral_loss, set_neutral_loss)

    def generate_neutral_losses(self, losses=None):
        if losses is None:
            losses = generic_neutral_losses_composition.keys()

        for loss in losses:
            frag = self.clone()
            frag.neutral_loss = NeutralLoss(loss)
            yield frag


class PeptideFragment(FragmentBase):
    """Glycopeptide Fragment"""

    # parser = re.compile(r"(?P<kind>[abcxyz])(?P<position>[0-9]+)(?P<modificaiton>\+.*)?")
    concerned_mods = ['HexNAc']

    def __init__(self, frag_type, position, modification_dict, mass, golden_pairs=None,
                 flanking_amino_acids=None, glycosylation=None, neutral_loss=None):
        if golden_pairs is None:
            golden_pairs = []
        self.type = frag_type
        # The mass value is the bare backbone's mass
        self.bare_mass = mass

        self.mass = mass

        self.flanking_amino_acids = flanking_amino_acids

        self.position = position
        self.modification_dict = modification_dict

        for key, value in self.modification_dict.items():
            self.mass += Modification(key).mass * value

        self.golden_pairs = golden_pairs
        self.glycosylation = glycosylation
        self.neutral_loss = neutral_loss

    def clone(self):
        corrected_mass = self.mass
        if self._neutral_loss is not None:
            corrected_mass - self.neutral_loss.mss
        return self.__class__(
            self.type, self.position, dict(self.modification_dict),
            corrected_mass, self.golden_pairs, tuple(self.flanking_amino_acids),
            self.glycosylation.clone() if self.glycosylation is not None else None,
            self._neutral_loss.clone() if self._neutral_loss is not None else None)

    def base_name(self):
        """Simply return string like b2, y3 with no modificaiton information."""
        fragment_name = []
        fragment_name.append(self.type)
        fragment_name.append(str(self.position))
        return ''.join(fragment_name)

    def get_fragment_name(self):
        """Connect the information into string."""
        fragment_name = []
        fragment_name.append(str(self.type))
        fragment_name.append(str(self.position))

        # Only concerned modifications are reported.
        for mod_name in self.concerned_mods:
            if mod_name in self.modification_dict:
                if self.modification_dict[mod_name] > 1:
                    fragment_name.extend(['+', str(self.modification_dict[mod_name]), mod_name])
                elif self.modification_dict[mod_name] == 1:
                    fragment_name.extend(['+', mod_name])
                else:
                    pass

        if self.neutral_loss is not None:
            fragment_name.append(str(self.neutral_loss))

        return ''.join(fragment_name)

    def partial_loss(self, modifications=None):
        if modifications is None:
            modifications = self.concerned_mods
        mods = dict(self.modification_dict)
        mods_of_interest = defaultdict(int, {k: v for k, v in mods.items() if k in modifications})
        mods_of_interest["HexNAc"] *= 2  # Allow partial destruction of N-glycan core

        other_mods = {k: v for k, v in mods.items() if k not in modifications}
        for mod in descending_combination_counter(mods_of_interest):
            other_mods.update({k: v for k, v in mod.items() if v != 0})
            yield PeptideFragment(
                self.type, self.position, dict(other_mods), self.bare_mass,
                golden_pairs=self.golden_pairs, flanking_amino_acids=self.flanking_amino_acids)

    def to_json(self):
        d = dict(self.__dict__)
        d['key'] = self.name

    @property
    def name(self):
        return self.get_fragment_name()

    def __repr__(self):
        return ("PeptideFragment(%(type)s %(position)s %(mass)s "
                "%(modification_dict)s %(flanking_amino_acids)s %(neutral_loss)r)") % {
            "type": self.type, "position": self.position, "mass": self.mass,
            "modification_dict": self.modification_dict, "flanking_amino_acids": self.flanking_amino_acids,
            "neutral_loss": self.neutral_loss
        }


class SimpleFragment(FragmentBase):
    def __init__(self, name, mass, kind, neutral_loss=None):
        self.name = name
        self.mass = mass
        self.kind = kind
        self.neutral_loss = neutral_loss

    def __repr__(self):
        return "SimpleFragment(name={self.name}, mass={self.mass:.04f}, kind={self.kind})".format(self=self)


class MemoizedIonSeriesMetaclass(type):
    def __call__(self, name=None, *args, **kwargs):
        if not hasattr(self, "_cache"):
            self._cache = dict()
        try:
            if name is not None:
                return self._cache[name]
            else:
                raise Exception("Must provide a name parameter")
        except KeyError:
            if name is not None:
                inst = type.__call__(self, name=name, *args, **kwargs)
                self._cache[inst.name] = inst
                return inst
            else:
                raise KeyError("Cannot find an IonSeries for %r" % (name))


class IonSeries(object):
    __metaclass__ = MemoizedIonSeriesMetaclass

    @classmethod
    def get(cls, name):
        return cls(name)

    def __init__(self, name, direction=None, includes_peptide=True, mass_shift=None, regex=None):
        if direction is None:
            if name in fragment_direction:
                direction = fragment_direction[name]
            else:
                direction = 0
        if mass_shift is None:
            if name in fragment_shift:
                mass_shift = fragment_shift[name]
            else:
                mass_shift = 0.
        self.name = name
        self.direction = direction
        self.includes_peptide = includes_peptide
        self.mass_shift = mass_shift
        self.regex = re.compile(regex) if regex is not None else regex

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        try:
            return self.name == other.name
        except AttributeError:
            return self.name == other

    def __ne__(self, other):
        return not self == other

    __repr__ = simple_repr

    def __str__(self):
        return str(self.name)

    def is_member(self, key):
        if self.regex is None:
            return key.startswith(self.name)
        else:
            return self.regex.search(key)

    __call__ = is_member

    def __radd__(self, other):
        return other + self.mass_shift

    def __rsub__(self, other):
        return other - self.mass_shift

    def __neg__(self):
        return -self.mass_shift

    def __pos__(self):
        return self.mass_shift

    def __add__(self, other):
        return self.mass_shift + other

    def __sub__(self, other):
        return self.mass_shift - other


IonSeries.b = IonSeries("b")
IonSeries.y = IonSeries("y")
IonSeries.oxonium_ion = IonSeries("oxonium_ion", includes_peptide=False)
IonSeries.stub_glycopeptide = IonSeries("stub_glycopeptide")
