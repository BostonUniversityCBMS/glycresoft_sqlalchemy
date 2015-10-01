from glypy import MonosaccharideResidue
from .residue import Residue as AminoAcidResidue
from .composition import Composition
from .modification import Modification
from glycresoft_sqlalchemy.utils.collectiontools import SqliteSet


class AminoAcidSequenceBuildingBlock(object):
    def __init__(self, residue_, modifications, neutral_mass=None):
        self.residue = residue_
        self.modifications = tuple(modifications)
        if neutral_mass is None:
            neutral_mass = residue_.mass + sum(m.mass for m in modifications)
        self.neutral_mass = neutral_mass

    def __iter__(self):
        yield self.residue
        yield self.modifications

    def __hash__(self):
        return hash((self.residue, self.modifications))

    def __eq__(self, other):
        return self.residue == other.residue and self.modifications == other.modifications

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "({}, {}):{:.2f}".format(self.residue.symbol, self.modifications, self.neutral_mass)

    def __str__(self):
        return "{}{}".format(
            self.residue.symbol,
            "({.name})".format(self.modifications[0]) if len(self.modifications) > 0 else "")


class SequenceSegmentBlock(object):
    def __init__(self, sequence, neutral_mass=None):
        self.sequence = sequence
        if neutral_mass is None:
            neutral_mass = 0.
            for residue, mods in sequence:
                neutral_mass += AminoAcidResidue.mass_by_name(residue)
                for mod in mods:
                    neutral_mass += Modification.mass_by_name(mod)

        self.neutral_mass = neutral_mass

    def __iter__(self):
        return iter(self.sequence)

    def __hash__(self):
        return hash(self.sequence)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __ne__(self, other):
        return not self == other


class MonosaccharideResidueAdapter(object):
    def __init__(self, residue):
        if isinstance(residue, basestring):
            residue = MonosaccharideResidue.from_iupac_lite(residue)
        self.residue = residue
        self.mass = residue.mass()
        self.symbol = "(%s)" % str(residue)

    def __hash__(self):
        return hash(self.residue)

    def __eq__(self, other):
        if isinstance(other, (MonosaccharideResidue, AminoAcidResidue)):
            return self.symbol == other.symbol
        else:
            return self.symbol == other

    def __ne__(self, other):
        return not self == other


class SequenceComposition(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        self._mass = None
        self._composition_offset = Composition("H2O")
        self.extend(*args)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self._mass = None

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    @property
    def mass(self):
        if self._mass is not None:
            return self._mass
        mass = self._composition_offset.mass
        for residue_type, count in self.items():
            mass += residue_type.neutral_mass * count
        self._mass = mass
        return mass

    def extend(self, *args):
        for residue in args:
            self[residue] += 1

    def __iadd__(self, other):
        for elem, cnt in (other.items()):
            self[elem] += cnt
        return self

    def __add__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        return self

    def __sub__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError(
                'Cannot multiply Composition by non-integer',
                other)
        prod = SequenceComposition()
        for k, v in self.items():
            prod[k] = v * other

        return prod

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not isinstance(other, dict):
            return False
        self_items = set([i for i in self.items() if i[1]])
        other_items = set([i for i in other.items() if i[1]])
        return self_items == other_items

    def __neg__(self):
        return -1 * self

    def __missing__(self, key):
        return 0

    @property
    def composition_offset(self):
        return self._composition_offset

    @composition_offset.setter
    def composition_offset(self, value):
        self._mass = None
        self._composition_offset = value

    def clone(self):
        return self.__class__(*list(self.flatten()))

    def flatten(self):
        for k in self:
            for i in range(0, self[k]):
                yield k

    def serialize(self):
        return "{%s}" % '; '.join("{}:{}".format(str(k), v) for k, v in sorted(
            self.items(), key=lambda x: x[0].neutral_mass))

    __str__ = serialize

    def __hash__(self):
        return hash(str(self))


def all_compositions(blocks, count=20):
    def extend_segment(base, blocks):
        for block in blocks:
            ext = base.clone()
            ext[block] += 1
            ext._mass = base.mass + block.neutral_mass
            yield ext

    candidates = [SequenceComposition()]
    next_round = SqliteSet()
    for i in range(count):
        for candidate in candidates:
            extents = (extend_segment(candidate, blocks))
            for case in extents:
                next_round.add(case)
        candidates = next_round
        next_round = SqliteSet()
    return set(candidates)
