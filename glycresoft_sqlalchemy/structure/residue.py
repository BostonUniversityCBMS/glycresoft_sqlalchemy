from . import ResidueBase
from .composition import Composition, composition_to_mass
from ..utils.memoize import memoize

symbol_to_residue = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val'
}


residue_to_symbol = {value: key for key, value in symbol_to_residue.items()}


residue_table = {
    'Ala': 'C3H5NO',
    'Arg': 'C6H12N4O1',
    'Asn': 'C4H6N2O2',
    'Asp': 'C4H5N1O3',
    'Cys': 'C3H5N1O1S1',
    'Glu': 'C5H7NO3',
    'Gln': 'C5H8N2O2',
    'Gly': 'C2H3N1O1',
    'His': 'C6H7N3O1',
    'Ile': 'C6H11N1O1',
    'Leu': 'C6H11N1O1',
    'Lys': 'C6H12N2O1',
    'Met': 'C5H9N1O1S1',
    'Phe': 'C9H9N1O1',
    'Pro': 'C5H7N1O1',
    'Ser': 'C3H5N1O2',
    'Thr': 'C4H7N1O2',
    'Trp': 'C11H10N2O1',
    'Tyr': 'C9H9N1O2',
    'Val': 'C5H9N1O1',
}

residue_chemical_property_group = {
    'Ala': 'hydrophobic',
    'Arg': 'positive_side_chain',
    'Asn': 'polar_uncharged',
    'Asp': 'negative_side_chain',
    'Cys': 'special_case',
    'Glu': 'negative_side_chain',
    'Gln': 'polar_uncharged',
    'Gly': 'special_case',
    'His': 'positive_side_chain',
    'Ile': 'hydrophobic',
    'Leu': 'hydrophobic',
    'Lys': 'positive_side_chain',
    'Met': 'hydrophobic',
    'Phe': 'hydrophobic',
    'Pro': 'special_case',
    'Ser': 'polar_uncharged',
    'Thr': 'polar_uncharged',
    'Trp': 'hydrophobic',
    'Tyr': 'hydrophobic',
    'Val': 'hydrophobic',
}


degeneracy_index = {
}


class MemoizedResidueMetaclass(type):
    '''
    A metaclass that memoizes Residues as they are constructed
    by overriding the class __call__ method. It will attempt to
    look up previously created residues by symbol, then by name.
    If a previous instance is not found, it will be created and
    saved.

    Attributes
    ----------
    _cache: dict

    '''
    def __call__(self, symbol=None, name=None, *args, **kwargs):
        if not hasattr(self, "_cache"):
            self._cache = dict()
        try:
            if symbol is not None:
                return self._cache[symbol]
            elif name is not None:
                return self._cache[name]
            else:
                raise Exception("Must provide a symbol or name parameter")
        except KeyError:
            if symbol is not None:
                inst = type.__call__(self, symbol=symbol, *args, **kwargs)
                self._cache[inst.symbol] = inst
                self._cache[inst.name] = inst
                return inst

            elif name is not None:
                inst = type.__call__(self, name=name, *args, **kwargs)
                self._cache[inst.symbol] = inst
                self._cache[inst.name] = inst
                return inst
            else:
                raise KeyError("Cannot find a Residue for %r" % ((symbol, name),))


class Residue(ResidueBase):
    '''
    Represent a single Amino Acid residue which compose peptide sequences. The
    structure itself is intended to be immutable.

    Attributes
    ----------
    name: str
    symbol: str
    mass: float
    composition: :class:`glypy.Composition`
    chemical_class: str
    is_degenerate: bool

    '''
    __slots__ = ["name", "symbol", "mass", "composition"]

    __metaclass__ = MemoizedResidueMetaclass

    @staticmethod
    @memoize()
    def mass_by_name(sym):
        name = symbol_to_residue.get(sym, sym)
        formula = residue_table.get(name)
        return composition_to_mass(formula)

    def __init__(self, symbol=None, name=None):
        self.symbol = symbol
        self.name = name
        self.mass = 0.0
        if symbol is not None:
            self.by_symbol(symbol)
        elif name is not None:
            self.by_name(name)

    def by_name(self, name):
        """Configure this instance by name information

        Parameters
        ----------
        name : str
            Amino Acid Name
        """
        self.composition = Composition(residue_table[name])
        self.name = name
        self.mass = self.composition.mass
        self.symbol = residue_to_symbol[name]

    def by_symbol(self, symbol):
        """Configure this instance by symbol information,
        by going from symbol to name, and from name to data

        Parameters
        ----------
        symbol : str
            Amino Acid symbol
        """
        try:
            name = symbol_to_residue[symbol]
            self.by_name(name)
        except KeyError:
            self.by_name(symbol)

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if self is other:
            return True
        if isinstance(other, Residue):
            other = other.name
        return self.name == other

    def __ne__(self, other):
        if self is other:
            return False
        if isinstance(other, Residue):
            other = other.name
        return self.name != other

    def __getstate__(self):
        return [self.name, self.symbol, self.mass]

    def __setstate__(self, state):
        self.name, self.symbol, self.mass = state

    @property
    def chemical_class(self):
        return residue_chemical_property_group[self.name]

    @property
    def is_degenerate(self):
        try:
            return degeneracy_index[self.name]
        except KeyError:
            return False


def register_residue(name, symbol, formula, chemical_class):
    assert symbol not in symbol_to_residue
    assert name not in residue_table
    residue_to_symbol[name] = symbol
    symbol_to_residue[symbol] = name
    residue_table[name] = formula
    residue_chemical_property_group[name] = chemical_class
    return Residue(symbol=symbol)


def register_degenerate(name, symbol, mappings):
    assert symbol not in symbol_to_residue
    assert name not in residue_table
    residue_to_symbol[name] = symbol
    symbol_to_residue[symbol] = name
    residue_table[name] = residue_table[mappings[0]]
    residue_chemical_property_group[name] = residue_chemical_property_group[mappings[0]]
    degeneracy_index[name] = frozenset(mappings)
    return Residue(symbol=symbol)


register_degenerate("Xle", "J", ["Leu", "Ile"])


def get_all_residues():
    symbols = [
        'A',
        'R',
        'N',
        'D',
        'C',
        'E',
        'Q',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V',
    ]
    return map(Residue, symbols)
