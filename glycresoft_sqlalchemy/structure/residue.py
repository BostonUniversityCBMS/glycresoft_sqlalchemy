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
    'V': 'Val'}


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
    'HexNAc': 'C8H13NO5',
    'Hex': 'C6H10O5',
    'dHex': 'C6H10O4',
    'NeuAc': 'C11H17NO8',
    'NeuGc': 'C11H17NO9',
    'Water': 'H2O'}


class Residue(ResidueBase):
    __slots__ = ["name", "symbol", "mass"]

    @staticmethod
    @memoize()
    def mass_by_name(sym):
        name = symbol_to_residue.get(sym, sym)
        formula = residue_table.get(name)
        return composition_to_mass(formula)

    """Basic mass values for peptide sequences"""
    def __init__(self, symbol=None, name=None):
        self.symbol = symbol
        self.name = name
        self.mass = 0.0
        if symbol is not None:
            self.by_symbol(symbol)
        elif name is not None:
            self.by_name(name)

    def by_name(self, name):
        self.compo = residue_table[name]
        self.name = name
        self.mass = composition_to_mass(self.compo)
        self.symbol = residue_to_symbol[name]

    def by_symbol(self, symbol):
        try:
            name = symbol_to_residue[symbol]
            self.by_name(name)
        except KeyError:
            self.by_name(symbol)

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, Residue):
            other = other.name
        return self.name == other

    def __ne__(self, other):
        if isinstance(other, Residue):
            other = other.name
        return self.name != other

    def __getstate__(self):
        return [self.name, self.symbol, self.mass]

    def __setstate__(self, state):
        self.name, self.symbol, self.mass = state
