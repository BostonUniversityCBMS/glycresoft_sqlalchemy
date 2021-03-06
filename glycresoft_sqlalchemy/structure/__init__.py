from constants import constants
from composition import Composition


__all__ = [
    "sequence",
    "modification",
    "composition",
    "fragment",
    "glycan",
    "residue",
    "sequence_space"
]


class MoleculeBase(object):
    mass = None

    def __copy__(self):
        return self.clone()


# A few base types for doing type-based behavior changes
class PeptideSequenceBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequences, with or without modifiations
    '''


class ModificationBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequence modifications
    '''


class ResidueBase(MoleculeBase):
    '''
    A base type for classes describing amino acid residues
    '''
