__all__ = [
    "sequence",
    "modification",
    "composition",
    "fragment",
    "glycans",
    "residue",
    "stub_glycopeptides",
    "sequence_space"
]

from constants import constants


class MoleculeBase(object):
    mass = None


# A few base types for doing type-based behavior changes
class PeptideSequenceBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequences, with or without modifiations
    '''
    def get_sequence(self, start=0):
        return self.getSequence(start)

    def getSequence(self, start=0):
        pass


class ModificationBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequence modifications
    '''
    pass


class ResidueBase(MoleculeBase):
    '''
    A base type for classes describing amino acid residues
    '''
