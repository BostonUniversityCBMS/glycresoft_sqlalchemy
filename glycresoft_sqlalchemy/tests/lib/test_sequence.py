import unittest

from glycresoft_sqlalchemy.structure import sequence, modification, residue
from glypy import GlycanComposition, Glycan


p1 = "PEPTIDE"
p2 = "YPVLN(HexNAc)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"


class TestPeptideSequence(unittest.TestCase):
    def test_parser(self):
        pass
