import unittest

from glycresoft_sqlalchemy.structure import sequence, modification, residue
from glypy import GlycanComposition, Glycan


R = residue.Residue


p1 = "PEPTIDE"
p2 = "YPVLN(HexNAc)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"


class TestPeptideSequence(unittest.TestCase):
    def test_parser(self):
        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p1)
        self.assertEqual(len(mods), 0)
        self.assertEqual(len(chunks), len(p1))
        self.assertEqual(glycan, "")

        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p2)
        self.assertEqual(GlycanComposition.parse("{Hex:9; HexNAc:2}"), glycan)
        self.assertEqual(len(mods), 2)
        self.assertEqual(len(chunks), 16)

    # def test_mass(self):
    #     case = sequence.PeptideSequence(p1)
    #     self.assertAlmostEqual(case.mass, )



if __name__ == '__main__':
    unittest.main()
