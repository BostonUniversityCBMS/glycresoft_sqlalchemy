import unittest

from glycresoft_sqlalchemy.structure import sequence, modification, residue
from glypy import GlycanComposition, Glycan, MonosaccharideResidue


R = residue.Residue


p1 = "PEPTIDE"
p2 = "YPVLN(HexNAc)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"
p3 = "NEEYN(HexNAc)K{Hex:5; HexNAc:4; NeuAc:2}"
hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
hexose_mass = MonosaccharideResidue.from_iupac_lite("Hex").mass()


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

    def test_stub_ions(self):
        peptide = sequence.parse(p3)
        stubs = sorted({f.mass for f in peptide.stub_fragments()})
        self.assertAlmostEqual(stubs[0], 795.3399, 3)
        self.assertAlmostEqual(stubs[1], 998.4193, 3)
        self.assertAlmostEqual(stubs[2], 1201.4986, 3)
        self.assertAlmostEqual(stubs[3], 1363.5515, 3)
        self.assertAlmostEqual(stubs[4], 1525.6043, 3)
        self.assertAlmostEqual(stubs[5], 1687.6571, 3)

    def test_glycan_fragments_stubs(self):
        peptide = sequence.parse(p3)
        stubs = {f.name: f.mass for f in peptide.glycan_fragments(all_series=True, allow_ambiguous=True)}
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:2}"], 1687.6571, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:3}"],
                               1687.6571 + hexnac_mass, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:4; HexNAc:3}"],
                               1687.6571 + hexnac_mass + hexose_mass, 3)

    # def test_mass(self):
    #     case = sequence.PeptideSequence(p1)
    #     self.assertAlmostEqual(case.mass, )


if __name__ == '__main__':
    unittest.main()
