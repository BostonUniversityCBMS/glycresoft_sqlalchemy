import unittest
from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms2 import make_decoys

sequences = [
    "QD(Dehydrated)QC(Carbamidomethyl)IYN(HexNAc)TTYLNVQR{Fuc:3; Hex:7; HexNAc:5; Neu5Ac:1}",
    "QDQC(Carbamidomethyl)IYN(HexNAc)TTYLNVQR{Hex:6; HexNAc:5; Neu5Ac:3}",
    "QDQC(Carbamidomethyl)IYN(HexNAc)TTYLNVQR{Fuc:5; Hex:4; HexNAc:3; Neu5Ac:3}",
    "(Gln->pyro-Glu)-QDQC(Carbamidomethyl)IYN(HexNAc)TTYLNVQR{Fuc:3; Hex:7; HexNAc:5; Neu5Ac:1}"
]


class TestDecoyProducer(unittest.TestCase):
    def test_make_decoy(self):
        for seq in sequences:
            target = sequence.parse(seq)
            decoy = make_decoys.reverse_preserve_sequon(seq)
            self.assertAlmostEqual(target.mass, decoy.mass, 4)

            target_fragments = make_decoys.fragments(target)
            decoy_fragments = make_decoys.fragments(decoy)

            target_stubs = target_fragments[-1]
            decoy_stubs = decoy_fragments[-1]

            for t, d in zip(target_stubs, decoy_stubs):
                self.assertAlmostEqual(t['mass'], d['mass'], 4)


if __name__ == '__main__':
    unittest.main()
