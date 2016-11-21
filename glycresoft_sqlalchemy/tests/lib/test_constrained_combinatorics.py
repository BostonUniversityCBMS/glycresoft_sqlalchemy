import unittest
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics


class TestConstrainedCombinatorics(unittest.TestCase):
    def test_constrained(self):
        rules_table = {
            "Hex": (3, 10),
            "HexNAc": (2, 10),
            "Fuc": (0, 5),
            "NeuAc": (0, 4)
        }

        constraints_list = [
            ["Fuc", "<", "HexNAc"],
            ["NeuAc", "<", "HexNAc - 1"]
        ]

        generator = constrained_combinatorics.CombinatoricCompositionGenerator(
            rules_table=rules_table, constraints=constraints_list)
        compositions = list(generator)
        for composit, classifications in compositions:
            self.assertTrue(composit['Fuc'] < composit['HexNAc'])
            self.assertTrue(composit['NeuAc'] < (composit['HexNAc'] - 1))


if __name__ == '__main__':
    unittest.main()
