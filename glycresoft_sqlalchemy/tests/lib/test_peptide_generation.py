import unittest

from glycresoft_sqlalchemy.data_model import NaivePeptide, Protein
from glycresoft_sqlalchemy.structure import sequence, modification, residue
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder import peptide_utilities


constant_modifications = ["Carbamidomethyl (C)"]
variable_modifications = ["Deamidation (N)", "Pyro-glu from Q (Q@N-term)"]

peptide_sequence = "SVQEIQATFFYFTPNK"
peptide_sequence2 = "QDQCIYNTTYLNVQR"

protein_sequence = "MALSWVLTVLSLLPLLEAQIPLCANLVPVPITNATLDQITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQDQCIYNTTYLNVQRENGTISRYVGGQEHFAHLLILRDTKTYMLAFDVNDEKNWGLSVYADKPETTKEQLGEFYEALDCLRIPKSDVVYTDWKKDKCEPLEKQHEKERKQEEGES"

protein = Protein(
    protein_sequence=protein_sequence,
    glycosylation_sites=sequence.parse(protein_sequence).n_glycan_sequon_sites)

peptide_obj = NaivePeptide(
    base_peptide_sequence=peptide_sequence,
    start_position=78-21,
    end_position=94-21)

peptide_obj2 = NaivePeptide(
    base_peptide_sequence=peptide_sequence2,
    start_position=107-21,
    end_position=122-21)


protein.naive_peptides.append(peptide_obj)
protein.naive_peptides.append(peptide_obj2)


class TestPeptidoformBuilder(unittest.TestCase):
    def test_unpositioned_isoforms(self):
        mt = modification.RestrictedModificationTable(
            None, constant_modifications=constant_modifications,
            variable_modifications=variable_modifications)
        solutions = list(peptide_utilities.unpositioned_isoforms(
            peptide_obj, constant_modifications, variable_modifications, mt))

        solution_1 = [('SVQEIQATFFYFTPNK', {}, 1918.94651437366, [14])]
        solution_2 = [('QDQC(Carbamidomethyl)IYNTTYLNVQR', {}, 1914.88941078313, [6]),
                      ('(Gln->pyro-Glu)-QDQC(Carbamidomethyl)IYNTTYLNVQR',
                       {},
                       1896.85503675106,
                       [6]),
                      ('QDQC(Carbamidomethyl)IYNTTYLNVQR',
                       {'Deamidation': 1},
                       1915.87342678313,
                       [6]),
                      ('(Gln->pyro-Glu)-QDQC(Carbamidomethyl)IYNTTYLNVQR',
                       {'Deamidation': 1},
                       1897.83905275106,
                       [6])]

        self.assertEqual(solutions, solution_1)

        solutions = list(peptide_utilities.unpositioned_isoforms(
            peptide_obj2, constant_modifications, variable_modifications, mt))

        self.assertEqual(solutions, solution_2)


if __name__ == '__main__':
    unittest.main()
