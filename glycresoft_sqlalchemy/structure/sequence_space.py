import copy
import itertools

from .modification import Modification, ModificationTable
from .sequence import Sequence


class SequenceSpace:
    """
    Generate all possible modified peptide sequences for a given
    base sequence and a collection of modifications and specific
    glycosylation_sites.

    Attributes
    ----------
    candidate_sites : list of int
        List of indices that may be glycosylated
    modifications : list of :class:`Modification`
        The modifications to place on the sequence
    sequence : str
        The base peptide sequence to build on
    """
    def __init__(self, sequence, glycan_sites, mod_list):
        self.sequence = Sequence(sequence)
        self.candidate_sites = glycan_sites
        self.modifications = mod_list

    def __repr__(self):
        rep = "SequenceSpace({sequence} {candidate_sites} {modifications})".format(**self.__dict__)
        return rep

    def get_modification_sites(self):
        modification_index_bound = []
        for mod in self.modifications:
            if mod.position != -1:
                # The modification position is pre-specified
                modification_index_bound.append([mod.position])
            else:
                valid_sites = mod.find_valid_sites(self.sequence)
                all_combinations = [position for position in itertools.combinations(valid_sites, mod.number)]
                if(len(all_combinations) == 0):
                    raise NoSitesFoundException("No Valid Sites Found: %r %r in\n %r" % (
                        mod, valid_sites, self.sequence))
                modification_index_bound.append(all_combinations)

        return modification_index_bound

    def compose_sequence(self, modification_sites, required_glycosylation_sites):
        seq_space = []
        num_modifications = len(self.modifications)
        combination_index_sites = itertools.product(*map(lambda x: range(0, len(x)), modification_sites))
        for modification_indices in combination_index_sites:
            mod_sites = [modification_sites[mod_i][pos_j] for mod_i, pos_j in enumerate(modification_indices)]

            common_sites = set(itertools.chain.from_iterable(mod_sites))
            glycosylation_sites = set(self.candidate_sites).difference(common_sites)

            if(len(common_sites) != sum(map(len, mod_sites)) or
               (required_glycosylation_sites > len(glycosylation_sites))):
                # Invalid Configuration, can't place all Glycans
                continue

            raw_sequence = (self.sequence).clone()

            for index in range(num_modifications):
                for mod_site in mod_sites[index]:
                    raw_sequence.add_modification(mod_site, self.modifications[index].rule)

            for sites in itertools.combinations(glycosylation_sites, required_glycosylation_sites):
                temp_seq = (raw_sequence).clone()
                for site in sites:
                    gly_mod = Modification(ModificationTable.other_modifications["HexNAc"], site, 1)
                    temp_seq.add_modification(modification_type=gly_mod)
                seq_space.append(temp_seq)

        return seq_space

    def get_theoretical_sequence(self, num_sites):
        modification_index_bound = self.get_modification_sites()
        sequence_space = self.compose_sequence(modification_index_bound, num_sites)
        return sequence_space


class NoSitesFoundException(Exception):
    pass


class UnqualifiedModifierException(Exception):
    pass
