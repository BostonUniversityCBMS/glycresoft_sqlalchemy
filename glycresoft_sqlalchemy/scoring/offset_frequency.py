import itertools

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatch, PipelineModule
from glycresoft_sqlalchemy.utils import collectiontools

from .simple_scoring_algorithm import split_ion_list

Sequence = sequence.sequence


def offset_frequency(gsms, kind='b'):
    total_sites = 0
    total_explained = 0
    total_unexplained = 0
    total_sparsity = 0
    for gsm in gsms:
        n_frag_sites = count_fragmentation_sites(gsm.glycopeptide_match.glycopeptide_sequence, kind)




def count_fragmentation_sites(sequence, kind='b'):
    sequence = Sequence(sequence)
    fragmentation_sites = len(collectiontools.flatten(sequence.get_fragments("b")))
    return fragmentation_sites
