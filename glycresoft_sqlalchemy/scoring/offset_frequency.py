from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatch, PipelineModule
from glycresoft_sqlalchemy.utils import collectiontools

Sequence = sequence.sequence


def offset_frequency(gsms, kind='b'):
    pass


def fragmentation_sites(sequence, kind='b'):
    sequence = Sequence(sequence)
    fragmentation_sites = collectiontools.flatten(sequence.get_fragments("b"))
    return fragmentation_sites
