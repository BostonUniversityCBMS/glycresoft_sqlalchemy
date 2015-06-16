import re
import itertools
import operator

from ..model import GlycopeptideMatch, TheoreticalGlycopeptide, DatabaseManager
from glycresoft_ms2_classification.utils import collectiontools

from glycresoft_ms2_classification.structure.sequence import Sequence
from glycresoft_ms2_classification.structure.stub_glycopeptides import StubGlycopeptide


imap = itertools.imap
chain = itertools.chain

key_getter = operator.itemgetter("key")

bare_ion_pattern = re.compile(r"[czby](\d+)")
glycosylated_ion_pattern = re.compile(r"[czby](\d+)\+(.+)")


def apply(structure):



def stream_backbone_coordinate(iterable):
    return imap(int, chain.from_iterable(imap(bare_ion_pattern.findall, imap(key_getter, iterable))))


def incvec_spread(vector, p, s=1):
    for i in range(s - 1, p, 1 if p > s else -1):
        vector[i] += 1
    return vector


def incvec_spot(vector, p, s=1):
    vector[p] = True
    return vector


def bitwise_or(a, b):
    return list(imap(max, zip(a, b)))


def mean_coverage(matched):
    sequence = Sequence(matched.glycopeptide_sequence)
    sequence_length = len(sequence)

    b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.bare_b_ions):
        incvec_spot(b_ions, b_ion)

    glycosylated_b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.glycosylated_b_ions):
        incvec_spot(glycosylated_b_ions, b_ion)

    y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.bare_y_ions):
        incvec_spot(y_ions, sequence_length - y_ion, sequence_length)

    glycosylated_y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.glycosylated_y_ions):
        incvec_spot(glycosylated_y_ions, sequence_length - y_ion, sequence_length)

    coverage = sum(bitwise_or(b_ions, glycosylated_b_ions)) + sum(bitwise_or(y_ions, glycosylated_y_ions))
    coverage /= 2. * sequence_length

    return coverage


def mean_hexnac_coverage(matched, theoretical):
    b_observed = set(stream_backbone_coordinate(matched.glycosylated_b_ions))
    y_observed = set(stream_backbone_coordinate(matched.glycosylated_y_ions))

    b_enumerated = set(stream_backbone_coordinate(theoretical.glycosylated_b_ions))
    y_enumerated = set(stream_backbone_coordinate(theoretical.glycosylated_y_ions))
    return (len(b_observed) + len(y_observed)) / float(len(b_enumerated) + len(y_enumerated))


def observed_vs_enumerated_stub(matched):
    obs = len({stub['key'].split(':')[0] for stub in matched.stub_ions})
    expected = len(StubGlycopeptide.from_sequence(
        Sequence(matched.glycopeptide_sequence)).get_stubs())
    return obs/float(expected)


def bad_oxonium(matched, theoretical):
    pass
