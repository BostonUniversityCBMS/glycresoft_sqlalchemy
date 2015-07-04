import re
import itertools
import operator

from ..data_model import GlycopeptideMatch, TheoreticalGlycopeptide, DatabaseManager
from glycresoft_ms2_classification.utils import collectiontools

from glycresoft_ms2_classification.structure.sequence import Sequence
from glycresoft_ms2_classification.structure.stub_glycopeptides import StubGlycopeptide


imap = itertools.imap
chain = itertools.chain

key_getter = operator.itemgetter("key")

bare_ion_pattern = re.compile(r"[czby](\d+)")
glycosylated_ion_pattern = re.compile(r"[czby](\d+)\+(.+)")


def apply(matched, theoretical, **parameters):
    matched.mean_coverage = mean_coverage(matched)
    matched.mean_hexnac_coverage = mean_hexnac_coverage(matched, theoretical)
    calculate_score(matched, **parameters)


def calculate_score(matched, backbone_weight=0.5, hexnac_weight=0.5):
    matched.ms2_score = (matched.mean_coverage * backbone_weight) +\
     (matched.mean_hexnac_coverage * hexnac_weight)
    return matched


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


def mean(vector):
    return sum(vector) / float(len(vector))


def vsum(a, b):
    return [x + y for x, y in zip(a, b)]


def compute_backbone_vectors(matched):
    sequence = Sequence(matched.glycopeptide_sequence)
    sequence_length = len(sequence)

    b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.bare_b_ions):
        incvec_spread(b_ions, b_ion)

    glycosylated_b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.glycosylated_b_ions):
        incvec_spread(glycosylated_b_ions, b_ion)

    y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.bare_y_ions):
        incvec_spread(y_ions, sequence_length - y_ion, sequence_length)

    glycosylated_y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.glycosylated_y_ions):
        incvec_spread(glycosylated_y_ions, sequence_length - y_ion, sequence_length)

    return b_ions, glycosylated_b_ions, y_ions, glycosylated_y_ions


def mean_coverage(matched):
    sequence = Sequence(matched.glycopeptide_sequence)
    sequence_length = len(sequence)

    b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.bare_b_ions):
        incvec_spread(b_ions, b_ion)

    glycosylated_b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.glycosylated_b_ions):
        incvec_spread(glycosylated_b_ions, b_ion)

    y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.bare_y_ions):
        incvec_spread(y_ions, sequence_length - y_ion, sequence_length)

    glycosylated_y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.glycosylated_y_ions):
        incvec_spread(glycosylated_y_ions, sequence_length - y_ion, sequence_length)

    coverage = mean(vsum(bitwise_or(b_ions, glycosylated_b_ions), bitwise_or(y_ions, glycosylated_y_ions)))
    coverage /= 1. * sequence_length

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


def compute_percent_uncovered(matched):
    bs, gbs, ys, gys = compute_backbone_vectors(matched)
    backbone = vsum(bitwise_or(bs, gbs), bitwise_or(ys, gys))
    is_zero = 0
    for i in backbone:
        if i == 0:
            is_zero += 1
    return is_zero / float(len(matched))
