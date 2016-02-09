import re
import itertools
import operator
import math

from ..utils import collectiontools

from ..structure.sequence import Sequence

from .base import ScorerBase, GlycopeptideSpectrumMatchScorer

imap = itertools.imap
chain = itertools.chain

key_getter = operator.itemgetter("key")

_normalizer = math.log(3, 2)
_expected_stub_ion_count = 3

bare_ion_pattern = re.compile(r"[czby](\d+)")
glycosylated_ion_pattern = re.compile(r"[czby](\d+)\+(.+)")


def split_ion_list(ion_list):
    """
    Split a list of ions into classed sub-lists
    Parameters
    ----------
    ion_list : list
        List of ions of mixed series/classes

    Returns
    -------
    dict
        Maps ion series names to ion series members from `ion_list`
    """
    ion_store = {
        "oxonium_ions": [],
        "stub_ions": [],
        "bare_b_ions": [],
        "bare_y_ions": [],
        "glycosylated_b_ions": [],
        "glycosylated_y_ions": [],
        "bare_c_ions": [],
        "bare_z_ions": [],
        "glycosylated_c_ions": [],
        "glycosylated_z_ions": []
    }
    for ion in ion_list:
        key = ion['key']
        # Pseudoswitch
        if 'pep' == key[:3]:
            ion_store['stub_ions'].append(ion)
            continue
        match = glycosylated_ion_pattern.match(key)
        if match:
            ion_store['glycosylated_%s_ions' % key[0]].append(ion)
            continue

        match = bare_ion_pattern.match(key)
        if match:
            ion_store['bare_%s_ions' % key[0]].append(ion)
            continue

        # Finally
        ion_store['oxonium_ions'].append(ion)
    return ion_store


def merge_ion_matches(matches):
    """Group multiple matches to the same fragment

    Parameters
    ----------
    matches : list of dict
        The list of ion matches to group

    Returns
    -------
    list of dict: Merged ion matches
    """
    groups = collectiontools.groupby(matches,
                                     key_getter)
    best_matches = []
    fabs = math.fabs
    for key, matched_key in groups.items():
        best_match = matched_key[0]
        best_ppm = fabs(best_match["ppm_error"])
        peak_map = {best_match["peak_id"]: best_match}
        for match in matched_key[1:]:
            peak_map[match["peak_id"]] = match
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
                best_ppm = fabs(match["ppm_error"])
        best_match = best_match.copy()
        best_match["peak_map"] = peak_map
        best_matches.append(best_match)
    return best_matches


def evaluate(matched, theoretical, **parameters):
    """
    Compute the score statistics of `matched`, the `mean_coverage` and
    `mean_hexnac_coverage`, and the `ms2_score` based on these by calling
    :func:`calculate_score`

    Parameters
    ----------
    matched : :class:`GlycopeptideMatch`
        Match to be scored
    theoretical : :class:`TheoreticalGlycopeptide`
        The theoretical reference of `matched`
        Description

    Returns
    -------
    float
    """
    matched.mean_coverage = mean_coverage2(matched)
    matched.mean_hexnac_coverage = mean_hexnac_coverage(matched, theoretical)
    return calculate_score(matched, **parameters)


def evaluate2(matched, theoretical, **parameters):
    matched.mean_coverage = mean_coverage2(matched)
    matched.mean_hexnac_coverage = mean_hexnac_coverage(matched, theoretical)
    return calculate_score2(matched, **parameters)


def calculate_score(matched, backbone_weight=0.5, hexnac_weight=0.5, stub_weight=0.2):
    """
    Calculates the `ms2_score` of `matched`, setting it and returning

    Parameters
    ----------
    matched : :class:`GlycopeptideMatch`
        Match to be scored
    backbone_weight : float, optional
        Weight to be given to bare backbone ions. Defaults to 0.5
    hexnac_weight : float, optional
        Weight to be given to glycosylated backbone ions. Defaults to 0.5
    stub_weight : float, optional
        Weight to be given to the stub ions. Defaults to 0.2

    Returns
    -------
    float
    """
    matched.ms2_score = (((matched.mean_coverage * backbone_weight) +
                          (matched.mean_hexnac_coverage * hexnac_weight)) * (1 - stub_weight)) +\
                        (observed_vs_enumerated_stub(matched) * stub_weight)
    return matched.ms2_score


def calculate_score2(matched, backbone_weight=0.5, hexnac_weight=0.5, stub_weight=0.2):
    score = ((matched.mean_coverage * backbone_weight) +
             (matched.mean_hexnac_coverage * hexnac_weight))
    scaler = (observed_vs_enumerated_stub(matched) * stub_weight)
    score *= 1 - scaler
    matched.ms2_score = score
    return score


def stream_backbone_coordinate(iterable):
    return set(imap(int, chain.from_iterable(imap(bare_ion_pattern.findall, imap(key_getter, iterable)))))


def incvec_spread(vector, p, s=1):
    for i in range(s - 1, p, 1 if p > s else -1):
        vector[i] += 1
    return vector


def incvec_spot(vector, p, s=1):
    vector[p] = True
    return vector


def bitwise_or(a, b):
    """
    Compute the max of each parallel position
    of two sequences

    Parameters
    ----------
    a : Iterable
    b : Iterable
        Iterable of numbers

    Returns
    -------
    list
        List of the maximum of each position between `a` and `b`
    """
    return list(imap(max, zip(a, b)))


def mean(vector):
    """
    Arithmetic Mean

    Parameters
    ----------
    vector : Iterable
        Iterable of numbers

    Returns
    -------
    float
        Description
    """
    return sum(vector) / float(len(vector))


def vsum(a, b):
    """
    Sum two iterables

    Parameters
    ----------
    a : Iterable
    b : Iterable
        Iterable of numbers

    Returns
    -------
    list
        The sum of the `a` and `b`
    """
    return [x + y for x, y in zip(a, b)]


def compute_backbone_vectors(matched, incvec=incvec_spread):
    sequence = Sequence(matched.glycopeptide_sequence)
    sequence_length = len(sequence)

    b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.bare_b_ions):
        incvec(b_ions, b_ion)

    glycosylated_b_ions = [0.] * sequence_length
    for b_ion in stream_backbone_coordinate(matched.glycosylated_b_ions):
        incvec(glycosylated_b_ions, b_ion)

    y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.bare_y_ions):
        incvec(y_ions, sequence_length - y_ion, sequence_length)

    glycosylated_y_ions = [0.] * sequence_length
    for y_ion in stream_backbone_coordinate(matched.glycosylated_y_ions):
        incvec(glycosylated_y_ions, sequence_length - y_ion, sequence_length)

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


def mean_coverage2(matched):
    b_ions, glycosylated_b_ions, y_ions, glycosylated_y_ions = compute_backbone_vectors(matched, incvec_spot)
    coverage = vsum(bitwise_or(b_ions, glycosylated_b_ions), bitwise_or(y_ions, glycosylated_y_ions))
    # This normalization constant log2(3, 2) keeps the sum very close to 1.0 in a perfect case.
    coverage = sum([math.log(1 + i, 2)/_normalizer for i in coverage]) / (1. * len(coverage))
    return coverage


def mean_hexnac_coverage(matched, theoretical):
    b_observed = (stream_backbone_coordinate(matched.glycosylated_b_ions))
    y_observed = (stream_backbone_coordinate(matched.glycosylated_y_ions))

    b_enumerated = (stream_backbone_coordinate(theoretical.glycosylated_b_ions))
    y_enumerated = (stream_backbone_coordinate(theoretical.glycosylated_y_ions))
    return (len(b_observed) + len(y_observed)) / float(len(b_enumerated) + len(y_enumerated))


def observed_vs_enumerated_stub(matched):
    obs = len({stub['key'].split(':')[0] for stub in matched.stub_ions})
    return min(obs / _expected_stub_ion_count, 1.0)


def compute_percent_uncovered(matched):
    bs, gbs, ys, gys = compute_backbone_vectors(matched)
    backbone = vsum(bitwise_or(bs, gbs), bitwise_or(ys, gys))
    is_zero = 0
    for i in backbone:
        if i == 0:
            is_zero += 1
    return is_zero / float(len(matched))


class SimpleSpectrumScorer(GlycopeptideSpectrumMatchScorer):
    def __init__(self):
        super(SimpleSpectrumScorer, self).__init__("simple_ms2_score", "ms2_score")

    def evaluate(self, *args, **kwargs):
        evaluate(*args, **kwargs)


class SimpleScorer(ScorerBase):
    def __init__(self):
        super(SimpleScorer, self).__init__("ms2_score")

    def evaluate(self, *args, **kwargs):
        evaluate(*args, **kwargs)


class SimpleSpectrumScorer2(GlycopeptideSpectrumMatchScorer):
    def __init__(self):
        super(SimpleSpectrumScorer2, self).__init__("simple2_ms2_score", "ms2_score")

    def evaluate(self, *args, **kwargs):
        evaluate2(*args, **kwargs)


class SimpleScorer2(ScorerBase):
    def __init__(self):
        super(SimpleScorer2, self).__init__("ms2_score")

    def evaluate(self, *args, **kwargs):
        evaluate2(*args, **kwargs)
