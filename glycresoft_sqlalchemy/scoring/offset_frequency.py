import itertools
import random
import operator
import numpy as np

from glycresoft_sqlalchemy.structure import sequence, composition, modification
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatch, PipelineModule
from glycresoft_sqlalchemy.utils import collectiontools
from glycresoft_sqlalchemy.utils import common_math

PODPeak = common_math.DPeak
ppm_error = common_math.ppm_error
mass_offset_match = common_math.mass_offset_match
MassOffsetFeature = common_math.MassOffsetFeature
search_spectrum = common_math.search_spectrum

Sequence = sequence.Sequence

chain_iterable = itertools.chain.from_iterable

get_intensity = operator.attrgetter("intensity")


# Lacking a reasonable definition of the "space between fragmentation sites"
SMALLEST_UNIT = 1000 * 2e-5

# BEGIN REFERENCE Estimation
# Simple one parameter estimator functions for learning the basic alpha, beta and p parameters


def offset_frequency(gsms, kind='b'):
    total_sites = 0
    total_explained = 0
    for gsm in gsms:
        n_frag_sites = count_fragmentation_sites(gsm.glycopeptide_match.glycopeptide_sequence, kind)
        kind_explained = sum([1 for i in chain_iterable(gsm.peak_match_map.values()) if i['key'][0] == kind])
        total_sites += n_frag_sites
        total_explained += kind_explained
    return (total_explained)/float(total_sites)


def unknown_peak_rate(gsms, kind='b'):
    total_sparsity = 0
    total_unexplained = 0

    for gsm in gsms:
        sequence = gsm.glycopeptide_match.glycopeptide_sequence
        n_frag_sites = count_fragmentation_sites(sequence, kind)
        kind_explained = sum([1 for i in chain_iterable(gsm.peak_match_map.values()) if i['key'][0] == kind])
        peaks_unexplained = gsm.peaks_unexplained + (gsm.peaks_explained - kind_explained)
        total_unexplained += peaks_unexplained
        total_sparsity += estimate_fragment_sparsity(sequence, kind) - n_frag_sites

    return total_unexplained / float(total_sparsity)


def count_fragmentation_sites(sequence, kind='b'):
    sequence = Sequence(sequence)
    fragmentation_sites = len(collectiontools.flatten(sequence.get_fragments(kind)))
    return fragmentation_sites


def prior_fragment_probability(gsms, kind='b'):
    hits = 0
    for gsm in gsms:
        sequence = Sequence(gsm.glycopeptide_match.glycopeptide_sequence)
        random_mass = np.random.uniform(0, sequence.mass)
        for fragment in collectiontools.flatten(sequence.get_fragments(kind)):
            if abs(ppm_error(fragment.mass, random_mass)) <= 2e-5:
                hits += 1
            elif fragment.mass - (random_mass + 230.) > 0:
                break
    return hits / float(len(gsms))

# END REFERENCE Estimation


def estimate_fragment_sparsity(sequence, kind='b'):
    return Sequence(sequence).mass / SMALLEST_UNIT


def estimate_offset_parameters(gsms, kind='b'):
    total_sites = 0
    total_explained = 0
    total_sparsity = 0
    total_unexplained = 0
    random_hits = 0

    i = 0.
    for gsm in gsms:
        sequence = Sequence(gsm.glycopeptide_match.glycopeptide_sequence)
        fragments = collectiontools.flatten(sequence.get_fragments(kind))

        n_frag_sites = len(fragments)

        kind_explained = sum([1 for d in chain_iterable(gsm.peak_match_map.values()) if d['key'][0] == kind])
        total_sites += n_frag_sites
        total_explained += kind_explained

        peaks_unexplained = gsm.peaks_unexplained + (gsm.peaks_explained - kind_explained)

        total_unexplained += peaks_unexplained
        total_sparsity += estimate_fragment_sparsity(sequence, kind) - n_frag_sites

        random_mass = np.random.uniform(56., sequence.mass + 38.)

        for fragment in fragments:
            if abs(ppm_error(fragment.mass, random_mass)) <= (2e-5):
                random_hits += 1
            elif fragment.mass - (random_mass + 230.) > 0:
                break
        i += 1

    alpha = total_explained / float(total_sites)
    beta = total_unexplained / float(total_sparsity)
    prior_fragment_probability = max(random_hits / i, 0.0001)
    return alpha, beta, prior_fragment_probability


def probability_of_peak_explained(offset_frequency, unknown_peak_rate, prior_fragment_probability):
    a = (prior_fragment_probability * offset_frequency)
    b = (1 - prior_fragment_probability) * unknown_peak_rate
    return a / (a + b)


def make_offset_function(offset=.0, tolerance=2e-5, name=None):
    return MassOffsetFeature(offset=offset, tolerance=tolerance, name=name)


def intensity_ratio_function(peak1, peak2):
    ratio = peak1.intensity / float(peak2.intensity)
    if ratio >= 5:
        return -4
    elif 2.5 <= ratio < 5:
        return -3
    elif 1.7 <= ratio < 2.5:
        return -2
    elif 1.3 <= ratio < 1.7:
        return -1
    elif 1.0 <= ratio < 1.3:
        return 0
    elif 0.8 <= ratio < 1.0:
        return 1
    elif 0.6 <= ratio < 0.8:
        return 2
    elif 0.4 <= ratio < 0.6:
        return 3
    elif 0.2 <= ratio < 0.4:
        return 4
    elif 0. <= ratio < 0.2:
        return 5


def intensity_rank(peak_list, minimum_intensity=100.):
    peak_list = sorted(peak_list, key=get_intensity, reverse=True)
    i = 0
    rank = 10
    tailing = 6
    for p in peak_list:
        if p.intensity < minimum_intensity:
            p.rank = 0
            continue
        i += 1
        if i == 10 and rank != 0:
            if rank == 1:
                if tailing != 0:
                    i = 0
                    tailing -= 1
                else:
                    i = 0
                    rank -= 1
            else:
                i = 0
                rank -= 1
        if rank == 0:
            break
        p.rank = rank


class PeakRelation(object):
    def __init__(self, from_peak, to_peak, feature, intensity_ratio=None, kind=None):
        if intensity_ratio is None:
            intensity_ratio = intensity_ratio_function(from_peak, to_peak)
        self.from_peak = from_peak
        self.to_peak = to_peak
        self.feature = feature
        self.intensity_ratio = intensity_ratio
        self.same_terminal = None
        self.from_charge = from_peak.charge
        self.to_charge = to_peak.charge
        self.kind = kind or "Noise"

    def __repr__(self):
        template = "<PeakRelation {s.from_peak.neutral_mass}({s.from_charge}) ->" +\
            " {s.to_peak.neutral_mass}({s.to_charge}) by {s.feature.name} on {s.kind}>"
        return template.format(s=self)


def feature_function_estimator(gsms, feature_function, kind='b'):
    total_on_kind_satisfied = 0.
    total_off_kind_satisfied = 0.
    total_on_kind = 0.
    total_off_kind = 0.
    peak_relations = []
    for gsm in gsms:
        spectrum = gsm.spectrum
        peaks = list(map(PODPeak, spectrum))
        intensity_rank(peaks)
        peaks = [p for p in peaks if p.rank > 0]
        related = []
        for peak in peaks:
            is_on_kind = any(k[0] == kind for k in gsm.peak_explained_by(peak.id))
            matches = search_spectrum(peak, peaks, feature_function)
            for match in matches:
                pr = PeakRelation(peak, match, feature_function, intensity_ratio_function(peak, match))
                related.append(pr)
                if is_on_kind:
                    print peak.id, match.id
                    total_on_kind_satisfied += 1
                    pr.kind = kind
                else:
                    total_off_kind_satisfied += 1
                    pr.kind = "Noise"
            if is_on_kind:
                total_on_kind += 1
            else:
                total_off_kind += 1
        peak_relations.append((gsm, related))

    return total_on_kind_satisfied / max(total_on_kind, 1), total_off_kind_satisfied / max(total_off_kind, 1), peak_relations


shifts = [
    make_offset_function(-composition.Composition("NH3").mass, name='Neutral-Loss-Ammonia'),
    make_offset_function(-composition.Composition("H2O").mass, name='Neutral-Loss-Water'),

]
