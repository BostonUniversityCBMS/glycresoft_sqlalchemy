import itertools
import random

import numpy as np

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.data_model import GlycopeptideSpectrumMatch, PipelineModule
from glycresoft_sqlalchemy.utils import collectiontools

from .simple_scoring_algorithm import split_ion_list

Sequence = sequence.Sequence

chain_iterable = itertools.chain.from_iterable


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


def ppm_error(x, y):
    return (x - y) / y


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


def feature_function_estimator(gsms, feature_function, kind='b'):
    total_on_kind_satisfied = 0
    total_off_kind_satisfied = 0
    total_on_kind = 0
    total_off_kind = 0
    for gsm in gsms:
        sequence = Sequence(gsm.glycopeptide_match.glycopeptide_sequence)
        fragments = collectiontools.flatten(sequence.get_fragments(kind))
        spectrum = gsm.spectrum
        

