from collections import Counter, defaultdict
import json
import itertools
import operator

import numpy as np

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.utils import simple_repr
from glycresoft_sqlalchemy.utils.collectiontools import flatten

from glycresoft_sqlalchemy.scoring.offset_frequency.peak_relations import MatchedSpectrum

from glycresoft_sqlalchemy.data_model import DatabaseManager, GlycopeptideMatch


def make_fragment_map(glycopeptide_sequence, include_y=True, include_b=True):
    fragment_map = {}
    seq = sequence.Sequence(glycopeptide_sequence)

    if include_y:
        for fragment in flatten(seq.get_fragments('y')):
            fragment_map[fragment.name] = fragment
    if include_b:
        for fragment in flatten(seq.get_fragments('b')):
            fragment_map[fragment.name] = fragment
    return fragment_map


class FrequencyCounter(object):
    def __init__(self):
        self.observed_pairs = Counter()
        self.observed_n_term = Counter()
        self.observed_c_term = Counter()
        self.observed_sequences = Counter()

        self.possible_pairs = Counter()
        self.possible_n_term = Counter()
        self.possible_c_term = Counter()

        self.include_b = True
        self.include_y = True

    def add_sequence(self, sequence):
        self.observed_sequences[sequence] += 1

    def add_fragment(self, fragment):
        n_term, c_term = fragment.flanking_amino_acids
        self.observed_pairs[n_term, c_term] += 1
        self.observed_n_term[n_term] += 1
        self.observed_c_term[c_term] += 1

    def total_possible_outcomes(self):
        possible_pairs = Counter()
        possible_n_term = Counter()
        possible_c_term = Counter()
        for seq, count in self.observed_sequences.items():
            seq = sequence.Sequence(seq)
            if self.include_y:
                for fragment in flatten(seq.get_fragments('y')):
                    n_term, c_term = fragment.flanking_amino_acids
                    possible_pairs[n_term, c_term] += count
                    possible_n_term[n_term] += count
                    possible_c_term[c_term] += count

            if self.include_b:
                for fragment in flatten(seq.get_fragments('b')):
                    n_term, c_term = fragment.flanking_amino_acids
                    possible_pairs[n_term, c_term] += count
                    possible_n_term[n_term] += count
                    possible_c_term[c_term] += count

        self.possible_pairs = possible_pairs
        self.possible_n_term = possible_n_term
        self.possible_c_term = possible_c_term

    def process_match(self, glycopeptide_match):
        key_seq = (glycopeptide_match.glycopeptide_sequence)
        fragment_map = make_fragment_map(key_seq, self.include_y, self.include_b)

        for spectrum_match in glycopeptide_match.spectrum_matches:
            observed = set()
            for case in itertools.chain.from_iterable(spectrum_match.peak_match_map.values()):
                fkey = case['key']
                if case['intensity'] < 250:
                    continue
                if fkey in fragment_map and fkey not in observed:
                    observed.add(fkey)
                    f = fragment_map[fkey]
                    self.add_fragment(f)
            self.add_sequence(key_seq)

    def n_term_probability(self, residue=None):
        if residue is not None:
            return self.observed_n_term[residue] / float(self.possible_n_term[residue])
        else:
            return {r: self.n_term_probability(r) for r in self.observed_n_term}

    def c_term_probability(self, residue=None):
        if residue is not None:
            return self.observed_c_term[residue] / float(self.possible_c_term[residue])
        else:
            return {r: self.c_term_probability(r) for r in self.observed_c_term}

    def pair_probability(self, pair=None):
        if pair is not None:
            return self.observed_pairs[pair] / float(self.possible_pairs[pair])
        else:
            return {r: self.pair_probability(r) for r in self.observed_pairs}

    def residue_probability(self, residue=None):
        if residue is not None:
            return (self.observed_n_term[residue] + self.observed_c_term[residue]) / float(
                self.possible_n_term[residue] + self.possible_c_term[residue])
        else:
            return {r: self.residue_probability(r) for r in (set(self.observed_c_term) | set(self.observed_n_term))}

    def save(self, stream):
        d = {
            "observed_pairs": self.observed_pairs,
            "observed_n_term": self.observed_n_term,
            "observed_c_term": self.observed_c_term,
            "observed_sequences": self.observed_sequences,

            "possible_pairs": self.possible_pairs,
            "possible_n_term": self.possible_n_term,
            "possible_c_term": self.possible_c_term,

            "include_b": self.include_b,
            "include_y": self.include_y
        }
        json.dump(d, stream)

    @classmethod
    def load(cls, stream):
        d = json.load(stream)
        inst = cls()
        inst.include_b = d['include_b']
        inst.include_y = d['include_y']
        inst.observed_pairs = Counter(d["observed_pairs"])
        inst.observed_n_term = Counter(d["observed_n_term"])
        inst.observed_c_term = Counter(d["observed_c_term"])
        inst.observed_sequences = Counter(d["observed_sequences"])

        inst.possible_pairs = Counter(d["possible_pairs"])
        inst.possible_n_term = Counter(d["possible_n_term"])
        inst.possible_c_term = Counter(d["possible_c_term"])
        return inst


_MAXIMA = 3e6
_NOISE_FACTOR = 2000 * 2e-9


def geometric_mean(values):
    return np.exp(1./len(values) * np.log(values).sum())


def fragment_transform(fragment_mean, min1=.99, max1=.99, min2=1, max2=1):
    lmin1, lmax1, lmin2, lmax2 = np.log([min1, max1, min2, max2])
    precursor_fragment_range_ratio = (min(lmax1 - lmin1, 0.1) / (lmax2 - lmin2))
    fragment_series_range = np.log(fragment_mean) - lmin2 + lmax1

    # print lmin1, lmax1, lmin2, lmax2
    # print "precursor_fragment_range_ratio", precursor_fragment_range_ratio
    # print "fragment_series_range", fragment_series_range, np.log(fragment_mean), lmin2, lmax1
    # print "(precursor_fragment_range_ratio * fragment_series_range)", (precursor_fragment_range_ratio * fragment_series_range)

    return np.exp(precursor_fragment_range_ratio * fragment_series_range)


class ScoreMatch(object):

    def __init__(self, match, frequency_counter, tolerance=2e-5, maxima=_MAXIMA, noise_factor=_NOISE_FACTOR):
        # Constants
        self.reference = match.glycopeptide_match.theoretical_reference
        match = MatchedSpectrum(match)
        match.reindex_peak_matches()
        self.spectrum = match.peak_list
        self.match = match
        self.frequency_counter = frequency_counter
        self.tolerance = tolerance
        self.maxima = maxima
        self.noise_factor = noise_factor

        # Containers
        self.possible_peak_scores = {}
        self.possible_peaks = []
        self.n_peaks_gen = None
        self.noise_peaks_total = None
        self.explainable_peaks_total = None
        self.score = None

    def calculate_score(self):
        MAXIMA = self.maxima
        NOISE_FACTOR = self.noise_factor
        tolerance = self.tolerance
        match = self.match
        spectrum = self.spectrum
        frequency_counter = self.frequency_counter

        possible_peak_scores = self.possible_peak_scores = {}
        reference = self.reference
        possible_peaks = [p for p in reference.fragments()]
        self.n_peaks_gen = n_peaks_gen = len(possible_peaks)
        self.noise_peaks_total = noise_peaks_total = NOISE_FACTOR * (MAXIMA - 2 * (n_peaks_gen * tolerance))
        self.explainable_peaks_total = explainable_peaks_total = 0.

        fragment_map = make_fragment_map(reference.glycopeptide_sequence)

        for peak in possible_peaks:
            fkey = peak['key']
            if fkey in fragment_map:
                n_term, c_term = fragment_map[fkey].flanking_amino_acids
                s = frequency_counter.residue_probability(n_term)
                s *= frequency_counter.residue_probability(c_term)
                s *= (2 * tolerance)
                explainable_peaks_total += s
                possible_peak_scores[fkey] = s
            else:
                # TODO:
                # No model for oxonium ions or stub ions. Future
                # work would build one.
                s = (2 * tolerance) * NOISE_FACTOR * 10
                explainable_peaks_total += s
                possible_peak_scores[fkey] = s

        total_area = noise_peaks_total + explainable_peaks_total

        self.accumulated_score = accumulated_score = defaultdict(float)
        peak_match_map = match.peak_match_map
        for peak in spectrum:
            if peak.id in peak_match_map:
                for peak_match in peak_match_map[peak.id]:
                    fkey = peak_match.key
                    accumulated_score[peak.id] = max(possible_peak_scores[fkey]/total_area, accumulated_score[peak.id])
            else:
                accumulated_score[peak.id] = NOISE_FACTOR / total_area

        self.score = geometric_mean(np.array(accumulated_score.values()))
        return self.score

    def min_fragment_score(self):
        return min(self.possible_peak_scores.values())

    def max_fragment_score(self):
        return max(self.possible_peak_scores.values())

    __repr__ = simple_repr
