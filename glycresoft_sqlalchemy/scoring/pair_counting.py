from collections import Counter, defaultdict
import json
import itertools
import operator
import pickle
import numpy as np

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.utils.collectiontools import flatten

from glycresoft_sqlalchemy.scoring import simple_scoring_algorithm


from .base import GlycopeptideSpectrumMatchScorer


def make_fragment_map(glycopeptide_sequence, fragment_series="by"):
    fragment_map = {}
    seq = sequence.Sequence(glycopeptide_sequence)

    for series in fragment_series:
        for fragment in flatten(seq.get_fragments(series)):
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

        self.series = "by"

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
            for series in self.series:
                for fragment in flatten(seq.get_fragments(series)):
                    n_term, c_term = fragment.flanking_amino_acids
                    possible_pairs[n_term, c_term] += count
                    possible_n_term[n_term] += count
                    possible_c_term[c_term] += count

        self.possible_pairs = possible_pairs
        self.possible_n_term = possible_n_term
        self.possible_c_term = possible_c_term

    def process_match(self, glycopeptide_match):
        key_seq = (glycopeptide_match.glycopeptide_sequence)
        fragment_map = make_fragment_map(key_seq, self.series)

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
            pair = tuple(pair)
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

            "series": self.series
        }
        pickle.dump(d, stream)

    @classmethod
    def load(cls, stream):
        d = pickle.load(stream)
        inst = cls()
        inst.series = d['series']
        inst.observed_pairs = Counter(d["observed_pairs"])
        inst.observed_n_term = Counter(d["observed_n_term"])
        inst.observed_c_term = Counter(d["observed_c_term"])
        inst.observed_sequences = Counter(d["observed_sequences"])

        inst.possible_pairs = Counter(d["possible_pairs"])
        inst.possible_n_term = Counter(d["possible_n_term"])
        inst.possible_c_term = Counter(d["possible_c_term"])
        return inst

    def train(self, matches):
        map(self.process_match, matches)
        self.total_possible_outcomes()


def extract_matched_fragments(match):
    keys = set()
    for f in match.bare_b_ions:
        if f['key'] in keys:
            continue
        keys.add(f['key'])
        yield f
    for f in match.bare_y_ions:
        if f['key'] in keys:
            continue
        keys.add(f['key'])
        yield f
    for f in match.glycosylated_b_ions:
        if f['key'] in keys:
            continue
        keys.add(f['key'])
        yield f
    for f in match.glycosylated_y_ions:
        if f['key'] in keys:
            continue
        keys.add(f['key'])
        yield f


class FrequencyScorer(GlycopeptideSpectrumMatchScorer):
    def __init__(self, frequency_counter, use_interaction=False):
        super(FrequencyScorer, self).__init__("pair_counting_ms2_score", "ms2_score")
        self.use_interaction = use_interaction
        self.frequency_counter = frequency_counter

    def evaluate(self, matched, theoretical, **parameters):
        matched.mean_coverage = simple_scoring_algorithm.mean_coverage2(matched)
        matched.mean_hexnac_coverage = simple_scoring_algorithm.mean_hexnac_coverage(matched, theoretical)
        ms2_score = matched.ms2_score = self.calculate_score(matched, theoretical=theoretical, **parameters)
        return ms2_score

    def backbone_score(self, matched, **parameters):
        fragment_map = make_fragment_map(matched.glycopeptide_sequence)
        total = self._maximum_score(fragment_map)
        observed = 0
        if self.use_interaction:
            for frag in extract_matched_fragments(matched):
                observed += self.frequency_counter.pair_probability(fragment_map[frag['key']].flanking_amino_acids)
        else:
            track_site = set()
            for frag in extract_matched_fragments(matched):
                f = fragment_map[frag['key']]
                position = f.position
                n_term, c_term = f.flanking_amino_acids
                score = self.frequency_counter.n_term_probability(
                    n_term) * self.frequency_counter.c_term_probability(c_term)
                weight = 0.6 if position not in track_site else 0.4
                track_site.add(position)
                observed += score * weight

        return observed / total

    def calculate_score(self, matched, stub_weight=0.2, **parameters):
        backbone_score = self.backbone_score(matched, **parameters)
        stub_ion_score = simple_scoring_algorithm.observed_vs_enumerated_stub(matched)
        backbone_weight = 1 - stub_weight
        return (backbone_weight * backbone_score) + (stub_weight * stub_ion_score)

    def _maximum_score(self, fragment_map):
        total = 0
        frags = fragment_map.values()
        if self.use_interaction:
            for frag in frags:
                total += self.frequency_counter.pair_probability(frag.flanking_amino_acids)
        else:
            for frag in frags:
                n_term, c_term = fragment_map[frag.name].flanking_amino_acids
                score = self.frequency_counter.n_term_probability(
                    n_term) * self.frequency_counter.c_term_probability(c_term)
                total += score * 0.5

        return total
