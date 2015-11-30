from collections import Counter, defaultdict
import json
import itertools
import operator

import numpy as np

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.utils.collectiontools import flatten

from glycresoft_sqlalchemy.scoring import simple_scoring_algorithm


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


class FrequencyScorer(object):
    def __init__(self, frequency_counter):
        self.frequency_counter = frequency_counter

    def evaluate(self, matched, theoretical, **parameters):
        matched.mean_coverage = simple_scoring_algorithm.mean_coverage2(matched)
        matched.mean_hexnac_coverage = simple_scoring_algorithm.mean_hexnac_coverage(matched, theoretical)
        matched.ms2_score = self.calculate_score(matched, theoretical, **parameters)

    def calculate_score(self, matched, theoretical, **parameters):
        return 0
