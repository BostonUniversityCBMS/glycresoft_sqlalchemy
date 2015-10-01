import operator

import numpy as np
import networkx as nx
import matplotlib
from matplotlib import pyplot
from glypy import MonosaccharideResidue

from collections import defaultdict
from ..data_model import TheoreticalGlycopeptide
from ..utils.common_math import ppm_error, tol_ppm_error, search_spectrum_by_mass


from ..structure import sequence, residue, fragment
from glycresoft_sqlalchemy.structure.sequence_composition import (
    AminoAcidResidue, MonosaccharideResidueAdapter, SequenceComposition,
    AminoAcidSequenceBuildingBlock as Block)

fragment_shift = fragment.fragment_shift
neutral_mass_getter = operator.attrgetter('neutral_mass')
intensity_getter = operator.attrgetter("intensity")
default_match_tolerance = 2e-5
Sequence = sequence.Sequence


def find_tags(peak_list, blocks, tolerance=default_match_tolerance):
    peak_list = sorted(peak_list, key=neutral_mass_getter)
    spectrum_graph = nx.DiGraph()
    for peak in peak_list:
        spectrum_graph.add_node(
            peak.id, neutral_mass=peak.neutral_mass,
            charge=peak.charge, intensity=peak.intensity)

    for a_peak in peak_list:
        a_mass = a_peak.neutral_mass
        for block in blocks:
            a_mass += block.neutral_mass
            for b_peak in peak_list:
                if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                    spectrum_graph.add_edge(a_peak.id, b_peak.id, residue=block)
            a_mass -= block.neutral_mass
    return spectrum_graph


def longest_paths(G):
    dist = {}  # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0] + 1, v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    pathes = sorted(dist.items(), key=lambda x: x[1], reverse=True)
    last = set()
    for node, length_ in pathes:
        length, _ = length_
        if length < 1:
            break
        path = []
        while length > 0:
            path.append(node)
            length, node = dist[node]
        node_set = set(path)
        if node_set.issubset(last):
            continue
        else:
            last = node_set
        yield list(reversed(path))


def label_path(G, path):
    last = path[0]
    sequences = [[]]
    for i in path[1:]:
        new_sequences = []
        for edge in G.edge[last][i].values():
            for sequence in sequences:
                new_seq = sequence[:]
                new_seq.append(edge['residue'])
                new_sequences.append(new_seq)
        sequences = new_sequences
        last = i
    return sequences


def layout_graph(graph):
    nodes = graph.node.items()
    coords = {}
    for ix_node in (nodes):
        ix, node = ix_node
        coords[ix] = (node['neutral_mass'], node['intensity'])
    return coords


def edge_labels(graph):
    labels = defaultdict(list)
    for ij in graph.edges():
        i, j = ij
        edge = graph.edge[i][j]
        for edge_case, edge_data in edge.items():
            labels[ij].append("{residue}".format(**edge_data))
    labels = {k: u'/\n'.join(v) for k, v in labels.items()}
    return labels


def draw_graph(graph, **kwargs):
    position = layout_graph(graph)
    kwargs.setdefault("node_size", 5)
    kwargs.setdefault("with_labels", False)
    edge_label_dict = edge_labels(graph)
    fig, ax = pyplot.subplots(1, 1)
    nx.draw_networkx(graph, pos=position, ax=ax, arrows=False, **kwargs)
    nx.draw_networkx_edge_labels(
        graph,
        pos=position,
        edge_labels=edge_label_dict,
        arrows=False,
        ax=ax)
    return ax


def delta_peak_list(peak_list):
    p = peak_list[0]
    deltas = []
    for next_peak in peak_list[1:]:
        deltas.append(next_peak.neutral_mass - p.neutral_mass)
        p = next_peak
    return deltas


class FeatureOffsetBuildingBlockAdapter(object):
    def __init__(self, feature):
        self.neutral_mass = feature.offset
        self.symbol = "F(%s)" % str(feature)
        self.feature = feature

    def __hash__(self):
        return hash(self.symbol)

    def __eq__(self, other):
        if hasattr(other, 'symbol'):
            return self.symbol == other.symbol
        else:
            return self.symbol == other

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return self.symbol


class SequencePath(object):
    def __init__(self, peaks, sequence_=None):
        self.peaks = peaks
        try:
            self.sequence = Sequence(sequence_)
        except:
            self.sequence = sequence_
        self.weight = 0.
        for peak in peaks:
            self.weight += peak.intensity
        self.low_mass = self.peaks[0].neutral_mass
        self.high_mass = self.peaks[-1].neutral_mass

    def __repr__(self):
        return "SequencePath(%s %d %0.3f %0.3f-%0.3f)" % (
            self.sequence, len(self.peaks), self.weight,
            self.peaks[0].neutral_mass, self.peaks[-1].neutral_mass)

    def __iter__(self):
        return iter(self.peaks)


class TagFinder(object):
    def __init__(self, peak_list, blocks=None):
        peak_list = [p for p in peak_list if p.intensity > 150.]
        self.peak_list = sorted(peak_list, key=neutral_mass_getter)
        self.spectrum_graph = nx.MultiDiGraph()
        self.blocks = blocks

    def find_tags(self, blocks=None, tolerance=default_match_tolerance):
        peak_list = self.peak_list
        spectrum_graph = self.spectrum_graph

        if blocks is None:
            blocks = self.blocks

        for peak in peak_list:
            spectrum_graph.add_node(
                peak.id, neutral_mass=peak.neutral_mass,
                charge=peak.charge, intensity=peak.intensity)

        for i, a_peak in enumerate(peak_list):
            a_mass = a_peak.neutral_mass
            for block in blocks:
                a_mass += block.neutral_mass
                for b_peak in peak_list[:]:
                    if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                        spectrum_graph.add_edge(a_peak.id, b_peak.id, residue=block)
                a_mass -= block.neutral_mass
        return spectrum_graph

    def __getitem__(self, ids):
        try:
            iter(ids)
        except:
            ids = [ids]
        result = []
        for peak in self.peak_list:
            if peak.id in ids:
                result.append(peak)
        return result

    def tag_paths(self):
        return longest_paths(self.spectrum_graph)

    def label_paths(self, paths=None):
        if paths is None:
            paths = self.tag_paths()
        seen = set()
        results = []
        for path in paths:
            g = (''.join(map(str, seq)) for seq in label_path(self.spectrum_graph, path))
            for seq in g:
                seen.add(tuple(path[:-1]))
                if tuple(path) in seen:
                    continue
                seq = SequencePath(self[path], seq)
                results.append(seq)
        results = sorted(results, key=operator.attrgetter("weight"), reverse=True)
        return results
    sequences = label_paths

    def join_paths(self, blocks=None):
        if blocks is None:
            blocks = self.blocks
        paths = sorted(list(self.sequencs()), lambda x: (x.low_mass, x.high_mass))
        for path in paths:
            for other_path in paths:
                if other_path.high_mass < path.low_mass:
                    gap = path.low_mass - other_path.high_mass
                    if gap < 2000:
                        segment = resolve_sequence_segment(gap, blocks)
                        print other_path, "->", segment, "->", path
                elif path.high_mass < other_path.low_mass:
                    gap = other_path.low_mass - path.high_mass
                    if gap < 2000:
                        segment = resolve_sequence_segment(gap, blocks)
                        print path, "->", segment, "->", other_path
            else:
                if path.low_mass < 2000:
                    y_start = resolve_sequence_segment(gap, blocks, shift=fragment_shift['y'])
                    b_start = resolve_sequence_segment(gap, blocks, shift=fragment_shift['b'])
                    print y_start, b_start


class DegenerateSymbol(object):
    def __init__(self, symbols):
        self.symbols = symbols

    def __str__(self):
        return "[%s]" % '|'.join(map(str, self.symbols))

    def __eq__(self, other):
        for sym in self.symbols:
            if sym == other:
                return True
        return False

    def __ne__(self, other):
        return not self == other


def resolve_sequence_segment(total_mass, blocks, shift=0, tolerance=default_match_tolerance):
    solutions = set()

    def extend_segment(base, blocks):
        for block in blocks:
            ext = base.clone()
            ext[block] += 1
            ext._mass = base.mass + block.neutral_mass
            yield ext

    candidates = [SequenceComposition()]
    next_round = set()
    while True:
        for candidate in candidates:
            extents = set(extend_segment(candidate, blocks))
            for case in extents:
                case_mass = case.mass + shift
                if abs(ppm_error(case_mass, total_mass)) <= tolerance:
                    solutions.add(case)
                elif case_mass < total_mass + 1:
                    next_round.add(case)
        if len(next_round) == 0:
            break
        else:
            candidates = next_round
            next_round = set()

    return solutions
