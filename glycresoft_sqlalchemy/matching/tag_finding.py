import operator

import numpy as np
import networkx as nx
import matplotlib
from matplotlib import pyplot

from collections import defaultdict
from ..data_model import TheoreticalGlycopeptide
from ..utils.common_math import ppm_error


from ..structure import sequence, residue, fragment
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms2 import residue_counter

fragment_shift = fragment.fragment_shift
neutral_mass_getter = operator.attrgetter('neutral_mass')
intensity_getter = operator.attrgetter("intensity")


Block = residue_counter.Block
default_match_tolerance = 2e-5


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
    nx.draw_networkx(graph, pos=position, ax=ax, **kwargs)
    nx.draw_networkx_edge_labels(
        graph,
        pos=position,
        edge_labels=edge_label_dict,
        ax=ax)
    return ax


def delta_peak_list(peak_list):
    p = peak_list[0]
    deltas = []
    for next_peak in peak_list[1:]:
        deltas.append(next_peak.neutral_mass - p.neutral_mass)
        p = next_peak
    return deltas


class TagFinder(object):
    def __init__(self, peak_list):
        peak_list = [p for p in peak_list if p.intensity > 150.]
        self.peak_list = sorted(peak_list, key=neutral_mass_getter)
        self.spectrum_graph = nx.MultiGraph()

    def find_tags(self, blocks, tolerance=default_match_tolerance):
        peak_list = self.peak_list
        spectrum_graph = self.spectrum_graph

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

    def label_paths(self, paths=None, peak_ids=False):
        if paths is None:
            paths = self.tag_paths()
        seen = set()
        for path in paths:
            g = (''.join(map(str, seq)) for seq in label_path(self.spectrum_graph, path))
            for seq in g:
                seen.add(seq[:-1])
                if seq in seen:
                    continue
                yield seq if not peak_ids else seq, path
    sequences = label_paths
