import operator
import itertools
from collections import defaultdict

import numpy as np

from ..data_model import Peak, TheoreticalGlycanCombination
from ..utils.common_math import ppm_error, DPeak
from ..utils.tree import SuffixTree


from ..structure import sequence, fragment, composition

Composition = composition.Composition
IonSeries = fragment.IonSeries
neutral_mass_getter = operator.attrgetter('neutral_mass')
intensity_getter = operator.attrgetter("intensity")
default_match_tolerance = 2e-5
Sequence = sequence.Sequence


def delta_peak_list(peak_list):
    p = peak_list[0]
    deltas = []
    for next_peak in peak_list[1:]:
        deltas.append(next_peak.neutral_mass - p.neutral_mass)
        p = next_peak
    return deltas


class DegenerateSegmentBlock(object):
    def __init__(self, composition):
        self.composition = composition
        self.neutral_mass = composition.mass
        self._string = str(composition)
        self._hash = hash(self.composition.keys()[0])

    def __eq__(self, other):
        try:
            return self.composition == other.composition
        except:
            return str(self) == str(other)

    def __ne__(self, other):
        try:
            return self.composition != other.composition
        except:
            return str(self) != str(other)

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return self._string

    def orderings(self):
        for seq in self.composition.orderings():
            yield str(seq)

    def __contains__(self, item):
        return item in self.composition


class SpectrumGraph(object):
    def __init__(self, peaks, precursor_mass, glycosylation_resolver):
        self.peaks = sorted(peaks, key=neutral_mass_getter)
        self.source_node = SpectrumNode(DPeak(Peak(id=-1, neutral_mass=0, charge=0, intensity=0, scan_peak_index=-1)))
        self.nodes = [SpectrumNode(p) for p in self.peaks]
        self.edges = defaultdict(list)
        self.sink_nodes = []
        self.glycosylation_resolver = glycosylation_resolver

    def iteredges(self):
        for nodes, edges in self.edges.items():
            for edge in edges:
                yield edge

    def find_starting_points(self, blocks, offset, tolerance=2e-5):
        if isinstance(offset, IonSeries):
            offset_mass = offset.mass_shift
        else:
            offset_mass = offset
        a_mass = offset_mass
        a_peak = self.source_node
        for block in blocks:
            a_mass += block.neutral_mass
            for b_peak in self.nodes:
                if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                    self.edges[a_peak, b_peak].append(SourceSpectrumEdge(
                        a_peak, b_peak, block.neutral_mass, annotation=block, series=offset))
                elif b_peak.neutral_mass + (b_peak.neutral_mass * tolerance) > a_mass:
                    break
            a_mass -= block.neutral_mass

    def find_edges(self, blocks, tolerance=2e-5):
        for i, a_peak in enumerate(self.nodes):
            a_mass = a_peak.neutral_mass
            for block in blocks:
                a_mass += block.neutral_mass
                for b_peak in self.nodes[i:]:
                    if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                        self.edges[a_peak, b_peak].append(SpectrumEdge(
                            a_peak, b_peak, block.neutral_mass, annotation=block))
                a_mass -= block.neutral_mass

    def as_adjacency_matrix(self):
        matrix = np.zeros((len(self.nodes), len(self.nodes)))
        index_of = {node: i for i, node in enumerate(self.nodes)}
        for node, i in index_of.items():
            for edge in node.out_edges:
                j = index_of[edge.child]
                matrix[i][j] = 1
            for edge in node.in_edges:
                j = index_of[edge.parent]
                matrix[i][j] = 1
        return matrix

    def generate_contiguous(self, denovo_type=None):
        if denovo_type is None:
            denovo_type = DeNovoSequence
        suffix_tree = SuffixTree()

        for node in [self.source_node] + sorted(self.nodes, key=neutral_mass_getter):
            for contig in self._link(node, denovo_type=denovo_type):
                is_suffix = contig in suffix_tree

                if is_suffix:
                    continue

                suffix_tree.add_ngram(contig)
        return suffix_tree

    def _link(self, node, current_sequence=None, denovo_type=None):
        if denovo_type is None:
            denovo_type = DeNovoSequence

        if current_sequence is None:
            current_sequence = []
        streams = []
        if len(node.out_edges) > 0:
            for edge in node.out_edges:
                seq = list(current_sequence)
                seq.append(edge)
                streams.append(self._link(edge.child, seq, denovo_type=denovo_type))
            for stream in streams:
                for seq in stream:
                    yield seq
        else:
            yield denovo_type(current_sequence, graph=self)

    def search_from(self, node, blocks, tolerance=2e-5):
        a_mass = node.neutral_mass
        for block in blocks:
            a_mass += block.neutral_mass
            for b_peak in self.nodes:
                if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                    self.edges[node, b_peak].append(SpectrumEdge(
                        node, b_peak, block.neutral_mass, annotation=block))
                if a_mass < b_peak.neutral_mass + 10:
                    break
            a_mass -= block.neutral_mass

    def search_below(self, node, blocks, tolerance=2e-5):
        for a_peak in [self.source_peak] + self.peaks:
            a_mass = a_peak.neutral_mass
            for block in blocks:
                a_mass += block.neutral_mass
                if abs(ppm_error(a_mass, node.neutral_mass)) <= tolerance:
                    self.edges[node, node].append(SpectrumEdge(
                        node, node, block.neutral_mass, annotation=block))
                if a_mass > node.neutral_mass + 10:
                    break

    def _do_join_internal_segments(self, start, internal_segments, blocks, tolerance=2e-5):
        a_peak = start.last_peak
        for internal in internal_segments:
            b_peak = internal.first_peak
            if a_peak.neutral_mass > b_peak.neutral_mass:
                continue
            a_mass = a_peak.neutral_mass
            for block in blocks:
                if abs(ppm_error(a_mass + block.neutral_mass, b_peak.neutral_mass)) <= tolerance:
                    e = SpectrumEdge(a_peak, b_peak, block.neutral_mass, annotation=block)
                    self.edges[a_peak, b_peak].append(e)

    def join_internal_segments(self, denovo_sequences, blocks, tolerance=2e-5):
        sourced_sequences, internal_segments = partition_sequences(denovo_sequences)
        for sourced_sequence in sourced_sequences:
            self._do_join_internal_segments(
                sourced_sequence, internal_segments,
                blocks, tolerance=2e-5)


def partition_sequences(denovo_sequences):
    sourced_sequences = []
    internal_segments = []

    for seq in denovo_sequences:
        if seq.has_source():
            sourced_sequences.append(seq)
        else:
            internal_segments.append(seq)

    return sourced_sequences, internal_segments


class GlycanCombinationResolverBase(object):
    def __init__(self, total_mass, **kwargs):
        self.total_mass = total_mass

    @property
    def total_mass(self):
        return self._total_mass

    @total_mass.setter
    def total_mass(self, value):
        self._total_mass = value
        self.possible_solutions = self.possible_peptide_masses(value)

    def possible_peptide_masses(self):
        raise NotImplementedError()

    def solutions(self, peptide, tolerance=2e-5):
        raise NotImplementedError()


class GlycanCombinationResolverLinearSearcher(GlycanCombinationResolverBase):
    def __init__(self, total_mass, glycan_combinations, **kwargs):
        self.glycan_combinations = glycan_combinations
        super(GlycanCombinationResolverLinearSearcher, self).__init__(total_mass)

    def possible_peptide_masses(self, total_mass=None):
        if total_mass is None:
            total_mass = self.total_mass
        possible_solutions = [(total_mass - combn.dehydrated_mass(), combn) for combn in self.glycan_combinations]
        return possible_solutions

    def solutions(self, peptide, tolerance=2e-5):
        peptide_mass = peptide.sequence_mass()
        solutions = []
        for solution_mass, glycan_combination in self.possible_solutions:
            error = ppm_error(peptide_mass, solution_mass)
            # print error, glycan_combination, peptide_mass, solution_mass
            if abs(error) < tolerance:
                solutions.append(glycan_combination)
        return solutions


class GlycanCombinationResolverQuerySearcher(GlycanCombinationResolverBase):
    def __init__(self, total_mass, session, glycan_combination_filters,
                 glycan_combination_type=TheoreticalGlycanCombination, **kwargs):
        self.glycan_combination_filters = glycan_combination_filters
        self.session = session
        self.glycan_combination_type = glycan_combination_type
        self.glycan_combinations = session.query(glycan_combination_type).filter(*self.glycan_combination_filters)
        super(GlycanCombinationResolverQuerySearcher, self).__init__(total_mass)

    def possible_peptide_masses(self, total_mass=None):
        if total_mass is None:
            total_mass = self.total_mass

        expr = total_mass - self.glycan_combination_type.dehydrated_mass()
        expr.type.asdecimal = False

        possible_solutions = self.session.query(
            expr, self.glycan_combination_type).filter(*self.glycan_combination_filters).filter()
        return possible_solutions

    def solutions(self, peptide, tolerance=2e-5):
        peptide_mass = peptide.sequence_mass()
        gct = self.glycan_combination_type
        adjusted_mass = self.total_mass - gct.dehydrated_mass()
        mass_ppm_error = (peptide_mass - adjusted_mass) / adjusted_mass

        solutions = self.session.query(
            gct).filter(
            *self.glycan_combination_filters).filter(
            mass_ppm_error.between(-tolerance, tolerance)).all()
        return solutions


class SpectrumNode(object):
    def __init__(self, peak):
        self.in_edges = set()
        self.out_edges = set()
        self.peak = peak

        self.neutral_mass = peak.neutral_mass

    def __repr__(self):
        return "SpectrumNode({} {} {})".format(self.peak.id, self.peak.neutral_mass, self.peak.intensity)

    def __eq__(self, other):
        try:
            return self.peak == other.peak
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(int(self.neutral_mass))

    def iterout(self):
        for edge in self.out_edges:
            yield edge, [edge.child.iterout()]


class SpectrumEdge(object):
    def __init__(self, parent, child, mass_shift, annotation):
        self.parent = parent
        self.child = child
        self._hash = None
        self.mass_shift = mass_shift
        self.annotation = annotation
        parent.out_edges.add(self)
        child.in_edges.add(self)

    def __repr__(self):
        return "SpectrumEdge({}, {}, {}, {})".format(
            self.parent, self.child, self.mass_shift, self.annotation)

    def __eq__(self, other):
        try:
            return (self.parent == other.parent and self.child == other.child and
                    self.mass_shift == other.mass_shift and self.annotation == other.annotation)
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.annotation)
        return self._hash


class SourceSpectrumEdge(SpectrumEdge):
    def __init__(self, parent, child, mass_shift, annotation, series):
        super(SourceSpectrumEdge, self).__init__(parent, child, mass_shift, annotation)
        self.series = series

    def __repr__(self):
        return "SourceSpectrumEdge({}, {}, {}, {}|{})".format(
            self.parent, self.child, self.mass_shift, self.annotation, self.series)


class DeNovoSequence(object):

    def __init__(self, edges, graph):
        self._edges = tuple(edges)
        self.sequence = tuple(e.annotation for e in self._edges)
        self._string = None
        self._mass = None
        self._sequence = None
        self.graph = graph

    def reset(self):
        self._string = None

    def __str__(self):
        if self._string is None:
            self._string = ''.join(map(str, self.sequence))
        return self._string

    def __repr__(self):
        if self._string is None:
            self._string = ''.join(map(str, self.sequence))
        rep = "%s:%0.3f-%0.3f:%0.3e" % (
            self._string, self._edges[0].parent.neutral_mass,
            self._edges[-1].child.neutral_mass, self.total_energy())
        series = self.has_source()
        if series is not None:
            rep += ' [%s]' % series
        return rep

    def __hash__(self):
        return hash(self._edges)

    def total_energy(self):
        total = 0
        for edge in self._edges:
            total += edge.child.peak.intensity
        total += self._edges[0].parent.peak.intensity
        return total

    def __eq__(self, other):
        try:
            return self._edges == other._edges
        except:
            return str(self) == other

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return len(self._edges)

    def __iter__(self):
        return iter(self._edges)

    def __getitem__(self, i):
        return self._edges[i]

    def __contains__(self, seq):
        if isinstance(seq, str):
            return seq in str(self)
        else:
            return seq in self.sequence

    def mass_range(self):
        return self._edges[0].parent.peak.neutral_mass, self._edges[-1].child.peak.neutral_mass

    def reverse(self):
        return self.__class__(reversed(self._edges))

    def isin(self, collection):
        return self in collection

    def has_source(self):
        if isinstance(self._edges[0], SourceSpectrumEdge):
            return self._edges[0].series
        elif isinstance(self._edges[-1], SourceSpectrumEdge):
            return self._edges[-1].series
        return None

    @property
    def last_peak(self):
        return self._edges[-1].child

    @property
    def first_peak(self):
        return self._edges[0].parent

    def search_from_end(self, blocks, tolerance=2e-5):
        self.graph.search_from(self.last_peak, blocks=blocks, tolerance=tolerance)

    def search_from_start(self, blocks, tolerance=2e-5):
        self.graph.search_below(self.first_peak, blocks=blocks, tolerance=tolerance)

    def as_peptide(self):
        if self._sequence is None:
            self._sequence = Sequence(str(self))
        return self._sequence

    def sequence_mass(self, n_term=Composition("H"), c_term=Composition("OH")):
        if self._mass is None:
            self._mass = sum(a.annotation.neutral_mass for a in self) + n_term.mass + c_term.mass
        return self._mass
