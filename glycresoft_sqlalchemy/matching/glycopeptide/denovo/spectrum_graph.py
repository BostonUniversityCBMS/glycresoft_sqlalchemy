import operator
import itertools
from collections import defaultdict

import numpy as np

from glycresoft_sqlalchemy.data_model import Peak, TheoreticalGlycanCombination, AminoAcidComposition
from glycresoft_sqlalchemy.utils.common_math import ppm_error, DPeak
from glycresoft_sqlalchemy.utils.tree import SuffixTree
from glycresoft_sqlalchemy.utils.collectiontools import groupby

from glycresoft_sqlalchemy.structure import fragment


from .graph_components import (
    SpectrumNode, SpectrumEdge, SourceSpectrumEdge,
    DeNovoSequence, partition_sequences, is_valid_n_glycopeptide)

from .building_blocks import DegenerateSegmentBlock


IonSeries = fragment.IonSeries
neutral_mass_getter = operator.attrgetter('neutral_mass')
intensity_getter = operator.attrgetter("intensity")
default_match_tolerance = 2e-5


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
            elif error > tolerance:
                break
        return solutions


class GlycanCombinationResolverQuerySearcher(GlycanCombinationResolverBase):
    def __init__(self, total_mass, session, hypothesis_id,
                 glycan_combination_type=TheoreticalGlycanCombination, **kwargs):
        self.hypothesis_id = hypothesis_id
        self.session = session
        self.glycan_combination_type = glycan_combination_type
        self.glycan_combinations = session.query(glycan_combination_type).filter(
            glycan_combination_type.hypothesis_id == self.hypothesis_id)
        super(GlycanCombinationResolverQuerySearcher, self).__init__(total_mass)

    def possible_peptide_masses(self, total_mass=None):
        if total_mass is None:
            total_mass = self.total_mass

        expr = total_mass - self.glycan_combination_type.dehydrated_mass()
        expr.type.asdecimal = False

        possible_solutions = self.session.query(
            expr, self.glycan_combination_type).filter(
            self.glycan_combination_type.hypothesis_id == self.hypothesis_id)
        return possible_solutions

    def solutions(self, peptide=None, mass=None, tolerance=2e-5):
        if peptide is not None:
            peptide_mass = peptide.sequence_mass()
        elif mass is not None:
            peptide_mass = mass
        delta_mass = self.total_mass - peptide_mass

        solutions = self.glycan_combination_type.ppm_error_tolerance_search(
            self.session, delta_mass, tolerance, self.hypothesis_id, dehydrate=True).all()

        return solutions


class AminoAcidCompositionResolverBase(object):
    def __init__(self, enzyme_termini, **kwargs):
        self.enzyme_termini = enzyme_termini

    def has_cleavage_site(self, sequence):
        for terminal in self.enzyme_termini:
            if terminal in sequence:
                return True
        return False

    def solutions(self, mass, tolerance=2e-5):
        raise NotImplementedError()


class AminoAcidCompositionResolverLinearSearcher(AminoAcidCompositionResolverBase):
    def __init__(self, enzyme_termini, amino_acid_compositions, **kwargs):
        super(AminoAcidCompositionResolverLinearSearcher, self).__init__(enzyme_termini, **kwargs)
        self.amino_acid_compositions = amino_acid_compositions
        self.amino_acid_compositions.sort(key=lambda x: x.neutral_mass)

    def solutions(self, mass, tolerance=2e-5):
        solutions = []
        for composition in self.amino_acid_compositions:
            error = ppm_error(mass, composition.neutral_mass)
            if abs(error) < tolerance:
                solutions.append(composition)
            if composition.neutral_mass > mass:
                break
        return solutions


class AminoAcidCompositionResolverQuerySearcher(AminoAcidCompositionResolverBase):
    def __init__(self, enzyme_termini, session, hypothesis_id,
                 amino_acid_composition_type=AminoAcidComposition,
                 **kwargs):
        super(AminoAcidCompositionResolverQuerySearcher, self).__init__(enzyme_termini)

        self.session = session
        self.amino_acid_composition_type = amino_acid_composition_type
        self.hypothesis_id = hypothesis_id

    def solutions(self, mass, tolerance=2e-5):
        solutions = self.amino_acid_composition_type.ppm_error_tolerance_search(
            self.session, mass, tolerance  # , self.hypothesis_id
            )
        solutions = [DegenerateSegmentBlock(a.to_sequence_composition()) for a in solutions]
        return solutions


class SpectrumGraph(object):
    def __init__(self, peaks, precursor_mass, glycosylation_resolver):
        self.peaks = sorted(peaks, key=neutral_mass_getter)
        self.source_node = SpectrumNode(DPeak(Peak(
            id=-1, neutral_mass=0, charge=0, intensity=0, scan_peak_index=-1)))
        self.nodes = [SpectrumNode(p) for p in self.peaks]
        self.edges = defaultdict(list)
        self.precursor_mass = precursor_mass
        self.sink_node = SpectrumNode(DPeak(Peak(
            id=-2, neutral_mass=precursor_mass, charge=0, intensity=0, scan_peak_index=-2)))
        self.glycosylation_resolver = glycosylation_resolver
        self.solutions = {}

    def iteredges(self):
        for nodes, edges in self.edges.items():
            for edge in edges:
                yield edge

    def find_starting_points(self, blocks, offset, terminal_modifications=None, tolerance=2e-5):
        if isinstance(offset, IonSeries):
            offset_mass = offset.mass_shift
        else:
            offset_mass = offset
        if terminal_modifications is None:
            terminal_modifications = []
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

        for modification in terminal_modifications:
            a_mass = offset_mass + modification.mass
            a_peak = self.source_node
            for block in blocks:
                a_mass += block.neutral_mass
                for b_peak in self.nodes:
                    if abs(ppm_error(a_mass, b_peak.neutral_mass)) <= tolerance:
                        self.edges[a_peak, b_peak].append(SourceSpectrumEdge(
                            a_peak, b_peak, block.neutral_mass,
                            annotation=block, series=offset, modification=modification))
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

    def find_sinks(self, sequences, blocks, terminal_modifications=None, tolerance=2e-5):
        # b_peak = self.sink_node
        if terminal_modifications is None:
            terminal_modifications = []
        had_solutions = list()
        for seq in sequences:
            # a_peak = seq.last_peak
            a_mass = seq.sequence_mass()
            series = seq.has_source()
            is_n_glycopeptide = is_valid_n_glycopeptide(seq.annotations, series)
            if is_n_glycopeptide:
                solutions = self.glycosylation_resolver.solutions(mass=a_mass)
                if len(solutions) > 0:
                    had_solutions.append(seq)
            for block in blocks:
                extended = seq._append(block)
                if not is_n_glycopeptide and not is_valid_n_glycopeptide(extended.annotations, series):
                    continue
                a_mass += block.neutral_mass

                solutions = self.glycosylation_resolver.solutions(mass=a_mass)
                if len(solutions) > 0:
                    extended.glycoforms.extend(solutions)
                    had_solutions.append(extended)

                a_mass -= block.neutral_mass
            for mod in terminal_modifications:
                a_mass += mod.mass
                for block in blocks:
                    extended = seq._append(block, mod)
                    if not is_n_glycopeptide and not is_valid_n_glycopeptide(extended.annotations, series):
                        continue
                    a_mass += block.neutral_mass

                    solutions = self.glycosylation_resolver.solutions(mass=a_mass)
                    if len(solutions) > 0:
                        extended.glycoforms.extend(solutions)
                        had_solutions.append(extended)

                    a_mass -= block.neutral_mass
                a_mass -= mod.mass
        return had_solutions

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

    def run(self, single_blocks, blocks_list=None, ion_series=(IonSeries.b, IonSeries.y), tolerance=2e-5):
        if blocks_list is None:
            blocks_list = []
        for offset in ion_series:
            self.find_starting_points(single_blocks, offset, tolerance)
            for blocks_opt in blocks_list:
                self.find_starting_points(blocks_opt, offset, tolerance)
        self.find_edges(single_blocks, tolerance)
        for blocks_opt in blocks_list:
            self.find_edges(blocks_opt, tolerance)
        contigs = self.generate_contiguous()

        self.join_internal_segments(contigs, single_blocks, tolerance)
        for blocks_opt in blocks_list:
            self.join_internal_segments(contigs, blocks_opt, tolerance)

        contigs = self.generate_contiguous()
        by_source = groupby(contigs, lambda x: x.has_source())
        by_source.pop(None, None)  # remove those contigs without a source

        forward_sequences = []
        reverse_sequences = []
        for k, v in by_source.items():
            if k.direction > 0:
                forward_sequences.extend(v)
            else:
                reverse_sequences.extend(v)

        forward_solutions = []
        reverse_solutions = []
        for blocks in blocks_list:
            forward_solutions.extend(self.find_sinks(forward_sequences, blocks, tolerance=tolerance))
            reverse_solutions.extend(self.find_sinks(forward_sequences, blocks, tolerance=tolerance))

        solutions = set(itertools.chain.from_iterable(
            [sol.search_for_opposite_ladder(tolerance=tolerance) for sol in
             itertools.chain(forward_solutions, reverse_solutions)]))
        return solutions


class QuerySpectrumGraph(object):
    def __init__(self, peaks, precursor_mass, glycosylation_resolver, amino_acid_resolver):
        self.peaks = sorted(peaks, key=neutral_mass_getter)
        self.source_node = SpectrumNode(DPeak(Peak(
            id=-1, neutral_mass=0, charge=0, intensity=0, scan_peak_index=-1)))
        self.nodes = [SpectrumNode(p) for p in self.peaks]
        self.edges = defaultdict(list)
        self.precursor_mass = precursor_mass
        self.sink_node = SpectrumNode(DPeak(Peak(
            id=-2, neutral_mass=precursor_mass, charge=0, intensity=0, scan_peak_index=-2)))
        self.glycosylation_resolver = glycosylation_resolver
        self.amino_acid_resolver = amino_acid_resolver
        self.solutions = {}

    def find_starting_points(self, ion_series, modifications=None, tolerance=2e-5):
        if modifications is None:
            modifications = tuple()
        a_peak = self.source_node

        for b_peak in self.nodes:
            delta_mass = b_peak.neutral_mass + ion_series
            for solution in self.amino_acid_resolver.solutions(delta_mass, tolerance=tolerance):
                self.edges[a_peak, b_peak].append(
                    SourceSpectrumEdge(
                        a_peak, b_peak, delta_mass, annotation=solution, series=ion_series))
            for modification in modifications:
                delta_mass_mod = delta_mass + modification.mass
                for solution in self.amino_acid_resolver.solutions(delta_mass_mod, tolerance=tolerance):
                    self.edges[a_peak, b_peak].append(
                        SourceSpectrumEdge(
                            a_peak, b_peak, delta_mass, annotation=solution, series=ion_series,
                            modification=modification))

    def find_edges(self, tolerance=2e-5):
        for i, a_peak in enumerate(self.nodes):
            a_mass = a_peak.neutral_mass
            for b_peak in self.nodes[i:]:
                delta_mass = b_peak.neutral_mass - a_mass
                for solution in self.amino_acid_resolver.solutions(delta_mass, tolerance=tolerance):
                    self.edges[a_peak, b_peak].append(
                        SpectrumEdge(
                            a_peak, b_peak, delta_mass, annotation=solution))

    def _has_solution(self, mass, tolerance=2e-5):
        for solution in self.glycosylation_resolver.solutions(mass=mass, tolerance=tolerance):
            yield solution

    def _link(self, node, current_sequence=None, glycosylation_sites=0,
              cleavage_sites=0, current_mass=0, ion_series=None,
              glycosite_resolver=(lambda sequence, series: False)):
        if current_sequence is None:
            current_sequence = []
        streams = []
        if len(node.out_edges) > 0:
            for edge in node.out_edges:
                new_mass = edge.mass() + current_mass
                seq = list(current_sequence)
                if len(seq) == 1:
                    try:
                        ion_series = edge.series
                    except AttributeError:
                        ion_series = None
                _cleavage_count = cleavage_sites
                _glycosylation_sites = glycosylation_sites + glycosite_resolver(seq, ion_series)
                has_cleavage_site = self.amino_acid_resolver.has_cleavage_site(edge.annotation)
                if has_cleavage_site:
                    _cleavage_count + 1
                seq.append(edge)
                if has_cleavage_site and _glycosylation_sites > 0:
                    for solution in self._has_solution(new_mass):
                        self.solutions.add((seq, solution))
                streams.append(self._link(
                    edge.child, seq, _glycosylation_sites, _cleavage_count, new_mass,
                    ion_series, glycosite_resolver))

            for stream in streams:
                for seq in stream:
                    yield seq
        else:
            yield DeNovoSequence(current_sequence)
