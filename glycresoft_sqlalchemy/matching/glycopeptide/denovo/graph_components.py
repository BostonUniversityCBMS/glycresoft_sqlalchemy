from itertools import product, chain, combinations

from glycresoft_sqlalchemy.utils.collectiontools import groupby
from glycresoft_sqlalchemy.structure import sequence, composition, modification, fragment

from glycresoft_sqlalchemy.data_model import Peak
from glycresoft_sqlalchemy.utils.common_math import DPeak, ppm_error

NullPeak = DPeak(Peak(id=-3, neutral_mass=0, charge=0, intensity=0, scan_peak_index=-3))

Sequence = sequence.Sequence
Residue = sequence.Residue
Composition = composition.Composition
Modification = modification.Modification


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


NullNode = SpectrumNode(NullPeak)


class SpectrumEdge(object):
    series = None

    def __init__(self, parent, child, mass_shift, annotation, modification=None):
        self.parent = parent
        self.child = child
        self._hash = hash(annotation)
        self.mass_shift = mass_shift
        self.annotation = annotation
        parent.out_edges.add(self)
        child.in_edges.add(self)
        self.modification = modification

    def _reset(self):
        self._hash = hash(self.annotation)

    def __repr__(self):
        if self.modification is not None:
            return "SpectrumEdge({}, {}, {}, {}|{})".format(
                self.parent, self.child, self.mass_shift, self.annotation,
                self.modification)
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
        return self._hash

    def components(self):
        annote = self.annotation.clone()
        for ordering in annote.orderings():
            yield ordering

    def mass(self, include_modification=True):
        mass = self.annotation.neutral_mass
        if self.modification is not None and include_modification:
            mass += self.modification.mass
        return mass


class SourceSpectrumEdge(SpectrumEdge):
    def __init__(self, parent, child, mass_shift, annotation, series, modification=None):
        super(SourceSpectrumEdge, self).__init__(parent, child, mass_shift, annotation, modification)
        self.series = series

    def __repr__(self):
        if self.modification is not None:
            return "SourceSpectrumEdge({}, {}, {}, {}|{}|{})".format(
                self.parent, self.child, self.mass_shift, self.annotation, self.series,
                self.modification)
        return "SourceSpectrumEdge({}, {}, {}, {}|{})".format(
            self.parent, self.child, self.mass_shift, self.annotation, self.series)


class ModifiedSourceSpectrumEdge(SourceSpectrumEdge):
    def __init__(self, parent, child, mass_shift, annotation, series, modification=None):
        super(ModifiedSourceSpectrumEdge, self).__init__(parent, child, mass_shift, annotation, series)
        self.modification = modification

    def __repr__(self):
        return "ModifiedSourceSpectrumEdge({}, {}, {}, {}|{}|{})".format(
            self.parent, self.child, self.mass_shift, self.annotation, self.series,
            self.modification)

    def components(self):
        annote = self.annotation.clone()
        # annote = annote.add_modification(self.modification)
        for ordering in annote.orderings():
            yield ordering


class UnlinkedSpectrumEdge(SpectrumEdge):
    def __init__(self, parent, child, annotation, modification=None):
        self.parent = parent
        self.child = child
        self._hash = hash(annotation)
        self.annotation = annotation
        self.mass_shift = annotation.neutral_mass
        self.modification = modification

    def __repr__(self):
        if self.modification is not None:
            return "UnlinkedSpectrumEdge({}, {}, {}, {}|{})".format(
                self.parent, self.child, self.mass_shift, self.annotation,
                self.modification)
        return "UnlinkedSpectrumEdge({}, {}, {}, {})".format(
            self.parent, self.child, self.mass_shift, self.annotation)


class DeNovoSequence(object):

    def __init__(self, edges, graph, glycoforms=None):
        if glycoforms is None:
            glycoforms = []
        self._edges = tuple(edges)
        self.sequence = tuple(e.annotation for e in self._edges)
        self._string = None
        self._mass = None
        self._sequence = None
        self.graph = graph
        self.glycoforms = glycoforms

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
        return self.__class__(reversed(self._edges), self.graph)

    def isin(self, collection):
        return self in collection

    def has_source(self):
        series = self._edges[0].series
        if series is None:
            series = self._edges[-1].series
        return series
        # if isinstance(self._edges[0], SourceSpectrumEdge):
        #     return self._edges[0].series
        # elif isinstance(self._edges[-1], SourceSpectrumEdge):
        #     return self._edges[-1].series
        # return None

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

    def sequence_mass(self):
        n_term, c_term = self.terminal_modifications()
        if self._mass is None:
            self._mass = sum(a.annotation.neutral_mass for a in self) + n_term.mass + c_term.mass
        return self._mass

    @property
    def annotations(self):
        return [e.annotation for e in self]

    def terminal_modifications(self, default_n_term=Modification("H"), default_c_term=Modification("OH")):
        series = self.has_source()
        if series.direction > 0:
            n_term = self[0].modification or default_n_term
            c_term = self[-1].modification or default_c_term
        else:
            n_term = self[-1].modification or default_n_term
            c_term = self[0].modification or default_c_term
        return n_term, c_term

    def possible_peptides(self, as_sequence=False, include_glycan=False):
        n_term, c_term = self.terminal_modifications()

        converter = Sequence.from_iterable if as_sequence else tuple

        series = self.has_source()
        possibilities = []
        if series.direction > 0 or series is None:
            possibilities = [
                converter(s) for s in map(
                    flatlist, product(*[list(e.components()) for e in self]))]
        else:
            possibilities = [
                converter(s) for s in map(
                    flatlist, product(*[list(e.components()) for e in self[::-1]]))]

        if as_sequence:
            for peptide in possibilities:
                peptide.n_term = n_term
                peptide.c_term = c_term
            if include_glycan:
                glycosylated = []
                for peptide in possibilities:
                    for site in combinations(peptide.n_glycan_sequon_sites, 1):
                        for glycan in self.glycoforms:
                            glycoform = peptide.clone()
                            glycoform.add_modification(site[0], Modification("HexNAc"))
                            glycoform.glycan = glycan.as_composition()
                            glycosylated.append(glycoform)
                possibilities = glycosylated

        return possibilities

    def append(self, edge):
        return self.__class__(self._edges + (edge,), self.graph)

    def _append(self, annotation, modification=None):
        return self.__class__(self._edges + (UnlinkedSpectrumEdge(
            self.last_peak, NullNode, annotation, modification),), self.graph)

    def is_valid_n_glycopeptide(self):
        for possible in self.possible_peptides(as_sequence=True):
            if len(possible.n_glycan_sequon_sites):
                return True
        return False

    def search_for_opposite_ladder(self, graph=None, tolerance=2e-5):
        if graph is None:
            graph = self.graph
        series = self.has_source()
        if series == "b":
            opp_series = fragment.IonSeries.y
        elif series == 'y':
            opp_series = fragment.IonSeries.b
        else:
            raise Exception("Unsupported Ion Serise %r" % series)

        used_nodes = set()
        used_nodes.add(self.first_peak.peak.scan_peak_index)
        for edge in self:
            used_nodes.add(edge.child.peak.scan_peak_index)

        solutions = {}
        for peptide in self.possible_peptides(as_sequence=True, include_glycan=True):
            fragment_matches = {}

            fragments = [y for x in peptide.get_fragments(opp_series) for y in x]
            fragments.sort(key=lambda x: x.mass)
            i = 0
            for frag in fragments:
                while i < len(graph.peaks):
                    peak = graph.peaks[i]
                    if abs(ppm_error(frag.mass, peak.neutral_mass)) < tolerance:
                        if peak.scan_peak_index in used_nodes:
                            i += 1
                            continue
                        fragment_matches[frag.name] = peak
                        break
                    elif peak.neutral_mass > frag.mass:
                        i -= 1
                        break
                    i += 1
            total_energy = sum(p.intensity for p in fragment_matches.values())
            solutions[str(peptide)] = DeNovoSolution(peptide, self, {
                "total_energy": total_energy,
                "fragment_matches": fragment_matches
            })

        return solutions.values()


class DeNovoSolution(object):
    def __init__(self, peptide, path, reverse):
        self.peptide = peptide
        self.path = path
        self.reverse = reverse
        self.total_energy = path.total_energy() + reverse['total_energy']

    def __repr__(self):
        return "%s:%0.3e" % (self.peptide, self.total_energy)


def partition_sequences(denovo_sequences):
    sourced_sequences = []
    internal_segments = []

    for seq in denovo_sequences:
        if seq.has_source():
            sourced_sequences.append(seq)
        else:
            internal_segments.append(seq)

    return sourced_sequences, internal_segments


def unique_sequences(denovo_sequences):
    denovo_sequences = sorted(denovo_sequences, key=DeNovoSequence.total_energy, reverse=True)
    return [v[0] for v in groupby(denovo_sequences, str).values()]


def flatlist(iterable):
    result = []
    for i in iterable:
        if isinstance(i, (list, tuple)):
            result.extend(i)
        else:
            result.append(i)
    return result


def possible_peptides(annotations, ion_series, n_term=Modification("H"), c_term=Modification("OH")):
    annotations = annotations[::ion_series.direction]
    possibilities = [
        (s) for s in map(
            flatlist, product(*[list(a.orderings()) for a in annotations]))]
    return possibilities


def is_valid_n_glycopeptide(annotations, ion_series):
        for possible in possible_peptides(annotations, ion_series):
            if find_n_glycosylation_sequons(possible, True):
                return True
        return False


__asn = Residue("Asn")
__pro = Residue("Pro")
__ser = Residue("Ser")
__thr = Residue("Thr")
__hexnac = Modification("HexNAc")


def find_n_glycosylation_sequons(sequence, justcheck=False):
    hexnac = __hexnac
    state = "seek"  # [seek, n, ^p, st]
    asn = __asn
    pro = __pro
    ser = __ser
    thr = __thr
    i = 0
    positions = []
    n_pos = None
    while(i < len(sequence)):
        next_pos = tuple(sequence[i])
        if state == "seek":
            # A sequon starts with an Asn residue without modifications, or for counting
            # purposes one that has already been glycosylated
            if next_pos[0] is asn:
                if ((len(next_pos[1]) == 0) or (next_pos[1][0].rule is hexnac.rule)):
                    n_pos = i
                    state = "n"
        elif state == "n":
            if next_pos[0] is not pro:
                state = "^p"
            else:
                state = "seek"
                i = n_pos
                n_pos = None
        elif state == "^p":
            if next_pos[0] is ser or next_pos is thr:
                positions.append(n_pos)
                if justcheck:
                    return positions
            i = n_pos
            n_pos = None
            state = "seek"
        i += 1
    return positions
