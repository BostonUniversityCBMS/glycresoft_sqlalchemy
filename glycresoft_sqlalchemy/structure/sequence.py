import re
import copy
import itertools
import operator
from collections import defaultdict, deque, Counter, namedtuple

from . import PeptideSequenceBase, MoleculeBase
from . import constants as structure_constants
from .composition import Composition
from .fragment import PeptideFragment, fragment_shift, fragment_shift_composition, SimpleFragment, IonSeries
from .modification import Modification, SequenceLocation, ModificationCategory
from .residue import Residue
from glypy import GlycanComposition, Glycan, MonosaccharideResidue, Substituent, ReducedEnd
from glypy.composition.glycan_composition import FrozenGlycanComposition, FrozenMonosaccharideResidue

from .parser import sequence_tokenizer, sequence_length, strip_modifications
from .glycan import GlycosylationType, GlycosylationManager, GlycosylationSite, glycosylation_site_detectors

from ..utils.iterators import peekable
from ..utils.memoize import memoize
from ..utils.collectiontools import descending_combination_counter


b_series = IonSeries.b
y_series = IonSeries.y
oxonium_ion_series = IonSeries.oxonium_ion
stub_glycopeptide_series = IonSeries.stub_glycopeptide


def list_to_sequence(seq_list, wrap=True):
    flat_chunks = []
    for chunk in seq_list:
        if(isinstance(chunk[0], list)):
            flat_chunks.extend(chunk)
        else:
            flat_chunks.append(chunk)
    seq = Sequence.from_iterable(flat_chunks) if wrap else flat_chunks
    return seq


@memoize()
def sequence_to_mass(sequence):
    mass = 0.0
    chunks, modifications, glycan, n_term, c_term = sequence_tokenizer(sequence)
    for residue, mods in chunks:
        mass += Residue.mass_by_name(residue)
        for mod in mods:
            mass += Modification.mass_by_name(mod)
    if n_term is not None:
        mass += Modification.mass_by_name(n_term)
    else:
        mass += Composition("H").mass
    if c_term is not None:
        mass += Modification.mass_by_name(c_term)
    else:
        mass += Composition("OH").mass
    return mass


@glycosylation_site_detectors(GlycosylationType.o_linked)
def find_o_glycosylation_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification("HexNAc")}
    except:
        allow_modified = (Modification("HexNAc"),)

    positions = []
    ser = Residue("S")
    thr = Residue("T")

    site_set = (ser, thr)

    if isinstance(sequence, basestring):
        sequence = parse(sequence)

    for i, position in enumerate(sequence):
        if position[0] in site_set:
            if ((len(position[1]) == 0) or position[1][0] in allow_modified):
                positions.append(i)
    return positions


@glycosylation_site_detectors(GlycosylationType.n_linked)
def find_n_glycosylation_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification("HexNAc")}
    except:
        allow_modified = (Modification("HexNAc"),)
    state = "seek"  # [seek, n, ^p, st]
    if isinstance(sequence, basestring):
        sequence = PeptideSequence(sequence)

    asn = Residue("Asn")
    pro = Residue("Pro")
    ser = Residue("Ser")
    thr = Residue("Thr")

    i = 0
    positions = []
    n_pos = None
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            # A sequon starts with an Asn residue without modifications, or for counting
            # purposes one that has already been glycosylated
            if next_pos[0] == asn:
                if ((len(next_pos[1]) == 0) or next_pos[1][0] in allow_modified):
                    n_pos = i
                    state = "n"
        elif state == "n":
            if next_pos[0] != pro:
                state = "^p"
            else:
                state = "seek"
                i = n_pos
                n_pos = None
        elif state == "^p":
            if next_pos[0] in {ser, thr}:
                positions.append(n_pos)
            i = n_pos
            n_pos = None
            state = "seek"
        i += 1
    return positions


@glycosylation_site_detectors(GlycosylationType.glycosaminoglycan)
def find_glycosaminoglycan_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification("HexNAc")}
    except:
        allow_modified = (Modification("HexNAc"),)
    state = "seek"  # [seek, s, g1, x, g2]
    ser = Residue("Ser")
    gly = Residue("Gly")

    i = 0
    positions = []

    s_position = None
    if isinstance(sequence, basestring):
        sequence = PeptideSequence(sequence)
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            # A sequon starts with an Asn residue without modifications, or for counting
            # purposes one that has already been glycosylated
            if next_pos[0] == ser:
                if ((len(next_pos[1]) == 0) or next_pos[1][0] in allow_modified):
                    s_position = i
                    state = "s"
        elif state == "s":
            if next_pos[0] == gly:
                state = 'g1'
            else:
                state = "seek"
                i = s_position
                s_position = None
        elif state == 'g1':
            # Anything will satisfy this position pattern
            state = 'x'
        elif state == "x":
            if next_pos[0] == gly:
                positions.append(s_position)
            state = "seek"
            i = s_position
            s_position = None
        i += 1
    return positions


def total_composition(sequence):
    if isinstance(sequence, basestring):
        sequence = parse(sequence)
    return _total_composition(sequence)


def _total_composition(sequence):
    glycan = sequence.glycan
    total = Composition()
    if glycan is not None:
        for position in sequence:
            total += position[0].composition
            for mod in position[1]:
                if mod.name == "HexNAc":
                    continue
                total += mod.composition
        total += sequence.n_term.composition
        total += sequence.c_term.composition
        total += glycan.total_composition()
    else:
        for position in sequence:
            total += position[0].composition
            for mod in position[1]:
                total += mod.composition
        total += sequence.n_term.composition
        total += sequence.c_term.composition

    return total


def _calculate_mass(sequence):
    glycan = sequence.glycan
    total = 0
    if glycan is not None:
        for position in sequence:
            total += position[0].mass
            for mod in position[1]:
                if mod.name == "HexNAc":
                    continue
                total += mod.mass
        total += sequence.n_term.mass
        total += sequence.c_term.mass
        total += glycan.mass()
    else:
        for position in sequence:
            total += position[0].mass
            for mod in position[1]:
                total += mod.mass
        total += sequence.n_term.mass
        total += sequence.c_term.mass

    return total


class PeptideSequence(PeptideSequenceBase):
    '''
    Represents a peptide that may have post-translational modifications
    including glycosylation.

    Attributes
    ----------
    seq: list
        The underlying container for positions in the amino acid sequence
    modification_index: defaultdict(int)
        A count of different modifications attached to the amino acid sequence
    n_term: Modification or Composition
    c_term: Modification or Composition
        Terminal modifications (N-terminus and C-terminus respectively) which
        default to H and OH respectively.
    glycan: Glycan or GlycanComposition
        The total glycan moiety attached to the molecule. The current semantics
        do not cleanly support more than one glycosylation per sequence for generating
        glycan sequence fragments.
    mass: float
        The pre-calculated monoisotopic mass of the molecule. This quantity is
        assumes that the glycan's glycosidic bonds have been broken, leaving only
        the amide-bound HexNAc as a modification attached to the amino acid backbone
    '''
    position_class = list

    @classmethod
    def from_iterable(cls, iterable):
        seq = cls("")
        n_term = structure_constants.N_TERM_DEFAULT
        c_term = structure_constants.C_TERM_DEFAULT
        i = 0
        for pos, next_pos in peekable(iterable):
            i += 1
            try:
                resid, mods = pos
            except ValueError:
                if i == 0:
                    n_term = pos
                elif next_pos is peekable.sentinel:
                    c_term = pos
                else:
                    raise
            if not isinstance(resid, Residue):
                resid = Residue(symbol=resid)
            seq.mass += resid.mass
            mod_list = []
            for mod in mods:
                if mod == "":
                    continue
                if not isinstance(mod, Modification):
                    mod = Modification(mod)
                mod_list.append(mod)
                seq.mass += mod.mass
                seq.modification_index[mod.name] += 1
            seq.sequence.append(cls.position_class([resid, mod_list]))
        if not isinstance(n_term, MoleculeBase):
            n_term = Modification(n_term)
        if not isinstance(c_term, MoleculeBase):
            c_term = Modification(c_term)

        seq.n_term = n_term
        seq.c_term = c_term
        return seq

    def __init__(self, sequence=None, **kwargs):
        self.mass = 0.0
        self.sequence = []
        self.modification_index = defaultdict(int)

        self._glycan = None

        self._n_term = None
        self._c_term = None
        if sequence == "" or sequence is None:
            pass
        else:
            seq_list, modifications, glycan, n_term, c_term = sequence_tokenizer(sequence)
            for item in seq_list:
                res = Residue(item[0])
                self.mass += res.mass
                mods = []
                for mod in item[1]:
                    if mod != '':
                        mod = Modification(mod)
                        mods.append(mod)
                        self.modification_index[mod.name] += 1
                        self.mass += mod.mass
                self.sequence.append(self.position_class([res, mods]))

            self.glycan = glycan if glycan != "" else None
            self.n_term = Modification(n_term) if isinstance(n_term, basestring) else n_term
            self.c_term = Modification(c_term) if isinstance(c_term, basestring) else c_term
        self._fragments_map = {}

    def _patch_glycan_composition(self):
        occupied_sites = self.modification_index['HexNAc']
        offset = Composition({"H": 2, "O": 1}) * occupied_sites
        self.glycan.composition_offset -= offset

    def __repr__(self):
        n_term = ""
        if self.n_term is not None:
            n_term = "({0})-".format(self.n_term)
        c_term = ""
        if self.c_term is not None:
            c_term = "-({0})".format(self.c_term)
        rep = "{n_term}{sequence}{c_term}{glycan}[{mass}]".format(
            n_term=n_term, c_term=c_term,
            glycan=self._glycan if self._glycan is not None else "",
            **self.__dict__)
        return rep

    def __len__(self):
        return len(self.sequence)

    @property
    def total_mass(self):
        return self.total_composition().mass

    @property
    def glycan(self):
        return self._glycan

    @glycan.setter
    def glycan(self, value):
        self._glycan = value
        if isinstance(value, GlycanComposition):
            self._patch_glycan_composition()
        elif isinstance(value, Glycan):
            value.reducing_end = ReducedEnd(-Composition("H2O"), valence=0)

    @property
    def n_term(self):
        return self._n_term

    @n_term.setter
    def n_term(self, value):
        reset_mass = 0
        try:
            reset_mass = self._n_term.mass
        except:
            pass
        self._n_term = value
        new_mass = 0
        try:
            new_mass = value.mass
        except:
            pass
        self.mass += new_mass - reset_mass

    @property
    def c_term(self):
        return self._c_term

    @c_term.setter
    def c_term(self, value):
        reset_mass = 0
        try:
            reset_mass = self._c_term.mass
        except:
            pass
        self._c_term = value
        new_mass = 0
        try:
            new_mass = value.mass
        except:
            pass
        self.mass += new_mass - reset_mass

    def __iter__(self):
        '''Makes the sequence iterable'''
        for i in self.sequence:
            yield (i)

    def __getitem__(self, index):
        sub = self.sequence[index]
        return sub

    def __setitem__(self, index, value):
        self.sequence[index] = value

    def subseq(self, slice_obj):
        sub = self[slice_obj]
        subseq = Sequence.from_iterable(sub)
        subseq.n_term = self.n_term
        return subseq

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return str(self) != str(other)

    def deglycosylate(self):
        _glycosylation_enum = ModificationCategory.glycosylation
        for i, pos in enumerate(self):
            mods = [mod.name for mod in pos[1] if mod.is_a(
                _glycosylation_enum)]
            for mod in mods:
                self.drop_modification(i, mod)

    def break_at(self, idx, neutral_losses=None):
        b_shift = fragment_shift['b']
        if self.n_term != structure_constants.N_TERM_DEFAULT:
            b_shift = b_shift + self.n_term.mass
        y_shift = fragment_shift['y']
        if self.c_term != structure_constants.C_TERM_DEFAULT:
            y_shift = y_shift + self.c_term.mass

        mod_b = defaultdict(int)
        mass_b = 0

        mod_y = defaultdict(int)
        mass_y = 0

        pos = 0
        residues_in_b = []
        for pos in range(idx):
            for mod in self.sequence[pos][1]:
                mod_b[mod.serialize()] += 1
            residues_in_b.append(self.sequence[pos][0].symbol)
            mass_b += self.sequence[pos][0].mass

        flanking_residues = [self.sequence[pos][0], self.sequence[pos + 1][0]]
        b_frag = PeptideFragment(
            b_series, pos + structure_constants.FRAG_OFFSET, mod_b, mass_b + b_shift,
            flanking_amino_acids=flanking_residues)

        break_point = pos + 1
        residues_in_y = []
        for pos in range(break_point, len(self)):
            for mod in self.sequence[pos][1]:
                mod_y[mod.serialize()] += 1
            residues_in_y.append(self.sequence[pos][0].symbol)
            mass_y += self.sequence[pos][0].mass

        y_frag = PeptideFragment(
            y_series, len(self) - (break_point - 1 + structure_constants.FRAG_OFFSET),
            mod_y, mass_y + y_shift, flanking_amino_acids=flanking_residues)
        if structure_constants.PARTIAL_HEXNAC_LOSS:
            b_frag.golden_pairs = [frag.name for frag in y_frag.partial_loss()]
            y_frag.golden_pairs = [frag.name for frag in b_frag.partial_loss()]
            b_frag = list(b_frag.partial_loss())
            y_frag = list(y_frag.partial_loss())
        else:
            b_frag.golden_pairs = [y_frag.name]
            y_frag.golden_pairs = [b_frag.name]

        if neutral_losses is not None:
            b_frag = list(b_frag.generate_neutral_losses(neutral_losses))
            y_frag = list(y_frag.generate_neutral_losses(neutral_losses))

        return b_frag, y_frag

    def get_fragments(self, kind, neutral_losses=None, **kwargs):
        """Return a list of mass values for each fragment of `kind`"""

        mass_shift = 0.0

        # The set of modification names.
        mod_dict = {}

        # The key is the position, the value is an array of fragments.
        # And the first element is always bare fragment.
        # The total number of HexNAc on the fragment should be recorded.
        kind = IonSeries(kind)
        running_composition = Composition()

        if kind in (b_series, "a", "c"):
            seq_list = self.sequence

            mass_shift = fragment_shift[kind]
            running_composition += fragment_shift_composition[kind]
            # if self.n_term != structure_constants.N_TERM_DEFAULT:
                # mass_shift = mass_shift + self.n_term.mass
            mass_shift = mass_shift + self.n_term.mass
            running_composition += self.n_term.composition

        elif kind in (y_series, 'x', 'z'):
            mass_shift = fragment_shift[kind]
            running_composition += fragment_shift_composition[kind]
            # if self.c_term != structure_constants.C_TERM_DEFAULT:
            #     mass_shift = mass_shift + self.c_term.mass
            mass_shift = mass_shift + self.c_term.mass
            seq_list = list(reversed(self.sequence))
            running_composition += self.c_term.composition

        else:
            raise Exception("Can't recognize ion series %r" % kind)

        def composition_of_position(position):
            residue, modifications = position[0], position[1]
            composition = Composition(residue.composition)
            for mod in modifications:
                composition += mod.composition
            return composition

        hexnac_composition = Modification("HexNAc").composition

        current_mass = mass_shift

        for idx in range(len(seq_list) - 1):
            for mod in seq_list[idx][1]:
                mod_serial = mod.serialize()
                if mod_serial in mod_dict:
                    mod_dict[mod_serial] += 1
                else:
                    mod_dict[mod_serial] = 1

            current_mass += seq_list[idx][0].mass

            running_composition += composition_of_position(seq_list[idx])

            fragments_from_site = []
            flanking_residues = [seq_list[idx][0], seq_list[idx + 1][0]]
            if kind is y_series:
                flanking_residues = flanking_residues[::-1]
            # If incremental loss of HexNAc is not allowed, only one fragment of a given type is generated
            if not structure_constants.PARTIAL_HEXNAC_LOSS:
                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, dict(mod_dict), current_mass,
                    flanking_amino_acids=flanking_residues, composition=running_composition)
                fragments_from_site.append(frag)
                bare_dict = dict(mod_dict)

                lost_composition = running_composition - mod_dict['HexNAc'] * hexnac_composition

                bare_dict["HexNAc"] = 0

                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, dict(bare_dict), current_mass,
                    flanking_amino_acids=flanking_residues, composition=lost_composition)
                fragments_from_site.append(frag)
            # Else a fragment for each incremental loss of HexNAc must be generated
            else:
                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, dict(mod_dict), current_mass,
                    flanking_amino_acids=flanking_residues, composition=running_composition)
                fragments_from_site.extend(frag.partial_loss())

            if neutral_losses is not None:
                all_frags = []
                for frag in fragments_from_site:
                    all_frags.append(frag)
                    all_frags.extend(frag.generate_neutral_losses(neutral_losses))
                fragments_from_site = all_frags

            yield fragments_from_site

    def drop_modification(self, position, modification_type):
        '''
        Drop a modification by name from a specific residue. If the
        position is the N-term or the C-term, the terminal modification will
        be reset to the default.

        Parameters
        ----------
        position: int
            The position of the modification to drop
        modification_type: str or Modification
            The modification to drop
        '''
        dropped_index = None

        if position is SequenceLocation.n_term:
            self.n_term = Modification("H")
            return
        elif position is SequenceLocation.c_term:
            self.c_term = Modification("OH")
            return

        for i, mod in enumerate(self.sequence[position][1]):
            if modification_type == mod.name:
                dropped_index = i
                break
        try:
            drop_mod = self.sequence[position][1].pop(dropped_index)
            self.mass -= drop_mod.mass
            self.modification_index[drop_mod.name] -= 1
        except:
            raise ValueError("Modification not found! %s @ %s" % (modification_type, position))

    def add_modification(self, position=None, modification_type=None):
        if position is None and isinstance(modification_type, Modification):
            position = modification_type.position

        if isinstance(modification_type, Modification):
            mod = modification_type
        else:
            mod = Modification(rule=modification_type, mod_pos=position)

        if position is SequenceLocation.n_term:
            self.n_term = mod
        elif position is SequenceLocation.c_term:
            self.c_term = mod
        else:
            if (position == -1) or (position >= len(self.sequence)):
                raise IndexError(
                    "Invalid modification position. %s, %s, %s" %
                    (position, str(self.sequence), modification_type))

            self.sequence[position][1].append(mod)
            self.mass += mod.mass
            self.modification_index[mod.name] += 1

    def fragment(self, key):
        try:
            return self._fragments_map[key]
        except KeyError:
            for group in self.get_fragments(key[0]):
                for frag in group:
                    self._fragments_map[frag.name] = frag
            return self._fragments_map[key]

    def get_sequence(self, start=0, include_glycan=True, include_termini=True,
                     implicit_n_term=None, implicit_c_term=None):
        """
        Generate human readable sequence string. Called by :meth:`__str__`

        Parameters
        ----------
        start: int
            The position to start from
        include_glycan: bool
            Whether to include the glycan in the resulting string. Defaults to `True`
        include_termini: bool
            Whether to include the N- and C-termini. Make sure this is `True` if you want non-standard
            termini to be properly propagated.


        Returns
        -------
        str
        """
        if implicit_n_term is None:
            implicit_n_term = structure_constants.N_TERM_DEFAULT
        if implicit_c_term is None:
            implicit_c_term = structure_constants.C_TERM_DEFAULT

        seq_list = []
        for x, y in self.sequence[start:]:
            mod_str = ''
            if y != []:
                mod_str = '|'.join(mod.serialize() for mod in y)
                mod_str = ''.join(['(', mod_str, ')'])
            seq_list.append(''.join([x.symbol, mod_str]))
        rep = ''.join(seq_list)
        if include_termini:
            n_term = ""
            if self.n_term is not None and self.n_term != implicit_n_term:
                n_term = "({0})-".format(self.n_term.serialize())
            c_term = ""
            if self.c_term is not None and self.c_term != implicit_c_term:
                c_term = "-({0})".format(self.c_term.serialize())
            rep = "{0}{1}{2}".format(n_term, rep, c_term)
        if include_glycan:
            if self._glycan is not None:
                rep += str(self._glycan)
        return rep

    __str__ = get_sequence

    def clone(self):
        inst = self.__class__()
        inst.n_term = self.n_term.clone()
        inst.c_term = self.c_term.clone()
        inst.sequence = [self.position_class([o[0], list(o[1])]) for o in self]
        inst.mass = self.mass
        if self.glycan is not None:
            inst.glycan = self.glycan.clone()
        return inst

    def insert(self, position, residue, modifications=None):
        if modifications is None:
            modifications = []
        self.sequence.insert(position, [residue, modifications])
        self.mass += residue.mass
        for mod in modifications:
            self.mass += mod.mass

    def delete(self, position):
        residue, mods = self.sequence.pop(position)
        self.mass -= residue.mass
        for mod in mods:
            self.mass -= mod.mass

    def append(self, residue, modification=None):
        self.mass += residue.mass
        next_pos = [residue]
        if modification is None:
            next_pos.append([])
        else:
            next_pos.append([modification])
            self.mass += modification.mass
            self.modification_index[modification.name] += 1
        self.sequence.append(self.position_class(next_pos))

    def extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence.extend(sequence.seq)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count

    def leading_extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence = sequence.seq + self.sequence
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count

    @property
    def n_glycan_sequon_sites(self):
        return find_n_glycosylation_sequons(self, structure_constants.ALLOW_MODIFIED_ASPARAGINE)

    @property
    def o_glycan_sequon_sites(self):
        return find_o_glycosylation_sequons(self)

    @property
    def glycosaminoglycan_sequon_sites(self):
        return find_glycosaminoglycan_sequons(self)

    def stub_fragments(self):
        if isinstance(self.glycan, Glycan):
            glycan = GlycanComposition.from_glycan(self.glycan)
        elif isinstance(self.glycan, GlycanComposition):
            glycan = self.glycan
        else:
            raise TypeError((
                "Cannot infer monosaccharides from non-Glycan"
                "or GlycanComposition {}").format(self.glycan))
        fucose_count = glycan['Fuc'] or glycan['dHex']
        core_count = self.modification_index['HexNAc']

        per_site_shifts = []
        hexose_mass = FrozenMonosaccharideResidue.from_iupac_lite("Hex").mass()
        hexnac_mass = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc").mass()
        fucose_mass = FrozenMonosaccharideResidue.from_iupac_lite("Fuc").mass()
        base_mass = self.mass - (hexnac_mass * core_count)
        for i in range(core_count):
            core_shifts = []
            for hexnac_count in range(3):
                if hexnac_count == 0:
                    shift = {
                        "mass": 0,
                        "key": ""
                    }
                    core_shifts.append(shift)
                elif hexnac_count == 1:
                    shift = {
                        "mass": (hexnac_count * hexnac_mass),
                        "key": {"HexNAc": hexnac_count}
                    }
                    core_shifts.append(shift)
                    if i < fucose_count:
                        fucosylated = shift.copy()
                        fucosylated['key'] = fucosylated['key'].copy()
                        fucosylated['mass'] += fucose_mass
                        fucosylated['key']["Fuc"] = 1
                        core_shifts.append(fucosylated)
                elif hexnac_count == 2:
                    shift = {
                        "mass": (hexnac_count * hexnac_mass),
                        "key": {"HexNAc": hexnac_count}
                    }
                    core_shifts.append(shift)

                    if i < fucose_count:
                        fucosylated = shift.copy()
                        fucosylated['key'] = fucosylated['key'].copy()
                        fucosylated['mass'] += fucose_mass
                        fucosylated['key']["Fuc"] = 1
                        core_shifts.append(fucosylated)

                    for hexose_count in range(1, 4):
                        shift = {
                            "mass": (hexnac_count * hexnac_mass) + (hexose_count * hexose_mass),
                            "key": {"HexNAc": hexnac_count, "Hex": hexose_count}
                        }
                        core_shifts.append(shift)
                        if i < fucose_count:
                            fucosylated = shift.copy()
                            fucosylated['key'] = fucosylated['key'].copy()
                            fucosylated['mass'] += fucose_mass
                            fucosylated['key']["Fuc"] = 1
                            core_shifts.append(fucosylated)
            per_site_shifts.append(core_shifts)
        for positions in itertools.product(*per_site_shifts):
            key_base = 'peptide'
            names = Counter()
            mass = base_mass
            for site in positions:
                mass += site['mass']
                names += Counter(site['key'])
            extended_key = ''.join("%s%d" % kv for kv in names.items())
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            yield SimpleFragment(name=key_base, mass=mass, kind='stub_glycopeptide')

    def glycan_fragments(self, oxonium=True, all_series=False, allow_ambiguous=False):
        r'''
        Generate all oxonium ions for the attached glycan, and
        if `all_series` is `True`, then include the B/Y glycan
        ion ladder, with the peptide attached to the Y ion ladder.

        Parameters
        ----------
        all_series: bool
            Generate the B/Y+peptide ion ladder, otherwise just 2-residue
            pairs for all monosaccharides in :attr:`self.glycan`

        Yields
        ------
        SimpleFragment
        '''
        WATER = Composition("H2O").mass
        TAIL = Composition("CH2O").mass
        if oxonium:
            glycan = None
            if isinstance(self.glycan, Glycan):
                glycan = FrozenGlycanComposition.from_glycan(self.glycan)
            elif isinstance(self.glycan, GlycanComposition):
                glycan = FrozenGlycanComposition(self.glycan)
            else:
                raise TypeError(
                    "Cannot infer monosaccharides from non-Glycan or GlycanComposition {}".format(
                        self.glycan))
            for k in glycan:
                key = str(k)
                mass = k.mass()
                yield SimpleFragment(name=key, mass=mass, kind=oxonium_ion_series)
                yield SimpleFragment(name=key + "-H2O", mass=mass - WATER, kind=oxonium_ion_series)

                yield SimpleFragment(
                    name=key + "-2H2O", mass=mass - 2 * WATER, kind=oxonium_ion_series)
                yield SimpleFragment(
                    name=key + "-2H2O-CH2O", mass=mass - (2 * WATER) - TAIL, kind=oxonium_ion_series)
            for i in range(2, 4):
                for kk in itertools.combinations_with_replacement(glycan, i):
                    for k, v in Counter(kk).items():
                        if glycan[k] < v:
                            continue
                    key = ''.join(map(str, kk))
                    mass = sum(k.mass() for k in kk)
                    yield SimpleFragment(
                        name=key, mass=mass, kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=key + "-H2O", mass=mass - WATER, kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=key + "-2H2O", mass=mass - 2 * WATER, kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=key + "-2H2O-CH2O", mass=mass - (2 * WATER) - TAIL, kind=oxonium_ion_series)

        if isinstance(self.glycan, Glycan) and all_series:
            glycan = self.glycan
            hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
            base_mass = self.mass - (hexnac_mass)
            for fragment in glycan.fragments("BY"):
                if fragment.is_reducing():
                    # TODO:
                    # When self.glycan is adjusted for the attachment cost with the anchoring
                    # amino acid, this WATER penalty can be removed
                    yield SimpleFragment(
                        name="peptide+" + fragment.name, mass=base_mass + fragment.mass - WATER,
                        kind=stub_glycopeptide_series)
                else:
                    yield SimpleFragment(
                        name=fragment.name, mass=fragment.mass, kind=oxonium_ion_series)
        elif allow_ambiguous and all_series:
            _offset = Composition()
            total = FrozenGlycanComposition(self.glycan)
            base = FrozenGlycanComposition(Hex=3, HexNAc=2)
            remainder = total - base
            # GlycanComposition's clone semantics do not propagate the
            # composition_offset attribute yet. Should it?
            remainder.composition_offset = _offset
            stub_mass = self.mass + base.mass() - WATER - FrozenMonosaccharideResidue.from_iupac_lite("HexNAc").mass()
            for composition in descending_combination_counter(remainder):
                composition = FrozenGlycanComposition(composition)
                composition.composition_offset = _offset
                if sum(composition.values()) > 2:
                    yield SimpleFragment(
                        name=composition.serialize(), mass=composition.mass(),
                        kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=composition.serialize() + "-H2O", mass=composition.mass() - (WATER),
                        kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=composition.serialize() + "-H2O-H2O", mass=composition.mass() - (WATER * 2),
                        kind=oxonium_ion_series)
                f = SimpleFragment(
                    name="peptide+" + str(total - composition),
                    mass=stub_mass + remainder.mass() - composition.mass(),
                    kind=stub_glycopeptide_series)
                yield f
        elif all_series:
            raise TypeError("Cannot generate B/Y fragments from non-Glycan {}".format(self.glycan))

    total_composition = _total_composition


Sequence = PeptideSequence
parse = Sequence

_get1 = operator.itemgetter(1)


def cleave(sequence, rule, missed_cleavages=0, min_length=0, **kwargs):
    '''A reimplementation of pyteomics.parser.cleave which produces leaky cleavages
    of a peptide sequence by a regex rule. Includes the cut indices, not included in
    pyteomics.'''
    peptides = []
    if isinstance(sequence, Sequence):
        sequence = str(sequence)
    cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites) - 1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or sequence_length(seq) >= min_length:
                    peptides.append((seq, cleavage_sites[j], cleavage_sites[-1] if cleavage_sites[-1]
                                     is not None else sequence_length(sequence)))
    return sorted(set(peptides), key=_get1)


def itercleave(sequence, rule, missed_cleavages=0, min_length=0, **kwargs):
    if isinstance(sequence, Sequence):
        sequence = str(sequence)
    seen = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites) - 1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or sequence_length(seq) >= min_length:
                    case = ((seq, cleavage_sites[j], cleavage_sites[-1] if cleavage_sites[-1]
                             is not None else sequence_length(sequence)))
                    if case in seen:
                        continue
                    seen.add(case)
                    yield case


class FragmentationStrategyBase(object):
    def __init__(self, name, *args, **kwargs):
        self.name = name

    def handle_modifications(self, fragment_state, modifications):
        raise NotImplementedError()

    def handle_peptide_backbone(self, fragment_state, position):
        raise NotImplementedError()


class FragmentationState(object):
    def __init__(self, peptide):
        self.peptide = peptide

        self.index = 0
        self.current_fragment_mass = 0.
        self.current_fragment_modifications = dict()
        self.current_fragment_composition = Composition()

    def _reset(self):
        self.current_fragment_mass = 0.
        self.current_fragment_modifications = dict()
        self.current_fragment_composition = Composition()

    def __repr__(self):
        return "FragmentationState(%s)" % self.peptide

    def partition_forward(self, position):
        if position > len(self.peptide) - 1 or position < 1:
            raise IndexError(position)
        return self.peptide[:position], self.peptide[position:]

    def partition_reverse(self, position):
        if position > len(self.peptide) - 1 or position < 1:
            raise IndexError(position)
        return self.peptide[-position:], self.peptide[:-position]

    def configure_from_sequence_segment(self, segment):
        mass = 0.
        composition = Composition()
        modification_dict = Counter()
        for residue, modifications in segment:
            mass += residue.mass
            composition += residue.composition
            for mod in modifications:
                mass += mod.mass
                composition += mod.composition
                modification_dict[mod.rule] += 1
        self.current_fragment_composition = composition
        self.current_fragment_mass = mass
        self.current_fragment_modifications = dict(modification_dict)

    def get_flanking_amino_acids(self, front, back):
        return front[-1][0], back[0][0]

    def adjust_for_series(self, series):
        self.current_fragment_mass += series.mass_shift
        self.current_fragment_composition += series.composition_shift
        if series.direction > 0:
            self.current_fragment_mass += self.peptide.n_term.mass
            self.current_fragment_composition += self.peptide.n_term.composition
        elif series.direction < 0:
            self.current_fragment_mass += self.peptide.c_term.mass
            self.current_fragment_composition += self.peptide.c_term.composition
        else:
            self.current_fragment_mass += self.peptide.n_term.mass
            self.current_fragment_composition += self.peptide.n_term.composition

            self.current_fragment_mass += self.peptide.c_term.mass
            self.current_fragment_composition += self.peptide.c_term.composition

    def build_position(self, position, series):
        self._reset()

        if series.direction > 0:
            forward, backward = self.partition_forward(position)
            self.configure_from_sequence_segment(forward)
            self.adjust_for_series(series)
        elif series.direction < 0:
            backward, forward = self.partition_reverse(position)
            self.configure_from_sequence_segment(backward)
            self.adjust_for_series(series)

        flanking_residues = self.get_flanking_amino_acids(forward, backward)

        return PeptideFragment(
            series, position, self.current_fragment_modifications,
            self.current_fragment_mass,
            flanking_amino_acids=flanking_residues,
            composition=self.current_fragment_composition)


class HCDFragmentationStrategy(FragmentationStrategyBase):
    _name = "HCDFragmentationStrategy"

    def __init__(self):
        super(HCDFragmentationStrategy, self).__init__(self._name)

    def handle_peptide_backbone(self, fragment_state, position):
        pass


class ExDFragmentationStrategy(FragmentationStrategyBase):

    def __init__(self):
        super(ExDFragmentationStrategy, self).__init__("ExDFragmentationStrategy")

    def handle_peptide_backbone(self, fragment_state, position):
        pass
