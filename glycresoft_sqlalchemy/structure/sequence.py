import re
import copy
import warnings
import itertools
import operator
from collections import defaultdict, deque

from . import PeptideSequenceBase, MoleculeBase
from . import constants as structure_constants
from .composition import Composition, composition_to_mass, std_ion_comp
from .fragment import Fragment, fragment_pairing, fragment_direction, fragment_shift
from .modification import Modification
from .residue import Residue

from .parser import *

from ..utils.iterators import peekable
from ..utils.memoize import memoize


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


def find_n_glycosylation_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {"HexNAc"}
    except:
        allow_modified = {"HexNAc"}
    state = "seek"  # [seek, n, ^p, st]
    if(isinstance(sequence, Sequence)):
        asn = "Asn"
        pro = "Pro"
        ser = "Ser"
        thr = "Thr"
    else:
        sequence, mods, glycan, n_term, c_term = sequence_tokenizer(sequence)
        asn = "N"
        pro = "P"
        ser = "S"
        thr = "T"
    i = 0
    positions = []
    n_pos = None
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            # A sequon starts with an Asn residue without modifications, or for counting
            # purposes one that has already been glycosylated
            if next_pos[0] == asn:
                if ((len(next_pos[1]) == 0) or (next_pos[1][0] == '')
                   or next_pos[1][0] in allow_modified):
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


def golden_pair_map(sequence):
    seq_obj = Sequence(sequence)
    key_to_golden_pairs = {}
    fragments = map(seq_obj.break_at, range(1, len(seq_obj)))
    for pair in fragments:
        for frag in pair:
            if isinstance(frag, list):
                for item in frag:
                    key_to_golden_pairs[item.name] = item.golden_pairs
            else:
                key_to_golden_pairs[frag.name] = frag.golden_pairs
    return key_to_golden_pairs


class Sequence(PeptideSequenceBase):
    '''
    Represents a peptide that may have post-translational modifications
    including glycosylation.
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
                elif next_pos == StopIteration:
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
                seq.mod_index[mod.name] += 1
            seq.seq.append(cls.position_class([resid, mod_list]))
        if not isinstance(n_term, MoleculeBase):
            n_term = Modification(n_term)
        if not isinstance(c_term, MoleculeBase):
            c_term = Modification(c_term)

        seq.n_term = n_term
        seq.c_term = c_term
        return seq

    def __init__(self, sequence, **kwargs):
        seq_list, modifications, glycan, n_term, c_term = sequence_tokenizer(sequence)
        self.mass = 0.0
        self.seq = []
        self.mod_index = defaultdict(int)
        for item in seq_list:
            try:
                res = Residue()
                res.by_symbol(item[0])
                self.mass += res.mass
                mods = []
                for mod in item[1]:
                    if mod != '':
                        mod = Modification(mod)
                        mods.append(mod)
                        self.mod_index[mod.name] += 1
                        self.mass += mod.mass
                self.seq.append(self.position_class([res, mods]))
            except:
                print(sequence)
                print(item)
                raise
        # Termini
        # self.mass += Composition("H2O").mass

        # For compatibility with Glycopeptide sequences
        self.glycan = glycan

        # For compatibility with terminal modifications
        self._n_term = None
        self._c_term = None

        self.n_term = Modification(n_term) if isinstance(n_term, basestring) else n_term
        self.c_term = Modification(c_term) if isinstance(c_term, basestring) else c_term

    def __repr__(self):
        n_term = ""
        if self.n_term is not None:
            n_term = "({0})-".format(self.n_term)
        c_term = ""
        if self.c_term is not None:
            c_term = "-({0})".format(self.c_term)
        rep = "{n_term}{seq}{c_term}{glycan}[{mass}]".format(n_term=n_term, c_term=c_term, **self.__dict__)
        return rep

    def __len__(self):
        return self.length

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

    @property
    def length(self):
        return len(self.seq)

    def __iter__(self):
        '''Makes the sequence iterable'''
        for i in self.seq:
            yield(i)

    def __getitem__(self, index):
        '''A Pythonic way to access an index in the sequence'''
        sub = self.seq[index]
        return sub

    # Backwards compatibility
    at = __getitem__

    def __setitem__(self, index, value):
        '''A Pythonic way to set the value at an index, but does not
        try to validate the result.'''
        self.seq[index] = value

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
        for i, pos in enumerate(self):
            mods = [mod.name for mod in pos[1] if mod.name == "HexNAc" or "Glycan" in mod.name]
            for mod in mods:
                self.drop_modification(i, mod)

    def break_at(self, idx):
        b_shift = fragment_shift['b']
        y_shift = fragment_shift['y']

        mod_b = defaultdict(int)
        mass_b = 0

        mod_y = defaultdict(int)
        mass_y = 0

        pos = 0
        residues_in_b = []
        for pos in range(idx):
            for mod in self.seq[pos][1]:
                mod_b[mod.serialize()] += 1
            residues_in_b.append(self.seq[pos][0].symbol)
            mass_b += self.seq[pos][0].mass

        b_frag = Fragment("b", pos + structure_constants.FRAG_OFFSET, mod_b, mass_b + b_shift)

        break_point = pos + 1
        residues_in_y = []
        for pos in range(break_point, len(self)):
            for mod in self.seq[pos][1]:
                mod_y[mod.serialize()] += 1
            residues_in_y.append(self.seq[pos][0].symbol)
            mass_y += self.seq[pos][0].mass

        y_frag = Fragment(
            "y", len(self) - (break_point - 1 + structure_constants.FRAG_OFFSET), mod_y, mass_y + y_shift)
        if structure_constants.PARTIAL_HEXNAC_LOSS:
            b_frag.golden_pairs = [frag.name for frag in y_frag.partial_loss()]
            y_frag.golden_pairs = [frag.name for frag in b_frag.partial_loss()]
            b_frag = list(b_frag.partial_loss())
            y_frag = list(y_frag.partial_loss())
        else:
            b_frag.golden_pairs = [y_frag.name]
            y_frag.golden_pairs = [b_frag.name]

        return b_frag, y_frag

    def get_fragments(self, kind, include_golden_pairs=False):
        """Return a list of mass values for each fragment of `kind`"""

        mass_shift = 0.0

        # The set of modification names.
        mod_dict = {}

        # The key is the position, the value is an array of fragments.
        # And the first element is always bare fragment.
        # The total number of HexNAc on the fragment should be recorded.

        if kind == "b":
            seq_list = self.seq
            # Hydrogen ionized is from terminal modification
            mass_shift = fragment_shift['b']

        elif kind == "y":
            # y ions abstract a proton from the precursor
            mass_shift = fragment_shift['y']
            seq_list = list(reversed(self.seq))

        current_mass = 0
        for idx in range(len(seq_list) - 1):
            for mod in seq_list[idx][1]:
                mod_serial = mod.serialize()
                if mod_serial in mod_dict:
                    mod_dict[mod_serial] += 1
                else:
                    mod_dict[mod_serial] = 1

            if idx == 0:
                current_mass = seq_list[0][0].mass + mass_shift
            else:
                current_mass = current_mass + seq_list[idx][0].mass

            frag_dri = []
            # If incremental loss of HexNAc is not allowed, only one fragment of a given type is generated
            if not structure_constants.PARTIAL_HEXNAC_LOSS:
                frag = Fragment(kind, idx + structure_constants.FRAG_OFFSET, copy.copy(mod_dict), current_mass)
                frag_dri.append(frag)
                bare_dict = copy.copy(mod_dict)
                bare_dict["HexNAc"] = 0
                frag = Fragment(kind, idx + structure_constants.FRAG_OFFSET, copy.copy(bare_dict), current_mass)
                frag_dri.append(frag)
            # Else a fragment for each incremental loss of HexNAc must be generated
            else:
                frag = Fragment(kind, idx + structure_constants.FRAG_OFFSET, copy.copy(mod_dict), current_mass)
                frag_dri.extend(frag.partial_loss())
            yield frag_dri

    def drop_modification(self, pos, mod_type):
        '''
        Drop a modification by name from a specific residue
        '''
        dropped_index = None
        for i, mod in enumerate(self.seq[pos][1]):
            if mod_type == mod.name:
                dropped_index = i
                break
        try:
            drop_mod = self.seq[pos][1].pop(dropped_index)
            self.mass -= drop_mod.mass
            self.mod_index[drop_mod.name] -= 1
        except:
            raise ValueError("Modification not found! %s @ %s" % (mod_type, pos))

    def add_modification(self, pos=None, mod_type=None):
        if pos is None and isinstance(mod_type, Modification):
            pos = mod_type.position

        if (pos == -1) or (pos >= len(self.seq)):
            raise IndexError(
                "Invalid modification position. %s, %s, %s" %
                (pos, str(self.seq), mod_type))
            return -1
        if isinstance(mod_type, Modification):
            mod = mod_type
        else:
            mod = Modification(rule=mod_type, mod_pos=pos)
        self.seq[pos][1].append(mod)
        self.mass += mod.mass
        self.mod_index[mod.name] += 1

    def get_sequence(self, start=0, include_glycan=True, include_termini=True,
                     implicit_n_term=None, implicit_c_term=None):
        """
        Generate human readable sequence string.
        DQYELLC(Carbamidomethyl)LDN(HexNAc)TR
        """
        if implicit_n_term is None:
            implicit_n_term = structure_constants.N_TERM_DEFAULT
        if implicit_c_term is None:
            implicit_c_term = structure_constants.C_TERM_DEFAULT

        seq_list = []
        for x, y in self.seq[start:]:
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
        rep += str(self.glycan) if include_glycan else ""
        return rep

    __str__ = get_sequence

    def get_sequence_list(self, start=0):
        return self.seq[start:]

    def clone(self):
        return copy.deepcopy(self)

    def append(self, residue, modification=None):
        self.mass += residue.mass
        next_pos = [residue]
        if modification is None:
            next_pos.append([])
        else:
            next_pos.append([modification])
            self.mass += modification.mass
            self.mod_index[modification.name] += 1
        self.seq.append(self.position_class(next_pos))

    def extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = Sequence(sequence)
        self.seq.extend(sequence.seq)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.mod_index.items():
            self.mod_index[mod] += count

    def leading_extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = Sequence(sequence)
        self.seq = sequence.seq + self.seq
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.mod_index.items():
            self.mod_index[mod] += count

    @property
    def n_glycan_sequon_sites(self):
        return find_n_glycosylation_sequons(self, structure_constants.ALLOW_MODIFIED_ASPARAGINE)


class GrowingSequence(Sequence):

    def __init__(self, sequence, cleavage_pattern):
        super(GrowingSequence, self).__init__(sequence)
        self.cleavage_pattern = cleavage_pattern
        self.missed_cleavages = cleavage_pattern.count(sequence)

    def extend(self, sequence):
        super(GrowingSequence, self).extend(sequence)
        self.missed_cleavages += sequence.missed_cleavages

    def leading_extend(self, sequence):
        super(GrowingSequence, self).leading_extend(sequence)
        self.missed_cleavages += sequence.missed_cleavages

    def pad(self):
        for padded in self.cleavage_pattern.pad(self):
            yield GrowingSequence(padded, self.cleavage_pattern)


class Protease(object):

    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def count(self, sequence):
        cnt = 0
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = Sequence(sequence)
        cnt += sum([resid.symbol in self.prefix for resid,
                    mods in sequence[1:]])
        cnt += sum([resid.symbol in self.suffix for resid,
                    mods in sequence[:-1]])
        return cnt

    def pad(self, sequence):
        for pref in self.prefix:
            for suf in self.suffix:
                yield "{0}{1}{2}".format(pref, str(sequence), suf)

get1 = operator.itemgetter(1)


def cleave(sequence, rule, missed_cleavages=0, min_length=0, **kwargs):
    '''A reimplementation of pyteomics.parser.cleave which produces leaky cleavages
    of a peptide sequence by a regex rule. Includes the cut indices, not included in
    pyteomics.'''
    peptides = []
    if isinstance(sequence, Sequence):
        sequence = str(sequence)
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites)-1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or sequence_length(seq) >= min_length:
                    peptides.append((seq, cleavage_sites[j], cleavage_sites[-1] if cleavage_sites[-1]
                                     is not None else sequence_length(sequence)))
    return sorted(set(peptides), key=get1)


def itercleave(sequence, rule, missed_cleavages=0, min_length=0, **kwargs):
    if isinstance(sequence, Sequence):
        sequence = str(sequence)
    seen = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites)-1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or sequence_length(seq) >= min_length:
                    case = ((seq, cleavage_sites[j], cleavage_sites[-1] if cleavage_sites[-1]
                             is not None else sequence_length(sequence)))
                    if case in seen:
                        continue
                    seen.add(case)
                    yield case


@memoize()
def sequence_tokens_to_mass(tokens):
    mass = 0.0
    for residue, mods in tokens:
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
