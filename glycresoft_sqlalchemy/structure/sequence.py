import re
import copy
import itertools
import operator
from collections import defaultdict, deque, Counter, namedtuple

from . import PeptideSequenceBase, MoleculeBase
from . import constants as structure_constants
from .composition import Composition
from .fragment import PeptideFragment, fragment_shift, SimpleFragment
from .modification import Modification
from .residue import Residue
from glypy import GlycanComposition, Glycan, MonosaccharideResidue

from .parser import sequence_tokenizer, sequence_length, strip_modifications

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
        self.modification_index = defaultdict(int)
        for item in seq_list:
            try:
                res = Residue(item[0])
                self.mass += res.mass
                mods = []
                for mod in item[1]:
                    if mod != '':
                        mod = Modification(mod)
                        mods.append(mod)
                        self.modification_index[mod.name] += 1
                        self.mass += mod.mass
                self.seq.append(self.position_class([res, mods]))
            except:
                print(sequence)
                print(item)
                raise

        self._glycan = None
        self.glycan = glycan if glycan != "" else None

        self._n_term = None
        self._c_term = None

        self.n_term = Modification(n_term) if isinstance(n_term, basestring) else n_term
        self.c_term = Modification(c_term) if isinstance(c_term, basestring) else c_term

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
        rep = "{n_term}{seq}{c_term}{glycan}[{mass}]".format(
            n_term=n_term, c_term=c_term,
            glycan=self._glycan if self._glycan is not None else "",
            **self.__dict__)
        return rep

    def __len__(self):
        return len(self.seq)

    @property
    def glycan(self):
        return self._glycan

    @glycan.setter
    def glycan(self, value):
        self._glycan = value
        if isinstance(value, GlycanComposition):
            self._patch_glycan_composition()
        elif isinstance(value, Glycan):
            pass
            # TODO: Make this attach the reducing end to a substituent

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

        flanking_residues = [self.seq[pos][0]]
        b_frag = PeptideFragment(
            "b", pos + structure_constants.FRAG_OFFSET, mod_b, mass_b + b_shift,
            flanking_amino_acids=flanking_residues)

        break_point = pos + 1
        residues_in_y = []
        for pos in range(break_point, len(self)):
            for mod in self.seq[pos][1]:
                mod_y[mod.serialize()] += 1
            residues_in_y.append(self.seq[pos][0].symbol)
            mass_y += self.seq[pos][0].mass
        flanking_residues.append(self.seq[pos][0])

        y_frag = PeptideFragment(
            "y", len(self) - (break_point - 1 + structure_constants.FRAG_OFFSET),
            mod_y, mass_y + y_shift, flanking_amino_acids=flanking_residues)
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

        current_mass = mass_shift
        for idx in range(len(seq_list) - 1):
            for mod in seq_list[idx][1]:
                mod_serial = mod.serialize()
                if mod_serial in mod_dict:
                    mod_dict[mod_serial] += 1
                else:
                    mod_dict[mod_serial] = 1

            current_mass += seq_list[idx][0].mass

            fragments_from_site = []
            flanking_residues = [seq_list[idx][0], seq_list[idx + 1][0]]
            if kind == 'y':
                flanking_residues = flanking_residues[::-1]
            # If incremental loss of HexNAc is not allowed, only one fragment of a given type is generated
            if not structure_constants.PARTIAL_HEXNAC_LOSS:
                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, copy.copy(mod_dict), current_mass,
                    flanking_amino_acids=flanking_residues)
                fragments_from_site.append(frag)
                bare_dict = copy.copy(mod_dict)
                bare_dict["HexNAc"] = 0
                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, copy.copy(bare_dict), current_mass,
                    flanking_amino_acids=flanking_residues)
                fragments_from_site.append(frag)
            # Else a fragment for each incremental loss of HexNAc must be generated
            else:
                frag = PeptideFragment(
                    kind, idx + structure_constants.FRAG_OFFSET, copy.copy(mod_dict), current_mass,
                    flanking_amino_acids=flanking_residues)
                fragments_from_site.extend(frag.partial_loss())
            yield fragments_from_site

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
            self.modification_index[drop_mod.name] -= 1
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
        self.modification_index[mod.name] += 1

    def get_sequence(self, start=0, include_glycan=True, include_termini=True,
                     implicit_n_term=None, implicit_c_term=None):
        """
        Generate human readable sequence string.
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
        if include_glycan:
            if self._glycan is not None:
                rep += str(self._glycan)
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
            self.modification_index[modification.name] += 1
        self.seq.append(self.position_class(next_pos))

    def extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.seq.extend(sequence.seq)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count

    def leading_extend(self, sequence):
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.seq = sequence.seq + self.seq
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count

    @property
    def n_glycan_sequon_sites(self):
        return find_n_glycosylation_sequons(self, structure_constants.ALLOW_MODIFIED_ASPARAGINE)

    def stub_fragments(self):
        if isinstance(self.glycan, Glycan):
            glycan = GlycanComposition.from_glycan(self.glycan)
        elif isinstance(self.glycan, GlycanComposition):
            glycan = self.glycan
        else:
            raise TypeError("Cannot infer monosaccharides from non-Glycan\
             or GlycanComposition {}".format(self.glycan))
        fucose_count = glycan['Fuc'] or glycan['dHex']
        core_count = self.modification_index['HexNAc']

        per_site_shifts = []
        hexose_mass = MonosaccharideResidue.from_iupac_lite("Hex").mass()
        hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
        fucose_mass = MonosaccharideResidue.from_iupac_lite("Fuc").mass()
        base_mass = self.mass - (hexnac_mass * core_count) + Composition("H+").mass
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

                    for hexose_count in range(4):
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

    def glycan_fragments(self, all_series=False, include_constants=True):
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
        PROTON = Composition("H+").mass
        WATER = Composition("H2O").mass
        TAIL = Composition("CH2O").mass
        if not all_series:
            glycan = None
            if isinstance(self.glycan, Glycan):
                glycan = GlycanComposition.from_glycan(self.glycan)
            elif isinstance(self.glycan, GlycanComposition):
                glycan = self.glycan
            else:
                raise TypeError("Cannot infer monosaccharides from non-Glycan or\
                 GlycanComposition {}".format(self.glycan))
            for k in glycan:
                key = str(k)
                mass = k.mass()
                yield SimpleFragment(name=key, mass=mass + PROTON, kind="oxonium_ion")
                yield SimpleFragment(name=key + "-H2O", mass=mass + PROTON - WATER, kind="oxonium_ion")

                yield SimpleFragment(name=key + "-2H2O", mass=mass + PROTON - 2 * WATER, kind="oxonium_ion")
                yield SimpleFragment(name=key + "-2H2O-CH2O", mass=mass + PROTON - (2 * WATER) - TAIL, kind="oxonium_ion")
            for kk in itertools.combinations(glycan, 2):
                key = ''.join(map(str, kk))
                mass = sum(k.mass() for k in kk)
                yield SimpleFragment(name=key, mass=mass + PROTON, kind="oxonium_ion")
                yield SimpleFragment(name=key + "-H2O", mass=mass + PROTON - WATER, kind="oxonium_ion")
                yield SimpleFragment(name=key + "-2H2O", mass=mass + PROTON - 2 * WATER, kind="oxonium_ion")
                yield SimpleFragment(name=key + "-2H2O-CH2O", mass=mass + PROTON - (2 * WATER) - TAIL, kind="oxonium_ion")

        else:
            glycan = None
            if not isinstance(self.glycan, Glycan):
                raise TypeError("Cannot generate B/Y fragments from non-Glycan {}".format(self.glycan))
            glycan = self.glycan
            hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
            base_mass = self.mass - (hexnac_mass) + PROTON
            for fragment in glycan.fragments("BY"):
                if fragment.is_reducing():
                    # TODO:
                    # When self.glycan is adjusted for the attachment cost with the anchoring
                    # amino acid, this WATER penalty can be removed
                    yield SimpleFragment(name="peptide+" + fragment.name, mass=base_mass + fragment.mass - WATER, kind="stub_glycopeptide")
                else:
                    yield SimpleFragment(name=fragment.name, mass=fragment.mass + PROTON, kind="oxonium_ion")

Sequence = PeptideSequence
parse = Sequence

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
