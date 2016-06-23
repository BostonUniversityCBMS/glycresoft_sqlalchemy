import csv
import json
from copy import deepcopy
import re
from pkg_resources import resource_stream

from collections import defaultdict
from collections import Iterable
from glypy.composition.glycan_composition import (
    FrozenMonosaccharideResidue, FrozenGlycanComposition)

from .residue import residue_to_symbol, Residue
from .composition import composition_to_mass, Composition
from .parser import prefix_to_postfix_modifications
from ..utils import Enum
from . import PeptideSequenceBase
from . import ModificationBase
from . import ResidueBase


target_string_pattern = re.compile(
    r'(?P<amino_acid>[A-Z]+)?([@ ]?(?P<n_term>[Nn][-_][tT]erm)|(?P<c_term>[Cc][-_][tT]erm))?')
title_cleaner = re.compile(
    r'(?P<name>.*)\s\((?P<target>.+)\)$')

is_mass_delta = re.compile(r"^Delta:")


class SequenceLocation(Enum):
    anywhere = None
    n_term = -1
    c_term = -2
    protein_n_term = -3
    protein_c_term = -4


SequenceLocation.n_term.add_name("N-term")
SequenceLocation.n_term.add_name("N-Term")
SequenceLocation.n_term.add_name("n-term")
SequenceLocation.n_term.add_name("N_term")
SequenceLocation.c_term.add_name("C-term")
SequenceLocation.c_term.add_name("C-Term")
SequenceLocation.c_term.add_name("c-term")
SequenceLocation.c_term.add_name("C_term")
SequenceLocation.anywhere.add_name("Anywhere")
SequenceLocation.protein_n_term.add_name("Protein N-term")
SequenceLocation.protein_c_term.add_name("Protein C-term")
SequenceLocation.protein_n_term.add_name("Protein N-Term")
SequenceLocation.protein_c_term.add_name("Protein C-Term")


class ModificationCategory(Enum):
    glycosylation = 1
    artefact = 2


def composition_delta_parser(formula):
    '''Parses a Unimod composition offset formula where
    the isotope specification is the 'opposite' of the
    specification used by this library.

    Warning - Currently unused internally.

    Parameters
    ----------
    formula: str
        The formula parsed

    Return
    ------
    dict
    '''
    elemental_composition = {}
    components = re.findall(r"([A-Za-z0-9]+)\(([\-0-9]+)\)", formula)
    for element, amount in components:
        parts = re.search(r"(?P<isotope>\d*)(?P<element>.*)", element).groupdict()
        if parts['isotope'] != '':
            element = "{element}[{isotope}]".format(**parts)
        elemental_composition[element] = amount
    return elemental_composition


def extract_targets_from_string(target_string):
    '''Parses the Protein Prospector modification name string to
    extract the modification target specificity

    Format:

    `"Modification Name (Residues Targeted)"`

    Parameters
    ----------
    target_string: str

    Return
    ------
    :class:`ModificationTarget`

    '''
    position_only = re.search(r"^([NnCc][-_][tT]erm)", target_string)
    if position_only:
        return ModificationTarget(None, SequenceLocation[position_only.groups()[0]])
    else:
        amino_acid_and_position = re.search(
            r'(?P<amino_acid>[A-Z]+)([@ ]?(?P<n_term>[Nn][-_][tT]erm)|(?P<c_term>[Cc][-_][tT]erm))?', target_string)
        params = amino_acid_and_position.groupdict()
        aa = list(map(Residue, params["amino_acid"]))
        if params['n_term']:
            return ModificationTarget(aa, SequenceLocation.n_term)
        elif params['c_term']:
            return ModificationTarget(aa, SequenceLocation.c_term)
        else:
            return ModificationTarget(aa, SequenceLocation.anywhere)


def get_position_modifier_rules_dict(sequence):
    '''Labels the start and end indices of the sequence

    A convenience dictionary initializer for ease-of-access in
    :class:`ModificationTarget` site validation functions

    Parameters
    ----------
    sequence: PeptideSequence
        The sequence for which to project labels

    Return
    ------
    defaultdict
    '''
    return defaultdict(lambda: SequenceLocation.anywhere, **{
        0: SequenceLocation.n_term,
        (len(sequence) - 1): SequenceLocation.c_term
    })


class ModificationTarget(object):
    '''Specifies the range of targets that a given modification can be found at.

    A single instance of :class:`ModificationTarget` describes one or more position/residue
    rules for a single type of modification, bundled together underneath a single :class:`ModificationRule`
    instance. :class:`ModificationTarget` may be extracted from Unimod data dictionaries using the
    :meth:`from_unimod_specificity` class method.

    Attributes
    ----------
    amino_acid_targets: list
        The list of amino acid residues that this rule
        permits.
    position_modifier: SequenceLocation
        The sequence-position specific filter
        this rule enforces.
    '''

    @classmethod
    def from_unimod_specificity(cls, specificity):
        amino_acid = (specificity["site"])
        if amino_acid in set(["N-term", "C-term"]):
            amino_acid = None
        else:
            amino_acid = {Residue(amino_acid)}
        position_modifier = specificity["position"]
        if position_modifier == "Anywhere":
            position_modifier = SequenceLocation.anywhere
        elif position_modifier == "Protein N-term" or position_modifier == "Any N-term":
            position_modifier = SequenceLocation.n_term
        elif position_modifier == "Protein C-term" or position_modifier == "Any C-term":
            position_modifier = SequenceLocation.c_term
        else:
            raise Exception("Undefined Position, " + position_modifier)
        return cls(amino_acid, position_modifier)

    def __init__(self, site_targets=None, position_modifier=None):
        self.amino_acid_targets = site_targets
        self.position_modifier = position_modifier

    def valid_site(self, amino_acid=None, position_modifier=SequenceLocation.anywhere):
        '''Return if a given residue at a sequence position (N-term, C-term, None)
        is valid'''
        valid = False

        # Validate amino acid target target
        valid = (self.amino_acid_targets is None) or (
            amino_acid in self.amino_acid_targets)
        valid = valid and ((self.position_modifier is SequenceLocation.anywhere) or
                           (position_modifier == self.position_modifier))

        if valid and (position_modifier is not SequenceLocation.anywhere) and (
                self.position_modifier is not SequenceLocation.anywhere):
            return position_modifier

        return valid

    def valid_site_seq(self, sequence, position, position_modifier=SequenceLocation.anywhere):
        aa_mod_pair = sequence[position]
        if len(aa_mod_pair[1]) > 0:
            return False
        amino_acid = aa_mod_pair[0]

        return self.valid_site(amino_acid, position_modifier)

    def __repr__(self):
        rep = "{amino_acid_targets}@{position_modifier}".format(
            **self.__dict__)
        return rep

    def __eq__(self, other):
        return repr(self) == repr(other)

    def __hash__(self):
        return hash(repr(self))

    def __len__(self):
        position_modifiers_count = 0
        amino_acid_targets_count = 0
        if self.position_modifiers is not None:
            position_modifiers_count = len(self.position_modifiers)
        if self.amino_acid_targets is not None:
            amino_acid_targets_count = len(self.amino_acid_targets)
        return amino_acid_targets_count + position_modifiers_count


class ModificationRule(object):
    '''Represent the name, mass, and position specifities associated with a given modification.
    Additionally stores information on metadata about this rule computed from these values and
    categorical information from expert sources.

    Attributes
    ----------
    categories : list
        The list of named categories of PTMs that this modification
        falls under. Examples are "artefact", "glycosylation", or "biological process"
    mass : float
        The mass of the modification defined by this rule
    name : str
        A common name for the modification defined by this rule
    names : set
        The set of all names used to refer to this modification, beyond those
        commonly used and those for display purposes.
    options : dict
        Arbitrary information from external sources about this modification
        rule and its phenomena.
    preferred_name : str
        The name for this modification preferred for using in user-facing
        displays
    title : str
        A formal or "full" name for this modification that may be its chemical name
        or some other more verbose representation which is used to identify it but
        not convenient to read or write.
    targets : set
        All applicable :class:`.ModificationTarget` instances which apply to this
        rule. This set may differ from instance to instance of the same rule by
        configuration as per the :class:`RestrictedModificationTable`, for example
    '''

    @staticmethod
    def reduce(rules):
        collector = defaultdict(list)
        for rule in rules:
            collector[rule.name].append(rule)
        reductions = []
        for rule_name, rules in collector.items():
            accum = rules[0]
            for other in rules[1:]:
                accum = accum + other
            reductions.append(accum)
        return reductions

    @classmethod
    def from_protein_prospector(cls, amino_acid_specificity, modification_name,
                                title=None, monoisotopic_mass=None):
        # If the specificity is a string, parse it into rules
        targets = None
        if isinstance(amino_acid_specificity, basestring):
            try:
                targets = set(
                    [extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                print(modification_name, amino_acid_specificity)
                raise TypeError("Could not extract targets")

        # If the specificity is already a rule, store it as a list of one rules
        elif isinstance(amino_acid_specificity, ModificationTarget):
            targets = set([amino_acid_specificity])

        # If the specificity is a collection of rules, store it as is
        elif isinstance(amino_acid_specificity, Iterable) and\
                isinstance(iter(amino_acid_specificity).next(), ModificationTarget):
            targets = set(amino_acid_specificity)

        return cls(targets, modification_name, title, monoisotopic_mass)

    @classmethod
    def from_unimod(cls, unimod_entry):
        specificity = unimod_entry["specificity"]
        monoisotopic_mass = unimod_entry["mono_mass"]
        full_name = unimod_entry["full_name"]
        title = unimod_entry["title"]
        specificity = list(map(ModificationTarget.from_unimod_specificity, specificity))
        return cls(
            specificity, full_name, title, monoisotopic_mass,
            composition=Composition({str(k): int(v) for k, v in unimod_entry['composition'].items()}))

    def __init__(self, amino_acid_specificity, modification_name,
                 title=None, monoisotopic_mass=None, composition=None,
                 categories=None, **kwargs):
        # Attempt to parse the protein prospector name which contains the
        # target in
        try:
            self.name = title_cleaner.search(title).groupdict()['name']
        except:
            self.name = modification_name

        self.unimod_name = modification_name
        self.mass = float(monoisotopic_mass)
        self.title = title if title is not None else modification_name
        self.composition = composition
        self.names = {self.unimod_name, self.name, self.title}
        self.categories = categories or []
        self.options = kwargs
        self._n_term_target = None
        self._c_term_target = None

        self.preferred_name = min(self.unimod_name, self.title, self.name, key=len)

        # The type of the parameter passed for amino_acid_specificity is variable
        # so select the method correct for the passed type

        # If the specificity is a string, parse it into rules
        if isinstance(amino_acid_specificity, basestring):
            try:
                self.targets = set([
                    extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                raise TypeError("Could not extract targets from %s" % amino_acid_specificity)

        # If the specificity is already a rule, store it as a list of one rules
        elif isinstance(amino_acid_specificity, ModificationTarget):
            self.targets = set([amino_acid_specificity])

        # If the specificity is a collection of rules, store it as is
        elif isinstance(amino_acid_specificity, Iterable):
            self.targets = set(amino_acid_specificity)
        else:
            print(amino_acid_specificity, type(amino_acid_specificity))
            raise Exception("Could not interpret target specificity")

    def clone(self, propagated_targets=None):
        if propagated_targets is None:
            propagated_targets = (self.targets)

        dup = self.__class__(
            set(propagated_targets), self.name,
            self.title, self.mass, self.composition,
            self.categories, **self.options)
        dup.names.update(self.names)
        return dup

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return max([target.valid_site(amino_acid, position_modifiers) for
                    target in self.targets])

    def why_valid(self, amino_acid=None, position_modifiers=None):
        possible_targets = [target for target in self.targets if target.valid_site(
            amino_acid, position_modifiers)]
        if len(possible_targets) == 0:
            raise Exception("No valid targets. %s for %s" %
                            (str(amino_acid) + str(position_modifiers), self))
        minimized_target = min(possible_targets, key=len)
        return minimized_target

    def valid_site_seq(self, sequence, position, position_modifiers):
        return any([target.valid_site_seq(sequence, position, position_modifiers) for
                    target in self.targets])

    def find_valid_sites(self, sequence):
        valid_indices = []
        position_modifier_rules = get_position_modifier_rules_dict(sequence)
        for index in range(len(sequence)):
            position = sequence[index]
            if len(position[1]) > 0:
                continue
            amino_acid = position[0].name
            position_modifiers = position_modifier_rules.get(index)
            is_valid = self.valid_site(amino_acid, position_modifiers)
            if(is_valid):
                if (is_valid is SequenceLocation.n_term) or (is_valid is SequenceLocation.c_term):
                    valid_indices.append(is_valid)
                else:
                    valid_indices.append(index)

        return valid_indices

    def serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.preferred_name

    def __repr__(self):
        rep = "{preferred_name}:{mass}".format(
            **self.__dict__)
        return rep

    def __add__(self, other):
        if(isinstance(other, ModificationRule)):
            if (other.targets) is (self.targets):
                return self
            dup = self.clone()
            dup.targets = (set(self.targets) | set(other.targets))
            dup.composition = self.composition or other.composition
            dup.names = self.names | other.names
            dup.preferred_name = min(dup.names, key=len)
            dup.options.update(other.options)
            return dup
        else:
            raise TypeError(
                "Unsupported types. Can only add two ModificationRules together.")

    def __sub__(self, other):
        if isinstance(other, ModificationRule):
            dup = self.clone()
            dup.targets = set(self.targets) - set(other.targets)
            return dup
        else:
            raise TypeError(
                "Unsupported types. Can only subtract two ModificationRules together.")

    def __call__(self, **kwargs):
        return Modification(self, **kwargs)

    def is_a(self, category):
        return category in self.categories

    @property
    def n_term_targets(self):
        if self._n_term_target is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.n_term:
                    solutions.append(target)
            self._n_term_target = solutions
        return self._n_term_target

    @property
    def c_term_targets(self):
        if self._c_term_target is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.c_term:
                    solutions.append(target)
            self._c_term_target = solutions
        return self._c_term_target

    def __eq__(self, other):
        if self is other:
            return True
        idents = self.names
        try:
            other_idents = other.names
        except AttributeError:
            other_idents = {other}
        return len(idents & other_idents) > 0

    def __ne__(self, other):
        return not self == other


class AnonymousModificationRule(ModificationRule):
    '''
    Construct a modification rule from a string that can be placed
    in a peptide sequence and be reconstructed without being stored
    in a :class:`ModificationTable` instance and backed by a modification
    database.

    "...(@{name}-{mass})..."

    '''
    parser = re.compile(r"@(?P<name>.+?)-(?P<mass>[0-9\.]+)")

    @classmethod
    def try_parse(cls, rule_string):
        anon = cls.parser.search(rule_string)
        if anon is not None:
            anon_match_data = anon.groupdict()
            name = anon_match_data['name']
            mass = float(anon_match_data['mass'])
            mod = AnonymousModificationRule(name, mass)
            return mod
        else:
            return None

    def __init__(self, name, mass):
        try:
            super(AnonymousModificationRule, self).__init__(
                "", name, name, mass)
        except AttributeError:
            pass

    def valid_site(self, *args, **kwargs):
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def find_valid_sites(self, *args, **kwargs):
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def serialize(self):
        return "@" + self.name + "-" + str(self.mass)

    def clone(self):
        dup = self.__class__(self.name, self.mass, **self.options)
        return dup

    def __repr__(self):
        rep = "{name}:{mass}".format(
            **self.__dict__)
        return rep


def _simple_sequence_mass(seq_list):
    '''Because AminoAcidSubstitution may contain sequence-like entries that
    need massing but we cannot import `sequence` here, this defines a simple
    sequence list to total mass without any terminal groups, but handles residues
    and modifications.

    For a full sequence_to_mass converter, see the sequence module'''
    mass = 0
    for res, mods in seq_list:
        mass += Residue.mass_by_name(res)
        for mod in mods:
            if mod != "":
                mass += Modification.mass_by_name(mod)
    return mass


class AminoAcidSubstitution(AnonymousModificationRule):
    parser = re.compile(r"@(?P<original_residue>\S+)->(?P<substitution_residue>\S+)")

    @classmethod
    def try_parse(cls, rule_string):
        aa_sub = cls.parser.search(rule_string)
        if aa_sub is not None:
            match_data = aa_sub.groupdict()
            mod = cls(match_data['original_residue'], match_data['substitution_residue'])
            return mod
        else:
            return None

    def __init__(self, original_residue, substitution_residue, **kwargs):
        self.name = "{0}->{1}".format(original_residue, substitution_residue)
        self.mass = 0.0
        self.original = prefix_to_postfix_modifications(original_residue)
        self.substitution = prefix_to_postfix_modifications(substitution_residue)
        self.delta = _simple_sequence_mass(self.substitution) - _simple_sequence_mass(self.original)
        self.preferred_name = self.name
        self.names = {self.name}
        self.options = kwargs
        self.categories = ['substitution']

    def serialize(self):
        return "@" + self.name

    def __repr__(self):
        return "{name}:{delta}".format(**self.__dict__)


_hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
_hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")


class Glycosylation(ModificationRule):
    """
    Incubator Idea - Represent occupied glycosylation sites
    using the Modification interface.

    Attributes
    ----------
    categories : list
        Description
    mass : TYPE
        Description
    name : TYPE
        Description
    names : TYPE
        Description
    options : dict
        Description
    parser : TYPE
        Description
    preferred_name : TYPE
        Description
    title : TYPE
        Description
    """
    parser = re.compile(r"@(?P<glycan_composition>\{[^\}]+\})")

    @classmethod
    def try_parse(cls, rule_string):
        glycosylation = cls.parser.search(rule_string)
        glycosylation = FrozenGlycanComposition.parse(glycosylation.group(0))
        return cls(glycosylation)

    def __init__(self, glycan_composition):
        if isinstance(glycan_composition, basestring):
            glycan_composition = FrozenGlycanComposition.parse(glycan_composition)

        glycan_composition.composition_offset = Composition()

        self.name = "@" + str(glycan_composition)
        self.mass = glycan_composition.mass()
        self.title = self.name
        self.preferred_name = self.name
        self.categories = [ModificationCategory.glycosylation]
        self.names = {self.name, self.title}
        self.options = {}

    def losses(self):
        for i in []:
            yield i

    def clone(self):
        return self.__class__(FrozenGlycanComposition.parse(self.name[1:]))


class NGlycanCoreGlycosylation(Glycosylation):
    mass_ladder = {
        "{HexNAc:1}": _hexnac.mass(),
        "{HexNAc:2}": _hexnac.mass() * 2,
        "{HexNAc:2; Hex:1}": _hexnac.mass() * 2 + _hexose.mass() * 1,
        "{HexNAc:2; Hex:2}": _hexnac.mass() * 2 + _hexose.mass() * 2,
        "{HexNAc:2; Hex:3}": _hexnac.mass() * 2 + _hexose.mass() * 3
    }

    def __init__(self, base_mass=_hexnac.mass()):
        self.name = "HexNAc"
        self.mass = base_mass
        self.title = "N-Glycan Core Glycosylation"
        self.unimod_name = "HexNAc"
        self.targets = [(ModificationTarget("N"))]
        self.composition = _hexnac.total_composition().clone()
        self.title = self.name
        self.names = {self.name}
        self.preferred_name = self.name
        self.options = {}
        self.categories = [ModificationCategory.glycosylation]

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


class OGlcNAcylation(Glycosylation):
    mass_ladder = {
        "{GlcNAc:1}": _hexnac.mass()
    }

    def __init__(self, base_mass=_hexnac.mass()):
        self.name = "O-GlcNAc"
        self.mass = base_mass
        self.title = "O-GlcNAc"
        self.unimod_name = "O-GlcNAc"
        self.targets = [ModificationTarget("S"), ModificationTarget("T")]
        self.composition = _hexnac.total_composition().clone()
        self.title = self.name
        self.preferred_name = self.name
        self.categories = [ModificationCategory.glycosylation]
        self.options = {}

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


def load_from_csv(stream):
    modification_definitions = []
    parser = csv.DictReader(stream)
    modification_definitions = [rec for rec in parser]
    return modification_definitions


def load_from_json(stream):
    modification_definitions = []
    modification_definitions = map(ModificationRule.from_unimod, json.load(stream))
    return modification_definitions


class ModificationTable(dict):

    '''Represents the set of modification rules to apply when generating
    the theoretical sequence space.'''

    # Class Constants
    # Path to the bootstrapping data file for Protein Prospector-defined
    # modifications
    _table_definition_file = staticmethod(lambda: resource_stream(__name__,
                                          "data/ProteinProspectorModifications-for_gly2.csv"))
    _unimod_definitions = staticmethod(lambda: resource_stream(__name__, "data/unimod.json"))

    use_protein_prospector = True

    other_modifications = {
        "HexNAc": NGlycanCoreGlycosylation(),
        "O-GlcNAc": OGlcNAcylation(),
        "H": ModificationRule(
            "N-term", "H", "H",
            composition_to_mass("H"), composition=Composition("H")),
        "OH": ModificationRule(
            "C-term", "OH", "OH",
            composition_to_mass("OH"), composition=Composition("OH")),
    }

    # Class Methods
    @classmethod
    def _definitions_from_stream(cls, stream, format="csv"):
        if format == "csv":
            definitions = load_from_csv(stream)
        elif format == "json":
            definitions = load_from_json(stream)
        return definitions

    @classmethod
    def _definitions_from_stream_default(cls):
        defs = []
        defs += load_from_json(cls._unimod_definitions())
        if cls.use_protein_prospector:
            defs += load_from_csv(cls._table_definition_file())
        return defs

    @classmethod
    def load_from_file(cls, stream=None, format="csv"):
        '''Load the rules definitions from a CSV or JSON file and instantiate a ModificationTable from it'''
        if stream is None:
            return cls.load_from_file_default()
        if format == "csv":
            definitions = load_from_csv(stream)
        elif format == "json":
            definitions = load_from_json(stream)
        return ModificationTable(definitions)

    @classmethod
    def load_from_file_default(cls):
        '''Load the rules definitions from the default package files'''
        defs = []
        defs += load_from_json(cls._unimod_definitions())
        if cls.use_protein_prospector:
            defs += load_from_csv(cls._table_definition_file())
        return ModificationTable(defs)

    bootstrapped = None

    @classmethod
    def bootstrap(cls, reuse=True):
        '''Load the bootstrapping file's definitions and instantiate a ModificationTable from it'''
        if cls.bootstrapped is not None and reuse:
            return cls.bootstrapped
        instance = cls.load_from_file()
        cls.bootstrapped = instance
        return instance

    # Instance Methods
    def __init__(self, modification_definitions=None):
        '''Creates a table from a list of ModificationRules. Assumes that every possibility is available.
        Is an instance of a dictionary object, and therefore can be used as an iterable.

        modification_definitions - A list of ModificationRule objects which define all possible rules from
                                   external sources. They will be filtered down.
        '''
        super(ModificationTable, self).__init__()
        self.modification_rules = []

        if modification_definitions is None:
            modification_definitions = self._definitions_from_stream_default()

        for definition in modification_definitions:
            if not isinstance(definition, ModificationRule):
                definition = ModificationRule(**definition)
            self.modification_rules.append(definition)

        self.modification_rules = ModificationRule.reduce(self.modification_rules)

        self.modification_rules.extend(
            ModificationTable.other_modifications.values())

        self.by_title = {
            modification.title: modification for modification in self.modification_rules
        }

        for name, mod in ModificationTable.other_modifications.items():
            self.add_modification(name, mod)
        for mod in self.modification_rules:
            self.add_modification(mod)

    def add_modification(self, key, value=None):
        if isinstance(key, ModificationRule):
            if key.name in self:
                self[key.name] += key
            else:
                self[key.name] = key
            if key.title is not None:
                if key.title in self:
                    self[key.title] += key
                else:
                    self[key.title] = key
            if key.preferred_name is not None:
                if key.preferred_name in self:
                    self[key.preferred_name] += key
                else:
                    self[key.preferred_name] = key
        elif value is not None:
            if key in self:
                self[key] += value
            else:
                self[key] = value
        else:
            try:
                mod = self.by_title[key]
                if mod.name in self:
                    self[mod.name] += mod
                else:
                    self[mod.name] = mod
                if mod.title is not None:
                    if mod.title in self:
                        self[mod.title] += mod
                    else:
                        self[mod.title] = mod
                if mod.preferred_name is not None:
                    if mod.preferred_name in self:
                        self[mod.preferred_name] += mod
                    else:
                        self[mod.preferred_name] = mod
            except KeyError:
                try:
                    matches = [
                        mod for mod in self.modification_rules if mod.name == key]
                    mod = sum(matches[1:], matches[0])
                    if mod.name in self:
                        self[mod.name] += mod
                    else:
                        self[mod.name] = mod
                except KeyError, e:
                    raise KeyError(e + "Key: %s, Value: %s" % (key, value))

    def remove_modification(self, name, target=None):
        if target is not None:
            self[name].targets.discard(target)
        else:
            self.pop(name)

    def update_modifications(self, *args, **kwargs):
        '''Update with multiple modification rules'''
        for mod in args:
            self.add_modification(mod)
        for name, mod in kwargs.items():
            self.add_modification(name, mod)

    def __getitem__(self, key):
        '''Try to get the value directly, then try the alternative name
        and finally try a cleaned version of the alternative name'''
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            try:
                return self.by_title[key]
            except KeyError:
                try:
                    name = title_cleaner.search(
                        key).groupdict()["name"]
                except KeyError:
                    raise ModificationStringParseError("Could not parse {0}".format(key))
                except AttributeError:
                    raise ModificationStringParseError("Could not parse {0}".format(key))
                try:
                    return dict.__getitem__(self, name)
                except:
                    return self.by_title[name]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def keys(self):
        return dict.keys(self) + self.by_title.keys()

    def values(self):
        return dict.values(self) + self.by_title.values()

    def items(self):
        return dict.items(self) + self.by_title.items()

    def get_modification(self, name, position=-1, number=1):
        '''Build a Modification instance from the ModificationRule given by `name`'''
        try:
            mod_rule = self[name]
        except:
            mod_rule = AnonymousModificationRule.try_parse(name)
            if mod_rule is None:
                raise
        modification_instance = Modification(
            rule=mod_rule, mod_pos=position, mod_num=number)
        return modification_instance


class ModificationStringParseError(Exception):
    pass


class ModificationNameResolutionError(Exception):
    pass


class RestrictedModificationTable(ModificationTable):

    '''Creates a table from a list of ModificationRules. Assumes that we begin with a set of all rules,
    and then selects a subset of them for including in determining valid modifications.

    Is an instance of a dictionary object, and therefore can be used as an iterable.

    :param modification_definitions - A list of ModificationRule objects which define all possible rules from
                               external sources. They will be filtered down.
    :param constant_modifications   - A list of modification names present in the table to include in the
                              constant_modifications sub-table
    :param variable_modifications   - A list of modification names present in the table to include in the
                               variable_modifications sub-table
    '''
    # Class Methods
    @classmethod
    def load_from_file(cls, path=None, constant_modifications=None, variable_modifications=None):
        '''Load the rules definitions from a CSV file and instantiate a RestrictedModificationTable from it'''
        definitions = super(RestrictedModificationTable, cls).load_from_file(path).values()
        inst = RestrictedModificationTable(definitions, constant_modifications, variable_modifications)
        return inst

    bootstrapped = None

    @classmethod
    def bootstrap(cls, constant_modifications=None, variable_modifications=None, reuse=True):
        '''Load the bootstrapping file's definitions and instantiate a RestrictedModificationTable from it'''
        instance = cls.load_from_file(None, constant_modifications, variable_modifications)
        return instance

    def __init__(self, modification_rules=None, constant_modifications=None, variable_modifications=None):
        super(RestrictedModificationTable, self).__init__(modification_rules)
        if constant_modifications is None:
            constant_modifications = []
        self.constant_modifications = constant_modifications
        if variable_modifications is None:
            variable_modifications = []
        self.variable_modifications = variable_modifications

        self.collect_available_modifications()

    def collect_available_modifications(self):
        '''Clears the rules stored in the core dictionary (on self, but not its data members),
        and copies all information in the data member sub-tables into the core dictionary'''
        modifications = {}

        for name in self.constant_modifications + self.variable_modifications:
            mod = deepcopy(self[name])
            try:
                target = extract_targets_from_string(title_cleaner.search(name).groupdict()["target"])
                mod.targets = {target}
            except:
                pass
            if mod.name in modifications:
                modifications[mod.name] += mod
            else:
                modifications[mod.name] = mod

        self.clear()
        self.update(modifications)

        self.by_title = {}
        self.by_title = {
            modification.title: modification for modification in self.values()
        }


class Modification(ModificationBase):

    """Represents a molecular modification, which may be bound at a given position,
    or may be present at multiple locations. This class pulls double duty,
    and when describing the total number of modifications of a type on a molecule, its
    position attributes are unused."""

    _table = ModificationTable()

    __slots__ = ["name", "mass", "position", "number", "rule", "composition"]

    @classmethod
    def mass_by_name(cls, name, mod_num=1):
        if len(name) == 0:
            mass = 0.0
        elif "@" == name[0]:
            mass = float(
                AnonymousModificationRule.parser.search(name).groupdict()['mass'])
        else:
            mass = cls._table[name].mass
        return mass

    def __init__(self, rule, mod_pos=-1, mod_num=1):
        if(isinstance(rule, basestring)):
            try:
                rule = Modification._table[rule]
            except ModificationStringParseError:
                anon = AnonymousModificationRule.try_parse(rule)
                if anon is None:
                    aa_sub = AminoAcidSubstitution.try_parse(rule)
                    if aa_sub is None:
                        glycosylation = Glycosylation.try_parse(rule)
                        if glycosylation is None:
                            raise ModificationNameResolutionError(rule)
                        else:
                            rule = glycosylation
                    else:
                        rule = aa_sub
                else:
                    rule = anon

        self.name = rule.preferred_name
        self.mass = rule.mass
        self.position = mod_pos
        self.number = mod_num
        self.rule = rule
        try:
            self.composition = rule.composition
        except:
            self.composition = None

    def serialize(self):
        rep = self.rule.serialize()
        if self.number > 1:
            # Numbers large than 1 will cause errors in lookup. Need
            # to incorporate a method for inferring number from string
            rep = "{0}{{{1}}}".format(rep, self.number)
        return rep

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return self.rule.valid_site(amino_acid, position_modifiers)

    def why_valid(self, amino_acid=None, position_modifiers=None):
        return self.rule.why_valid(amino_acid, position_modifiers)

    def find_valid_sites(self, sequence):
        return self.rule.find_valid_sites(sequence)

    __repr__ = serialize

    __str__ = serialize

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        try:
            return (self.name == other.name) and (self.number == other.number)
        except AttributeError:
            return self.rule == other

    def __ne__(self, other):
        return not self == other

    def __getstate__(self):
        return [self.name, self.mass, self.position, self.number, self.rule, self.composition]

    def __setstate__(self, state):
        self.name, self.mass, self.position, self.number, self.rule, self.composition = state

    def clone(self):
        return self.__class__(self.rule)

    def is_a(self, category):
        return self.rule.is_a(category)
