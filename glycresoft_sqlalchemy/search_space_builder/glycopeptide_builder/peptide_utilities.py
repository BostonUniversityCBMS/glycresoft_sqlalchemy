import re
import itertools
import logging
try:
    logger = logging.getLogger("peptide_utilities")
except:
    logger = logging
from collections import Counter

from glycresoft_sqlalchemy.data_model import Protein, NaivePeptide
from glycresoft_sqlalchemy.structure import sequence, modification

from glycresoft_sqlalchemy.utils.collectiontools import SqliteSet, descending_combination_counter

Sequence = sequence.Sequence
RestrictedModificationTable = modification.RestrictedModificationTable
combinations = itertools.combinations
product = itertools.product
chain_iterable = itertools.chain.from_iterable


expasy_rules = {'arg-c': 'R',
                'asp-n': '\\w(?=D)',
                'bnps-skatole': 'W',
                'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
                'caspase 10': '(?<=IEA)D',
                'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
                'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
                'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
                'caspase 5': '(?<=[LW]EH)D',
                'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
                'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
                'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
                'caspase 9': '(?<=LEH)D',
                'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
                'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
                'clostripain': 'R',
                'cnbr': 'M',
                'enterokinase': '(?<=[DE]{3})K',
                'factor xa': '(?<=[AFGILTVM][DE]G)R',
                'formic acid': 'D',
                'glutamyl endopeptidase': 'E',
                'granzyme b': '(?<=IEP)D',
                'hydroxylamine': 'N(?=G)',
                'iodosobenzoic acid': 'W',
                'lysc': 'K',
                'ntcb': '\\w(?=C)',
                'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
                'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
                'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
                'proteinase k': '[AEFILTVWY]',
                'staphylococcal peptidase i': '(?<=[^E])E',
                'thermolysin': '[^DE](?=[AFILMV])',
                'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
                'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'}


def merge_enzyme_rules(enzyme_names):
    rules = ["(" + expasy_rules[name] + ")" for name in enzyme_names]
    return "|".join(rules)


def n_glycan_sequon_sites(peptide):
    sites = set(sequence.find_n_glycosylation_sequons(peptide.base_peptide_sequence))
    try:
        if peptide.protein is not None:
            protein = peptide.protein
            sites |= set(site - peptide.start_position for site in protein.glycosylation_sites
                         if peptide.start_position <= site < peptide.end_position)
    except AttributeError:
        pass
    return tuple(sites)


def growing_combinations(iterable):
    data = list(iterable)
    return chain_iterable(combinations(data, i) for i in range(len(data) + 1))


class SiteGenerator(object):
    '''
    A simple generator-generator
    '''
    def __init__(self, iterable):
        self.data = list(iterable)

    def __iter__(self):
        return growing_combinations(self.data)


class SiteCombinator(object):
    '''
    A more memory efficient version of itertools.product() to leverage
    the fact that the iterables being combined can be re-generated instead
    of completely stored in memory.

    This is still a brute force solution. We're more concerned with how many
    (length of) each position as opposed to the exact values in them. A better
    solution would look at how many ways each set can intersect. This gets tricky
    when there are more than two sets and intersection in one blocks the other. A
    similar stack-based algorithm would be effecient for this problem too.
    '''
    def __init__(self, *iterables):
        self.combinators = list(map(SiteGenerator, iterables))
        self.stack_index = 0
        self.generator_stack = None
        self.current = None

    def __iter__(self):
        self.current = [None] * len(self.combinators)
        self.generator_stack = [iter(g) for g in self.combinators]
        self.stack_index = -1
        for g in self.generator_stack:
            self.current[self.stack_index] = g.next()
            self.stack_index += 1

        running = True
        while running:
            try:
                self.current[self.stack_index] = self.generator_stack[self.stack_index].next()
                yield tuple(self.current)
            except StopIteration:
                self.reset_generator()

    def reset_generator(self):
        try:
            self.stack_index -= 1
            self.current[self.stack_index] = self.generator_stack[self.stack_index].next()
            self.stack_index += 1
            self.generator_stack[self.stack_index] = iter(self.combinators[self.stack_index])
        except StopIteration:
            if self.stack_index == -1:
                raise StopIteration("Really stop")
            else:
                self.reset_generator()


def all_combinations(site_assignments):
    all_positions = reduce(set().union, site_assignments.values(), set())
    if len(all_positions) > 10:
        logger.info('Begin Position Assignment Combinations (%d)\n%r', len(all_positions), all_positions)
    intersects = {}
    for i in all_positions:
        intersects[i] = []
        for t, t_sites in site_assignments.items():
            if i in t_sites:
                intersects[i].append(t)

    combinations = [Counter()]
    for i, options in intersects.items():
        if len(options) == 1:
            for c in combinations:
                c[options[0]] += 1
        else:
            new_combinations = []
            for opt in options:
                for c in combinations:
                    new_c = Counter(c)
                    new_c[opt] += 1
                    new_combinations.append(new_c)
            combinations = new_combinations
    if len(all_positions) > 10:
        logger.info("%d combinations computed. Extrapolating", len(combinations))
    unique_combinations = set()
    for combn in combinations:
        for combo in descending_combination_counter(combn):
            unique_combinations.add(frozenset((k, v) for k, v in combo.items() if v != 0))
    combinations = map(dict, unique_combinations)
    if len(all_positions) > 10:
        logger.info("%d combinations extrapolated. Uniquifying", len(combinations))
    unique_combinations = set()
    for combn in combinations:
        for i in range(1, len(combn) + 1):
            for key_set in itertools.combinations(combn, i):
                unique_combinations.add(frozenset((k, combn[k]) for k in key_set))
    result = map(dict, unique_combinations)
    return result


def unpositioned_isoforms(
        theoretical_peptide, constant_modifications, variable_modifications, modification_table, max_modifications=4):
    if variable_modifications is None:
        variable_modifications = []
    if constant_modifications is None:
        constant_modifications = []
    try:
        sequence = Sequence(theoretical_peptide.base_peptide_sequence)
    except KeyError:
        raise StopIteration()
    for mod in {modification_table[const_mod]
                for const_mod in constant_modifications}:
        for site in mod.find_valid_sites(sequence):
            sequence.add_modification(site, mod.name)
    try:
        sequons = theoretical_peptide.n_glycan_sequon_sites
    except:
        sequons = n_glycan_sequon_sites(theoretical_peptide)
    variable_modifications = {
        modification_table[mod] for mod in variable_modifications}
    variable_sites = {
        mod.name: set(
            mod.find_valid_sites(sequence)) for mod in variable_modifications}
    counter = 0
    kept = 0
    strseq = str(sequence)
    for i in range(len(sequons)):
        for sequons_occupied in (combinations(sequons, i + 1)):
            solutions = SqliteSet()
            sequons_occupied = set(sequons_occupied)
            _sequons_occupied = list(sequons_occupied)
            yield strseq, {}, sequence.mass, _sequons_occupied
            counter += 1
            kept += 1
            avail_sites = {
                mod: sites -
                sequons_occupied for mod,
                sites in variable_sites.items()}
            for modifications in all_combinations(avail_sites):
                if counter % 1000 == 0:
                    logger.info("%d modification configurations computed, %d unique (%s)", counter, kept, strseq)
                counter += 1
                hashable = frozenset(modifications.items())
                if hashable in solutions or sum(modifications.values()) > max_modifications:
                    continue
                else:
                    solutions.add(hashable)
                mass = sequence.mass
                for name, count in modifications.items():
                    mass += modification_table[name].mass * count
                if len(modifications) > 0:
                    yield strseq, dict(modifications), mass, _sequons_occupied
                    kept += 1


def generate_peptidoforms(reference_protein, constant_modifications,
                          variable_modifications, enzyme, missed_cleavages=1,
                          max_modifications=4):
    modtable = modification.RestrictedModificationTable.bootstrap(
        constant_modifications,
        variable_modifications, reuse=False)
    enzyme = expasy_rules.get(enzyme, enzyme)
    for peptide, start, end in sequence.cleave(
            reference_protein.protein_sequence, enzyme, missed_cleavages=missed_cleavages):
        if len(peptide) < 5:
            continue
        missed = len(re.findall(enzyme, peptide))
        if missed > missed_cleavages:
            continue

        ref_peptide = NaivePeptide(
            base_peptide_sequence=peptide,
            protein=reference_protein,
            protein_id=reference_protein.id,
            start_position=start,
            end_position=end)
        ref_peptide.protein = reference_protein
        for modseq, modifications, mass, sequons_occupied in unpositioned_isoforms(ref_peptide, constant_modifications,
                                                                                   variable_modifications,
                                                                                   modtable):
            peptidoform = NaivePeptide(
                base_peptide_sequence=peptide,
                modified_peptide_sequence=modseq,
                protein=reference_protein,
                protein_id=reference_protein.id,
                start_position=start,
                end_position=end,
                peptide_modifications='|'.join(str(v) + str(k) for k, v in modifications.items()),
                calculated_mass=mass,
                count_missed_cleavages=missed,
                count_glycosylation_sites=len(sequons_occupied),
                glycosylation_sites=sequons_occupied,
                sequence_length=len(peptide)
            )
            yield peptidoform


class FastaFileParser(object):
    def __init__(self, path):
        self.state = "defline"
        self.handle = open(path)
        self.defline = None
        self.sequence_chunks = []

    def process_result(self, d):
        return d

    def __iter__(self):
        for line in self.handle:
            if self.state == 'defline':
                if line[0] == ">":
                    self.defline = re.sub(r"[\n\r]", "", line[1:])
                    self.state = "sequence"
                else:
                    continue
            else:
                if not re.match(r"^(\s+|>)", line):
                    self.sequence_chunks.append(re.sub(r"[\n\r]", "", line))
                else:
                    if self.defline is not None:
                        yield self.process_result({
                            "name": self.defline,
                            "protein_sequence": ''.join(self.sequence_chunks)
                            })
                    self.sequence_chunks = []
                    self.defline = None
                    self.state = 'defline'
                    if line[0] == '>':
                        self.defline = re.sub(r"[\n\r]", "", line[1:])
                        self.state = "sequence"

        if len(self.sequence_chunks) > 0:
            yield self.process_result({"name": self.defline, "protein_sequence": ''.join(self.sequence_chunks)})


class ProteinFastaFileParser(FastaFileParser):
    def __init__(self, path):
        super(ProteinFastaFileParser, self).__init__(path)

    def process_result(self, d):
        p = Protein(**d)
        p.glycosylation_sites = sequence.find_n_glycosylation_sequons(p.protein_sequence)
        return p


class SiteListFastaFileParser(FastaFileParser):
    def __init__(self, path):
        super(SiteListFastaFileParser, self).__init__(path)

    def process_result(self, d):
        v = d.pop("protein_sequence")
        d['glycosylation_sites'] = set(map(int, v.split(" ")))
        return d
