import re
import itertools
import logging
import functools
import multiprocessing
from collections import Counter
import textwrap

from glycresoft_sqlalchemy.data_model import Protein, NaivePeptide
from glycresoft_sqlalchemy.structure import sequence, modification
from glycresoft_sqlalchemy.proteomics.enzyme import expasy_rules, merge_enzyme_rules

from glycresoft_sqlalchemy.utils.collectiontools import SqliteSet, descending_combination_counter

logger = logging.getLogger("peptide_utilities")

Sequence = sequence.Sequence
RestrictedModificationTable = modification.RestrictedModificationTable
combinations = itertools.combinations
product = itertools.product
chain_iterable = itertools.chain.from_iterable

SequenceLocation = modification.SequenceLocation


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


def split_terminal_modifications(modifications):
    n_terminal = []
    c_terminal = []
    internal = []

    for mod in modifications:
        n_term = mod.n_term_targets
        if n_term:
            n_term_rule = mod.clone(n_term)
            mod = mod - n_term_rule
            n_terminal.append(n_term_rule)
        c_term = mod.c_term_targets
        if c_term:
            c_term_rule = mod.clone(c_term)
            mod = mod - c_term_rule
            c_terminal.append(c_term_rule)
        if (mod.targets):
            internal.append(mod)

    return n_terminal, c_terminal, internal


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
        if i == SequenceLocation.n_term or i == SequenceLocation.c_term:
            for c in combinations:
                c[i] = True

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

    sequence = Sequence(theoretical_peptide.base_peptide_sequence)

    has_fixed_n_term = False
    has_fixed_c_term = False

    for mod in {modification_table[const_mod]
                for const_mod in constant_modifications}:
        for site in mod.find_valid_sites(sequence):
            if site == SequenceLocation.n_term:
                has_fixed_n_term = True
            elif site == SequenceLocation.c_term:
                has_fixed_c_term = True
            sequence.add_modification(site, mod.name)
    try:
        sequons = theoretical_peptide.n_glycan_sequon_sites
    except:
        sequons = n_glycan_sequon_sites(theoretical_peptide)

    variable_modifications = {
        modification_table[mod] for mod in variable_modifications}

    (n_term_modifications, c_term_modifications,
     variable_modifications) = split_terminal_modifications(variable_modifications)

    n_term_modifications = [mod for mod in n_term_modifications if mod.find_valid_sites(sequence)]
    c_term_modifications = [mod for mod in c_term_modifications if mod.find_valid_sites(sequence)]

    n_term_modifications.append(None)
    c_term_modifications.append(None)

    variable_sites = {
        mod.name: set(
            mod.find_valid_sites(sequence)) for mod in variable_modifications}
    counter = 0
    kept = 0
    strseq = str(sequence)
    seq_map = {(None, None): sequence}
    strseq_map = {(None, None): strseq}

    for n_term, c_term in itertools.product(n_term_modifications, c_term_modifications):
        seq = sequence.clone()
        if has_fixed_n_term:
            n_term = None

        if n_term is not None:
            seq.n_term = n_term()

        if has_fixed_c_term:
            c_term = None

        if c_term is not None:
            seq.c_term = c_term()

        seq_map[n_term, c_term] = seq
        strseq_map[n_term, c_term] = str(seq)

    solutions = SqliteSet()
    for i in range(len(sequons)):
        for sequons_occupied in (combinations(sequons, i + 1)):
            sequons_occupied = set(sequons_occupied)
            _sequons_occupied = list(sequons_occupied)
            yield strseq, {}, sequence.mass, _sequons_occupied
            for n_term, c_term in itertools.product(n_term_modifications, c_term_modifications):
                        if n_term is c_term is None:
                            continue
                        yield (strseq_map[n_term, c_term], {},
                               seq_map[n_term, c_term].mass, _sequons_occupied)
            counter += 1
            kept += 1
            avail_sites = {
                mod: sites - sequons_occupied
                for mod, sites in variable_sites.items()}
            for modifications in all_combinations(avail_sites):
                if counter % 1000 == 0:
                    logger.info("%d modification configurations computed, %d unique (%s)", counter, kept, strseq)
                counter += 1
                if sum(modifications.values()) > max_modifications:
                    continue
                hashable = frozenset(modifications.items())
                if hashable in solutions:
                    continue
                else:
                    solutions.add(hashable)
                mass = sequence.mass

                modifications = dict(modifications)

                mass_delta = 0
                for name, count in modifications.items():
                    mass_delta += modification_table[name].mass * count

                if len(modifications) > 0:
                    yield strseq, modifications, mass + mass_delta, _sequons_occupied
                    kept += 1

                    for n_term, c_term in itertools.product(n_term_modifications, c_term_modifications):
                        if n_term is c_term is None:
                            continue
                        yield (strseq_map[n_term, c_term], modifications,
                               seq_map[n_term, c_term].mass + mass_delta, _sequons_occupied)


def generate_peptidoforms(reference_protein, constant_modifications,
                          variable_modifications, enzyme, missed_cleavages=1,
                          max_modifications=4, peptide_class=NaivePeptide, **peptide_kwargs):
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

        ref_peptide = peptide_class(
            base_peptide_sequence=peptide,
            protein=reference_protein,
            protein_id=reference_protein.id,
            start_position=start,
            end_position=end,
            **peptide_kwargs)
        ref_peptide.protein = reference_protein
        for modseq, modifications, mass, sequons_occupied in unpositioned_isoforms(ref_peptide, constant_modifications,
                                                                                   variable_modifications,
                                                                                   modtable):
            peptidoform = peptide_class(
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
                sequence_length=len(peptide),
                **peptide_kwargs
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


class FastaFileWriter(object):
    def __init__(self, handle):
        self.handle = handle

    def write(self, defline, sequence):
        self.handle.write(defline)
        self.handle.write("\n")
        self.handle.write(sequence)
        self.handle.write("\n\n")

    def writelines(self, iterable):
        for defline, seq in iterable:
            self.write(defline, seq)


class ProteinFastFileWriter(FastaFileWriter):
    def write(self, protein):
        defline = ''.join([">", protein.name, " ", str(protein.glycosylation_sites)])
        seq = '\n'.join(textwrap.wrap(protein.protein_sequence, 80))
        super(ProteinFastFileWriter, self).write(defline, seq)

    def writelines(self, iterable):
        for protein in iterable:
            self.write(protein)


class ProteomeDigestor(object):
    def __init__(self, database_path, hypothesis_id, protein_ids, constant_modifications,
                 variable_modifications, enzyme, max_missed_cleavages, max_modifications=4,
                 peptide_class=NaivePeptide, peptide_kwargs=None, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications
        self.enzyme = enzyme
        self.max_modifications = max_modifications
        self.max_missed_cleavages = max_missed_cleavages
        self.hypothesis_id = hypothesis_id
        self.protein_ids = protein_ids
        self.n_processes = n_processes
        self.peptide_class = peptide_class
        self.peptide_kwargs = peptide_kwargs

    def stream_proteins(self):
        for i in self.protein_ids:
            yield i

    def prepare_task_fn(self):
        return functools.partial(
            generate_peptidoforms,
            constant_modifications=self.constant_modifications,
            variable_modifications=self.variable_modifications,
            enzyme=self.enzyme,
            max_missed_cleavages=self.max_missed_cleavages,
            max_modifications=self.max_modifications,
            peptide_class=self.peptide_class,
            **self.peptide_kwargs)
