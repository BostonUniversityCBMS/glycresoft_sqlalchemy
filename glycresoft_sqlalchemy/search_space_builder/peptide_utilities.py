import re
import itertools
from collections import Counter

from glycresoft_sqlalchemy.data_model import Hypothesis, Protein, NaivePeptide
from glycresoft_ms2_classification.structure import sequence, modification, parser

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


def n_glycan_sequon_sites(peptide):
    sites = set(sequence.find_n_glycosylation_sequons(peptide.base_peptide_sequence))
    try:
        if peptide.protein is not None:
            sites |= set(site - peptide.start_position for site in peptide.parent.glycosylation_sites
                         if peptide.start_position <= site < peptide.end_position)
    except AttributeError:
        pass
    return tuple(sites)


def growing_combinations(iterable):
    data = list(iterable)
    return chain_iterable(combinations(data, i) for i in range(len(data) + 1))


def unpositioned_isoforms(
        theoretical_peptide, constant_modifications, variable_modifications, modification_table):
    sequence = Sequence(theoretical_peptide.base_peptide_sequence)
    for mod in {modification_table[const_mod]
                for const_mod in constant_modifications}:
        for site in mod.find_valid_sites(sequence):
            sequence.add_modification(site, mod.name)
    try:
        sequons = theoretical_peptide.n_glycan_sequon_sites
    except:
        print 'ding', theoretical_peptide
        sequons = n_glycan_sequon_sites(theoretical_peptide)
    variable_modifications = {
        modification_table[mod] for mod in variable_modifications}
    variable_sites = {
        mod.name: set(
            mod.find_valid_sites(sequence)) for mod in variable_modifications}
    solutions = set()
    strseq = str(sequence)
    for i in range(len(sequons)):
        for sequons_occupied in (combinations(sequons, i + 1)):
            sequons_occupied = set(sequons_occupied)
            _sequons_occupied = list(sequons_occupied)
            yield strseq, {}, sequence.mass, _sequons_occupied

            avail_sites = {
                mod: sites -
                sequons_occupied for mod,
                sites in variable_sites.items()}
            mods, sites = zip(*avail_sites.items())
            site_combinations = product(*map(growing_combinations, sites))
            for comb in site_combinations:
                assignments = {}
                for i, taken_sites in enumerate(comb):
                    mod = mods[i]
                    for site in (taken_sites):
                        assignments[site] = mod
                modifications = Counter()
                for site, mod in assignments.items():
                    modifications[mod] += 1
                hashable = frozenset(modifications.items())
                if hashable in solutions:
                    continue
                else:
                    solutions.add(hashable)
                mass = sequence.mass
                for name, count in modifications.items():
                    mass += modification_table[name].mass * count
                if len(modifications) > 0:
                    yield strseq, dict(modifications), mass, _sequons_occupied


def generate_peptidoforms(reference_protein, constant_modifications,
                          variable_modifications, enzyme, missed_cleavages=1):
    modtable = modification.RestrictedModificationTable.bootstrap(
        constant_modifications,
        variable_modifications)
    for peptide, start, end in sequence.cleave(
            reference_protein.protein_sequence, expasy_rules[enzyme], missed_cleavages=missed_cleavages):
        if len(peptide) < 5:
            continue
        ref_peptide = NaivePeptide(
            base_peptide_sequence=peptide,
            protein=reference_protein,
            protein_id=reference_protein.id,
            start_position=start,
            end_position=end)
        ref_peptide.parent = reference_protein
        missed = len(re.findall(expasy_rules[enzyme], peptide))
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
