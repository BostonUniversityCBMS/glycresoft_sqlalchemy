import re
import multiprocessing
import functools
import itertools
import logging

from glycresoft_ms2_classification.structure import sequence, constants
from glycresoft_ms2_classification.structure import stub_glycopeptides
from glycresoft_ms2_classification.structure.parser import sequence_tokenizer_respect_sequons

from ..data_model import TheoreticalGlycopeptide, Experiment, Protein, DatabaseManager

Sequence = sequence.Sequence
StubGlycopeptide = stub_glycopeptides.StubGlycopeptide
strip_modifications = sequence.strip_modifications
list_to_sequence = sequence.list_to_sequence
logger = logging.getLogger(__name__)


def pair_rotate(sequence):
    chunks = []
    gen = iter(sequence)
    while True:
        chunk = []
        try:
            chunk = [gen.next()]
            chunk.append(gen.next())
            chunks.append(chunk)
        except StopIteration:
            if len(chunk) > 0:
                chunks.append(chunk)
            break
    rev_seq = []
    for chunk in reversed(chunks):
        rev_seq.extend(chunk)
    return rev_seq


def make_decoy(theoretical_sequence, prefix_len=0, suffix_len=1, protein_decoy_map=None, database_manager=None):
    session = database_manager.session()

    theoretical_sequence = session.query(TheoreticalGlycopeptide).filter(
        TheoreticalGlycopeptide.id == theoretical_sequence).first()

    if protein_decoy_map is None:
        protein_decoy_map = {}
    seq = sequence_tokenizer_respect_sequons(theoretical_sequence.glycopeptide_sequence)
    pref = seq[:prefix_len]
    if suffix_len == 0:
        suf = ""
        body = seq[prefix_len:]
    else:
        suf = seq[-suffix_len:]
        body = seq[prefix_len:-suffix_len]
    body = body[::-1]
    rev_seq = (list_to_sequence(pref + list(body) + suf))
    if str(list_to_sequence(seq)) == str(rev_seq):
        rot_body = pair_rotate(body)
        rev_seq = (list_to_sequence(pref + list(rot_body) + suf))

    oxonium_ions, bare_b_ions, bare_y_ions, glycosylated_b_ions, glycosylated_y_ions, stub_ions = fragments(rev_seq)
    decoy = TheoreticalGlycopeptide(
        ms1_score=theoretical_sequence.ms1_score,
        observed_mass=theoretical_sequence.observed_mass,
        calculated_mass=theoretical_sequence.calculated_mass,
        ppm_error=theoretical_sequence.ppm_error,
        volume=theoretical_sequence.volume,
        count_glycosylation_sites=theoretical_sequence.count_glycosylation_sites,
        count_missed_cleavages=theoretical_sequence.count_missed_cleavages,
        start_position=theoretical_sequence.start_position,
        end_position=theoretical_sequence.end_position,
        base_peptide_sequence=strip_modifications(str(rev_seq)),
        modified_peptide_sequence=str(rev_seq),
        peptide_modifications=theoretical_sequence.peptide_modifications,
        glycopeptide_sequence=str(rev_seq) + theoretical_sequence.glycan_composition_str,
        sequence_length=len(rev_seq),
        glycan_composition_str=theoretical_sequence.glycan_composition_str,
        bare_b_ions=bare_b_ions,
        bare_y_ions=bare_y_ions,
        oxonium_ions=oxonium_ions,
        stub_ions=stub_ions,
        glycosylated_b_ions=glycosylated_b_ions,
        glycosylated_y_ions=glycosylated_y_ions,
        protein_id=protein_decoy_map[theoretical_sequence.protein_id]
    )
    session.add(decoy)
    session.commit()
    return 1


def fragments(sequence):
    fragments = zip(*map(sequence.break_at, range(1, len(sequence))))
    b_type = fragments[0]
    b_ions = []
    b_ions_hexnac = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if key == ("b1" or re.search(r'b1\+', key)) and constants.EXCLUDE_B1:
                # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                b_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]  # seq.get_fragments('Y')
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    pep_stubs = StubGlycopeptide.from_sequence(sequence)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()
    return (oxonium_ions, b_ions, b_ions,
            b_ions_hexnac, y_ions_hexnac,
            stub_ions)


class DecoySearchSpaceBuilder(object):
    def __init__(self, database_path, prefix_len=0, suffix_len=1, experiment_ids=None, n_processes=4):
        self.manager = DatabaseManager(database_path)
        self.session = self.manager.session()
        if experiment_ids is None:
            experiment_ids = [eid for experiment in self.session.query(Experiment.id)
                              for eid in experiment]
        self.experiment_ids = experiment_ids
        self.n_processes = n_processes
        self.prefix_len = prefix_len
        self.suffix_len = suffix_len
        self.protein_decoy_map = {}
        for experiment_id in self.experiment_ids:
            reference_experiment = self.session.query(Experiment).filter(
                Experiment.id == experiment_id).first()
            decoy_experiment = Experiment(name=reference_experiment.name.replace("target", "decoy"))
            self.session.add(decoy_experiment)
            for protein in reference_experiment.proteins.values():
                decoy_protein = self.protein_decoy_map[protein.id] = Protein(name='decoy-' + protein.name,
                                                                             experiment_id=decoy_experiment.id)
                self.session.add(decoy_protein)
            self.session.commit()
        for k, v in list(self.protein_decoy_map.items()):
            self.protein_decoy_map[k] = v.id

    def stream_theoretical_glycopeptides(self):
        session = self.manager.session()
        for exp_id in self.experiment_ids:
            for name, protein_id in session.query(
                    Protein.name, Protein.id).filter(Protein.experiment_id == exp_id):
                print(name)
                theoretical_glycopeptide_ids = (session.query(
                       TheoreticalGlycopeptide.id).filter(TheoreticalGlycopeptide.protein_id == protein_id))
                for theoretical_id in itertools.chain.from_iterable(theoretical_glycopeptide_ids):
                    yield theoretical_id
        session.close()

    def prepare_task_fn(self):
        return functools.partial(make_decoy, prefix_len=self.prefix_len, suffix_len=self.suffix_len,
                                 protein_decoy_map=self.protein_decoy_map, database_manager=self.manager)

    def run(self):
        task_fn = self.prepare_task_fn()
        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap(task_fn, self.stream_theoretical_glycopeptides(), chunksize=500):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Decoys Complete." % cntr)
            pool.terminate()
        else:
            for res in itertools.imap(task_fn, self.stream_theoretical_glycopeptides()):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Decoys Complete." % cntr)
