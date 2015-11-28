import datetime
import multiprocessing
import functools
import itertools
import logging

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.structure.parser import sequence_tokenizer_respect_sequons, sequence_tokenizer

from glycresoft_sqlalchemy.data_model import TheoreticalGlycopeptide, Hypothesis, MS2GlycopeptideHypothesis, Protein
from glycresoft_sqlalchemy.data_model import PipelineModule

from ..utils import fragments

Sequence = sequence.Sequence

strip_modifications = sequence.strip_modifications
list_to_sequence = sequence.list_to_sequence
logger = logging.getLogger("make_decoys")


def pair_rotate(sequence):
    """Invert each token pair.

    ABCD -> BADC

    Parameters
    ----------
    sequence : iterable

    Returns
    -------
    list
    """
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


def reverse_preserve_sequon(sequence, prefix_len=0, suffix_len=1):
    original = Sequence(sequence)
    sequence_tokens = sequence_tokenizer_respect_sequons(sequence)
    pref = sequence_tokens[:prefix_len]
    if suffix_len == 0:
        suf = ""
        body = sequence_tokens[prefix_len:]
    else:
        suf = sequence_tokens[-suffix_len:]
        body = sequence_tokens[prefix_len:-suffix_len]
    body = body[::-1]
    rev_sequence = (list_to_sequence(pref + list(body) + suf))
    if str(list_to_sequence(sequence_tokens)) == str(rev_sequence):
        rot_body = pair_rotate(body)
        rev_sequence = (list_to_sequence(pref + list(rot_body) + suf))
    rev_sequence.n_term = original.n_term
    rev_sequence.c_term = original.c_term
    rev_sequence.glycan = original.glycan
    return rev_sequence


def reverse_sequence(sequence, prefix_len=0, suffix_len=1):
    sequence_tokens, mods, glycan, n_term, c_term = sequence_tokenizer(sequence)
    pref = sequence_tokens[:prefix_len]
    if suffix_len == 0:
        suf = ""
        body = sequence_tokens[prefix_len:]
    else:
        suf = sequence_tokens[-suffix_len:]
        body = sequence_tokens[prefix_len:-suffix_len]
    body = body[::-1]
    rev_sequence = (list_to_sequence(pref + list(body) + suf))
    if str(list_to_sequence(sequence_tokens)) == str(rev_sequence):
        rot_body = pair_rotate(body)
        rev_sequence = (list_to_sequence(pref + list(rot_body) + suf))
    rev_sequence.n_term = n_term
    rev_sequence.c_term = c_term
    rev_sequence.glycan = glycan
    return rev_sequence


def make_decoy(theoretical_sequence, prefix_len=0, suffix_len=1,
               protein_decoy_map=None, database_manager=None,
               permute_fn=reverse_preserve_sequon):
    try:
        session = database_manager.session()

        theoretical_sequence = session.query(TheoreticalGlycopeptide).filter(
            TheoreticalGlycopeptide.id == theoretical_sequence).first()

        if protein_decoy_map is None:
            protein_decoy_map = {}

        permuted_sequence = permute_fn(theoretical_sequence.glycopeptide_sequence,
                                       prefix_len=prefix_len, suffix_len=suffix_len)

        (oxonium_ions, bare_b_ions, bare_y_ions, glycosylated_b_ions,
            glycosylated_y_ions, stub_ions) = fragments(permuted_sequence)

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
            base_peptide_sequence=strip_modifications(str(permuted_sequence)),
            modified_peptide_sequence=str(permuted_sequence),
            peptide_modifications=theoretical_sequence.peptide_modifications,
            glycopeptide_sequence=str(permuted_sequence),
            sequence_length=len(permuted_sequence),
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
    except Exception, e:
        logger.exception("%r", locals(), exc_info=e)
        raise e
    finally:
        session.close()
    return 1


decoy_type_map = {
    0: reverse_preserve_sequon,
    1: reverse_sequence
}


class DecoySearchSpaceBuilder(PipelineModule):
    '''
    A pipeline step that builds 
    '''

    HypothesisType = MS2GlycopeptideHypothesis

    def __init__(self, database_path, prefix_len=0, suffix_len=1,
                 hypothesis_ids=None, n_processes=4, decoy_type=0):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()
        HypothesisType = self.HypothesisType
        if hypothesis_ids is None:
            hypothesis_ids = [eid for hypothesis in self.session.query(Hypothesis.id)
                              for eid in hypothesis]
        self.hypothesis_ids = hypothesis_ids
        self.decoy_type = decoy_type
        self.n_processes = n_processes
        self.prefix_len = prefix_len
        self.suffix_len = suffix_len
        self.protein_decoy_map = {}
        self.decoy_hypothesis_ids = []
        for hypothesis_id in self.hypothesis_ids:
            reference_hypothesis = self.session.query(Hypothesis).get(hypothesis_id)
            logger.info("Making decoys for %r", reference_hypothesis)
            # Build Decoy Hypothesis object
            name = reference_hypothesis.name
            if name is None:
                name = 'decoy-{}'.format(datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S"))

            decoy_hypothesis = HypothesisType(name=name.replace("target", "decoy"), is_decoy=True)
            decoy_hypothesis.parameters = {"mirror": reference_hypothesis.id}
            self.session.add(decoy_hypothesis)
            self.session.commit()
            # Create Protein objects to mirror the Reference Hypothesis associated with the Decoy Hypothesis
            for protein in reference_hypothesis.proteins.values():
                decoy_protein = self.protein_decoy_map[protein.id] = Protein(name='decoy-' + protein.name,
                                                                             hypothesis_id=decoy_hypothesis.id)
                self.session.add(decoy_protein)
            self.session.commit()

            parameter = reference_hypothesis.parameters.get("decoys")
            if parameter is None:
                reference_hypothesis.parameters['decoys'] = []
            reference_hypothesis.parameters['decoys'].append({"hypothesis_id": decoy_hypothesis.id, "type": self.decoy_type})
            self.session.add(reference_hypothesis)
            self.session.commit()
            self.decoy_hypothesis_ids.append(decoy_hypothesis.id)
        for k, v in list(self.protein_decoy_map.items()):
            self.protein_decoy_map[k] = v.id

    def stream_theoretical_glycopeptides(self):
        session = self.manager.session()
        i = 0
        for hypothesis_id in self.hypothesis_ids:
            for name, protein_id in session.query(
                    Protein.name, Protein.id).filter(Protein.hypothesis_id == hypothesis_id):
                logger.info("Streaming %s (%d)", name, protein_id)
                theoretical_glycopeptide_ids = (session.query(
                       TheoreticalGlycopeptide.id).filter(TheoreticalGlycopeptide.protein_id == protein_id))
                for theoretical_id in itertools.chain.from_iterable(theoretical_glycopeptide_ids):
                    yield theoretical_id
                    i += 1
        session.close()
        if i == 0:
            raise ValueError("No theoretical peptides streamed")

    def prepare_task_fn(self):
        return functools.partial(make_decoy, prefix_len=self.prefix_len, suffix_len=self.suffix_len,
                                 protein_decoy_map=self.protein_decoy_map, database_manager=self.manager,
                                 permute_fn=decoy_type_map[self.decoy_type])

    def run(self):
        task_fn = self.prepare_task_fn()
        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_theoretical_glycopeptides(), chunksize=500):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Decoys Complete." % cntr)
        else:
            for res in itertools.imap(task_fn, self.stream_theoretical_glycopeptides()):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Decoys Complete." % cntr)
        return self.decoy_hypothesis_ids
