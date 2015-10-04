import multiprocessing
import logging
import functools
from collections import Counter

from glycresoft_sqlalchemy.structure import sequence, residue, fragment
from glycresoft_sqlalchemy.structure import sequence_composition
from glycresoft_sqlalchemy.data_model import (
    PipelineModule, Hypothesis, MS2GlycopeptideHypothesis,
    HypothesisSampleMatch, PeakGroupMatch, Protein,
    TheoreticalGlycopeptideGlycanAssociation,
    TheoreticalGlycopeptide, GlycopeptideMatch)

Sequence = sequence.Sequence
Block = AminoAcidSequenceBuildingBlock = sequence_composition.AminoAcidSequenceBuildingBlock


def melt_sequence(sequence_string, counter=None):
    if counter is None:
        counter = Counter()
    for position in Sequence(sequence_string):
        bb = AminoAcidSequenceBuildingBlock(*position)
        counter[bb] += 1
    return counter


def get_residues_from_sequences(sequence_ids, manager, source=GlycopeptideMatch):
    session = manager.session()
    counter = Counter()
    for sid in sequence_ids:
        s, = session.query(source.glycopeptide_sequence).filter(source.id == sid.id).first()
        melt_sequence(s, counter)
    return counter


def yield_ids(session, hypothesis_id, chunk_size=1000, filter=lambda q: q, source=GlycopeptideMatch):
    base_query = filter(session.query(source.id).filter(
        source.protein_id == Protein.id,
        Protein.hypothesis_id == hypothesis_id))
    chunk = []

    for item in base_query:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
    yield chunk


class ResidueCounter(PipelineModule):
    def __init__(self, database_path, hypothesis_id, filter=lambda q: q, n_processes=4, source=GlycopeptideMatch):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.source = source
        self.filter = filter
        self.n_processes = n_processes

    def stream_id_batches(self):
        session = self.manager.session()
        for chunk in yield_ids(session, self.hypothesis_id, filter=self.filter, source=self.source):
            yield chunk
        session.close()

    def prepare_task_fn(self):
        return functools.partial(
            get_residues_from_sequences,
            manager=self.manager,
            source=self.source)

    def run(self):
        counter = Counter()
        task_fn = self.prepare_task_fn()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for result in pool.imap_unordered(task_fn, self.stream_id_batches()):
                counter += result
        else:
            for chunk in self.stream_id_batches():
                counter += task_fn(chunk)
        return counter
