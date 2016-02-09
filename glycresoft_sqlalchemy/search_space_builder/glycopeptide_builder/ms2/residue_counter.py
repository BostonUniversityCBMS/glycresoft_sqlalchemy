import multiprocessing
import functools
from collections import Counter, defaultdict

from scipy.stats import poisson
import numpy as np

from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.structure import sequence_composition
from glycresoft_sqlalchemy.data_model import (
    func, PipelineModule, Protein, InformedPeptide,
    GlycopeptideMatch)

Sequence = sequence.Sequence
Block = AminoAcidSequenceBuildingBlock = sequence_composition.AminoAcidSequenceBuildingBlock


def melt_sequence(sequence_string, counter=None):
    if counter is None:
        counter = Counter()
    for position in Sequence(sequence_string):
        bb = AminoAcidSequenceBuildingBlock(*position)
        counter[bb] += 1
    return counter


def get_residues_from_sequences(sequence_ids, manager, source=GlycopeptideMatch,
                                sequence_attr="glycopeptide_sequence"):
    session = manager.session()
    counter = Counter()
    for sid in sequence_ids:
        s, = session.query(getattr(source, sequence_attr)).filter(source.id == sid.id).first()
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
    def __init__(self, database_path, hypothesis_id, filter=lambda q: q, n_processes=4, source=GlycopeptideMatch,
                 sequence_attr="glycopeptide_sequence"):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.source = source
        self.sequence_attr = sequence_attr
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
            source=self.source,
            sequence_attr=self.sequence_attr)

    def run(self):
        counter = Counter()
        task_fn = self.prepare_task_fn()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for result in pool.imap_unordered(task_fn, self.stream_id_batches()):
                counter += result
            pool.terminate()
        else:
            for chunk in self.stream_id_batches():
                counter += task_fn(chunk)
        return counter


class PTMFrequencyFilter(ResidueCounter):
    def __init__(self, database_path, hypothesis_id, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.source = InformedPeptide
        self.sequence_attr = "modified_peptide_sequence"
        self.filter = lambda q: q
        self.n_processes = n_processes

    def run(self):
        counter = super(PTMFrequencyFilter, self).run()
        mat = make_sparse_modification_matrix(counter)
        session = self.manager()
        modiform_counts = session.query(
            func.count(InformedPeptide.count_variable_modifications),
            InformedPeptide.count_variable_modifications).filter(
            InformedPeptide.from_hypothesis(self.hypothesis_id)).all()
        total = sum(x*y for x, y in modiform_counts)
        average = float(len(mat) * len(mat.values()[0]))
        model = poisson(average)

        probability_dist = model.pmf(np.arange)
        ix_max = probability_dist.argmax()
        effectively_zero_probability = (np.gradient(probability_dist[ix_max:]) < 1e-6).argmax()

        transpose = _dict_transpose(mat)
        filtered = {}
        for modification, counts in transpose.items():
            filtered_entry = {}
            acc = 0
            for residue, count in counts.items():
                if count < effectively_zero_probability:
                    count = 0
                acc += count
                filtered_entry[residue] = count
            if acc > 0:
                filtered[modification] = filtered_entry
        result = _dict_transpose(filtered)
        return result



def _dict_transpose(dd):
    transpose = defaultdict(dict)
    for outer, inner in dd.items():
        for inner_key, inner_value in inner.items():
            transpose[inner_key][outer] = inner_value
    return transpose


def get_residues_from_counter(counter):
    for bb in counter:
        yield bb.residue


def get_modifications_from_counter(counter):
    for bb in counter:
        yield bb.modifications


def make_sparse_modification_matrix(counter, include_unmodified=True):
    matrix = defaultdict(lambda: defaultdict(int))
    for key, count in counter.items():
        if not key.modifications and not include_unmodified:
            continue
        matrix[key.residue][key.modifications] += count
    for residue_ in get_residues_from_counter(counter):
        for modifications in get_modifications_from_counter(counter):
            if not modifications and not include_unmodified:
                continue
            matrix[residue_][modifications]
    return matrix


def _pandas_frequency_matrix(counter):
    import pandas as pd
    mat = make_sparse_modification_matrix(counter)
    mat = pd.DataFrame.from_dict(mat, orient='index')
    return mat
