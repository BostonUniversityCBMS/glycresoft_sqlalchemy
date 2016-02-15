import re
import itertools
import functools
import multiprocessing
import operator
import logging

from glycresoft_sqlalchemy.structure import sequence, constants
from glycresoft_sqlalchemy.utils import collectiontools
from glycresoft_sqlalchemy.data_model import (
    TheoreticalPeptideProductIon, TheoreticalGlycopeptideStubIon,
    PipelineModule, TheoreticalGlycopeptide)

Sequence = sequence.Sequence
kind_getter = operator.attrgetter("kind")

logger = logging.getLogger("fragment_generation")


def oxonium_ions_and_stub_ions(seq, full=False):
    if full:
        generators = [seq.stub_fragments(), seq.glycan_fragments(True, all_series=True, allow_ambiguous=True)]
    else:
        generators = [seq.stub_fragments(), seq.glycan_fragments(True, all_series=False, allow_ambiguous=False)]

    def dictify(f):
        return {"key": f.name, "mass": f.mass}

    g = collectiontools.groupby(itertools.chain.from_iterable(generators), kind_getter, transform_fn=dictify)

    r = g.get(sequence.oxonium_ion_series, []), g.get(sequence.stub_glycopeptide_series, [])
    return r


def dbfragments(sequence, id=None, full=True):
    fragments = zip(*map(sequence.break_at, range(1, len(sequence))))
    b_ions = list(map(TheoreticalPeptideProductIon.from_fragment, itertools.chain.from_iterable(fragments[0])))
    y_ions = list(map(TheoreticalPeptideProductIon.from_fragment, itertools.chain.from_iterable(fragments[1])))
    if full:
        generators = [sequence.stub_fragments(), sequence.glycan_fragments(
            True, all_series=True, allow_ambiguous=True)]
    else:
        generators = [sequence.stub_fragments(), sequence.glycan_fragments(
            True, all_series=False, allow_ambiguous=False)]
    stub_ions = list(map(
        TheoreticalGlycopeptideStubIon.from_simple_fragment, itertools.chain.from_iterable(generators)))
    if id is not None:
        for b in b_ions:
            b.theoretical_glycopeptide_id = id
        for y in y_ions:
            y.theoretical_glycopeptide_id = id
        for stub in stub_ions:
            stub.theoretical_glycopeptide_id = id
    return b_ions, y_ions, stub_ions


def fragments(sequence):
    """Generate characteristic 'high energy' HCD fragments for a given glycopeptide sequence

    Parameters
    ----------
    sequence : Sequence

    Returns
    -------
    oxonium_ions : list
    b_ions : list
    y_ions : list
    b_ions_hexnac : list
    y_ions_hexnac : list
    stub_ions : list
    """
    fragments = zip(*map(sequence.break_at, range(1, len(sequence))))
    b_type = fragments[0]
    b_ions = []
    b_ions_hexnac = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if re.search(r'b1\+', key) and constants.EXCLUDE_B1:
                continue
            mass = fm.mass
            if "HexNAc" in key:
                b_ions_hexnac.append({"key": key, "mass": mass})
            else:
                b_ions.append({"key": key, "mass": mass})

    y_type = fragments[1]
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.mass
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass})
            else:
                y_ions.append({"key": key, "mass": mass})

    oxonium_ions, stub_ions = oxonium_ions_and_stub_ions(sequence)
    return (oxonium_ions, b_ions, y_ions,
            b_ions_hexnac, y_ions_hexnac,
            stub_ions)


class WorkItemCollection(object):
    def __init__(self, session):
        self.session = session
        self.accumulator = []

    def add(self, glycopeptide_record):
        self.accumulator.append(glycopeptide_record)

    def reset(self):
        self.accumulator = []

    def commit(self):
        session = self.session
        session.add_all(self.accumulator)
        session.commit()

        self.reset()


class WorkItemCollectionFlat(object):
    def __init__(self, session):
        self.session = session
        self.accumulator = []
        self.glycan_accumulator = []

    def add(self, glycopeptide_record):
        self.accumulator.append(glycopeptide_record)

    def reset(self):
        self.session.expunge_all()
        self.accumulator = []

    def commit(self):
        session = self.session
        n = len(self.accumulator)
        session.bulk_save_objects(self.accumulator)
        session.commit()

        self.reset()
        return n


def flatten(iterable):
    return tuple(itertools.chain.from_iterable(iterable))


def subslurp(ids, session):
    total = len(ids)
    last = 0
    size = 150
    while last <= total:
        chunk = ids[last:(last + size)]
        batch = session.query(
            TheoreticalGlycopeptide.glycopeptide_sequence,
            TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.id.in_(chunk)).all()
        last += size
        for glycopeptide_id_pair in batch:
            yield glycopeptide_id_pair


def generate_fragments_task(ids, database_manager, full=False, **kwargs):
    session = database_manager()
    items = list(subslurp(ids, session))

    backbone_product_ions = []
    stub_ion_container = []

    for glycopeptide, id in items:
        b_ions, y_ions, stub_ions = dbfragments(Sequence(glycopeptide), id=id, full=full)
        backbone_product_ions.extend(b_ions)
        backbone_product_ions.extend(y_ions)
        stub_ion_container.extend(stub_ions)

    session.bulk_save_objects(backbone_product_ions)
    session.bulk_save_objects(stub_ion_container)
    session.commit()
    return len(ids)


class GlycopeptideFragmentGenerator(PipelineModule):
    def __init__(self, database_path, hypothesis_id, full=False, n_processes=4, **kwargs):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.full = full
        self.n_processes = n_processes
        self.options = kwargs

    def clear_existing_fragments(self, session):

        logger.info("Clearing existing fragments")
        ids = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.from_hypothesis(self.hypothesis_id))

        x = session.query(TheoreticalPeptideProductIon).filter(
            TheoreticalPeptideProductIon.theoretical_glycopeptide_id.in_(
                ids)).first()
        if x is None:
            return
        session.query(TheoreticalPeptideProductIon).filter(
            TheoreticalPeptideProductIon.theoretical_glycopeptide_id.in_(
                ids)).delete(synchronize_session=False)
        session.query(TheoreticalGlycopeptideStubIon).filter(
            TheoreticalGlycopeptideStubIon.theoretical_glycopeptide_id.in_(
                ids)).delete(synchronize_session=False)

        session.commit()

    def stream_ids(self, chunk_size=500):
        session = self.manager()
        ids = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.from_hypothesis(
                self.hypothesis_id)).all()
        total = len(ids)
        last = 0
        while last <= total:
            yield flatten(ids[last:(last + chunk_size)])
            last += chunk_size
        session.close()

    def prepare_taskfn(self):
        task_fn = functools.partial(
            generate_fragments_task,
            database_manager=self.manager,
            full=self.full)
        return task_fn

    def run(self):
        session = self.manager()
        self.clear_existing_fragments(session)
        session.close()
        task_fn = self.prepare_taskfn()

        count = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for c in pool.imap_unordered(task_fn, self.stream_ids()):
                count += c
                logger.info("%d sequences fragmented")
        else:
            for c in itertools.imap(task_fn, self.stream_ids()):
                count += c
                logger.info("%d sequences fragmented")
