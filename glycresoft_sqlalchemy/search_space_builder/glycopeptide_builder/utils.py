import re
import itertools
import operator

from glycresoft_sqlalchemy.structure import sequence, constants
from glycresoft_sqlalchemy.utils import collectiontools
from glycresoft_sqlalchemy.structure import stub_glycopeptides
from glycresoft_sqlalchemy.data_model import (
    TheoreticalPeptideProductIon, TheoreticalGlycopeptideStubIon)

Sequence = sequence.Sequence
StubGlycopeptide = stub_glycopeptides.StubGlycopeptide
kind_getter = operator.attrgetter("kind")


def oxonium_ions_and_stub_ions(seq, full=False):
    if full:
        generators = [seq.stub_fragments(), seq.glycan_fragments(True, all_series=True, allow_ambiguous=True)]
    else:
        generators = [seq.stub_fragments(), seq.glycan_fragments(True, all_series=False, allow_ambiguous=False)]

    def dictify(f):
        return {"key": f.name, "mass": f.mass}

    g = collectiontools.groupby(itertools.chain.from_iterable(generators), kind_getter, transform_fn=dictify)

    r = g.get(sequence.oxonium_ion_series, []), g.get(sequence.stub_glycopeptide_series, [])
    # print map(len, r)
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
    """Generate characteristic 'high energy' CID fragments for a given glycopeptide sequence

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
                # b1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
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

    # pep_stubs = StubGlycopeptide.from_sequence(sequence)
    # stub_ions = pep_stubs.get_stubs()
    # # oxonium_ions = pep_stubs.get_oxonium_ions()
    # oxonium_ions = [{"key": f.name, "mass": f.mass} for f in sequence.glycan_fragments()]
    oxonium_ions, stub_ions = oxonium_ions_and_stub_ions(sequence)
    return (oxonium_ions, b_ions, y_ions,
            b_ions_hexnac, y_ions_hexnac,
            stub_ions)


class WorkItemCollection(object):
    def __init__(self, session):
        self.session = session
        self.accumulator = []
        self.glycan_accumulator = []

    def add(self, glycopeptide_record):
        self.accumulator.append(glycopeptide_record)

    def reset(self):
        self.session.expunge_all()
        self.accumulator = []
        self.glycan_accumulator = []

    def commit(self):
        session = self.session
        session.add_all(self.accumulator)
        session.flush()
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
