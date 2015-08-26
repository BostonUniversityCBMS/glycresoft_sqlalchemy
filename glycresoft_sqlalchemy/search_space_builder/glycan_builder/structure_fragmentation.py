import functools
import multiprocessing
import logging

from glypy.io import glycoct

from glycresoft_sqlalchemy.data_model import (
    Hypothesis,
    PipelineModule, TheoreticalGlycanComposition, TheoreticalGlycanStructure, MS2GlycanHypothesis,
    )


logger = logging.getLogger("glycan_structure_fragmentation")


def fragment_glycan(glycan_structure, database_manager, kind="ABCXYZ", max_cleavages=1):
    session = database_manager.session()
    try:
        with session.no_autoflush:
            glycan_structure = session.query(TheoreticalGlycanStructure).get(glycan_structure)
            structure = glycoct.loads(glycan_structure.glycoct)
            fragments = {}
            for f in structure.fragments(kind=kind, max_cleavages=max_cleavages):
                fragments[f.name] = f
            glycan_structure._fragments = fragments
            return glycan_structure
    except KeyboardInterrupt, e:
        logger.exception("An error occurred, %r", glycan_structure, exc_info=e)
    finally:
        session.close()


class GlycanStructureFragmenter(PipelineModule):

    def __init__(self, database_path, hypothesis_id, kind="ABCXYZ", max_cleavages=1, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.kind = kind
        self.max_cleavages = max_cleavages

        session = self.manager.session()
        hypothesis = session.query(Hypothesis).get(hypothesis_id)
        hypothesis.parameters = {
        } if hypothesis.parameters is None else hypothesis.parameters
        hypothesis.parameters["fragment_kind"] = kind
        hypothesis.parameters["fragment_max_cleavages"] = max_cleavages
        session.add(hypothesis)
        session.commit()
        session.close()
        self.n_processes = n_processes

    def stream_glycan_structure_ids(self):
        session = self.manager.session()
        try:
            for id, in session.query(TheoreticalGlycanStructure.id).filter(
                    TheoreticalGlycanStructure.hypothesis_id == self.hypothesis_id):
                yield id
        finally:
            session.close()

    def prepare_task_fn(self):
        fn = functools.partial(
            fragment_glycan, database_manager=self.manager, kind=self.kind,
            max_cleavages=self.max_cleavages)
        return fn

    def run(self):
        session = self.manager.session()
        task_fn = self.prepare_task_fn()

        acc = []
        i = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for structure in pool.imap_unordered(task_fn, self.stream_glycan_structure_ids()):
                acc.append(structure)
                i += 1

                if i % 100 == 0:
                    self.inform("%d records handled" % i)
                    map(session.merge, acc)
                    session.commit()
                    acc = []
        else:
            for structure in self.stream_glycan_structure_ids():
                acc.append(task_fn(structure))
                i += 1
                if i % 100 == 0:
                    self.inform("%d records handled" % i)
                    map(session.merge, acc)
                    session.commit()
                    acc = []
        map(session.merge, acc)
        session.commit()
        acc = []

        session.close()
        return self.hypothesis_id
