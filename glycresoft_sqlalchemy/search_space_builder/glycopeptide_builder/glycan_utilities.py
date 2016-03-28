import operator
import logging
import itertools

from collections import Counter


from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Hypothesis,
    TheoreticalGlycanComposition, TheoreticalGlycanCombination,
    TheoreticalGlycanCombinationTheoreticalGlycanComposition)

from glypy.composition.glycan_composition import FrozenGlycanComposition

logger = logging.getLogger("glycan_utilities")


composition_getter = operator.attrgetter("composition")
id_getter = operator.attrgetter("id")

GlycanTable = TheoreticalGlycanComposition.__table__
composition_c = GlycanTable.c.composition
mass_c = GlycanTable.c.calculated_mass
id_c = GlycanTable.c.id
get_id_c = operator.itemgetter(id_c)
get_composition_c = operator.itemgetter(composition_c)


def mirror_glycan_combinations(source_session, target_path, hypothesis_id):
    manager = DatabaseManager(target_path)
    manager.clear(reinitialize=True)
    target_session = manager()

    target_session.add(Hypothesis(id=hypothesis_id))
    target_session.flush()

    # migrate over compositions
    results = source_session.execute(
        TheoreticalGlycanComposition.__table__.select().where(
            TheoreticalGlycanComposition.hypothesis_id == hypothesis_id))
    while True:
        bunch = results.fetchmany(15000)
        if len(bunch) == 0:
            break
        target_session.execute(TheoreticalGlycanComposition.__table__.insert(), bunch)

    # migrate over combinations
    results = source_session.execute(
        TheoreticalGlycanCombination.__table__.select().where(
            TheoreticalGlycanCombination.hypothesis_id == hypothesis_id))

    while True:
        bunch = results.fetchmany(15000)
        if len(bunch) == 0:
            break
        target_session.execute(TheoreticalGlycanCombination.__table__.insert(), bunch)

    junction = TheoreticalGlycanCombinationTheoreticalGlycanComposition.join(
        TheoreticalGlycanCombination,
        TheoreticalGlycanCombinationTheoreticalGlycanComposition.c.combination_id == TheoreticalGlycanCombination.id)
    results = source_session.execute(
        TheoreticalGlycanCombinationTheoreticalGlycanComposition.select().select_from(
            junction).where(TheoreticalGlycanCombination.hypothesis_id == hypothesis_id))

    while True:
        bunch = results.fetchmany(15000)
        if len(bunch) == 0:
            break
        target_session.execute(TheoreticalGlycanCombinationTheoreticalGlycanComposition.insert(), bunch)

    target_session.commit()

    return manager


# def merge_compositions(composition_list):
#     composition_list = list(composition_list)
#     first = GlycanComposition()
#     first.update(**composition_list[0])
#     for comp in composition_list[1:]:
#         first += comp
#     return first.serialize()


def merge_compositions_frozen(composition_list):
    """Given a list of monosaccharide packed strings,
    sum the monosaccharide counts across lists and return
    the merged string

    Parameters
    ----------
    composition_list : list of str

    Returns
    -------
    str
    """
    composition_list = list(composition_list)
    first = FrozenGlycanComposition()
    first.update(composition_list[0])
    for comp in composition_list[1:]:
        first += comp
    return first.serialize()


def create_combinations(session, n, hypothesis_id, unique_unordered=True):
    compositions = list((glycan_id, FrozenGlycanComposition.parse(
        composition)) for glycan_id, composition in session.query(
        TheoreticalGlycanComposition.id, TheoreticalGlycanComposition.composition).filter(
            TheoreticalGlycanComposition.hypothesis_id == hypothesis_id).all())
    join_table_accumulator = []
    j = 0
    for i in range(1, n + 1):
        if i > 1:
            logger.info("Building combinations of size %d", i)
        for comb_compositions in itertools.combinations_with_replacement(compositions, i):
            j += 1
            counts = Counter(g[0] for g in comb_compositions)
            tgc = TheoreticalGlycanCombination(count=i, hypothesis_id=hypothesis_id)
            tgc.composition = merge_compositions_frozen(g[1] for g in comb_compositions)
            tgc.calculated_mass = sum(g[1].mass() for g in comb_compositions)
            session.add(tgc)
            session.flush()
            pk = tgc.id
            for glyc in comb_compositions:
                join_table_accumulator.append({"glycan_id": glyc[0], "combination_id": pk, "count": counts[glyc[0]]})
            if len(join_table_accumulator) % 100000 == 0:
                session.execute(
                    TheoreticalGlycanCombinationTheoreticalGlycanComposition.insert(),
                    join_table_accumulator)
                join_table_accumulator = []

    session.execute(
        TheoreticalGlycanCombinationTheoreticalGlycanComposition.insert(),
        join_table_accumulator)
    join_table_accumulator = []
    session.commit()
    logger.info("%d combinations created", j)
    return j


def query_chunker(query, chunk_size=20000):
    last = 0
    while 1:
        part = query.slice(last, last + chunk_size).all()
        last += chunk_size
        for item in part:
            yield item
        if len(part) < chunk_size:
            break


def get_glycan_combinations(session, n, hypothesis_id):
    return query_chunker(session.query(TheoreticalGlycanCombination).filter(
        TheoreticalGlycanCombination.count == n,
        TheoreticalGlycanCombination.hypothesis_id == hypothesis_id))


class GlycanCombinationProvider(object):
    """
    A simple class wrapping a DatabaseManager + lazy Session to make
    providing TheoreticalGlycanCombination values with state easier.

    This class should be pickle-able before it is first invoked, as
    the _session attribute won't be populated yet

    Attributes
    ----------
    hypothesis_id : int
        The Hypothesis.id value for the source hypothesis to draw combinations from
    manager : DatabaseManager
        The connection source to request data from
    """
    def __init__(self, manager, hypothesis_id):
        self.manager = manager
        self.hypothesis_id = hypothesis_id
        self._session = None

    def combinations(self, n):
        if self._session is None:
            self._session = self.manager()
        return get_glycan_combinations(self._session, n, self.hypothesis_id)

    __call__ = combinations
