import operator
import logging
import itertools
from collections import deque, Counter, OrderedDict

from sqlalchemy.sql import select

from glycresoft_sqlalchemy.data_model import (
    TheoreticalGlycanComposition, TheoreticalGlycanCombination, TheoreticalGlycanCombinationTheoreticalGlycanComposition)

from glycresoft_sqlalchemy.utils.collectiontools import flatten

from glypy import GlycanComposition
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


def _glycan_product(session, group, n, hypothesis_id, uniqueness_cache=None):
    for glycan in session.query(TheoreticalGlycanComposition).filter(
            TheoreticalGlycanComposition.hypothesis_id == hypothesis_id):
        if uniqueness_cache is not None:
            id_bunch = tuple(sorted(map(id_getter, group) + [glycan.id]))
            if id_bunch in uniqueness_cache:
                continue
            uniqueness_cache.add(id_bunch)
        group.append(glycan)
        if n == 1:
            yield tuple(group)
        else:
            for prod in _glycan_product(session, group, n - 1, hypothesis_id, uniqueness_cache):
                yield prod
        group.pop()


def get_glycan_combinations_in_place(session, n, hypothesis_id, unique_unordered=True):
    if n > 1:
        logger.info("Combining glycans %d at a time", n)
    group = deque(maxlen=n)
    if unique_unordered:
        seen_set = set()
        for group in _glycan_product(session, group, n, hypothesis_id, set()):
            composition = tuple(sorted(map(composition_getter, group)))
            if composition in seen_set:
                continue
            seen_set.add(composition)
            yield group
    else:
        for group in _glycan_product(session, group, n, hypothesis_id):
            if n == 1:
                group = [group]
            yield group


def merge_compositions(composition_list):
    composition_list = list(composition_list)
    first = GlycanComposition()
    first.update(**composition_list[0])
    for comp in composition_list[1:]:
        first += comp
    return first.serialize()


def merge_compositions2(composition_list):
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
        for comb_compositions in itertools.combinations_with_replacement(compositions, i):
            j += 1
            counts = Counter(g[0] for g in comb_compositions)
            tgc = TheoreticalGlycanCombination(count=i, hypothesis_id=hypothesis_id)
            tgc.composition = merge_compositions2(g[1] for g in comb_compositions)
            tgc.calculated_mass = sum(g[1].mass() for g in comb_compositions)
            session.add(tgc)
            session.flush()
            pk = tgc.id
            for glyc in comb_compositions:
                join_table_accumulator.append({"glycan_id": glyc[0], "combination_id": pk, "count": counts[glyc[0]]})
            if len(join_table_accumulator) % 1000 == 0:
                session.execute(
                    TheoreticalGlycanCombinationTheoreticalGlycanComposition.insert(),
                    join_table_accumulator)
                join_table_accumulator = []

    session.execute(
        TheoreticalGlycanCombinationTheoreticalGlycanComposition.insert(),
        join_table_accumulator)
    join_table_accumulator = []
    session.commit()


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
