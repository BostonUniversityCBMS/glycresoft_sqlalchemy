import operator
import logging
from collections import deque

from sqlalchemy.sql import select

from glycresoft_sqlalchemy.data_model import Glycan

logger = logging.getLogger("glycan_utilities")


composition_getter = operator.attrgetter("composition")
id_getter = operator.attrgetter("id")

GlycanTable = Glycan.__table__
composition_c = GlycanTable.c.composition
mass_c = GlycanTable.c.mass
id_c = GlycanTable.c.id
get_id_c = operator.itemgetter(id_c)
get_composition_c = operator.itemgetter(composition_c)


def _glycan_product(session, group, n, hypothesis_id, uniqueness_cache=None):
    for glycan in session.query(Glycan).filter(Glycan.hypothesis_id == hypothesis_id):
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


def get_glycan_combinations(session, n, hypothesis_id, unique_unordered=True):
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
    first = composition_list[0].clone()
    for comp in composition_list[1:]:
        first += comp
    return first.serialize()


def _glycan_product_coreblock(session, group, n, hypothesis_id, uniqueness_cache=None):
    for glycan in session.execute(select(
            [GlycanTable.c.id, GlycanTable.c.mass, GlycanTable.c.composition]).where(
            GlycanTable.c.hypothesis_id == hypothesis_id)):
        if uniqueness_cache is not None:
            id_bunch = tuple(sorted(map(get_id_c, group) + [glycan[id_c]]))
            if id_bunch in uniqueness_cache:
                continue
            uniqueness_cache.add(id_bunch)
        group.append(glycan)
        if n == 1:
            yield tuple(group)
        else:
            for prod in _glycan_product_coreblock(session, group, n - 1, hypothesis_id, uniqueness_cache):
                yield prod
        group.pop()


def get_glycan_combinations_coreblock(session, n, hypothesis_id, unique_unordered=True):
    if n > 1:
        logger.info("Combining glycans %d at a time", n)
    group = deque(maxlen=n)
    if unique_unordered:
        seen_set = set()
        for group in _glycan_product_coreblock(session, group, n, hypothesis_id, set()):
            composition = tuple(sorted(map(get_composition_c, group)))
            if composition in seen_set:
                continue
            seen_set.add(composition)
            yield group
    else:
        for group in _glycan_product_coreblock(session, group, n, hypothesis_id):
            if n == 1:
                group = [group]
            yield group
