import time
from sqlalchemy import Table, Column, Integer, ForeignKey, alias, func

from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, GlycopeptideMatch, Base)

from glycresoft_sqlalchemy.data_model.glycomics import make_composition_table


def make_distance_table(model=GlycopeptideMatch):
    return Table(
        "GlycanCombinationDistance_%d" % time.time(), Base.metadata,
        Column("from_id", Integer, ForeignKey(model.id), index=True),
        Column("to_id", Integer, ForeignKey(model.id), index=True),
        Column("distance", Integer, index=True))


def populate_glycan_combination_space(session, selector, model=GlycopeptideMatch):
    conn = session.connection()
    residues = [r[0] for r in model.glycan_composition_extents(session)]
    table = make_composition_table(residues)
    table.create(conn)

    q = session.query(model.GlycanCompositionAssociation).join(model).filter(
        selector).order_by(model.GlycanCompositionAssociation.referent).yield_per(1000)

    group = []
    last_ref = None
    for monosaccharide in q:
        if last_ref is None:
            last_ref = monosaccharide.referent
        elif last_ref != monosaccharide.referent:
            data = {
                'id': last_ref,
            }
            for item in group:
                data[item.base_type] = item.count
            session.execute(table.insert(), [data])
            group = [monosaccharide]
            last_ref = monosaccharide.referent
        else:
            group.append(monosaccharide)
    session.commit()

    distance_table = compute_distance(session, table, residues, model)

    return table, distance_table


def compute_distance(session, table, monosaccharide_names, model=GlycopeptideMatch):
    distance_table = make_distance_table(model)
    distance_table.create(session.connection())
    from_entity = alias(table)
    to_entity = alias(table)

    distances = [getattr(from_entity.c, name) - getattr(to_entity.c, name) for name in monosaccharide_names]
    selected = [from_entity.c.id, to_entity.c.id] + distances
    q = session.query(*selected).join(to_entity, from_entity.c.id != to_entity.c.id).yield_per(1000)
    for fields in q:
        from_id, to_id = fields[:2]
        distance = sum(fields[2:])
        # print from_id, to_id, distance
        session.execute(distance_table.insert(), [{'from_id': from_id, "to_id": to_id, "distance": distance}])
    session.commit()
    return distance_table
