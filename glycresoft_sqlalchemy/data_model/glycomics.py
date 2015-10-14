import re
import sys
import operator
import functools

from sqlalchemy.orm import relationship, backref
from sqlalchemy import alias, event, func
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey, Table, PickleType, Boolean
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.ext.associationproxy import association_proxy, _AssociationDict
from sqlalchemy.orm.collections import attribute_mapped_collection

import glypy
from glypy.composition import glycan_composition
from glypy.io import glycoct
from glypy.algorithms import subtree_search

from .base import Base, Namespace
from .generic import (
    MutableDict, MutableList, HasTaxonomy, Taxon, ReferenceDatabase,
    ReferenceAccessionNumber, HasReferenceAccessionNumber)

from ..utils.database_utils import get_or_create
from ..utils.memoize import memoclone


crossring_pattern = re.compile(r"\d,\d")
glycoct_parser = memoclone(100)(glycoct.loads)


class MassShift(Base):
    __tablename__ = "MassShift"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), index=True)
    mass = Column(Numeric(10, 6, asdecimal=False))

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return (self.name == other.name) and (self.mass == other.mass)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "<MassShift {self.name} {self.mass}>".format(self=self)

    @classmethod
    def get(cls, session, name, mass=None):
        obj, made = get_or_create(session, cls, name=name, mass=mass)
        return obj


class AssociationComposition(_AssociationDict):
    def __init__(self, lazy_collection, *args, **kwargs):
        _AssociationDict.__init__(self, lazy_collection, *args, **kwargs)
        self.lazy_collection = lazy_collection

    def mass(self, *args, **kwargs):
        return glycan_composition.GlycanComposition(**self).mass(*args, **kwargs)

    def total_composition(self, *args, **kwargs):
        return glycan_composition.GlycanComposition(**self).total_composition(*args, **kwargs)

    def serialize(self):
        return glycan_composition.GlycanComposition(**self).serialize()

    def clone(self):
        return glycan_composition.GlycanComposition(**self)

    def __getitem__(self, key):
        try:
            return super(AssociationComposition, self).__getitem__(key)
        except KeyError:
            return 0

    def __iadd__(self, other):
        for elem, cnt in (other.items()):
            self[elem] += cnt
        return self

    def __add__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        return self

    def __sub__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError(
                'Cannot multiply Composition by non-integer',
                other)
        prod = glycan_composition.GlycanComposition()
        for k, v in self.items():
            prod[k] = v * other

        return (prod)

    def copy(self):
        return glycan_composition.GlycanComposition(**self)


def association_composition_creator(k, v, reference_table):
    return reference_table(base_type=k, count=v)


def composition_association_factory(lazy_collection, creator, value_attr, assoc_prox):
    _getter = operator.attrgetter(value_attr)

    def getter(target):
        return _getter(target) if target is not None else None

    def setter(o, k, v):
        return setattr(o, value_attr, v)
    return AssociationComposition(lazy_collection, creator, getter, setter, assoc_prox)


def has_glycan_composition(model, composition_attr):
    if hasattr(model, "GlycanCompositionAssociation"):
        has_glycan_composition_listener(getattr(model, composition_attr))
        return model
    else:
        class MonosaccharideBaseCounter(Base):
            __tablename__ = "%s_GlycanCompositionAssociation" % model.__name__
            referent = Column(
                Integer, ForeignKey(model.id, ondelete="CASCADE"), primary_key=True)
            base_type = Column(Unicode(30), primary_key=True)
            count = Column(Integer)

            composition = relationship(model, backref=backref(
                "_glycan_composition",
                collection_class=attribute_mapped_collection("base_type"),
                cascade="all, delete-orphan"))

            def __repr__(self):
                return "<{} {}>".format(self.base_type, self.count)

        @event.listens_for(getattr(model, composition_attr), "set")
        def convert_composition(target, value, oldvalue, initiator):
            if value == "{}" or value is None:
                return
            for k, v in glycan_composition.parse(value).items():
                target.glycan_composition[k.name()] = v

        @event.listens_for(model, "load")
        def convert_composition_load(target, context):
            value = getattr(target, composition_attr)
            if value == "{}" or value is None:
                return
            for k, v in glycan_composition.parse(value).items():
                target.glycan_composition[k.name()] = v

        creator = functools.partial(
            association_composition_creator, reference_table=MonosaccharideBaseCounter)
        MonosaccharideBaseCounter.__name__ = "%s_MonosaccharideBaseCounter" % model.__name__
        model.GlycanCompositionAssociation = MonosaccharideBaseCounter

        model.glycan_composition = association_proxy(
            '_glycan_composition', 'count',
            creator=creator,
            proxy_factory=composition_association_factory
            )

        def qmonosaccharide(cls, monosaccharide_name):
            if monosaccharide_name in cls._qmonosaccharide_cache:
                return cls._qmonosaccharide_cache[monosaccharide_name]
            symbol = alias(cls.GlycanCompositionAssociation.__table__.select().where(
                cls.GlycanCompositionAssociation.__table__.c.base_type == monosaccharide_name),
                monosaccharide_name)
            cls._qmonosaccharide_cache[monosaccharide_name] = symbol
            return symbol

        model.qmonosaccharide = classmethod(qmonosaccharide)
        model._qmonosaccharide_cache = {}

        def with_monosaccharide(cls, monosaccharide_name):
            return cls.glycan_composition.any(
                cls.GlycanCompositionAssociation.base_type == monosaccharide_name)

        model.with_monosaccharide = classmethod(with_monosaccharide)

        def glycan_composition_extents(cls, session, filter_fn=lambda q: q):
            q = session.query(
                cls.GlycanCompositionAssociation.base_type,
                func.min(cls.GlycanCompositionAssociation.count),
                func.max(cls.GlycanCompositionAssociation.count),
                ).group_by(cls.GlycanCompositionAssociation.base_type)
            return filter_fn(q)

        model.glycan_composition_extents = classmethod(glycan_composition_extents)

        # This hack is necessary to make inner class locatable to the pickle
        # machinery. A possible alternative solution is to define a single class
        # once and use multiple tables that are mapped to it.
        setattr(
            sys.modules[__name__],
            MonosaccharideBaseCounter.__name__,
            MonosaccharideBaseCounter)
    return model


def with_glycan_composition(attr_name):
    def decorator(model):
        has_glycan_composition(model, attr_name)
        return model
    return decorator


# For cases where event handlers must be set without creating a new
# table such as with inheritance.
def has_glycan_composition_listener(attr):
    @event.listens_for(attr, "set")
    def convert_composition(target, value, oldvalue, initiator):
        if value == "{}" or value is None:
            return
        for k, v in glycan_composition.parse(value).items():
            target.glycan_composition[k.name()] = v

    @event.listens_for(attr.class_, "load")
    def convert_composition_load(target, context):
        value = getattr(target, attr.prop.key)
        if value == "{}" or value is None:
            return
        for k, v in glycan_composition.parse(value).items():
            target.glycan_composition[k.name()] = v


class GlycanBase(object):
    id = Column(Integer, primary_key=True, autoincrement=True)
    calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    derivatization = Column(Unicode(64), index=True)
    reduction = Column(Unicode(64))
    name = Column(Unicode(64))
    composition = Column(Unicode(128), index=True)

    @declared_attr
    def hypothesis_id(self):
        return Column(Integer, ForeignKey("Hypothesis.id"))

    @hybrid_method
    def from_hypothesis(self, hypothesis_id):
        return self.hypothesis_id == hypothesis_id

    def __repr__(self):
        rep = "<{self.__class__.__name__} {self.composition}>".format(self=self)
        return rep


@with_glycan_composition("composition")
class StructureMotif(GlycanBase, HasReferenceAccessionNumber, Base):
    __tablename__ = "StructureMotif"
    canonical_sequence = Column(Unicode(256), index=True)
    motif_class = Column(Unicode(64), index=True)
    is_core_motif = Column(Boolean)
    _structure = None

    def structure(self):
        if self._structure is not None:
            return self._structure
        self._structure = glycoct_parser(self.canonical_sequence)
        return self._structure

    def matches(self, other):
        motif = self.structure()
        target = other.structure()
        ix = subtree_search.subtree_of(motif, target, exact=True)
        if self.is_core_motif:
            return ix == 1
        else:
            return ix is not None

    @classmethod
    def get(cls, session, name, canonical_sequence, motif_class=None, is_core_motif=None):
        obj, made = get_or_create(
            session, cls, name=name, canonical_sequence=canonical_sequence,
            motif_class=motif_class, is_core_motif=is_core_motif)
        return obj

    def __repr__(self):
        rep = "<{self.__class__.__name__} {self.name}\n{self.canonical_sequence}>".format(self=self)
        return rep

    @classmethod
    def initialize(cls, session):
        for name, motif in glypy.motifs.items():
            inst = cls.get(
                session, name=name, canonical_sequence=str(motif),
                motif_class=motif.motif_class, is_core_motif=motif.is_core_motif)
            session.add(inst)
        session.commit()

Namespace.initialization_list.append(StructureMotif.initialize)


@with_glycan_composition("composition")
class TheoreticalGlycanComposition(GlycanBase, HasTaxonomy, HasReferenceAccessionNumber, Base):
    __tablename__ = "TheoreticalGlycanComposition"

    structures = relationship("TheoreticalGlycanStructure", lazy='dynamic')

    __mapper_args__ = {
        'polymorphic_identity': u'TheoreticalGlycanComposition',
        "concrete": True
    }


@with_glycan_composition("composition")
class TheoreticalGlycanStructure(GlycanBase, HasTaxonomy, HasReferenceAccessionNumber, Base):
    __tablename__ = "TheoreticalGlycanStructure"

    composition_reference_id = Column(Integer, ForeignKey(TheoreticalGlycanComposition.id), index=True)

    glycoct = Column(Unicode(256), index=True)
    _fragments = Column(MutableDict.as_mutable(PickleType))

    ms1_score = Column(Numeric(7, 6, asdecimal=False), index=True)
    volume = Column(Numeric(12, 6, asdecimal=False))

    peak_group_matches = relationship("PeakGroupMatch", secondary=lambda: TheoreticalGlycanStructureToPeakGroupMatch, lazy="dynamic")


    def fragments(self, kind='BY'):
        if self._fragments is None:
            self._fragments = {}
        ion_types = map(''.join, map(sorted, kind))
        return (f for f in self._fragments.values()
                if ''.join(
                    sorted(crossring_pattern.sub("", f.kind))) in (ion_types))

    def structure(self):
        return glycoct_parser(self.glycoct)

    motifs = relationship(StructureMotif, secondary=lambda: TheoreticalGlycanStructureToMotifTable)

    @classmethod
    def with_motif(cls, name):
        return cls.motifs.any(StructureMotif.name.like(name))

    __mapper_args__ = {
        'polymorphic_identity': u'TheoreticalGlycanStructure',
        "concrete": True
    }

    def __repr__(self):
        rep = "<{self.__class__.__name__}\n{self.glycoct}>".format(self=self)
        return rep


TheoreticalGlycanStructureToPeakGroupMatch = Table(
    "TheoreticalGlycanStructureToPeakGroupMatch", Base.metadata,
    Column("glycan_id", Integer, ForeignKey("TheoreticalGlycanStructure.id"), index=True),
    Column("peak_group_match_id", Integer, ForeignKey("PeakGroupMatch.id"), index=True)
    )

TheoreticalGlycanStructureToMotifTable = Table(
    "TheoreticalGlycanStructureToMotifTable", Base.metadata,
    Column("glycan_id", Integer, ForeignKey("TheoreticalGlycanStructure.id"), index=True),
    Column("motif_id", Integer, ForeignKey("StructureMotif.id"), index=True)
    )
