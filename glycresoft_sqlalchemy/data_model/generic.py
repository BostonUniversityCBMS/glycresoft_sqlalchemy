from contextlib import contextmanager
import operator

from sqlalchemy.orm.session import object_session
from sqlalchemy.orm import validates

from sqlalchemy.ext.mutable import Mutable, MutableDict
from sqlalchemy import (
    Table, Column, Integer, ForeignKey, Unicode, ForeignKeyConstraint,
    PickleType)

from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declared_attr

from .base import Base
from glycresoft_sqlalchemy.utils.database_utils import get_or_create


class MutableList(Mutable, list):
    @classmethod
    def coerce(cls, key, value):
        if not isinstance(value, MutableList):
            if isinstance(value, list):
                return MutableList(value)
            value = Mutable.coerce(key, value)

        return value

    def __setitem__(self, key, value):
        old_value = list.__getitem__(self, key)
        for obj, key in self._parents.items():
            old_value._parents.pop(obj, None)

        list.__setitem__(self, key, value)
        for obj, key in self._parents.items():
            value._parents[obj] = key

        self.changed()

    def __getstate__(self):
        return list(self)

    def __setstate__(self, state):
        self[:] = state


class ParameterStore(object):
    parameters = Column(MutableDict.as_mutable(PickleType))

    def cache(self, obj, callpath, *args, **kwargs):
        self.parameters.setdefault("__cache__", {})
        cache = self.parameters['__cache__']
        getter = operator.attrgetter(callpath)
        fn = getter(obj)
        source_id, source_type = obj.id, obj.__class__.__name__

        try:
            value = cache[source_id, source_type, callpath, tuple(args), frozenset(kwargs.items())]
            return value
        except KeyError:
            value = fn(*args, **kwargs)
            cache[source_id, source_type, callpath, tuple(args), frozenset(kwargs.items())] = value
            return value

    @contextmanager
    def transaction(self):
        session = object_session(self)
        yield
        session.commit()


class Taxon(Base):
    __tablename__ = "Taxon"
    id = Column(Integer, primary_key=True)
    name = Column(Unicode(128), index=True)

    @classmethod
    def get(cls, session, id, name=None):
        obj, made = get_or_create(session, cls, id=id, name=name)
        return obj

    def __repr__(self):
        return "<Taxon {} {}>".format(self.id, self.name)


class HasTaxonomy(object):

    @declared_attr
    def taxa(cls):
        taxon_association = Table(
            "%s_Taxa" % cls.__tablename__,
            cls.metadata,
            Column("taxon_id", Integer, ForeignKey(Taxon.id, ondelete="CASCADE"), primary_key=True),
            Column("entity_id", Integer, ForeignKey(
                "%s.id" % cls.__tablename__, ondelete="CASCADE"), primary_key=True))
        cls.TaxonomyAssociationTable = taxon_association
        return relationship(Taxon, secondary=taxon_association)

    @classmethod
    def with_taxa(cls, ids):
        try:
            iter(ids)
            return cls.taxa.any(Taxon.id.in_(tuple(ids)))
        except:
            return cls.taxa.any(Taxon.id == ids)


class ReferenceDatabase(Base):
    __tablename__ = "ReferenceDatabase"
    id = Column(Integer, primary_key=True)
    name = Column(Unicode(128))
    url = Column(Unicode(128))

    @classmethod
    def get(cls, session, id=None, name=None, url=None):
        obj, made = get_or_create(session, cls, id=id, name=name, url=url)
        return obj

    def __repr__(self):
        return "<ReferenceDatabase {} {}>".format(self.id, self.name)


class ReferenceAccessionNumber(Base):
    __tablename__ = "ReferenceAccessionNumber"
    id = Column(Unicode(64), primary_key=True)
    database_id = Column(Integer, ForeignKey(ReferenceDatabase.id), primary_key=True)
    database = relationship(ReferenceDatabase)

    @classmethod
    def get(cls, session, id, database_id):
        obj, made = get_or_create(session, cls, id=id, database_id=database_id)
        return obj

    def __repr__(self):
        return "<ReferenceAccessionNumber {} {}>".format(self.id, self.database.name)


class HasReferenceAccessionNumber(object):
    @declared_attr
    def references(cls):
        reference_number_association = Table(
            "%s_ReferenceAccessionNumber" % cls.__tablename__,
            cls.metadata,
            Column("accession_code", Unicode(64), primary_key=True),
            Column("database_id", Integer, primary_key=True),
            Column(
                "entity_id", Integer, ForeignKey("%s.id" % cls.__tablename__, ondelete="CASCADE"), primary_key=True),
            ForeignKeyConstraint(
                ["accession_code", "database_id"],
                ["ReferenceAccessionNumber.id", "ReferenceAccessionNumber.database_id"]))
        cls.ReferenceAccessionAssocationTable = reference_number_association
        return relationship(ReferenceAccessionNumber, secondary=reference_number_association)


def find_by_name(session, model_class, name):
    return session.query(model_class).filter(model_class.name == name).first()


def make_unique_name(session, model_class, name):
    marked_name = name
    i = 1
    while find_by_name(session, model_class, marked_name) is not None:
        marked_name = "%s (%d)" % (name, i)
        i += 1
    return marked_name


class HasUniqueName(object):
    name = Column(Unicode(128), default=u"", unique=True)

    @classmethod
    def make_unique_name(cls, session, name):
        return make_unique_name(session, cls, name)

    @classmethod
    def find_by_name(cls, session, name):
        return find_by_name(session, cls, name)

    @validates("name")
    def ensure_unique_name(self, key, name):
        session = object_session(self)
        if session is not None:
            model_class = self.__class__
            name = make_unique_name(session, model_class, name)
            return name
        else:
            return name


_TemplateNumberStore = Table("_TemplateNumberStore", Base.metadata, Column("value", Integer))
