from sqlalchemy.ext.mutable import Mutable
from sqlalchemy import Table, Column, Integer, ForeignKey, Unicode, ForeignKeyConstraint
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declared_attr

from .base import Base
from glycresoft_sqlalchemy.utils.database_utils import get_or_create


class MutableDict(Mutable, dict):
    @classmethod
    def coerce(cls, key, value):
        if not isinstance(value, MutableDict):
            if isinstance(value, dict):
                return MutableDict(value)
            return Mutable.coerce(key, value)
        else:
            return value

    def __delitem(self, key):
        dict.__delitem__(self, key)
        self.changed()

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self.changed()

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, state):
        self.update(self)


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
            Column("entity_id", Integer, ForeignKey("%s.id" % cls.__tablename__, ondelete="CASCADE"), primary_key=True))
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
            Column("entity_id", Integer, ForeignKey("%s.id" % cls.__tablename__, ondelete="CASCADE"), primary_key=True),
            ForeignKeyConstraint(
                ["accession_code", "database_id"],
                ["ReferenceAccessionNumber.id", "ReferenceAccessionNumber.database_id"]))
        cls.ReferenceAccessionAssocationTable = reference_number_association
        return relationship(ReferenceAccessionNumber, secondary=reference_number_association)
