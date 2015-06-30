# -*- coding: utf-8 -*-

import os

from glycresoft_ms2_classification.structure import sequence


from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy import (PickleType, Numeric, Unicode, create_engine, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .generic import MutableDict, MutableList

Base = declarative_base()


class Experiment(Base):
    __tablename__ = "Experiment"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    # proteins = relationship("Protein", backref=backref("experiment", order_by=id),
    #                         collection_class=attribute_mapped_collection('name'))

    # glycans = relationship("Glycan", backref=backref("experiment", order_by=id),
    #                        collection_class=attribute_mapped_collection('name'))

    is_decoy = Column(Boolean, default=False)
    parameters = Column(MutableDict.as_mutable(PickleType), default={})

    def __repr__(self):
        return "<Experiment {0} {1} {2} proteins {3} glycans>".format(
            self.id, self.name, len(self.proteins), len(self.glycans))


class Hypothesis(Base):
    '''
    Represents a database of theoretical sequences to search against.
    '''
    __tablename__ = "Hypothesis"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    proteins = relationship("Protein", backref=backref("hypothesis", order_by=id),
                            collection_class=attribute_mapped_collection('name'))

    glycans = relationship("Glycan", backref=backref("hypothesis", order_by=id),
                           collection_class=attribute_mapped_collection('name'))

    is_decoy = Column(Boolean, default=False)
    parameters = Column(MutableDict.as_mutable(PickleType), default={})

    def __repr__(self):
        return "<Hypothesis {0} {1} {2} proteins {3} glycans>".format(
            self.id, self.name, len(self.proteins), len(self.glycans))


class Glycan(Base):
    __tablename__ = "Glycan"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    mass = Column(Numeric(10, 6, asdecimal=False), default=0)
    composition = Column(Unicode(128), default=u"")
    other = Column(MutableDict.as_mutable(PickleType))
    hypothesis_id = Column(Integer, ForeignKey("Hypothesis.id"))


class Protein(Base):
    __tablename__ = "Protein"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_sequence = Column(UnicodeText, default=u"")
    name = Column(Unicode(128), default=u"", index=True)
    other = Column(MutableDict.as_mutable(PickleType))
    hypothesis_id = Column(Integer, ForeignKey("Hypothesis.id"))
    glycosylation_sites = Column(MutableList.as_mutable(PickleType))

    theoretical_glycopeptides = relationship(
        "TheoreticalGlycopeptide", backref=backref('protein', order_by=id), lazy='dynamic')

    glycopeptide_matches = relationship(
        "GlycopeptideMatch", lazy='dynamic')

    def __repr__(self):
        return "<Protein {0} {1} {2} {3}...>".format(
            self.id, self.name, (self.glycopeptide_matches.count()),
            self.protein_sequence[:20] if self.protein_sequence is not None else "")


class PeptideBase(Base):
    __tablename__ = "PeptideBase"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_type = Column(Unicode(20))

    protein_id = Column(Integer, ForeignKey(Protein.id), index=True)
    calculated_mass = Column(Numeric(10, 6, asdecimal=False), index=True)

    count_glycosylation_sites = Column(Integer)
    count_missed_cleavages = Column(Integer)

    start_position = Column(Integer)
    end_position = Column(Integer)

    base_peptide_sequence = Column(Unicode(128), index=True)
    modified_peptide_sequence = Column(Unicode(128), index=True)

    sequence_length = Column(Integer, index=True)

    peptide_modifications = Column(Unicode(128))
    glycosylation_sites = Column(MutableList.as_mutable(PickleType))

    @hybrid_method
    def spans(self, point):
        return (self.start_position <= point) & (point < self.end_position)

    @property
    def n_glycan_sequon_sites(peptide):
        sites = set(sequence.find_n_glycosylation_sequons(peptide.base_peptide_sequence))
        try:
            if peptide.protein is not None:
                sites |= set(site - peptide.start_position for site in peptide.parent.glycosylation_sites
                             if peptide.start_position <= site < peptide.end_position)
        except AttributeError:
            pass
        return list(sites)

    __mapper_args__ = {
        'polymorphic_identity': u'PeptideBase',
        'polymorphic_on': sequence_type
    }


PeptideGlycanAssociation = Table(
    "PeptideGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("PeptideBase.id")),
    Column("glycan_id", Integer, ForeignKey("Glycan.id")))


class TheoreticalGlycopeptide(PeptideBase):
    __tablename__ = "TheoreticalGlycopeptide"
    id = Column(Integer, ForeignKey(PeptideBase.id), primary_key=True)
    glycans = relationship(Glycan, secondary=PeptideGlycanAssociation, backref='glycopeptides', lazy='dynamic')
    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)

    observed_mass = Column(Numeric(10, 6, asdecimal=False))
    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    ppm_error = Column(Numeric(10, 6, asdecimal=False))
    volume = Column(Numeric(10, 6, asdecimal=False))

    glycopeptide_sequence = Column(Unicode(128), index=True)
    glycan_composition_str = Column(Unicode(128), index=True)

    oxonium_ions = Column(MutableList.as_mutable(PickleType))
    stub_ions = Column(MutableList.as_mutable(PickleType))

    bare_b_ions = Column(MutableList.as_mutable(PickleType))
    glycosylated_b_ions = Column(MutableList.as_mutable(PickleType))

    bare_y_ions = Column(MutableList.as_mutable(PickleType))
    glycosylated_y_ions = Column(MutableList.as_mutable(PickleType))

    __mapper_args__ = {
        'polymorphic_identity': u'TheoreticalGlycopeptide',
    }

    def __repr__(self):
        rep = "<TheoreticalGlycopeptide {} {}>".format(self.glycopeptide_sequence, self.observed_mass)
        return rep


class ConnectionManager(object):
    echo = False

    def __init__(self, database_uri, database_uri_prefix, connect_args=None):
        self.database_uri = database_uri
        self.database_uri_prefix = database_uri_prefix
        self.connect_args = connect_args or {}

    def connect(self):
        return create_engine(
            "{}{}".format(self.database_uri_prefix, self.database_uri),
            echo=self.echo,
            connect_args=self.connect_args)

    def clear(self):
        pass


class SQLiteConnectionManager(ConnectionManager):
    connect_args = {"timeout": 30}
    database_uri_prefix = "sqlite:///"

    def __init__(self, path, connect_args=None):
        if connect_args is None:
            connect_args = self.connect_args
        super(SQLiteConnectionManager, self).__init__(path, self.database_uri_prefix, connect_args)

    def clear(self):
        try:
            os.remove(self.database_uri)
        except:
            pass


class DatabaseManager(object):
    connection_manager_type = SQLiteConnectionManager

    def __init__(self, path, clear=False):
        self.connection_manager = self.connection_manager_type(path)
        if clear:
            self.connection_manager.clear()
        self.path = path

    def connect(self):
        return self.connection_manager.connect()

    def initialize(self, conn=None):
        if conn is None:
            conn = self.connect()
        Base.metadata.create_all(conn)

    def session(self, connection=None):
        if connection is None:
            connection = self.connect()
        return sessionmaker(bind=connection)()


def initialize(database_path):
    manager = DatabaseManager(database_path, clear=True)
    manager.initialize()
    return manager


def session(database_path):
    manager = DatabaseManager(database_path)
    return manager.session()
