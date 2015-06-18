import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy import PickleType, Numeric, Unicode, create_engine, Column, Integer, ForeignKey, UnicodeText
from .json_type import JSONType
Base = declarative_base()


class Experiment(Base):
    __tablename__ = "Experiment"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    proteins = relationship("Protein", backref=backref("experiment", order_by=id),
                            collection_class=attribute_mapped_collection('name'))

    glycans = relationship("Glycan", backref=backref("experiment", order_by=id),
                            collection_class=attribute_mapped_collection('name'))
    parameters = relationship("ExperimentParameter", backref=backref("experiment", order_by=id),
                              collection_class=attribute_mapped_collection('name'),
                              cascade="all, delete-orphan")


    def __repr__(self):
        return "<Experiment {0} {1} {2} proteins {3} glycans>".format(
            self.id, self.name, len(self.proteins), len(self.glycans))


class ExperimentParameter(Base):
    __tablename__ = "ExperimentParameter"

    experiment_id = Column(Integer, ForeignKey("Experiment.id"))
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"", index=True)
    value = Column(PickleType)

    @classmethod
    def parameters_for(cls, param_dict, experiment_id, session):
        for k, v in param_dict.items():
            session.add(cls(name=k, value=v, experiment_id=experiment_id))
        session.commit()

    def __repr__(self):
        return "<{name} -> {value}>".format(**self.__dict__)


class Glycan(Base):
    __tablename__ = "Glycan"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    mass = Column(Numeric(10, 6, asdecimal=False), default=-1.0)
    composition = Column(Unicode(128), default=u"")
    experiment_id = Column(Integer, ForeignKey("Experiment.id"))


class Protein(Base):
    __tablename__ = "Protein"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_sequence = Column(UnicodeText, default=u"")
    name = Column(Unicode(128), default=u"", index=True)
    experiment_id = Column(Integer, ForeignKey("Experiment.id"))
    glycosylation_sites = Column(PickleType)
    theoretical_glycopeptides = relationship("TheoreticalGlycopeptide", backref=backref('protein', order_by=id), lazy='dynamic')
    glycopeptide_matches = relationship("GlycopeptideMatch", backref=backref('protein', order_by=id), lazy='dynamic')

    def __repr__(self):
        return "<Protein {0} {1} {2} {3}...>".format(
            self.id, self.name, (self.glycopeptide_matches.count()), self.protein_sequence[:20])


class TheoreticalGlycopeptide(Base):
    __tablename__ = "TheoreticalGlycopeptide"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("Protein.id"), index=True)
    # glycan_id = Column(Integer, ForeignKey("Glycan.id"), index=True)
    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)

    observed_mass = Column(Numeric(10, 6, asdecimal=False))
    calculated_mass = Column(Numeric(10, 6, asdecimal=False), index=True)
    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    ppm_error = Column(Numeric(10, 6, asdecimal=False))
    volume = Column(Numeric(10, 6, asdecimal=False))

    count_glycosylation_sites = Column(Integer)
    count_missed_cleavages = Column(Integer)

    start_position = Column(Integer)
    end_position = Column(Integer)

    base_peptide_sequence = Column(Unicode(128), index=True)
    modified_peptide_sequence = Column(Unicode(128))
    glycopeptide_sequence = Column(Unicode(128), index=True)
    glycan_composition_str = Column(Unicode(128), index=True)
    sequence_length = Column(Integer, index=True)

    peptide_modifications = Column(Unicode(128))

    glycosylation_sites = Column(PickleType)

    oxonium_ions = Column(PickleType)
    stub_ions = Column(PickleType)

    bare_b_ions = Column(PickleType)
    glycosylated_b_ions = Column(PickleType)

    bare_y_ions = Column(PickleType)
    glycosylated_y_ions = Column(PickleType)

    def __repr__(self):
        rep = "<TheoreticalGlycopeptide {} {}>".format(self.glycopeptide_sequence, self.observed_mass)
        return rep


class GlycopeptideMatch(Base):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    theoretical_glycopeptide = Column(Integer, ForeignKey("TheoreticalGlycopeptide.id"), index=True)
    protein_id = Column(Integer, ForeignKey("Protein.id"), index=True)
    # glycan_id = Column(Integer, ForeignKey("Glycan.id"), index=True)
    ms1_score = Column(Numeric(10, 6, asdecimal=False))
    ms2_score = Column(Numeric(10, 6, asdecimal=False))

    observed_mass = Column(Numeric(10, 6, asdecimal=False))
    calculated_mass = Column(Numeric(10, 6, asdecimal=False))
    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    ppm_error = Column(Numeric(10, 6, asdecimal=False))
    volume = Column(Numeric(10, 6, asdecimal=False))

    count_glycosylation_sites = Column(Integer)
    count_missed_cleavages = Column(Integer)

    total_bare_b_ions_possible = Column(Integer)
    total_glycosylated_b_ions_possible = Column(Integer)
    total_bare_y_ions_possible = Column(Integer)
    total_glycosylated_y_ions_possible = Column(Integer)
    percent_bare_b_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_bare_y_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_glycosylated_b_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_glycosylated_y_ion_coverage = Column(Numeric(10, 6, asdecimal=False))

    start_position = Column(Integer)
    end_position = Column(Integer)

    base_peptide_sequence = Column(Unicode(128))
    modified_peptide_sequence = Column(Unicode(128))
    glycopeptide_sequence = Column(Unicode(128))
    sequence_length = Column(Integer, index=True)

    peptide_modifications = Column(Unicode(128))
    glycan_composition_str = Column(Unicode(128))

    scan_id_range = Column(PickleType)
    first_scan = Column(Integer)
    last_scan = Column(Integer)

    glycosylation_sites = Column(PickleType)

    oxonium_ions = Column(PickleType)
    stub_ions = Column(PickleType)

    bare_b_ions = Column(PickleType)
    glycosylated_b_ions = Column(PickleType)

    bare_y_ions = Column(PickleType)
    glycosylated_y_ions = Column(PickleType)

    glycan_composition = Column(PickleType)

    mean_coverage = Column(Numeric(10, 6, asdecimal=False))
    mean_hexnac_coverage = Column(Numeric(10, 6, asdecimal=False))

    def __repr__(self):
        rep = "<GlycopeptideMatch {} {} {}>".format(self.glycopeptide_sequence, self.ms2_score, self.observed_mass)
        return rep


class SpectrumMatch(Base):
    __tablename__ = "SpectrumMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    spectrum_id = Column(Integer)
    theoretical_glycopeptide_id = Column(Integer, ForeignKey("TheoreticalGlycopeptide.id"), index=True)
    peak_match_map = Column(PickleType)


class DatabaseManager(object):
    connect_args = {"timeout": 30}
    database_uri_prefix = "sqlite:///"

    def __init__(self, path, clear=False):
        if clear:
            try:
                os.remove(path)
            except:
                pass
        self.path = path

    def connect(self):
        return create_engine("{}{}".format(self.database_uri_prefix, self.path), connect_args=self.connect_args)

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
