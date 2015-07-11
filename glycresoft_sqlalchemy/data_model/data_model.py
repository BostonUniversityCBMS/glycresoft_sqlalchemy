# -*- coding: utf-8 -*-

from glycresoft_ms2_classification.structure import sequence


from sqlalchemy.ext.declarative import declarative_base, AbstractConcreteBase, declared_attr
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .generic import MutableDict, MutableList
from .base import Base
from .glycomics import TheoreticalGlycanComposition as Glycan


class Hypothesis(Base):
    '''
    Represents a database of theoretical sequences to search against.
    '''
    __tablename__ = "Hypothesis"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128), default=u"")
    proteins = relationship("Protein", backref=backref("hypothesis", order_by=id),
                            collection_class=attribute_mapped_collection('name'))

    glycans = relationship(Glycan, backref=backref("hypothesis", order_by=id),
                           collection_class=attribute_mapped_collection('name'))

    is_decoy = Column(Boolean, default=False)
    parameters = Column(MutableDict.as_mutable(PickleType), default={})

    hypothesis_type = Column(Unicode(56), index=True)

    def __repr__(self):
        return "<Hypothesis {0} {1} {2} proteins {3} glycans>".format(
            self.id, self.name, len(self.proteins), len(self.glycans))

    __mapper_args__ = {
        'polymorphic_identity': u'Hypothesis',
        'polymorphic_on': hypothesis_type,
    }


class MS1GlycopeptideHypothesis(Hypothesis):
    __tablename__ = "MS1GlycopeptideHypothesis"

    id = Column(Integer, ForeignKey(Hypothesis.id), primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity': u'MS1GlycopeptideHypothesis',
    }


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
        "GlycopeptideMatch", lazy='dynamic', backref=backref('protein', order_by=id))

    def __repr__(self):
        return "<Protein {0} {1} {2} {3}...>".format(
            self.id, self.name, (self.glycopeptide_matches.count()),
            self.protein_sequence[:20] if self.protein_sequence is not None else "")


class PeptideBase(AbstractConcreteBase, Base):
    # __tablename__ = "PeptideBase"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_type = Column(Unicode(20))

    @declared_attr
    def protein_id(self):
        return Column(Integer, ForeignKey(Protein.id), index=True)

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

    @hybrid_method
    def from_hypothesis(self, hypothesis_id):
        return self.protein.hypothesis_id == hypothesis_id

    @from_hypothesis.expression
    def from_hypothesis(self, hypothesis_id):
        return self.protein_id == Protein.id & Protein.hypothesis_id == hypothesis_id

    def __len__(self):
        if self.sequence_length is not None:
            return self.sequence_length
        else:
            return sequence.sequence_length(self.base_peptide_sequence)

    @property
    def most_detailed_sequence(self):
        try:
            s = self.glycopeptide_sequence
            if s is not None and s != "":
                return s
        except:
            try:
                s = self.modified_peptide_sequence
                if s is not None and s != "":
                    return s
            except:
                return self.base_peptide_sequence

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
        'polymorphic_on': sequence_type,
    }


TheoreticalGlycopeptideGlycanAssociation = Table(
    "TheoreticalGlycopeptideGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("TheoreticalGlycopeptide.id")),
    Column("glycan_id", Integer, ForeignKey(Glycan.id)))


class TheoreticalGlycopeptide(PeptideBase):
    __tablename__ = "TheoreticalGlycopeptide"

    id = Column(Integer, primary_key=True)
    glycans = relationship(Glycan, secondary=TheoreticalGlycopeptideGlycanAssociation,
                           backref='glycopeptides', lazy='dynamic')
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

    def fragments(self, kind=('ox', 'b', 'y', 'gb', 'gy', 'stub')):
        '''
        A unified API for accessing lazy sequences of molecular fragments

        Parameters
        ----------
        kind: sequence of str
            A sequence of sigil strings for the different types of fragments
            to generate. These sigils follow the standard nomenclature of bc-yz for peptides
            and the BC-YZ nomenclature for glycans, as well as special combinations for oxonium
            ions and stub glycopeptide ions.

        Yields
        ------
        dict: A simple mapping object that defines a fragment
        '''
        kind = set(kind)
        if 'ox' in kind:
            for ox in self.oxonium_ions:
                yield ox
        if 'b' in kind:
            for b in self.bare_b_ions:
                yield b
        if 'y' in kind:
            for y in self.bare_y_ions:
                yield y
        if 'gb' in kind:
            for bg in self.glycosylated_b_ions:
                yield bg
        if 'gy' in kind:
            for yg in self.glycosylated_y_ions:
                yield yg
        if 'stub' in kind:
            for stub in self.stub_ions:
                yield stub

    __mapper_args__ = {
        'polymorphic_identity': u'TheoreticalGlycopeptide',
        "concrete": True
    }

    def __repr__(self):
        rep = "<TheoreticalGlycopeptide {} {}>".format(self.glycopeptide_sequence, self.observed_mass)
        return rep
