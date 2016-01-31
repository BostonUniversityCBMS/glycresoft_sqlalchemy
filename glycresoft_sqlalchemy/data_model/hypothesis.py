from collections import OrderedDict

from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection

from sqlalchemy import (
    Column, Integer, ForeignKey, PickleType, Unicode, Boolean)

from .generic import ParameterStore, HasUniqueName
from .base import Base

from glycresoft_sqlalchemy.utils import classproperty


class Hypothesis(Base, ParameterStore, HasUniqueName):
    '''
    Represents a database of theoretical sequences to search against.
    '''
    __tablename__ = "Hypothesis"

    id = Column(Integer, primary_key=True, autoincrement=True)
    proteins = relationship("Protein", backref=backref("hypothesis", order_by=id),
                            collection_class=attribute_mapped_collection('name'),
                            cascade="delete")
    proteins_query = relationship("Protein", lazy='dynamic')

    glycans = relationship("TheoreticalGlycanComposition", backref=backref("hypothesis"),
                           lazy="dynamic",
                           cascade="delete")

    @declared_attr
    def glycans(self):
        return relationship("TheoreticalGlycanComposition", backref=backref("hypothesis"),
                            lazy="dynamic",
                            cascade="delete")

    is_decoy = Column(Boolean, default=False)
    hypothesis_type = Column(Unicode(56), index=True)

    sample_matches = relationship(
        "HypothesisSampleMatch", foreign_keys='HypothesisSampleMatch.target_hypothesis_id')

    def __init__(self, **kwargs):
        kwargs.setdefault('parameters', {})
        super(Hypothesis, self).__init__(**kwargs)

    def to_json(self, *args, **kwargs):
        serialize_parameters = dict(self.parameters)
        if "decoys" in serialize_parameters:
            decoy_ops = serialize_parameters['decoys']
            for decoy in decoy_ops:
                decoy['type'] = str(decoy['type'])
            serialize_parameters['decoys'] = decoy_ops
        d = OrderedDict((
            ("id", self.id),
            ("name", str(self.name)),
            ("parameters", serialize_parameters),
            ("hypothesis_type", str(self.hypothesis_type)),
            ("is_decoy", self.is_decoy)
        ))
        return d

    def __repr__(self):
        return "<{} {} {} {} proteins {} glycans>".format(
            self.__class__.__name__,
            self.id, self.name, len(self.proteins), self.glycans.count())

    __mapper_args__ = {
        'polymorphic_identity': u'Hypothesis',
        'polymorphic_on': hypothesis_type,
    }

    def summarize(self):
        fields = self.to_json()
        fields['proteins'] = [
            p.to_json() for p in self.proteins_query
        ]
        return fields

    def decoy_hypothesis_id(self):
        try:
            return self.parameters['decoys'][0]['hypothesis_id']
        except:
            return None

    ms_level = 0


class MS1GlycanHypothesis(Hypothesis):
    __tablename__ = "MS1GlycanHypothesis"

    id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity': u'MS1GlycanHypothesis',
    }

    @classproperty
    def theoretical_structure_type(cls):
        from . import glycomics
        return glycomics.TheoreticalGlycanComposition

    ms_level = 1


class MS2GlycanHypothesis(Hypothesis):
    __tablename__ = "MS2GlycanHypothesis"

    id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)
    ms1_source_hypothesis_id = Column(Integer, ForeignKey("MS1GlycanHypothesisSampleMatch.id"))

    __mapper_args__ = {
        'polymorphic_identity': u'MS2GlycanHypothesis',
    }

    @classproperty
    def theoretical_structure_type(cls):
        from . import glycomics
        return glycomics.TheoreticalGlycanStructure

    ms_level = 2


class MS1GlycopeptideHypothesis(Hypothesis):

    __tablename__ = "MS1GlycopeptideHypothesis"

    id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity': u'MS1GlycopeptideHypothesis',
    }

    @classproperty
    def theoretical_structure_type(cls):
        from .sequence_model import peptide
        return peptide.TheoreticalGlycopeptideComposition

    ms_level = 1


class ExactMS1GlycopeptideHypothesis(MS1GlycopeptideHypothesis):
    __tablename__ = "ExactMS1GlycopeptideHypothesis"

    id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity': u'ExactMS1GlycopeptideHypothesis',
    }


class MS2GlycopeptideHypothesis(Hypothesis):
    __tablename__ = "MS2GlycopeptideHypothesis"

    id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)
    ms1_source_hypothesis_id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id))

    __mapper_args__ = {
        'polymorphic_identity': u'MS2GlycopeptideHypothesis',
    }

    @classproperty
    def theoretical_structure_type(cls):
        from .sequence_model import peptide
        return peptide.TheoreticalGlycopeptide

    ms_level = 2


class ExactMS2GlycopeptideHypothesis(MS2GlycopeptideHypothesis):
    __tablename__ = "ExactMS2GlycopeptideHypothesis"

    id = Column(Integer, ForeignKey(MS2GlycopeptideHypothesis.id, ondelete="CASCADE"), primary_key=True)
    ms1_source_hypothesis_id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id))

    __mapper_args__ = {
        'polymorphic_identity': u'ExactMS2GlycopeptideHypothesis',
    }

    @classproperty
    def theoretical_structure_type(cls):
        from .sequence_model import peptide
        return peptide.TheoreticalGlycopeptide


class HypothesisToGlycanHypothesis(Base):
    __tablename__ = "HypothesisToGlycanHypothesis"

    glycan_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id), index=True, primary_key=True)
    linked_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id), index=True, primary_key=True)
