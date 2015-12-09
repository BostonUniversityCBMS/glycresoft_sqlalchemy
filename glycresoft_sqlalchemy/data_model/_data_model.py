# # -*- coding: utf-8 -*-

# from collections import OrderedDict


# from sqlalchemy.ext.declarative import declared_attr
# from sqlalchemy import event
# from sqlalchemy.orm import relationship, backref, make_transient, Query
# from sqlalchemy.orm.collections import attribute_mapped_collection
# from sqlalchemy.ext.hybrid import hybrid_method
# from sqlalchemy import (PickleType, Numeric, Unicode, Table,
#                         Column, Integer, ForeignKey, UnicodeText, Boolean)
# from sqlalchemy.orm.exc import DetachedInstanceError
# from .generic import MutableDict, MutableList
# from .base import Base
# from .glycomics import TheoreticalGlycanComposition as Glycan, with_glycan_composition
# from . import glycomics

# from ..structure import sequence, fragment


# class Hypothesis(Base):
#     '''
#     Represents a database of theoretical sequences to search against.
#     '''
#     __tablename__ = "Hypothesis"

#     id = Column(Integer, primary_key=True, autoincrement=True)
#     name = Column(Unicode(128), default=u"")
#     proteins = relationship("Protein", backref=backref("hypothesis", order_by=id),
#                             collection_class=attribute_mapped_collection('name'),
#                             cascade="delete")
#     proteins_query = relationship("Protein", lazy='dynamic')

#     glycans = relationship(Glycan, backref=backref("hypothesis", order_by=id),
#                            lazy="dynamic",
#                            cascade="delete")

#     is_decoy = Column(Boolean, default=False)
#     parameters = Column(MutableDict.as_mutable(PickleType))

#     hypothesis_type = Column(Unicode(56), index=True)

#     sample_matches = relationship(
#         "HypothesisSampleMatch", foreign_keys='HypothesisSampleMatch.target_hypothesis_id')

#     def __init__(self, **kwargs):
#         kwargs.setdefault('parameters', {})
#         super(Hypothesis, self).__init__(**kwargs)

#     def to_json(self, *args, **kwargs):
#         serialize_parameters = dict(self.parameters)
#         if "decoys" in serialize_parameters:
#             decoy_ops = serialize_parameters['decoys']
#             for decoy in decoy_ops:
#                 decoy['type'] = str(decoy['type'])
#             serialize_parameters['decoys'] = decoy_ops
#         d = OrderedDict((
#             ("id", self.id),
#             ("name", str(self.name)),
#             ("parameters", serialize_parameters),
#             ("hypothesis_type", str(self.hypothesis_type)),
#             ("is_decoy", self.is_decoy)
#         ))
#         return d

#     def __repr__(self):
#         return "<{} {} {} {} proteins {} glycans>".format(
#             self.__class__.__name__,
#             self.id, self.name, len(self.proteins), self.glycans.count())

#     __mapper_args__ = {
#         'polymorphic_identity': u'Hypothesis',
#         'polymorphic_on': hypothesis_type,
#     }

#     def summarize(self):
#         fields = self.to_json()
#         fields['proteins'] = [
#             p.to_json() for p in self.proteins_query
#         ]
#         return fields


# class Protein(Base):
#     __tablename__ = "Protein"

#     id = Column(Integer, primary_key=True, autoincrement=True)
#     protein_sequence = Column(UnicodeText, default=u"")
#     name = Column(Unicode(128), default=u"", index=True)
#     other = Column(MutableDict.as_mutable(PickleType))
#     hypothesis_id = Column(Integer, ForeignKey("Hypothesis.id", ondelete="CASCADE"))
#     glycosylation_sites = Column(MutableList.as_mutable(PickleType))

#     theoretical_glycopeptides = relationship(
#         "TheoreticalGlycopeptide", backref=backref('protein', order_by=id), lazy='dynamic')

#     glycopeptide_matches = relationship(
#         "GlycopeptideMatch", lazy='dynamic', backref=backref('protein', order_by=id))

#     def __repr__(self):
#         return "<Protein {0} {1} {2} {3}...>".format(
#             self.id, self.name, (self.glycopeptide_matches.count()),
#             self.protein_sequence[:20] if self.protein_sequence is not None else "")

#     def to_json(self, full=False):
#         d = OrderedDict((
#             ('id', self.id),
#             ('name', self.name),
#             ("glycosylation_sites", list(self.glycosylation_sites)),
#             ('other', self.other)
#         ))
#         if full:
#             d.update({
#                 "protein_sequence": self.protein_sequence
#                 })
#             for k, v in self.__dict__.items():
#                 if isinstance(v, Query):
#                     d[k + '_count'] = v.count()
#         return d


# class PeptideBase(object):

#     id = Column(Integer, primary_key=True, autoincrement=True)
#     sequence_type = Column(Unicode(45))

#     @declared_attr
#     def protein_id(self):
#         return Column(Integer, ForeignKey(Protein.id, ondelete="CASCADE"), index=True)

#     calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

#     count_glycosylation_sites = Column(Integer)
#     count_missed_cleavages = Column(Integer)

#     start_position = Column(Integer)
#     end_position = Column(Integer)

#     base_peptide_sequence = Column(Unicode(128), index=True)
#     modified_peptide_sequence = Column(Unicode(128), index=True)

#     sequence_length = Column(Integer, index=True)

#     peptide_modifications = Column(Unicode(128))
#     glycosylation_sites = Column(MutableList.as_mutable(PickleType))

#     @hybrid_method
#     def spans(self, point):
#         return (self.start_position <= point) & (point < self.end_position)

#     @hybrid_method
#     def from_hypothesis(self, hypothesis_id):
#         return self.protein.hypothesis_id == hypothesis_id

#     @from_hypothesis.expression
#     def from_hypothesis(self, hypothesis_id):
#         return (self.protein_id == Protein.id) & (Protein.hypothesis_id == hypothesis_id)

#     def __len__(self):
#         if self.sequence_length is not None:
#             return self.sequence_length
#         else:
#             return sequence.sequence_length(self.base_peptide_sequence)

#     @property
#     def most_detailed_sequence(self):
#         try:
#             s = self.glycopeptide_sequence
#             if s is not None and s != "":
#                 return s
#         except:
#             try:
#                 s = self.modified_peptide_sequence
#                 if s is not None and s != "":
#                     return s
#             except:
#                 return self.base_peptide_sequence

#     @property
#     def n_glycan_sequon_sites(peptide):
#         sites = set(sequence.find_n_glycosylation_sequons(peptide.base_peptide_sequence))
#         try:
#             if peptide.protein is not None:
#                 sites |= set(site - peptide.start_position for site in peptide.protein.glycosylation_sites
#                              if peptide.start_position <= site < peptide.end_position)
#         except AttributeError:
#             pass
#         return list(sites)

#     def __hash__(self):
#         return hash((self.most_detailed_sequence, self.protein_id))

#     def __eq__(self, other):
#         try:
#             result = self.protein is None or other.protein is None or self.protein == other.protein
#         except DetachedInstanceError:
#             result = True
#         result &= self.most_detailed_sequence == other.most_detailed_sequence
#         return result

#     def __ne__(self, other):
#         return not (self == other)

#     def __repr__(self):
#         return "<{} {} {}>".format(self.__class__.__name__, self.id, self.most_detailed_sequence)


# TheoreticalGlycopeptideGlycanAssociation = Table(
#     "TheoreticalGlycopeptideGlycanAssociation", Base.metadata,
#     Column("peptide_id", Integer, ForeignKey("TheoreticalGlycopeptide.id", ondelete="CASCADE")),
#     Column("glycan_id", Integer, ForeignKey(Glycan.id, ondelete="CASCADE")))


# class GlycopeptideSequenceMS2Base(object):

#     ms1_score = Column(Numeric(12, 6, asdecimal=False), index=True)

#     observed_mass = Column(Numeric(12, 6, asdecimal=False))
#     glycan_mass = Column(Numeric(12, 6, asdecimal=False))
#     ppm_error = Column(Numeric(12, 6, asdecimal=False))
#     volume = Column(Numeric(12, 4, asdecimal=False))

#     glycopeptide_sequence = Column(Unicode(128), index=True)
#     glycan_composition_str = Column(Unicode(128), index=True)

#     oxonium_ions = Column(MutableList.as_mutable(PickleType))
#     stub_ions = Column(MutableList.as_mutable(PickleType))

#     bare_b_ions = Column(MutableList.as_mutable(PickleType))
#     glycosylated_b_ions = Column(MutableList.as_mutable(PickleType))

#     bare_y_ions = Column(MutableList.as_mutable(PickleType))
#     glycosylated_y_ions = Column(MutableList.as_mutable(PickleType))

#     @declared_attr
#     def glycan_combination_id(self):
#         return Column(Integer, ForeignKey("TheoreticalGlycanCombination.id", ondelete="CASCADE"), index=True)


# @with_glycan_composition("glycan_composition_str")
# class TheoreticalGlycopeptide(PeptideBase, Base, GlycopeptideSequenceMS2Base):
#     __tablename__ = "TheoreticalGlycopeptide"

#     id = Column(Integer, primary_key=True)
#     glycans = relationship(Glycan, secondary=TheoreticalGlycopeptideGlycanAssociation,
#                            backref='glycopeptides', lazy='dynamic',
#                            cascade="delete")

#     base_composition_id = Column(Integer, ForeignKey("PeakGroupMatch.id"), index=True)

#     def fragments(self, kind=('ox', 'b', 'y', 'gb', 'gy', 'stub')):
#         '''
#         A unified API for accessing lazy sequences of molecular fragments

#         Parameters
#         ----------
#         kind: sequence of str
#             A sequence of sigil strings for the different types of fragments
#             to generate. These sigils follow the standard nomenclature of bc-yz for peptides
#             and the BC-YZ nomenclature for glycans, as well as special combinations for oxonium
#             ions and stub glycopeptide ions.

#         Yields
#         ------
#         dict: A simple mapping object that defines a fragment
#         '''
#         kind = set(kind)
#         if 'ox' in kind:
#             for ox in self.oxonium_ions:
#                 yield ox
#         if 'b' in kind:
#             for b in self.bare_b_ions:
#                 yield b
#         if 'y' in kind:
#             for y in self.bare_y_ions:
#                 yield y
#         if 'gb' in kind:
#             for bg in self.glycosylated_b_ions:
#                 yield bg
#         if 'gy' in kind:
#             for yg in self.glycosylated_y_ions:
#                 yield yg
#         if 'stub' in kind:
#             for stub in self.stub_ions:
#                 yield stub

#     def __repr__(self):
#         rep = "<TheoreticalGlycopeptide {} {}>".format(self.glycopeptide_sequence, self.observed_mass)
#         return rep


# class MS1GlycanHypothesis(Hypothesis):
#     __tablename__ = "MS1GlycanHypothesis"

#     id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)

#     __mapper_args__ = {
#         'polymorphic_identity': u'MS1GlycanHypothesis',
#     }

#     theoretical_structure_type = glycomics.TheoreticalGlycanComposition


# class MS2GlycanHypothesis(Hypothesis):
#     __tablename__ = "MS2GlycanHypothesis"

#     id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)
#     ms1_source_hypothesis_id = Column(Integer, ForeignKey("MS1GlycanHypothesisSampleMatch.id"))

#     __mapper_args__ = {
#         'polymorphic_identity': u'MS2GlycanHypothesis',
#     }

#     theoretical_structure_type = glycomics.TheoreticalGlycanStructure


# class MS1GlycopeptideHypothesis(Hypothesis):
#     from . import proteomics

#     __tablename__ = "MS1GlycopeptideHypothesis"

#     id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)

#     __mapper_args__ = {
#         'polymorphic_identity': u'MS1GlycopeptideHypothesis',
#     }

#     theoretical_structure_type = proteomics.TheoreticalGlycopeptideComposition


# class ExactMS1GlycopeptideHypothesis(MS1GlycopeptideHypothesis):
#     __tablename__ = "ExactMS1GlycopeptideHypothesis"

#     id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id, ondelete="CASCADE"), primary_key=True)

#     __mapper_args__ = {
#         'polymorphic_identity': u'ExactMS1GlycopeptideHypothesis',
#     }


# class MS2GlycopeptideHypothesis(Hypothesis):
#     __tablename__ = "MS2GlycopeptideHypothesis"

#     id = Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), primary_key=True)
#     ms1_source_hypothesis_id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id))

#     __mapper_args__ = {
#         'polymorphic_identity': u'MS2GlycopeptideHypothesis',
#     }

#     theoretical_structure_type = TheoreticalGlycopeptide


# class ExactMS2GlycopeptideHypothesis(MS2GlycopeptideHypothesis):
#     __tablename__ = "ExactMS2GlycopeptideHypothesis"

#     id = Column(Integer, ForeignKey(MS2GlycopeptideHypothesis.id, ondelete="CASCADE"), primary_key=True)
#     ms1_source_hypothesis_id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesis.id))

#     __mapper_args__ = {
#         'polymorphic_identity': u'ExactMS2GlycopeptideHypothesis',
#     }

#     theoretical_structure_type = TheoreticalGlycopeptide


# class HypothesisToGlycanHypothesis(Base):
#     __tablename__ = "HypothesisToGlycanHypothesis"

#     glycan_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id), index=True, primary_key=True)
#     linked_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id), index=True, primary_key=True)


# class TheoreticalPeptideProductIon(Base):
#     __tablename__ = "TheoreticalPeptideProductIon"

#     id = Column(Integer, primary_key=True)
#     name = Column(Unicode(64), index=True)
#     calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

#     n_term = Column(Unicode(10))
#     c_term = Column(Unicode(10))
#     glycosylated = Column(Boolean)
#     sequence_index = Column(Integer, index=True)
#     ion_ladder = Column(Unicode(10), index=True)

#     theoretical_glycopeptide_id = Column(Integer, ForeignKey("TheoreticalGlycopeptide.id"), index=True)
#     theoretical_glycopeptide = relationship("TheoreticalGlycopeptide")

#     @classmethod
#     def from_fragment(cls, frag):
#         inst = cls(
#             name=frag.name,
#             calculated_mass=frag.mass,
#             n_term=frag.flanking_amino_acids[0].name,
#             c_term=frag.flanking_amino_acids[1].name,
#             glycosylated=('HexNAc' in frag.mod_dict),
#             sequence_index=frag.pos,
#             ion_ladder=frag.name[0]
#             )
#         return inst

#     def __repr__(self):
#         return "<TheoreticalPeptideProductIon %s %s-%s %0.3f>" % (
#             self.name, self.n_term, self.c_term, self.calculated_mass)


# class TheoreticalGlycopeptideStubIon(Base):
#     __tablename__ = "TheoreticalGlycopeptideStubIon"

#     id = Column(Integer, primary_key=True)
#     name = Column(Unicode(64), index=True)
#     calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
#     ion_ladder = Column(Unicode(10), index=True)
#     includes_peptide = Column(Boolean)

#     theoretical_glycopeptide_id = Column(Integer, ForeignKey("TheoreticalGlycopeptide.id"), index=True)
#     theoretical_glycopeptide = relationship("TheoreticalGlycopeptide")

#     def __repr__(self):
#         return "<TheoreticalGlycopeptideStubIon %s %s %0.3f>" % (
#             self.name, self.ion_ladder, self.calculated_mass)

#     @classmethod
#     def from_simple_fragment(cls, frag):
#         inst = cls(
#             name=frag.name,
#             calculated_mass=frag.mass,
#             ion_ladder=frag.kind,
#             includes_peptide=fragment.IonSeries(frag.kind).includes_peptide
#             )
#         return inst
