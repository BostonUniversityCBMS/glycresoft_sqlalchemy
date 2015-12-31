import re
from collections import OrderedDict

from sqlalchemy.ext.baked import bakery
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref, make_transient, Query
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy import (PickleType, Numeric, Unicode, Table, bindparam,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.exc import DetachedInstanceError
from ..generic import MutableDict, MutableList
from ..base import Base
from ..glycomics import (
    with_glycan_composition, TheoreticalGlycanCombination, has_glycan_composition_listener)

from ...structure import sequence


peptide_bakery = bakery()


class Protein(Base):
    __tablename__ = "Protein"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_sequence = Column(UnicodeText, default=u"")
    name = Column(Unicode(128), default=u"", index=True)
    other = Column(MutableDict.as_mutable(PickleType))
    hypothesis_id = Column(Integer, ForeignKey("Hypothesis.id", ondelete="CASCADE"))
    glycosylation_sites = Column(MutableList.as_mutable(PickleType))

    def __repr__(self):
        return "<Protein {0} {1} {2} {3}...>".format(
            self.id, self.name, self.glycosylation_sites,
            self.protein_sequence[:20] if self.protein_sequence is not None else "")

    def to_json(self, full=False):
        d = OrderedDict((
            ('id', self.id),
            ('name', self.name),
            ("glycosylation_sites", list(self.glycosylation_sites)),
            ('other', self.other)
        ))
        if full:
            d.update({
                "protein_sequence": self.protein_sequence
                })
            for k, v in self.__dict__.items():
                if isinstance(v, Query):
                    d[k + '_count'] = v.count()
        return d


def _convert_class_name_to_collection_name(name):
    parts = re.split(r"([A-Z]+[a-z]+)", name)
    parts = [p.lower() for p in parts if p]
    return '_'.join(parts) + 's'


class PeptideBase(object):

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_type = Column(Unicode(45))

    @declared_attr
    def protein_id(self):
        return Column(Integer, ForeignKey(Protein.id, ondelete="CASCADE"), index=True)

    @declared_attr
    def protein(self):
        if not hasattr(self, "__collection_name__"):
            name = _convert_class_name_to_collection_name(self.__name__)
        else:
            name = self.__collection_name__
        return relationship(Protein, backref=backref(name, lazy='dynamic'))

    calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

    count_glycosylation_sites = Column(Integer)
    count_missed_cleavages = Column(Integer)

    start_position = Column(Integer)
    end_position = Column(Integer)

    base_peptide_sequence = Column(Unicode(128), index=True)
    modified_peptide_sequence = Column(Unicode(128), index=True)

    sequence_length = Column(Integer, index=True)

    peptide_modifications = Column(Unicode(128))
    glycosylation_sites = Column(MutableList.as_mutable(PickleType))

    _query_ppm_tolerance_search_no_hypothesis_id = None
    _query_ppm_tolerance_search = None

    @classmethod
    def ppm_error_tolerance_search(cls, session, mass, tolerance, hypothesis_id=None):
        width = (mass * tolerance)
        lower = mass - width
        upper = mass + width
        if hypothesis_id is not None:
            if cls._query_ppm_tolerance_search is None:
                q = peptide_bakery(lambda session: session.query(cls).join(Protein))
                q += lambda q: q.filter(cls.calculated_mass.between(bindparam("lower"), bindparam("upper")))
                q += lambda q: q.filter(
                    cls.protein_id == Protein.id, Protein.hypothesis_id == bindparam("hypothesis_id"))
                cls._query_ppm_tolerance_search = q
            return cls._query_ppm_tolerance_search(session).params(
                lower=lower, upper=upper, hypothesis_id=hypothesis_id)
        else:
            if cls._query_ppm_tolerance_search_no_hypothesis_id is None:
                q = peptide_bakery(lambda session: session.query(cls))
                q += lambda q: q.filter(cls.calculated_mass.between(bindparam("lower"), bindparam("upper")))
                cls._query_ppm_tolerance_search_no_hypothesis_id = q
            return cls._query_ppm_tolerance_search_no_hypothesis_id(session).params(lower=lower, upper=upper)

    @hybrid_method
    def spans(self, point):
        return (self.start_position <= point) & (point < self.end_position)

    @hybrid_method
    def from_hypothesis(self, hypothesis_id):
        return self.protein.hypothesis_id == hypothesis_id

    @from_hypothesis.expression
    def from_hypothesis(self, hypothesis_id):
        return (self.protein_id == Protein.id) & (Protein.hypothesis_id == hypothesis_id)

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
                sites |= set(site - peptide.start_position for site in peptide.protein.glycosylation_sites
                             if peptide.start_position <= site < peptide.end_position)
        except AttributeError:
            pass
        return list(sites)

    def __hash__(self):
        return hash((self.most_detailed_sequence, self.protein_id))

    def __eq__(self, other):
        try:
            result = self.protein is None or other.protein is None or self.protein == other.protein
        except DetachedInstanceError:
            result = True
        result &= self.most_detailed_sequence == other.most_detailed_sequence
        return result

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "<{} {} {}>".format(self.__class__.__name__, self.id, self.most_detailed_sequence)


class GlycopeptideBase(PeptideBase):
    glycopeptide_sequence = Column(Unicode(128), index=True)

    glycan_mass = Column(Numeric(12, 6, asdecimal=False))
    glycan_composition_str = Column(Unicode(128), index=True)

    @declared_attr
    def glycan_combination_id(self):
        return Column(Integer, ForeignKey("TheoreticalGlycanCombination.id", ondelete="CASCADE"), index=True)

    @declared_attr
    def glycans(self):
        return relationship(TheoreticalGlycanCombination)


class NaivePeptide(PeptideBase, Base):
    __tablename__ = "NaivePeptide"

    id = Column(Integer, primary_key=True)

    def __repr__(self):
        return "<NaivePeptide {0} {1} {2} {3}...>".format(
            self.id, self.most_detailed_sequence, self.peptide_modifications,
            self.glycosylation_sites)


@with_glycan_composition("glycan_composition_str")
class TheoreticalGlycopeptideComposition(GlycopeptideBase, Base):
    __tablename__ = "TheoreticalGlycopeptideComposition"

    def __repr__(self):
        return "<{} {} {} {} {}>".format(
            self.__class__.__name__,
            self.id, self.glycopeptide_sequence, self.peptide_modifications, self.calculated_mass)

    __mapper_args__ = {
        'polymorphic_on': PeptideBase.sequence_type,
        'polymorphic_identity': u'TheoreticalGlycopeptideComposition',
    }


class InformedPeptide(PeptideBase, Base):
    __tablename__ = "InformedPeptide"
    id = Column(Integer, primary_key=True)
    peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
    peptide_score_type = Column(Unicode(56))
    theoretical_glycopeptides = relationship(
        "InformedTheoreticalGlycopeptideComposition",
        lazy="dynamic", backref="base_peptide")
    other = Column(PickleType)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedPeptide',
        "concrete": True
    }


class InformedTheoreticalGlycopeptideComposition(TheoreticalGlycopeptideComposition):
    __tablename__ = "InformedTheoreticalGlycopeptideComposition"

    id = Column(Integer, ForeignKey(TheoreticalGlycopeptideComposition.id, ondelete="CASCADE"), primary_key=True)

    peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
    other = Column(PickleType)
    base_peptide_id = Column(Integer, ForeignKey(InformedPeptide.id), index=True)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedTheoreticalGlycopeptideComposition',
    }

has_glycan_composition_listener(InformedTheoreticalGlycopeptideComposition.glycan_composition_str)


class HasMS1Information(object):
    ms1_score = Column(Numeric(12, 6, asdecimal=False), index=True)

    observed_mass = Column(Numeric(12, 6, asdecimal=False))
    ppm_error = Column(Numeric(12, 6, asdecimal=False))
    volume = Column(Numeric(12, 4, asdecimal=False))


class HasTheoreticalFragments(object):
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


@with_glycan_composition("glycan_composition_str")
class TheoreticalGlycopeptide(GlycopeptideBase, Base, HasMS1Information, HasTheoreticalFragments):
    __tablename__ = "TheoreticalGlycopeptide"

    base_composition_id = Column(Integer, ForeignKey("PeakGroupMatch.id"), index=True)

    def __repr__(self):
        rep = "<TheoreticalGlycopeptide {} {}>".format(self.glycopeptide_sequence, self.observed_mass)
        return rep
