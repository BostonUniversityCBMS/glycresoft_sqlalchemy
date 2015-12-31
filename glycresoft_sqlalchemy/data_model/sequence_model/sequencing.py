from sqlalchemy.orm import relationship, backref
from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey, Table, UnicodeText, bindparam

from ..base import Base
from sqlalchemy.ext.baked import bakery

from glycresoft_sqlalchemy.structure.sequence_composition import SequenceComposition, Composition


sequencing_bakery = bakery()


class SequenceBuildingBlock(Base):
    __tablename__ = "SequenceBuildingBlock"

    id = Column(Integer, primary_key=True)
    name = Column(Unicode(120), index=True)
    mass = Column(Numeric(10, 6), index=True)
    hypothesis_id = Column(Integer, ForeignKey("Hypothesis.id"), index=True)


class SequenceSegment(Base):
    __tablename__ = "SequenceSegment"
    id = Column(Integer, primary_key=True)
    sequence = Column(Unicode(128))
    mass = Column(Numeric(12, 6), index=True)
    count_n_glycosylation = Column(Integer)


class AminoAcidComposition(Base):
    __tablename__ = "AminoAcidComposition"
    id = Column(Integer, primary_key=True)
    composition = Column(Unicode(400), index=True)
    mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    size = Column(Integer, index=True)
    count_n_glycosylation = Column(Integer)

    def __repr__(self):
        return "AminoAcidComposition({} {} {} {})".format(
            self.id, self.composition, self.mass, self.count_n_glycosylation)

    _query_ppm_tolerance_search = None

    @classmethod
    def ppm_error_tolerance_search(cls, session, mass, tolerance):
        width = (mass * tolerance)
        lower = mass - width
        upper = mass + width
        if cls._query_ppm_tolerance_search is None:
            q = sequencing_bakery(lambda session: session.query(cls))
            q += lambda q: q.filter(cls.mass.between(bindparam("lower"), bindparam("upper")))
            cls._query_ppm_tolerance_search = q
        return cls._query_ppm_tolerance_search(session).params(
            lower=lower, upper=upper)

    @property
    def neutral_mass(self):
        return self.mass

    def to_sequence_composition(self):
        sc = SequenceComposition.parse(self.composition)
        sc.composition_offset = Composition()
        return sc
