from sqlalchemy.orm import relationship, backref
from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey, Table, Text

from .base import Base


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
