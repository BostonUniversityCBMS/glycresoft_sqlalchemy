from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy import PickleType, Numeric, Unicode, create_engine, Column, Integer, ForeignKey

from .data_model import Base, Experiment, Protein, TheoreticalGlycopeptide, GlycopeptideMatch


class InformedPeptide(Base):
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey(Protein.id), index=True)
    base_peptide_sequence = Column(Unicode(128), index=True)
    modified_peptide_sequence = Column(Unicode(128), index=True)
    glycosylation_sites = Column(PickleType)
    peptide_score = Column(Numeric(10, 6), index=True)
    start_position = Column(Integer, index=True)
    end_position = Column(Integer, index=True)


class InformedPeptideToTheoreticalGlycopeptide(Base):
    id = Column(Integer, primary_key=True, autoincrement=True)
    informed_peptide = Column(Integer, ForeignKey(InformedPeptide.id), index=True)
    theoretical_glycopeptide = Column(Integer, ForeignKey(TheoreticalGlycopeptide.id), index=True)
