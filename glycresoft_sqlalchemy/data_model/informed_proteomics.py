from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy import PickleType, Numeric, Unicode, create_engine, Column, Integer, ForeignKey

from .data_model import Base, Experiment, Protein, TheoreticalGlycopeptide, GlycopeptideMatch


class InformedPeptide(Base):
    id = Column(Integer, primary_key=True, autoincrement=True)
    base_peptide_sequence = Column(Unicode(128), index=True)
    modified_peptide_sequence = Column(Unicode(128), index=True)
    glycosylation_sites = Column(PickleType)

