from sqlalchemy.orm import relationship, backref
from sqlalchemy import PickleType, Numeric, Unicode, Column, Integer, ForeignKey, Table

from .data_model import Base, Glycan, TheoreticalGlycopeptide, PeptideBase
from .naive_proteomics import TheoreticalGlycopeptideComposition


class InformedPeptide(PeptideBase):
    __tablename__ = "InformedPeptide"
    id = Column(Integer, primary_key=True)
    peptide_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    other = Column(PickleType)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedPeptide',
        "concrete": True
    }


InformedPeptideToTheoreticalGlycopeptide = Table(
    "InformedPeptideToTheoreticalGlycopeptide", Base.metadata,
    Column("informed_peptide", Integer, ForeignKey(InformedPeptide.id)),
    Column("theoretical_glycopeptide", Integer, ForeignKey(TheoreticalGlycopeptide.id)))


class InformedTheoreticalGlycopeptideComposition(TheoreticalGlycopeptideComposition):
    __tablename__ = "InformedTheoreticalGlycopeptideComposition"
    id = Column(Integer, ForeignKey(TheoreticalGlycopeptideComposition.id), primary_key=True)

    peptide_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    other = Column(PickleType)
    base_peptide_id = Column(Integer, ForeignKey(InformedPeptide.id), index=True)
    base_peptide = relationship(InformedPeptide)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedTheoreticalGlycopeptideComposition',
    }
