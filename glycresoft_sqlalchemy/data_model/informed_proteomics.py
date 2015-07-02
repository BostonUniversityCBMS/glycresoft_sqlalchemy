from sqlalchemy.orm import relationship, backref
from sqlalchemy import PickleType, Numeric, Unicode, Column, Integer, ForeignKey, Table

from .data_model import Base, Glycan, TheoreticalGlycopeptide, PeptideBase


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


class InformedMS1Glycopeptide(InformedPeptide):
    __tablename__ = "InformedMS1Glycopeptide"
    id = Column(Integer, ForeignKey(InformedPeptide.id), primary_key=True)

    base_peptide = relationship(InformedPeptide, remote_side=[id])
    glycan_id = Column(Integer, ForeignKey(Glycan.id), index=True)

    glycopeptide_sequence = Column(Unicode(128), index=True)

    glycan_composition_str = Column(Unicode(128), index=True)

    glycan = relationship(
        Glycan, backref=backref('informed_ms1_glycopeptides', lazy='dynamic'))

    __mapper_args__ = {
        'polymorphic_identity': u'InformedMS1Glycopeptide',
    }
