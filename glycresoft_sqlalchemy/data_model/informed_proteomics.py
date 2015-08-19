from sqlalchemy.orm import relationship, backref
from sqlalchemy import PickleType, Numeric, Unicode, Column, Integer, ForeignKey, Table
from sqlalchemy import event

from .data_model import Base, PeptideBase, Protein
from .naive_proteomics import TheoreticalGlycopeptideComposition
from .glycomics import has_glycan_composition_listener


class InformedPeptide(PeptideBase, Base):
    __tablename__ = "InformedPeptide"
    id = Column(Integer, primary_key=True)
    peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
    peptide_score_type = Column(Unicode(56))
    protein = relationship(Protein, backref=backref('informed_peptides', lazy='dynamic'))
    other = Column(PickleType)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedPeptide',
        "concrete": True
    }


InformedPeptideToTheoreticalGlycopeptide = Table(
    "InformedPeptideToTheoreticalGlycopeptide", Base.metadata,
    Column("informed_peptide", Integer, ForeignKey(InformedPeptide.id, ondelete="CASCADE"), index=True),
    Column("theoretical_glycopeptide", Integer, ForeignKey(TheoreticalGlycopeptideComposition.id, ondelete="CASCADE")))


class InformedTheoreticalGlycopeptideComposition(TheoreticalGlycopeptideComposition):
    __tablename__ = "InformedTheoreticalGlycopeptideComposition"

    id = Column(Integer, ForeignKey(TheoreticalGlycopeptideComposition.id, ondelete="CASCADE"), primary_key=True)

    peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
    other = Column(PickleType)
    base_peptide_id = Column(Integer, ForeignKey(InformedPeptide.id), index=True)
    base_peptide = relationship(InformedPeptide)

    __mapper_args__ = {
        'polymorphic_identity': u'InformedTheoreticalGlycopeptideComposition',
    }

has_glycan_composition_listener(InformedTheoreticalGlycopeptideComposition.glycan_composition_str)
