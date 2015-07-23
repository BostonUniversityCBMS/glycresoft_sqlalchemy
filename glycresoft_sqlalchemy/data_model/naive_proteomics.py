from sqlalchemy.orm import relationship, backref
from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey, Table

from .data_model import Base, PeptideBase, Glycan


class NaivePeptide(PeptideBase):
    __tablename__ = "NaivePeptide"

    id = Column(Integer, primary_key=True)
    protein = relationship("Protein", backref=backref("naive_peptides", lazy='dynamic'))
    __mapper_args__ = {
        'polymorphic_identity': u'NaivePeptide',
        "concrete": True
    }

    def __repr__(self):
        return "<NaivePeptide {0} {1} {2} {3}...>".format(
            self.id, self.base_peptide_sequence, self.peptide_modifications,
            self.glycosylation_sites)


TheoreticalGlycopeptideCompositionGlycanAssociation = Table(
    "TheoreticalGlycopeptideCompositionGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("TheoreticalGlycopeptideComposition.id")),
    Column("glycan_id", Integer, ForeignKey(Glycan.id)))


class TheoreticalGlycopeptideComposition(PeptideBase):
    __tablename__ = "TheoreticalGlycopeptideComposition"

    id = Column(Integer, primary_key=True)

    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    glycopeptide_sequence = Column(Unicode(128), index=True)
    glycan_composition_str = Column(Unicode(128), index=True)
    glycans = relationship(Glycan, secondary=TheoreticalGlycopeptideCompositionGlycanAssociation)

    protein = relationship("Protein", backref=backref("theoretical_glycopeptide_compositions", lazy='dynamic'))

    def __repr__(self):
        return "<TheoreticalGlycopeptideComposition {} {} {}>".format(
            self.id, self.glycopeptide_sequence, self.calculated_mass)

    __mapper_args__ = {
        'polymorphic_on': "sequence_type",
        'polymorphic_identity': u'TheoreticalGlycopeptideComposition',
        'concrete': True
    }
