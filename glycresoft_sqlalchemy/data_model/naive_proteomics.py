from sqlalchemy.orm import relationship, backref
from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey

from .data_model import PeptideBase


class NaivePeptide(PeptideBase):
    __tablename__ = "NaivePeptide"

    id = Column(Integer, ForeignKey(PeptideBase.id), primary_key=True)
    protein = relationship("Protein", backref=backref("naive_peptides", lazy='dynamic'))
    __mapper_args__ = {
        'polymorphic_identity': u'NaivePeptide',
    }

    def __repr__(self):
        return "<NaivePeptide {0} {1} {2} {3}...>".format(
            self.id, self.base_peptide_sequence, self.peptide_modifications,
            self.glycosylation_sites)


class TheoreticalGlycopeptideComposition(NaivePeptide):
    __tablename__ = "TheoreticalGlycopeptideComposition"

    id = Column(Integer, ForeignKey(NaivePeptide.id), primary_key=True)
    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    glycopeptide_sequence = Column(Unicode(128), index=True)
    glycan_composition_str = Column(Unicode(128), index=True)

    def __repr__(self):
        return "<TheoreticalGlycopeptideComposition {} {} {}>"

    __mapper_args__ = {
        'polymorphic_identity': u'TheoreticalGlycopeptideComposition',
    }
