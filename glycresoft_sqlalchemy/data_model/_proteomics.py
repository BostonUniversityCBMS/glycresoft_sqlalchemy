# from sqlalchemy.orm import relationship, backref
# from sqlalchemy import Numeric, Unicode, Column, Integer, ForeignKey, Table

# from .data_model import Base, PeptideBase
# from .glycomics import TheoreticalGlycanCombination, TheoreticalGlycanComposition, with_glycan_composition


# class NaivePeptide(PeptideBase, Base):
#     __tablename__ = "NaivePeptide"

#     id = Column(Integer, primary_key=True)
#     protein = relationship("Protein", backref=backref("naive_peptides", lazy='dynamic'))

#     def __repr__(self):
#         return "<NaivePeptide {0} {1} {2} {3}...>".format(
#             self.id, self.most_detailed_sequence, self.peptide_modifications,
#             self.glycosylation_sites)


# TheoreticalGlycopeptideCompositionGlycanAssociation = Table(
#     "TheoreticalGlycopeptideCompositionGlycanAssociation", Base.metadata,
#     Column("peptide_id", Integer, ForeignKey("TheoreticalGlycopeptideComposition.id", ondelete="CASCADE"), index=True),
#     Column("glycan_id", Integer, ForeignKey(TheoreticalGlycanComposition.id, ondelete="CASCADE")))


# @with_glycan_composition("glycan_composition_str")
# class TheoreticalGlycopeptideComposition(PeptideBase, Base):
#     __tablename__ = "TheoreticalGlycopeptideComposition"

#     id = Column(Integer, primary_key=True)

#     glycan_mass = Column(Numeric(12, 6, asdecimal=False))
#     glycopeptide_sequence = Column(Unicode(128), index=True)
#     glycan_composition_str = Column(Unicode(128), index=True)

#     glycan_combination_id = Column(Integer, ForeignKey(TheoreticalGlycanCombination.id), index=True)
#     glycans = relationship(TheoreticalGlycanCombination)

#     protein = relationship("Protein", backref=backref("theoretical_glycopeptide_compositions", lazy='dynamic'))

#     def __repr__(self):
#         return "<{} {} {} {} {}>".format(
#             self.__class__.__name__,
#             self.id, self.glycopeptide_sequence, self.peptide_modifications, self.calculated_mass)

#     __mapper_args__ = {
#         'polymorphic_on': PeptideBase.sequence_type,
#         'polymorphic_identity': u'TheoreticalGlycopeptideComposition',
#     }