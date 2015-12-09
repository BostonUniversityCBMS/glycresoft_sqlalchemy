# from sqlalchemy.orm import relationship, backref
# from sqlalchemy import PickleType, Numeric, Unicode, Column, Integer, ForeignKey, Table
# from sqlalchemy import event

# from .data_model import Base, PeptideBase, Protein
# from .proteomics import TheoreticalGlycopeptideComposition
# from .glycomics import has_glycan_composition_listener


# class InformedPeptide(PeptideBase, Base):
#     __tablename__ = "InformedPeptide"
#     id = Column(Integer, primary_key=True)
#     peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
#     peptide_score_type = Column(Unicode(56))
#     protein = relationship(Protein, backref=backref('informed_peptides', lazy='dynamic'))
#     theoretical_glycopeptides = relationship(
#         "InformedTheoreticalGlycopeptideComposition",# secondary=lambda: InformedPeptideToTheoreticalGlycopeptide,
#         lazy="dynamic", backref="base_peptide")
#     other = Column(PickleType)

#     __mapper_args__ = {
#         'polymorphic_identity': u'InformedPeptide',
#         "concrete": True
#     }


# class InformedTheoreticalGlycopeptideComposition(TheoreticalGlycopeptideComposition):
#     __tablename__ = "InformedTheoreticalGlycopeptideComposition"

#     id = Column(Integer, ForeignKey(TheoreticalGlycopeptideComposition.id, ondelete="CASCADE"), primary_key=True)

#     peptide_score = Column(Numeric(12, 6, asdecimal=False), index=True)
#     other = Column(PickleType)
#     base_peptide_id = Column(Integer, ForeignKey(InformedPeptide.id), index=True)

#     __mapper_args__ = {
#         'polymorphic_identity': u'InformedTheoreticalGlycopeptideComposition',
#     }

# has_glycan_composition_listener(InformedTheoreticalGlycopeptideComposition.glycan_composition_str)
