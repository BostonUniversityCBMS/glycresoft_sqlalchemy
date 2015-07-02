# -*- coding: utf-8 -*-

from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .generic import MutableDict, MutableList
from .data_model import Base, Hypothesis, TheoreticalGlycopeptide, PeptideBase, Glycan


class HypothesisSampleMatch(Base):
    __tablename__ = "HypothesisSampleMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    sample_run_name = Column(Unicode(128))
    target_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))
    decoy_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))



GlycopeptideMatchGlycanAssociation = Table(
    "GlycopeptideMatchGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("GlycopeptideMatch.id")),
    Column("glycan_id", Integer, ForeignKey("Glycan.id")))


class GlycopeptideMatch(PeptideBase):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, primary_key=True)
    theoretical_glycopeptide_id = (Integer, ForeignKey(TheoreticalGlycopeptide.id))
    hypothesis_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id))
    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    ms2_score = Column(Numeric(10, 6, asdecimal=False))

    glycans = relationship(Glycan, secondary=GlycopeptideMatchGlycanAssociation, lazy='dynamic')

    observed_mass = Column(Numeric(10, 6, asdecimal=False))
    glycan_mass = Column(Numeric(10, 6, asdecimal=False))
    ppm_error = Column(Numeric(10, 6, asdecimal=False))
    volume = Column(Numeric(10, 6, asdecimal=False))

    glycopeptide_sequence = Column(Unicode(128), index=True)
    glycan_composition_str = Column(Unicode(128), index=True)

    oxonium_ions = Column(MutableList.as_mutable(PickleType))
    stub_ions = Column(MutableList.as_mutable(PickleType))

    bare_b_ions = Column(MutableList.as_mutable(PickleType))
    glycosylated_b_ions = Column(MutableList.as_mutable(PickleType))

    bare_y_ions = Column(MutableList.as_mutable(PickleType))
    glycosylated_y_ions = Column(MutableList.as_mutable(PickleType))

    # As defined in [1]
    p_value = Column(Numeric(10, 6, asdecimal=False))
    q_value = Column(Numeric(10, 6, asdecimal=False))

    scan_id_range = Column(MutableList.as_mutable(PickleType))
    first_scan = Column(Integer)
    last_scan = Column(Integer)

    mean_coverage = Column(Numeric(10, 6, asdecimal=False))
    mean_hexnac_coverage = Column(Numeric(10, 6, asdecimal=False))

    def __repr__(self):
        rep = "<GlycopeptideMatch {} {} {}>".format(self.glycopeptide_sequence, self.ms2_score, self.observed_mass)
        return rep

    __mapper_args__ = {
        'polymorphic_identity': u'GlycopeptideMatch',
        'concrete': True
    }


class SpectrumMatch(Base):
    __tablename__ = "SpectrumMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    spectrum_id = Column(Integer)
    glycopeptide_match_id = Column(Integer, ForeignKey(GlycopeptideMatch.id), index=True)
    peak_match_map = Column(PickleType)


'''
Citations
---------

[1] L. Käll, J. D. Storey, M. J. MacCoss, and W. S. Noble,
“Assigning significance to peptides identified by tandem mass spectrometry using decoy databases,”
J. Proteome Res., vol. 7, pp. 29–34, 2008.

'''
