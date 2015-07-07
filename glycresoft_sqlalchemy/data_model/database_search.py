# -*- coding: utf-8 -*-

from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.session import object_session

from .generic import MutableDict, MutableList
from .data_model import Base, Hypothesis, TheoreticalGlycopeptide, PeptideBase, Glycan
from .naive_proteomics import TheoreticalGlycopeptideComposition
from .glycomics import TheoreticalGlycanComposition, MassShift


class HypothesisSampleMatch(Base):
    __tablename__ = "HypothesisSampleMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    sample_run_name = Column(Unicode(128))
    target_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))
    decoy_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))

    def results(self):
        if self.glycopeptide_matches.count() > 0:
            yield GlycopeptideMatch, self.glycopeptide_matches
        if self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition").count() > 0:
            yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition")
        if self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition").count() > 0:
            yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition")

GlycopeptideMatchGlycanAssociation = Table(
    "GlycopeptideMatchGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("GlycopeptideMatch.id")),
    Column("glycan_id", Integer, ForeignKey(Glycan.id)))


class GlycopeptideMatch(PeptideBase):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, primary_key=True)
    theoretical_glycopeptide_id = (Integer, ForeignKey(TheoreticalGlycopeptide.id))
    hypothesis_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    hypothesis_sample_match = relationship(
        HypothesisSampleMatch, backref=backref("glycopeptide_matches", lazy='dynamic'))

    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    ms2_score = Column(Numeric(10, 6, asdecimal=False), index=True)

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

    spectrum_matches = relationship("SpectrumMatch", backref='glycopeptide_match', lazy='dynamic')

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

    def __repr__(self):
        return "<SpectrumMatch {} -> Spectrum {} | {} Peaks Matched>".format(
            self.glycopeptide_match.glycopeptide_sequence,
            self.spectrum_id, len(self.peak_match_map))


TheoreticalCompositionMap = {
    "TheoreticalGlycopeptideComposition": TheoreticalGlycopeptideComposition,
    "TheoreticalGlycanComposition": TheoreticalGlycanComposition
}


class PeakGroupMatch(Base):
    __tablename__ = "PeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    hypothesis_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    hypothesis_sample_match = relationship(HypothesisSampleMatch, backref=backref("peak_group_matches", lazy='dynamic'))

    theoretical_match_type = Column(Unicode(128), index=True)
    theoretical_match_id = Column(Integer)

    # ForeignKey across database boundaries
    peak_group_id = Column(Integer, index=True)

    @property
    def theoretical_match(self):
        matched_type = TheoreticalCompositionMap[self.theoretical_match_type]
        return object_session(self).query(matched_type).get(self.theoretical_match_id)

    charge_state_count = Column(Integer)
    scan_count = Column(Integer)
    first_scan_id = Column(Integer)
    last_scan_id = Column(Integer)

    ppm_error = Column(Numeric(10, 6, asdecimal=False))

    scan_density = Column(Numeric(10, 6, asdecimal=False))
    weighted_monoisotopic_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

    total_volume = Column(Numeric(12, 6, asdecimal=False))
    average_a_to_a_plus_2_ratio = Column(Numeric(12, 6, asdecimal=False))
    a_peak_intensity_error = Column(Numeric(10, 6, asdecimal=False))
    centroid_scan_estimate = Column(Numeric(12, 6, asdecimal=False))
    centroid_scan_error = Column(Numeric(10, 6, asdecimal=False))
    average_signal_to_noise = Column(Numeric(10, 6, asdecimal=False))

    peak_ids = Column(MutableList.as_mutable(PickleType))
    peak_data = Column(MutableDict.as_mutable(PickleType))

    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    matched = Column(Boolean, index=True)

    mass_shift_type = Column(Integer, ForeignKey(MassShift.id))

    mass_shift_count = Column(Integer)


'''
Citations
---------

[1] L. Käll, J. D. Storey, M. J. MacCoss, and W. S. Noble,
“Assigning significance to peptides identified by tandem mass spectrometry using decoy databases,”
J. Proteome Res., vol. 7, pp. 29–34, 2008.

'''
