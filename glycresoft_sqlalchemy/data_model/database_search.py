# -*- coding: utf-8 -*-

from sqlalchemy.ext.baked import bakery
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table, bindparam,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.session import object_session

from .generic import MutableDict, MutableList
from .data_model import Base, Hypothesis, TheoreticalGlycopeptide, PeptideBase, Glycan, Protein
from .naive_proteomics import TheoreticalGlycopeptideComposition
from .glycomics import TheoreticalGlycanComposition, MassShift
from .informed_proteomics import InformedTheoreticalGlycopeptideComposition
from .json_type import tryjson


class HypothesisSampleMatch(Base):
    __tablename__ = "HypothesisSampleMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    sample_run_name = Column(Unicode(128))
    target_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))
    decoy_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))
    target_hypothesis = relationship(Hypothesis, foreign_keys=[target_hypothesis_id])
    decoy_hypothesis = relationship(Hypothesis, foreign_keys=[decoy_hypothesis_id])

    hypothesis_sample_match_type = Column(Unicode(56))

    def __init__(self, **kwargs):
        kwargs.setdefault('parameters', {})
        super(HypothesisSampleMatch, self).__init__(**kwargs)

    def to_json(self):
        d = {
            "id": self.id,
            "parameters": self.parameters,
            "hypothesis_sample_match_type": self.hypothesis_sample_match_type,
            "name": self.name,
            "sample_run_name": self.sample_run_name,
            "target_hypothesis": tryjson(self.target_hypothesis),
            "decoy_hypothesis": tryjson(self.decoy_hypothesis)
        }
        return d

    def results(self):
        if self.glycopeptide_matches.first() is not None:
            yield GlycopeptideMatch, self.glycopeptide_matches.filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == Hypothesis.id,
                ~Hypothesis.is_decoy)
        if self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition").first() is not None:
            yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition") |
                (PeakGroupMatch.theoretical_match_type == None))
        if self.peak_group_matches.filter(
                PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition").first() is not None:
            yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition") |
                (PeakGroupMatch.theoretical_match_type == None))

    def __repr__(self):
        return "<{self.__class__.__name__} {self.id} {self.target_hypothesis_id} {self.sample_run_name} {matches}>".format(
            self=self, matches=','.join(res[0].__name__ for res in self.results()))

    __mapper_args__ = {
        "polymorphic_identity": u"HypothesisSampleMatch",
        "polymorphic_on": hypothesis_sample_match_type
    }


GlycopeptideMatchGlycanAssociation = Table(
    "GlycopeptideMatchGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("GlycopeptideMatch.id", ondelete="CASCADE")),
    Column("glycan_id", Integer, ForeignKey(Glycan.id, ondelete="CASCADE")))


class GlycopeptideMatch(PeptideBase, Base):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, primary_key=True)
    theoretical_glycopeptide_id = Column(Integer, ForeignKey(TheoreticalGlycopeptide.id))
    hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
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
    "TheoreticalGlycanComposition": TheoreticalGlycanComposition,
    "InformedTheoreticalGlycopeptideComposition": InformedTheoreticalGlycopeptideComposition
}


class PeakGroupMatch(Base):
    __tablename__ = "PeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    hypothesis_sample_match = relationship(HypothesisSampleMatch, backref=backref("peak_group_matches", lazy='dynamic'))

    theoretical_match_type = Column(Unicode(128), index=True)
    theoretical_match_id = Column(Integer)

    # ForeignKey across database boundaries
    peak_group_id = Column(Integer, index=True)

    @property
    def theoretical_match(self):
        try:
            matched_type = TheoreticalCompositionMap[self.theoretical_match_type]
            return object_session(self).query(matched_type).get(self.theoretical_match_id)
        except KeyError:
            return None

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

    def __repr__(self):
        rep = "<PeakGroupMatch {id} {ms1_score} {weighted_monoisotopic_mass} {theoretical_match_type}>"
        return rep.format(**self.__dict__)


class TempPeakGroupMatch(Base):
    __tablename__ = "TempPeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_run_id = Column(Integer, index=True)

    charge_state_count = Column(Integer)
    scan_count = Column(Integer)
    first_scan_id = Column(Integer)
    last_scan_id = Column(Integer)

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

'''
Citations
---------

[1] L. Käll, J. D. Storey, M. J. MacCoss, and W. S. Noble,
“Assigning significance to peptides identified by tandem mass spectrometry using decoy databases,”
J. Proteome Res., vol. 7, pp. 29–34, 2008.

'''
