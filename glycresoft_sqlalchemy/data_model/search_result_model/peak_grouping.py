from collections import defaultdict

from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.session import object_session

import numpy as np


from ..base import Hierarchy, Namespace, Base2 as Base

from ..generic import MutableDict, MutableList
from ..sequence_model.peptide import TheoreticalGlycopeptideComposition, Protein
from ..glycomics import TheoreticalGlycanComposition, FrozenGlycanComposition, MassShift
from ..observed_ions import HasPeakChromatogramData, PeakGroupBase


TheoreticalCompositionMap = {
    "TheoreticalGlycopeptideComposition": TheoreticalGlycopeptideComposition,
    "TheoreticalGlycanComposition": TheoreticalGlycanComposition,
}


class PeakGroupMatchBase(object):
    @declared_attr
    def hypothesis_sample_match_id(cls):
        return Column(Integer, ForeignKey("HypothesisSampleMatch.id"), index=True)

    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    matched = Column(Boolean, index=True)

    ppm_error = Column(Numeric(10, 8, asdecimal=False))

    theoretical_match_type = Column(Unicode(128), index=True)
    theoretical_match_id = Column(Integer, index=True)

    @property
    def theoretical_match(self):
        try:
            matched_type = TheoreticalCompositionMap[self.theoretical_match_type]
            return object_session(self).query(matched_type).get(self.theoretical_match_id)
        except KeyError:
            return None

    _glycan_composition = None

    @property
    def glycan_composition(self):
        if self._glycan_composition is None:
            tm = self.theoretical_match
            self._glycan_composition = FrozenGlycanComposition(dict(tm.glycan_composition))
        return self._glycan_composition


class PeakGroupMatch(Base, HasPeakChromatogramData, PeakGroupBase, PeakGroupMatchBase):
    __tablename__ = "PeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)

    hypothesis_sample_match = relationship(
        "HypothesisSampleMatch", backref=backref("individual_peak_group_matches", lazy='dynamic'))
    # ForeignKey across database boundaries
    peak_group_id = Column(Integer, index=True)

    mass_shift_type = Column(Integer, ForeignKey(MassShift.id))
    mass_shift = relationship(MassShift)
    mass_shift_count = Column(Integer)

    def __repr__(self):
        mass_shift_part = self.mass_shift.name if self.mass_shift is not None else ""
        if mass_shift_part != "":
            if self.mass_shift_count > 1:
                mass_shift_part += " (%d)" % self.mass_shift_count

        rep = "<PeakGroupMatch {id} {total_volume} {weighted_monoisotopic_mass:0.4f}" +\
              " {mass_shift_part} {theoretical_match_type}>"
        return rep.format(
            id=self.id, weighted_monoisotopic_mass=self.weighted_monoisotopic_mass,
            theoretical_match_type=self.theoretical_match_type,
            mass_shift_part=mass_shift_part,
            total_volume=self.total_volume if self.total_volume else 0.0)


def protein_peak_group_matches(protein):
    return object_session(protein).query(PeakGroupMatchType).join(
        TheoreticalGlycopeptideComposition,
        PeakGroupMatchType.theoretical_match_id == TheoreticalGlycopeptideComposition.id).filter(
        TheoreticalGlycopeptideComposition.protein_id == protein.id)

Protein.peak_group_matches = property(fget=protein_peak_group_matches)


class TempPeakGroupMatch(Base, HasPeakChromatogramData, PeakGroupBase):
    __tablename__ = "TempPeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_run_id = Column(Integer, index=True)

    peak_ids = Column(MutableList.as_mutable(PickleType))

    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    matched = Column(Boolean, index=True)


class JointPeakGroupMatch(Base, HasPeakChromatogramData, PeakGroupBase, PeakGroupMatchBase):
    __tablename__ = "JointPeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    modification_state_count = Column(Integer)

    fingerprint = Column(UnicodeText, index=True)

    hypothesis_sample_match = relationship(
        "HypothesisSampleMatch", backref=backref("peak_group_matches", lazy='dynamic'))

    subgroups = relationship(PeakGroupMatch, secondary=lambda: PeakGroupMatchToJointPeakGroupMatch, lazy='dynamic')

    def modification_states(self):
        states = []
        for group in self.subgroups:
            states.append((group.mass_shift, group.mass_shift_count))
        return states

    def as_feature_vector(self):
        names = [
            "charge_state_count",
            "scan_density",
            "modification_state_count",
            "total_volume",
            "a_peak_intensity_error",
            "centroid_scan_error",
            "scan_count",
            "average_signal_to_noise"
        ]
        return [getattr(self, n) for n in names]

    def __repr__(self):
        if self.ms1_score is None:
            ms1_score = ""
        else:
            ms1_score = "%0.4f" % self.ms1_score
        weighted_monoisotopic_mass = self.weighted_monoisotopic_mass
        id = self.id
        modification_state_count = self.modification_state_count

        theoretical_match_type = self.theoretical_match_type

        rep = "<JointPeakGroupMatch {id} {ms1_score} {weighted_monoisotopic_mass:0.4f}" +\
            " {modification_state_count} {theoretical_match_type}>"
        return rep.format(
            id=id, ms1_score=ms1_score, weighted_monoisotopic_mass=weighted_monoisotopic_mass,
            modification_state_count=modification_state_count, theoretical_match_type=theoretical_match_type)


PeakGroupMatchToJointPeakGroupMatch = Table(
    "PeakGroupMatchToJointPeakGroupMatch", Base.metadata,
    Column("peak_group_id", Integer, ForeignKey("PeakGroupMatch.id"), index=True),
    Column("joint_group_id", Integer, ForeignKey("JointPeakGroupMatch.id"), index=True)
)


PeakGroupMatchType = JointPeakGroupMatch
