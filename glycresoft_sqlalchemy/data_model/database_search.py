# -*- coding: utf-8 -*-
from sqlalchemy.ext.baked import bakery
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table, bindparam,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.session import object_session

from .base import Hierarchy
from .generic import MutableDict, MutableList
from .data_model import (
    Base, Hypothesis, TheoreticalGlycopeptide, TheoreticalPeptideProductIon,
    PeptideBase, Glycan, Protein, GlycopeptideSequenceMS2Base,
    MS1GlycopeptideHypothesis, MS1GlycanHypothesis,
    MS2GlycopeptideHypothesis, MS2GlycanHypothesis,
    ExactMS1GlycopeptideHypothesis, ExactMS2GlycopeptideHypothesis)

from .proteomics import TheoreticalGlycopeptideComposition
from .glycomics import TheoreticalGlycanComposition, MassShift, with_glycan_composition, TheoreticalGlycanStructure

from .informed_proteomics import InformedTheoreticalGlycopeptideComposition
from .json_type import tryjson, clean_dict
from .observed_ions import SampleRun, Peak, ScanBase, TandemScan, HasPeakChromatogramData, PeakGroupBase
from glycresoft_sqlalchemy.utils.data_migrator import Migrator


hypothesis_sample_match_root = Hierarchy()


@hypothesis_sample_match_root.references(Hypothesis)
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

    samples = relationship("SampleRun", secondary=lambda: HypothesisSampleMatchToSample.__table__)

    def __init__(self, **kwargs):
        kwargs.setdefault('parameters', {})
        super(HypothesisSampleMatch, self).__init__(**kwargs)

    def to_json(self):
        d = {
            "id": self.id,
            "parameters": clean_dict(self.parameters),
            "hypothesis_sample_match_type": self.hypothesis_sample_match_type,
            "name": self.name,
            "sample_run_name": self.sample_run_name,
            "target_hypothesis": tryjson(self.target_hypothesis),
            "decoy_hypothesis": tryjson(self.decoy_hypothesis)
        }
        return d

    def copy_tandem_sample_run(self, sample_run):

        target = object_session(self)
        source = object_session(sample_run)

        sid = sample_run.id
        migrator = Migrator(source, target)
        migrator.copy_model(SampleRun, lambda q: q.filter(SampleRun.id == sid))
        print SampleRun
        migrator.copy_model(TandemScan, lambda q: q.filter(TandemScan.sample_run_id == sid))
        print TandemScan
        migrator.copy_model(Peak, lambda q: q.join(ScanBase).filter(ScanBase.sample_run_id == sid))

        target.commit()
        link = HypothesisSampleMatchToSample(
            hypothesis_sample_match_id=self.id,
            sample_run_id=sample_run.id)
        target.add(link)
        target.commit()

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
        return ("<{self.__class__.__name__} {self.id} {self.target_hypothesis_id}" +
                " {self.sample_run_name} {matches}>").format(
            self=self, matches=self.results().next()[0].__name__)

    __mapper_args__ = {
        "polymorphic_identity": u"HypothesisSampleMatch",
        "polymorphic_on": hypothesis_sample_match_type
    }

    is_exact = False

    @property
    def results_type(self):
        return (GlycopeptideMatch, PeakGroupMatch)

    @property
    def basis_for(self):
        return NotImplemented

    hierarchy_root = hypothesis_sample_match_root

    def monosaccharides_in_results(self):
        return NotImplemented


class HypothesisSampleMatchToSample(Base):
    __tablename__ = "HypothesisSampleMatchToSample"

    id = Column(Integer, primary_key=True)
    hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)


@hypothesis_sample_match_root.references(MS1GlycanHypothesis)
class MS1GlycanHypothesisSampleMatch(HypothesisSampleMatch):
    __tablename__ = "MS1GlycanHypothesisSampleMatch"
    id = Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    def results(self):
        yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition") |
                (PeakGroupMatch.theoretical_match_type == None))

    __mapper_args__ = {
        "polymorphic_identity": u"MS1GlycanHypothesisSampleMatch",
    }

    ms_level = 1

    @property
    def results_type(self):
        return (PeakGroupMatch)

    @property
    def hypothesis_type(self):
        return MS1GlycanHypothesis

    @property
    def basis_for(self):
        return MS2GlycanHypothesis

    def monosaccharides_in_results(self):
        return NotImplemented


@hypothesis_sample_match_root.references(MS2GlycanHypothesis)
class MS2GlycanHypothesisSampleMatch(HypothesisSampleMatch):
    __tablename__ = "MS2GlycanHypothesisSampleMatch"
    id = Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    def results(self):
        if self.peak_group_matches.first() is not None:
            yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycanComposition") |
                (PeakGroupMatch.theoretical_match_type == None))
        yield TheoreticalGlycanStructure, ()

    __mapper_args__ = {
        "polymorphic_identity": u"MS2GlycanHypothesisSampleMatch",
    }

    ms_level = 1

    @property
    def results_type(self):
        return (PeakGroupMatch)

    @property
    def hypothesis_type(self):
        return MS2GlycanHypothesis

    @property
    def basis_for(self):
        return NotImplemented


@hypothesis_sample_match_root.references(MS1GlycopeptideHypothesis)
class MS1GlycopeptideHypothesisSampleMatch(HypothesisSampleMatch):
    __tablename__ = "MS1GlycopeptideHypothesisSampleMatch"

    id = Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        "polymorphic_identity": u"MS1GlycopeptideHypothesisSampleMatch",
    }

    def results(self):
        yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition") |
                (PeakGroupMatch.theoretical_match_type == None))

    is_exact = False
    ms_level = 1

    @property
    def hypothesis_type(self):
        return MS1GlycopeptideHypothesis

    @property
    def results_type(self):
        return (PeakGroupMatch)

    @property
    def basis_for(self):
        return MS2GlycopeptideHypothesis

    def monosaccharides_in_results(self):
        session = object_session(self)
        search_type, res_query = self.results().next()
        extents = search_type.glycan_composition_extents(
            session, lambda q: q.join(search_type).join(
                self.results_type, search_type.id == self.results_type.theoretical_match_id).filter(
                self.results_type.hypothesis_sample_match_id == self.id,
                self.results_type.ms1_score > 0.1))
        return [str(r[0]) for r in extents]


@hypothesis_sample_match_root.references(ExactMS1GlycopeptideHypothesis)
class ExactMS1GlycopeptideHypothesisSampleMatch(MS1GlycopeptideHypothesisSampleMatch):
    __tablename__ = "ExactMS1GlycopeptideHypothesisSampleMatch"

    id = Column(Integer, ForeignKey(MS1GlycopeptideHypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        "polymorphic_identity": u"ExactMS1GlycopeptideHypothesisSampleMatch",
    }

    is_exact = True

    def results(self):
        yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                (PeakGroupMatch.theoretical_match_type == "TheoreticalGlycopeptideComposition") |
                (PeakGroupMatch.theoretical_match_type == None))

    @property
    def results_type(self):
        return (PeakGroupMatch)

    @property
    def hypothesis_type(self):
        return ExactMS1GlycopeptideHypothesis

    @property
    def basis_for(self):
        return ExactMS2GlycopeptideHypothesis


@hypothesis_sample_match_root.references(MS2GlycopeptideHypothesis)
class MS2GlycopeptideHypothesisSampleMatch(HypothesisSampleMatch):
    __tablename__ = "MS2GlycopeptideHypothesisSampleMatch"

    id = Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    __mapper_args__ = {
        "polymorphic_identity": u"MS2GlycopeptideHypothesisSampleMatch",
    }

    is_exact = False
    ms_level = 2

    hypothesis_type = MS2GlycopeptideHypothesis

    def results(self):
        yield GlycopeptideMatch, self.glycopeptide_matches.filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == Hypothesis.id,
                ~Hypothesis.is_decoy)

    @property
    def hypothesis_type(self):
        return MS2GlycopeptideHypothesis

    @property
    def results_type(self):
        return (GlycopeptideMatch)

    def monosaccharides_in_results(self, threshold=0.1):
        session = object_session(self)
        search_type, res_query = self.results().next()
        extents = search_type.glycan_composition_extents(
            session, lambda q: q.join(
                search_type).filter(
                search_type.hypothesis_sample_match_id == self.id,
                search_type.ms2_score >= threshold, search_type.is_not_decoy())).all()
        return [str(r[0]) for r in extents]


@hypothesis_sample_match_root.references(ExactMS2GlycopeptideHypothesis)
class ExactMS2GlycopeptideHypothesisSampleMatch(MS2GlycopeptideHypothesisSampleMatch):
    __tablename__ = "ExactMS2GlycopeptideHypothesisSampleMatch"
    __mapper_args__ = {
        "polymorphic_identity": u"ExactMS2GlycopeptideHypothesisSampleMatch",
    }

    id = Column(Integer, ForeignKey(MS2GlycopeptideHypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    @property
    def hypothesis_type(self):
        return ExactMS2GlycopeptideHypothesis

    is_exact = True
    ms_level = 2


def hypothesis_sample_match_type(theoretical_type):
    if theoretical_type == TheoreticalGlycopeptideComposition or\
            issubclass(TheoreticalGlycopeptideComposition):
        return MS1GlycopeptideHypothesisSampleMatch
    elif theoretical_type == TheoreticalGlycanComposition or\
            issubclass(TheoreticalGlycanComposition):
        return MS1GlycanHypothesisSampleMatch
    elif theoretical_type == TheoreticalGlycopeptide or issubclass(TheoreticalGlycopeptide):
        return MS2GlycopeptideHypothesisSampleMatch
    else:
        raise KeyError(theoretical_type)


GlycopeptideMatchGlycanAssociation = Table(
    "GlycopeptideMatchGlycanAssociation", Base.metadata,
    Column("peptide_id", Integer, ForeignKey("GlycopeptideMatch.id", ondelete="CASCADE")),
    Column("glycan_id", Integer, ForeignKey(Glycan.id, ondelete="CASCADE")))


@with_glycan_composition("glycan_composition_str")
class GlycopeptideMatch(PeptideBase, Base, GlycopeptideSequenceMS2Base):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, primary_key=True)
    theoretical_glycopeptide_id = Column(Integer, ForeignKey(TheoreticalGlycopeptide.id), index=True)
    theoretical_reference = relationship(TheoreticalGlycopeptide)
    hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    hypothesis_sample_match = relationship(
        HypothesisSampleMatch, backref=backref("glycopeptide_matches", lazy='dynamic'))

    ms2_score = Column(Numeric(10, 6, asdecimal=False), index=True)

    glycans = relationship(Glycan, secondary=GlycopeptideMatchGlycanAssociation, lazy='dynamic')

    # As defined in [1]
    p_value = Column(Numeric(10, 6, asdecimal=False))
    q_value = Column(Numeric(10, 6, asdecimal=False))

    scan_id_range = Column(MutableList.as_mutable(PickleType))
    first_scan = Column(Integer)
    last_scan = Column(Integer)

    mean_coverage = Column(Numeric(10, 6, asdecimal=False))
    mean_hexnac_coverage = Column(Numeric(10, 6, asdecimal=False))

    spectrum_matches = relationship("GlycopeptideSpectrumMatch", backref='glycopeptide_match', lazy='dynamic')

    @classmethod
    def is_not_decoy(cls):
        return ((cls.protein_id == Protein.id) & (Protein.hypothesis_id == Hypothesis.id) & (~Hypothesis.is_decoy))

    def __repr__(self):
        rep = "<GlycopeptideMatch {} {} {}>".format(self.glycopeptide_sequence, self.ms2_score, self.observed_mass)
        return rep


class SpectrumMatchBase(object):
    scan_time = Column(Integer, index=True)
    peak_match_map = Column(MutableDict.as_mutable(PickleType))
    peaks_explained = Column(Integer, index=True)
    peaks_unexplained = Column(Integer)
    best_match = Column(Boolean, index=True)

    @declared_attr
    def hypothesis_sample_match_id(cls):
        return Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), index=True)

    @declared_attr
    def hypothesis_id(cls):
        return Column(Integer, ForeignKey(Hypothesis.id, ondelete="CASCADE"), index=True)

    def __iter__(self):
        return iter(self.spectrum)

    @property
    def spectrum(self):
        session = object_session(self)
        s = session.query(TandemScan).join(
            HypothesisSampleMatchToSample,
            (TandemScan.sample_run_id == HypothesisSampleMatchToSample.sample_run_id)).filter(
            HypothesisSampleMatch.id == self.hypothesis_sample_match_id).filter(
            TandemScan.sample_run_id == SampleRun.id,
            SampleRun.name == HypothesisSampleMatch.sample_run_name,
            TandemScan.time == self.scan_time)
        return s.first()

    def peak_explained_by(self, peak_id):
        explained = set()
        try:
            matches = self.peak_match_map[peak_id]
            for match in matches:
                explained.add(match["key"])
        except KeyError:
            pass
        return explained

    def add_peak_explaination(self, peak, match_record):
        recs = self.peak_match_map.get(peak.id, [])
        recs.append(match_record)
        self.peak_match_map[peak.id] = recs

    def total_signal_explained(self):
        total = 0.
        for k, v in self.peak_match_map.items():
            ion = v[0]
            total += ion['intensity']
        return total

    def matches_for(self, key):
        for peak_index, matches in self.peak_match_map.items():
            for match in matches:
                if match['key'] == key:
                    yield match


class GlycopeptideSpectrumMatch(Base, SpectrumMatchBase):
    __tablename__ = "GlycopeptideSpectrumMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    glycopeptide_match_id = Column(Integer, ForeignKey(GlycopeptideMatch.id), index=True)

    def __repr__(self):
        return "<GlycopeptideSpectrumMatch {} -> Spectrum {} | {} Peaks Matched>".format(
            self.glycopeptide_match.glycopeptide_sequence,
            self.scan_time, len(self.peak_match_map))

    @property
    def glycopeptide_sequence(self):
        return self.glycopeptide_match.glycopeptide_sequence


class GlycanStructureMatch(Base):
    __tablename__ = "GlycanStructureMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)

    theoretical_reference_id = Column(
        Integer, ForeignKey(TheoreticalGlycanStructure.id, ondelete="CASCADE"), index=True)
    theoretical_reference = relationship(TheoreticalGlycanStructure)

    hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    hypothesis_sample_match = relationship(
        HypothesisSampleMatch, backref=backref("glycan_structure_matches", lazy='dynamic'))

    ms1_score = Column(Numeric(10, 6, asdecimal=False), index=True)
    ms2_score = Column(Numeric(10, 6, asdecimal=False), index=True)

    observed_mass = Column(Numeric(12, 6, asdecimal=False))
    ppm_error = Column(Numeric(12, 6, asdecimal=False))
    volume = Column(Numeric(12, 4, asdecimal=False))

    fragment_matches = Column(MutableList.as_mutable(PickleType))

    scan_id_range = Column(MutableList.as_mutable(PickleType))
    first_scan = Column(Integer)
    last_scan = Column(Integer)

    spectrum_matches = relationship("GlycanSpectrumMatch", backref='glycan_structure_match', lazy='dynamic')

    def structure(self):
        return self.theoretical_reference.structure()

    @property
    def name(self):
        return self.theoretical_reference.name

    @property
    def glycoct(self):
        return self.theoretical_reference.glycoct

    @property
    def calculated_mass(self):
        return self.theoretical_reference.calculated_mass

    @property
    def composition(self):
        return self.theoretical_reference.composition

    @property
    def glycan_composition(self):
        return self.theoretical_reference.glycan_composition

    @property
    def reduction(self):
        return self.theoretical_reference.reduction

    @property
    def derivatization(self):
        return self.theoretical_reference.derivatization

    def __repr__(self):
        rep = "<{self.__class__.__name__} {self.ms2_score}\n{self.theoretical_reference.glycoct}>".format(self=self)
        return rep


class GlycanSpectrumMatch(Base, SpectrumMatchBase):
    __tablename__ = "GlycanSpectrumMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    glycan_structure_match_id = Column(Integer, ForeignKey(GlycanStructureMatch.id), index=True)

    def __repr__(self):
        return "<GlycanSpectrumMatch {} -> Spectrum {} | {} Peaks Matched>".format(
            self.glycan_structure_match.name or self.glycan_structure_match.id,
            self.scan_time, len(self.peak_match_map))


TheoreticalCompositionMap = {
    "TheoreticalGlycopeptideComposition": TheoreticalGlycopeptideComposition,
    "TheoreticalGlycanComposition": TheoreticalGlycanComposition,
    "InformedTheoreticalGlycopeptideComposition": InformedTheoreticalGlycopeptideComposition
}


class PeakGroupMatchBase(object):
    @declared_attr
    def hypothesis_sample_match_id(cls):
        return Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)

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


class PeakGroupMatch(Base, HasPeakChromatogramData, PeakGroupBase, PeakGroupMatchBase):
    __tablename__ = "PeakGroupMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)

    hypothesis_sample_match = relationship(HypothesisSampleMatch, backref=backref("peak_group_matches", lazy='dynamic'))
    # ForeignKey across database boundaries
    peak_group_id = Column(Integer, index=True)

    mass_shift_type = Column(Integer, ForeignKey(MassShift.id))
    mass_shift = relationship(MassShift)
    mass_shift_count = Column(Integer)

    def __repr__(self):
        rep = "<PeakGroupMatch {id} {ms1_score:0.4f} {weighted_monoisotopic_mass:0.4f}" +\
              " {mass_shift_name} {theoretical_match_type}>"
        return rep.format(mass_shift_name=self.mass_shift.name if self.mass_shift is not None else "", **self.__dict__)


def protein_peak_group_matches(protein):
    return object_session(protein).query(PeakGroupMatch).join(
        TheoreticalGlycopeptideComposition,
        PeakGroupMatch.theoretical_match_id == TheoreticalGlycopeptideComposition.id).filter(
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

    subgroups = relationship(PeakGroupMatch, secondary=lambda: PeakGroupMatchToJointPeakGroupMatch, lazy='dynamic')

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


'''
Citations
---------

[1] L. Käll, J. D. Storey, M. J. MacCoss, and W. S. Noble,
“Assigning significance to peptides identified by tandem mass spectrometry using decoy databases,”
J. Proteome Res., vol. 7, pp. 29–34, 2008.

'''
