# -*- coding: utf-8 -*-

from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm.session import object_session

from .hypothesis import (
    Hypothesis, MS1GlycanHypothesis, MS1GlycopeptideHypothesis,
    ExactMS1GlycopeptideHypothesis, MS2GlycanHypothesis, MS2GlycopeptideHypothesis,
    ExactMS2GlycopeptideHypothesis)

from .base import Hierarchy, Base2 as Base
from .generic import MutableDict
from .json_type import tryjson, clean_dict
from .observed_ions import SampleRun, TandemScan, ScanBase, Peak

from .sequence_model.peptide import (
    Protein, TheoreticalGlycopeptideComposition,
    TheoreticalGlycopeptide
    )

from .glycomics import TheoreticalGlycanComposition

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

    @declared_attr
    def samples(self):
        return relationship(
            SampleRun, secondary=lambda: HypothesisSampleMatchToSample.__table__)

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
        from .search_result_model.peak_grouping import PeakGroupMatchType
        from .search_result_model.sequence_identification import GlycopeptideMatch
        if self.glycopeptide_matches.first() is not None:
            yield GlycopeptideMatch, self.glycopeptide_matches.filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == Hypothesis.id,
                ~Hypothesis.is_decoy)
        if self.peak_group_matches.filter(
                PeakGroupMatchType.theoretical_match_type == "TheoreticalGlycanComposition").first() is not None:
            yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (PeakGroupMatchType.theoretical_match_type == "TheoreticalGlycanComposition") |
                (PeakGroupMatchType.theoretical_match_type == None))
        if self.peak_group_matches.filter(
                PeakGroupMatchType.theoretical_match_type == "TheoreticalGlycopeptideComposition").first() is not None:
            yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                (PeakGroupMatchType.theoretical_match_type == "TheoreticalGlycopeptideComposition") |
                (PeakGroupMatchType.theoretical_match_type == None))

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
        from .search_result_model import GlycopeptideMatch, PeakGroupMatchType
        return (GlycopeptideMatch, PeakGroupMatchType)

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
        results_type = self.results_type
        yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (results_type.theoretical_match_type == "TheoreticalGlycanComposition") |
                (results_type.theoretical_match_type == None))

    __mapper_args__ = {
        "polymorphic_identity": u"MS1GlycanHypothesisSampleMatch",
    }

    ms_level = 1

    @property
    def results_type(self):
        from .search_result_model import GlycopeptideMatch, PeakGroupMatchType
        return (PeakGroupMatchType)

    @property
    def hypothesis_type(self):
        return MS1GlycanHypothesis

    @property
    def basis_for(self):
        return MS2GlycanHypothesis

    def monosaccharides_in_results(self):
        session = object_session(self)
        search_type, res_query = self.results().next()
        extents = search_type.glycan_composition_extents(
            session, lambda q: q.join(search_type).join(
                self.results_type, search_type.id == self.results_type.theoretical_match_id).filter(
                self.results_type.hypothesis_sample_match_id == self.id,
                self.results_type.ms1_score > 0.1))
        return [str(r[0]) for r in extents]


@hypothesis_sample_match_root.references(MS2GlycanHypothesis)
class MS2GlycanHypothesisSampleMatch(HypothesisSampleMatch):
    __tablename__ = "MS2GlycanHypothesisSampleMatch"
    id = Column(Integer, ForeignKey(HypothesisSampleMatch.id, ondelete="CASCADE"), primary_key=True)

    def results(self):
        from .search_result_model import TheoreticalGlycanStructure, PeakGroupMatchType
        if self.peak_group_matches.first() is not None:
            yield TheoreticalGlycanComposition, self.peak_group_matches.filter(
                (PeakGroupMatchType.theoretical_match_type == "TheoreticalGlycanComposition") |
                (PeakGroupMatchType.theoretical_match_type == None))
        yield TheoreticalGlycanStructure, ()

    __mapper_args__ = {
        "polymorphic_identity": u"MS2GlycanHypothesisSampleMatch",
    }

    ms_level = 1

    @property
    def results_type(self):
        from .search_result_model import PeakGroupMatchType
        return (PeakGroupMatchType)

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
        results_type = self.results_type
        yield TheoreticalGlycopeptideComposition, self.peak_group_matches.filter(
                (results_type.theoretical_match_type == "TheoreticalGlycopeptideComposition") |
                (results_type.theoretical_match_type == None))

    is_exact = False
    ms_level = 1

    @property
    def hypothesis_type(self):
        return MS1GlycopeptideHypothesis

    @property
    def results_type(self):
        from .search_result_model import PeakGroupMatchType
        return (PeakGroupMatchType)

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
        from .search_result_model import GlycopeptideMatch

        yield GlycopeptideMatch, self.glycopeptide_matches.filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == Hypothesis.id,
                ~Hypothesis.is_decoy)

    @property
    def hypothesis_type(self):
        return MS2GlycopeptideHypothesis

    @property
    def results_type(self):
        from .search_result_model import GlycopeptideMatch
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
