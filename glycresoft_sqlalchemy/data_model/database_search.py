# -*- coding: utf-8 -*-

from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.orm import relationship, backref
from sqlalchemy import (PickleType, Numeric, Unicode, create_engine, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .generic import MutableDict, MutableList
from .data_model import Base, Hypothesis, TheoreticalGlycopeptide
from .observed_ions import SampleRun


class HypothesisSampleMatch(Base):
    __tablename__ = "HypothesisSampleMatch"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id))
    target_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))
    decoy_hypothesis_id = Column(Integer, ForeignKey(Hypothesis.id))


class GlycopeptideMatch(TheoreticalGlycopeptide):
    __tablename__ = "GlycopeptideMatch"

    id = Column(Integer, ForeignKey(TheoreticalGlycopeptide.id), primary_key=True)
    theoretical_glycopeptide = relationship(TheoreticalGlycopeptide, remote_side=[id])
    hypothesis_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id))
    ms2_score = Column(Numeric(10, 6, asdecimal=False))

    # As defined in [1]
    p_value = Column(Numeric(10, 6, asdecimal=False))
    q_value = Column(Numeric(10, 6, asdecimal=False))

    total_bare_b_ions_possible = Column(Integer)
    total_glycosylated_b_ions_possible = Column(Integer)
    total_bare_y_ions_possible = Column(Integer)
    total_glycosylated_y_ions_possible = Column(Integer)

    percent_bare_b_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_bare_y_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_glycosylated_b_ion_coverage = Column(Numeric(10, 6, asdecimal=False))
    percent_glycosylated_y_ion_coverage = Column(Numeric(10, 6, asdecimal=False))

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
