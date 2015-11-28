import operator

from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.baked import bakery
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy import (PickleType, Numeric, Unicode, Table, bindparam,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

import numpy as np

from .data_model import Base
from .connection import DatabaseManager
from .generic import MutableDict

from ..structure.composition import Composition
from ..utils.common_math import DPeak
from ..utils.collectiontools import groupby

PROTON = Composition("H+").mass

observed_ions_bakery = bakery()


class HasPeakChromatogramData(object):

    peak_data = Column(MutableDict.as_mutable(PickleType))

    def get_chromatogram(self):
        peak_data = self.peak_data
        scans = peak_data['scan_times']
        intensity = peak_data['intensities']

        scan_groups = groupby(zip(scans, intensity), key_fn=operator.itemgetter(0), transform_fn=operator.itemgetter(1))
        scans = []
        intensity = []
        for scan, peak_group in scan_groups.items():
            scans.append(scan)
            intensity.append(sum(peak_group))

        intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
        time = np.arange(0, scans.max() + 1)
        abundance_over_time = np.zeros_like(time)
        abundance_over_time[scans] = intensity

        return time, abundance_over_time



class SampleRun(Base):
    __tablename__ = "SampleRun"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(Unicode(32), index=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    ms_scans = relationship("MSScan", backref=backref("sample_run"), lazy='dynamic')
    tandem_scans = relationship("TandemScan", backref=backref("sample_run"), lazy='dynamic')
    sample_type = Column(Unicode(128))

    def __init__(self, **kwargs):
        kwargs.setdefault('parameters', {})
        super(SampleRun, self).__init__(**kwargs)

    def to_json(self):
        d = {
            "name": self.name,
            "id": self.id,
            "parameters": self.parameters,
            "sample_type": self.sample_type
        }
        return d

    def __repr__(self):
        return "<{} {} {} {} {}>".format(
            self.__class__.__name__, self.id, self.name,
            self.ms_scans.count(), self.tandem_scans.count())

    __mapper_args__ = {
        'polymorphic_identity': u'SampleRun',
        'polymorphic_on': sample_type,
    }


class BUPIDDeconvolutedLCMSMSSampleRun(SampleRun):
    __tablename__ = "BUPIDDeconvolutedLCMSMSSampleRun"
    id = Column(Integer, ForeignKey(SampleRun.id), primary_key=True, autoincrement=True)

    __mapper_args__ = {
        "polymorphic_identity": u"BUPIDDeconvolutedLCMSMSSampleRun"
    }


class Decon2LSLCMSSampleRun(SampleRun):
    __tablename__ = "Decon2LSLCMSSampleRun"
    id = Column(Integer, ForeignKey(SampleRun.id), primary_key=True, autoincrement=True)

    __mapper_args__ = {
        "polymorphic_identity": u"Decon2LSLCMSSampleRun"
    }


class ScanBase(Base):
    __tablename__ = "ScanBase"

    id = Column(Integer, primary_key=True, autoincrement=True)
    time = Column(Integer, index=True)
    scan_type = Column(Unicode(20), index=True)
    peaks = relationship("Peak", backref="scan", lazy="dynamic")
    decon2ls_peaks = relationship("Decon2LSPeak", backref="scan", lazy='dynamic')
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)

    __mapper_args__ = {
        'polymorphic_identity': u'Scan',
        'polymorphic_on': scan_type,
    }

    def __repr__(self):
        return "<{} {} {}>".format(self.scan_type, self.id, self.peaks.count())


class MSScan(ScanBase):
    __tablename__ = "MSScan"

    id = Column(Integer, ForeignKey(ScanBase.id), primary_key=True)

    __mapper_args__ = {
        "polymorphic_identity": "MSScan",
    }


class TandemScan(ScanBase):
    __tablename__ = "TandemScan"

    id = Column(Integer, ForeignKey(ScanBase.id), primary_key=True)
    precursor_id = Column(Integer, ForeignKey("Peak.id"), index=True)
    precursor_neutral_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    precursor_charge_state = Column(Integer)

    def __repr__(self):
        return "<{} {} {} {}>".format(self.scan_type, self.id, self.precursor_neutral_mass, self.peaks.count())

    @property
    def tandem_data(self):
        return list(map(DPeak, self.peaks))

    def __iter__(self):
        return iter(map(DPeak, self.peaks))

    @hybrid_method
    def ppm_match_tolerance_search(self, mass, tolerance, filter_callable=lambda q: q):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return filter_callable(self.peaks.filter(Peak.neutral_mass.between(lo, hi)))

    @ppm_match_tolerance_search.expression
    def ppm_match_tolerance_search(cls, mass, tolerance):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return cls.precursor_neutral_mass.between(lo, hi)

    __mapper_args__ = {
        "polymorphic_identity": "TandemScan",
    }


class Peak(Base):
    __tablename__ = "Peak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    charge = Column(Integer)
    neutral_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    intensity = Column(Numeric(12, 6, asdecimal=False))
     
    scan_peak_index = Column(Integer, index=True)

    scan_id = Column(Integer, ForeignKey("ScanBase.id"), index=True)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.join(ScanBase.peaks).filter(ScanBase.sample_run_id == sample_run_id)

    def __repr__(self):
        return "<Peak {} {} {} {}>".format(self.id, self.neutral_mass, self.intensity, self.charge)

    @hybrid_property
    def mz(self):
        return (self.neutral_mass + (self.charge * PROTON)) / self.charge


class Decon2LSPeak(Base):
    __tablename__ = "Decon2LSPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    scan_peak_index = Column(Integer, index=True)

    charge = Column(Integer)
    intensity = Column(Integer, index=True)
    monoisotopic_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    monoisotopic_intensity = Column(Integer)
    monoisotopic_plus_2_intensity = Column(Integer)
    average_mass = Column(Numeric(12, 6, asdecimal=False))
    most_abundant_mass = Column(Numeric(12, 6, asdecimal=False))
    full_width_half_max = Column(Numeric(12, 6, asdecimal=False))
    signal_to_noise = Column(Numeric(12, 6, asdecimal=False))
    isotopic_fit = Column(Numeric(7, 6, asdecimal=False))
    scan_id = Column(Integer, ForeignKey(ScanBase.id), index=True)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.join(ScanBase.decon2ls_peaks).filter(ScanBase.sample_run_id == sample_run_id)

    def __repr__(self):
        return "<Decon2LSPeak {} {} {} {}>".format(self.id, self.monoisotopic_mass, self.intensity, self.charge)


class PeakGroupBase(object):
    charge_state_count = Column(Integer)
    scan_count = Column(Integer)
    first_scan_id = Column(Integer)
    last_scan_id = Column(Integer)

    scan_density = Column(Numeric(10, 6, asdecimal=False))
    weighted_monoisotopic_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

    total_volume = Column(Numeric(12, 4, asdecimal=False))
    average_a_to_a_plus_2_ratio = Column(Numeric(12, 4, asdecimal=False))
    a_peak_intensity_error = Column(Numeric(10, 6, asdecimal=False))
    centroid_scan_estimate = Column(Numeric(12, 4, asdecimal=False))
    centroid_scan_error = Column(Numeric(10, 6, asdecimal=False))
    average_signal_to_noise = Column(Numeric(10, 6, asdecimal=False))



class Decon2LSPeakGroup(Base, HasPeakChromatogramData, PeakGroupBase):
    __tablename__ = "Decon2LSPeakGroup"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)

    matched = Column(Boolean, index=True)
    peaks = relationship(Decon2LSPeak, secondary="Decon2LSPeakToPeakGroupMap", lazy='dynamic')

    def __repr__(self):
        id = getattr(self, 'id', None)
        return "<Decon2LSPeakGroup {} {} {} {} {}>".format(
            id, self.weighted_monoisotopic_mass,
            self.total_volume, self.peaks.count(), self.charge_state_count)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.filter(self.sample_run_id == sample_run_id)


Decon2LSPeakToPeakGroupMap = Table(
    "Decon2LSPeakToPeakGroupMap", Base.metadata,
    Column("peak_id", ForeignKey(Decon2LSPeak.id), index=True),
    Column("group_id", ForeignKey(Decon2LSPeakGroup.id), index=True)
    )


class PeakGroupDatabase(DatabaseManager):
    def _baked_ppm_match_tolerance_query(self):
        build = observed_ions_bakery(lambda session: session.query(Decon2LSPeakGroup))
        build += lambda q: q.filter(Decon2LSPeakGroup.weighted_monoisotopic_mass.between(
            bindparam("lower"), bindparam("upper")))
        build += lambda q: q.filter(Decon2LSPeakGroup.sample_run_id == bindparam("sample_run_id"))
        return build

    def ppm_match_tolerance_search(self, mass, tolerance, sample_run_id=None):
        session = self.session()
        width = (mass * tolerance)
        lower = mass - width
        upper = mass + width
        return self._baked_ppm_match_tolerance_query()(session).params(
            lower=lower, upper=upper, sample_run_id=sample_run_id)


class MSMSSqlDB(DatabaseManager):
    def _baked_ppm_match_tolerance_query(self):
        build = observed_ions_bakery(lambda session: session.query(TandemScan))
        build += lambda q: q.filter(TandemScan.precursor_neutral_mass.between(
            bindparam("lower"), bindparam("upper")))
        build += lambda q: q.filter(TandemScan.sample_run_id == bindparam("sample_run_id"))
        return build

    def ppm_match_tolerance_search(self, mass, tolerance, sample_run_id=None):
        session = self.session()
        width = (mass * tolerance)
        lower = mass - width
        upper = mass + width
        return self._baked_ppm_match_tolerance_query()(session).params(
            lower=lower, upper=upper, sample_run_id=sample_run_id)
