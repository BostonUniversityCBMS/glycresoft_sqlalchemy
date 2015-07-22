from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .data_model import Base
from .connection import DatabaseManager
from .generic import MutableDict, MutableList


class SampleRun(Base):
    __tablename__ = "SampleRun"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(Unicode(128))
    parameters = Column(MutableDict.as_mutable(PickleType))
    ms_scans = relationship("MSScan", backref=backref("hypothesis_run"), lazy='dynamic')
    tandem_scans = relationship("TandemScan", backref=backref("hypothesis_run"), lazy='dynamic')

    def __repr__(self):
        return "<SampleRun {} {} {} {}>".format(self.id, self.name, self.ms_scans.count(), self.tandem_scans.count())


class ScanBase(Base):
    __tablename__ = "ScanBase"

    id = Column(Integer, primary_key=True, autoincrement=True)
    scan_type = Column(Unicode(20), index=True)
    peaks = relationship("Peak", backref="scan", lazy="dynamic")
    decon2ls_peaks = relationship("Decon2LSPeak", backref="scan", lazy='dynamic')
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)

    __mapper_args__ = {
        'polymorphic_identity': u'Scan',
        'polymorphic_on': scan_type,
        'concrete': True
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
        return self.peaks.all()

    @hybrid_method
    def ppm_match_tolerance_search(self, mass, tolerance):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return (lo <= self.precursor_neutral_mass) & (self.precursor_neutral_mass <= hi)

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

    scan_id = Column(Integer, ForeignKey("ScanBase.id"), index=True)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.join(ScanBase.peaks).filter(ScanBase.sample_run_id == sample_run_id)

    def __repr__(self):
        return "<Peak {} {} {} {}>".format(self.id, self.neutral_mass, self.intensity, self.charge)


class Decon2LSPeak(Base):
    __tablename__ = "Decon2LSPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    charge = Column(Integer)
    intensity = Column(Integer, index=True)
    monoisotopic_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    monoisotopic_intensity = Column(Integer)
    monoisotopic_plus_2_intensity = Column(Integer)
    average_mass = Column(Numeric(12, 6, asdecimal=False))
    most_abundant_mass = Column(Numeric(12, 6, asdecimal=False))
    full_width_half_max = Column(Numeric(12, 6, asdecimal=False))
    signal_to_noise = Column(Numeric(12, 6, asdecimal=False))

    scan_id = Column(Integer, ForeignKey(ScanBase.id), index=True)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.join(ScanBase.decon2ls_peaks).filter(ScanBase.sample_run_id == sample_run_id)

    def __repr__(self):
        return "<Decon2LSPeak {} {} {} {}>".format(self.id, self.monoisotopic_mass, self.intensity, self.charge)


class Decon2LSPeakGroup(Base):
    __tablename__ = "Decon2LSPeakGroup"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)

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

    matched = Column(Boolean, index=True)
    peaks = relationship(Decon2LSPeak, secondary="Decon2LSPeakToPeakGroupMap", lazy='dynamic')

    def __repr__(self):
        return "<Decon2LSPeakGroup {} {} {} {} {}>".format(
            self.id, self.weighted_monoisotopic_mass,
            self.total_volume, self.peaks.count(), self.charge_state_count)

    @hybrid_method
    def ppm_match_tolerance_search(self, mass, tolerance):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return (lo <= self.weighted_monoisotopic_mass) & (self.weighted_monoisotopic_mass <= hi)

    @ppm_match_tolerance_search.expression
    def ppm_match_tolerance_search(cls, mass, tolerance):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return cls.weighted_monoisotopic_mass.between(lo, hi)

    @classmethod
    def from_sample_run(self, query, sample_run_id):
        return query.filter(self.sample_run_id == sample_run_id)


Decon2LSPeakToPeakGroupMap = Table(
    "Decon2LSPeakToPeakGroupMap", Base.metadata,
    Column("peak_id", ForeignKey(Decon2LSPeak.id), index=True),
    Column("group_id", ForeignKey(Decon2LSPeakGroup.id), index=True)
    )


class PeakGroupDatabase(DatabaseManager):
    def ppm_match_tolerance_search(self, mass, tolerance, sample_run_id=None):
        session = self.session()
        query = session.query(Decon2LSPeakGroup).filter(
            Decon2LSPeakGroup.ppm_match_tolerance_search(mass, tolerance))
        if sample_run_id is not None:
            query = query.filter(Decon2LSPeakGroup.sample_run_id == sample_run_id)
        return query


class MSMSSqlDB(DatabaseManager):
    def ppm_match_tolerance_search(self, mass, tolerance, sample_run_id=None):
        session = self.session()
        query = session.query(TandemScan).filter(
            TandemScan.ppm_match_tolerance_search(mass, tolerance))
        if sample_run_id is not None:
            query = query.filter(TandemScan.sample_run_id == sample_run_id)
        return query
