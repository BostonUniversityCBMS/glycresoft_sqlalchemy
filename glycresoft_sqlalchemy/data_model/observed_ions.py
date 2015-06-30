from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy import (PickleType, Numeric, Unicode, Table,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from .data_model import Base, DatabaseManager
from .generic import MutableDict


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
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id), index=True)

    __mapper_args__ = {
        'polymorphic_identity': u'Scan',
        'polymorphic_on': scan_type
    }

    def __repr__(self):
        return "<{} {} {}>".format(self.scan_type, self.id, self.peaks.count())


class MSScan(ScanBase):
    __tablename__ = "MSScan"

    id = Column(Integer, ForeignKey(ScanBase.id), primary_key=True)

    __mapper_args__ = {
        "polymorphic_identity": "MSScan"
    }


class TandemScan(ScanBase):
    __tablename__ = "TandemScan"

    id = Column(Integer, ForeignKey(ScanBase.id), primary_key=True)
    precursor_id = Column(Integer, ForeignKey("Peak.id"), index=True)
    precursor_neutral_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

    def __repr__(self):
        return "<{} {} {} {}>".format(self.scan_type, self.id, self.precursor_neutral_mass, self.peaks.count())

    @property
    def tandem_data(self):
        return self.peaks.all()

    @hybrid_method
    def ppm_match_tolerance_search(self, mass, tolerance, mass_shift=0):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return (lo <= self.precursor_neutral_mass) & (self.precursor_neutral_mass <= hi)

    @ppm_match_tolerance_search.expression
    def ppm_match_tolerance_search(cls, mass, tolerance, mass_shift=0):
        spread = mass * tolerance
        hi = mass + spread
        lo = mass - spread
        return cls.precursor_neutral_mass.between(lo, hi)

    __mapper_args__ = {
        "polymorphic_identity": "TandemScan"
    }


class Peak(Base):
    __tablename__ = "Peak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    charge = Column(Integer)
    neutral_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    intensity = Column(Numeric(12, 6, asdecimal=False))

    scan_id = Column(Integer, ForeignKey("ScanBase.id"), index=True)

    def __repr__(self):
        return "<Peak {} {} {} {}>".format(self.id, self.neutral_mass, self.intensity, self.charge)


class MSMSSqlDB(DatabaseManager):

    def ppm_match_tolerance_search(self, mass, tolerance, mass_shift=0):
        session = self.session()
        return session.query(TandemScan).filter(TandemScan.ppm_match_tolerance_search(mass, tolerance, mass_shift))
