from collections import OrderedDict
from sqlalchemy.orm import relationship
from sqlalchemy import (
    Numeric, Unicode, Column, Integer, ForeignKey, PickleType)

import numpy as np

from ..base import Namespace, Base
from ..hypothesis_sample_match import HypothesisSampleMatch


GENERIC_MODEL_NAME = "Generic Scoring Model"


class PeakGroupScoringModel(Base):
    __tablename__ = "PeakGroupScoringModel"

    id = Column(Integer, primary_key=True)
    name = Column(Unicode(30), index=True)
    source_hypothesis_sample_match_id = Column(Integer, ForeignKey(HypothesisSampleMatch.id), index=True)
    source = relationship(HypothesisSampleMatch)

    intercept = Column(Numeric(8, asdecimal=False))
    modification_state_count = Column(Numeric(8, asdecimal=False))
    charge_state_count = Column(Numeric(8, asdecimal=False))
    total_volume = Column(Numeric(8, asdecimal=False))
    a_peak_intensity_error = Column(Numeric(8, asdecimal=False))
    centroid_scan_error = Column(Numeric(8, asdecimal=False))
    average_signal_to_noise = Column(Numeric(8, asdecimal=False))
    scan_density = Column(Numeric(8, asdecimal=False))
    scan_count = Column(Numeric(8, asdecimal=False))

    def to_parameter_vector(self):
        vector = [
            self.charge_state_count,
            self.scan_density,
            self.modification_state_count,
            self.total_volume,
            self.a_peak_intensity_error,
            self.centroid_scan_error,
            self.scan_count,
            self.average_signal_to_noise
        ]
        return np.array([vector])

    @classmethod
    def from_parameter_vector(cls, vector):
        vector = vector[0]
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
        return cls(**dict(zip(names, vector)))

    def to_parameter_dict(self):
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
        return OrderedDict(zip(names, self.to_parameter_vector()[0]))

    def __repr__(self):
        return "PeakGroupScoringModel(%r)" % self.to_parameter_dict()

    @classmethod
    def initialize(cls, session):
        legacy_model_coefficients = [
            [
                0.65784006, 1.345376317, 0.219899787,
                0.0000383503, -0.000349839, -0.001610525,
                -0.000947516, 0.011453828
            ]
        ]
        intercept = -3.69288417
        inst = cls.from_parameter_vector(legacy_model_coefficients)
        inst.intercept = intercept
        inst.name = GENERIC_MODEL_NAME
        session.add(inst)

        session.commit()

    GENERIC_MODEL_NAME = GENERIC_MODEL_NAME

Namespace.initializer(PeakGroupScoringModel.initialize)
