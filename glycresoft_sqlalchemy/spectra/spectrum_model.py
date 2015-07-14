import logging

from ..structure.composition import Composition

PROTON = Composition("H+").mass
db_logger = logging.getLogger(__name__)


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


def mass_charge_ratio(neutral_mass, z):
    return (neutral_mass + (z * PROTON)) / z
