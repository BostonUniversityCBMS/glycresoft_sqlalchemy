from collections import defaultdict

from ..structure.composition import Composition
from ..utils.common_math import DPeak, ppm_error

PROTON = Composition("H+").mass


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


def mass_charge_ratio(neutral_mass, z):
    return (neutral_mass + (z * PROTON)) / z


def aggregate_spectra(spectra_collection, match_tolerance=2e-5):
    peak_collection = defaultdict(list)

    for spectrum in spectra_collection:
        for peak in spectrum:
            extended = False
            for mz in list(peak_collection):
                if abs(ppm_error(peak.mass_charge_ratio, mz)) <= match_tolerance:
                    peak_collection[mz].append(peak)
                    extended = True
            if not extended:
                peak_collection[peak.mass_charge_ratio].append(peak)
    return peak_collection
