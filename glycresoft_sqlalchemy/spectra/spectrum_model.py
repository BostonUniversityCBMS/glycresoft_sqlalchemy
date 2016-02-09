from collections import defaultdict

from ..structure.composition import Composition
from ..utils.common_math import DPeak, ppm_error

PROTON = Composition("H+").mass


def neutral_mass(mz, z, charge_carrier=PROTON):
    return (mz * abs(z)) - (z * charge_carrier)


def mass_charge_ratio(neutral_mass, z, charge_carrier=PROTON):
    return (neutral_mass + (z * charge_carrier)) / abs(z)


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


def normalize_spectrum(spectrum):
    peaks = [p for p in list(spectrum)]
    max_intensity = max([p.intensity for p in peaks])
    cloned = [p.clone() for p in peaks]
    for p in cloned:
        p.intensity /= max_intensity
    return cloned
