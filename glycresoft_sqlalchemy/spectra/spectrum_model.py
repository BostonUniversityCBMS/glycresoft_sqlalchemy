from operator import attrgetter
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


class PeakSet(object):
    _charge = 0

    @classmethod
    def from_list(cls, peaks, charge=None):
        if charge is None:
            charge = peaks[0].charge
        inst = cls(peaks)
        inst.charge = charge
        return inst

    def __init__(self, peaks):
        self.peaks = list(peaks)
        self._charge = peaks[0].charge

    def __iter__(self):
        return iter(self.peaks)

    def __getitem__(self, i):
        return self.peaks[i]

    def __setitem__(self, i, v):
        self.peaks[i] = v

    def add(self, peak):
        self.peaks.append(peak)
        peak.mz = mass_charge_ratio(neutral_mass(peak.mz, peak.charge), self._charge)
        peak.charge = self._charge
        self.peaks.sort(key=attrgetter("mz"))

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = value
        self._adjust_charge(value)

    def _adjust_charge(self, charge):
        for peak in self:
            peak.mz = mass_charge_ratio(neutral_mass(peak.mz, peak.charge), charge)
            peak.charge = charge

    def scale(self, factor):
        for peak in self:
            peak.intensity *= factor
