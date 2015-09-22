from math import fabs
import operator
try:
    from ccommon_math import (
        ppm_error, tol_ppm_error, DPeak, mass_offset_match, MassOffsetFeature,
        search_spectrum, pintensity_ratio_function as intensity_ratio_function,
        pintensity_rank as intensity_rank)

except ImportError, e:
    print e

    get_intensity = operator.attrgetter("intensity")

    def ppm_error(x, y):
        return (x - y) / y

    def tol_ppm_error(x, y, tolerance):
        err = (x - y) / y
        if fabs(err) <= tolerance:
            return err
        else:
            return None

    def mass_offset_match(mass, peak, offset=0., tolerance=2e-5):
        return abs(ppm_error(mass + offset, peak.neutral_mass)) <= tolerance

    class DPeak(object):
        '''
        Defines a type for holding the same relevant information that the database model Peak does
        without the non-trivial overhead of descriptor access on the mapped object to check for
        updated data.
        '''
        def __init__(self, peak):
            self.neutral_mass = peak.neutral_mass
            self.id = peak.id
            self.charge = peak.charge
            self.intensity = peak.intensity
            self.rank = 0
            self.peak_relations = []

    class MassOffsetFeature(object):

        def __init__(self, offset, tolerance):
            self.offset = offset
            self.tolerance = tolerance

        def test(self, mass, peak):
            return abs(ppm_error(mass + self.offset, peak.neutral_mass)) <= self.tolerance

        __call__ = test

    def search_spectrum(peak, peak_list, feature, tolerance=2e-5):
        matches = []
        for peak in peak_list:
            if feature(peak, peak):
                matches.append(peak)
        return matches

    def intensity_ratio_function(peak1, peak2):
        ratio = peak1.intensity / float(peak2.intensity)
        if ratio >= 5:
            return -4
        elif 2.5 <= ratio < 5:
            return -3
        elif 1.7 <= ratio < 2.5:
            return -2
        elif 1.3 <= ratio < 1.7:
            return -1
        elif 1.0 <= ratio < 1.3:
            return 0
        elif 0.8 <= ratio < 1.0:
            return 1
        elif 0.6 <= ratio < 0.8:
            return 2
        elif 0.4 <= ratio < 0.6:
            return 3
        elif 0.2 <= ratio < 0.4:
            return 4
        elif 0. <= ratio < 0.2:
            return 5

    def intensity_rank(peak_list, minimum_intensity=100.):
        peak_list = sorted(peak_list, key=get_intensity, reverse=True)
        i = 0
        rank = 10
        tailing = 6
        for p in peak_list:
            if p.intensity < minimum_intensity:
                p.rank = 0
                continue
            i += 1
            if i == 10 and rank != 0:
                if rank == 1:
                    if tailing != 0:
                        i = 0
                        tailing -= 1
                    else:
                        i = 0
                        rank -= 1
                else:
                    i = 0
                    rank -= 1
            if rank == 0:
                break
            p.rank = rank
