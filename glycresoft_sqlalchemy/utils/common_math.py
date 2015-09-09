from math import fabs
try:
    from ccommon_math import ppm_error, tol_ppm_error, DPeak, mass_offset_match, MassOffsetFeature, search_spectrum

except ImportError, e:
    print e

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
