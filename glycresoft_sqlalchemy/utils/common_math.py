from math import fabs
import operator
try:
    from ccommon_math import (
        ppm_error, tol_ppm_error, DPeak, MassOffsetFeature,
        search_spectrum, pintensity_ratio_function as intensity_ratio_function,
        pintensity_rank as intensity_rank, MatchedSpectrum)

except ImportError, e:
    print "common_math", e

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

    OUT_OF_RANGE_INT = -999

    class MassOffsetFeature(object):

        def __init__(self, offset, tolerance, name=None, intensity_ratio=OUT_OF_RANGE_INT,
                     from_charge=OUT_OF_RANGE_INT, to_charge=OUT_OF_RANGE_INT, feature_type=''):
            if name is None:
                name = "F:" + str(offset)
            if intensity_ratio is not OUT_OF_RANGE_INT:
                name += ", %r" % (intensity_ratio if intensity_ratio > OUT_OF_RANGE_INT else '')

            self.name = name
            self.offset = offset
            self.tolerance = tolerance
            self.intensity_ratio = intensity_ratio
            self.from_charge = from_charge
            self.to_charge = to_charge
            self.feature_type = feature_type

        def test(self, peak1, peak2):
            if (self.intensity_ratio == OUT_OF_RANGE_INT or intensity_ratio_function(peak1, peak2) == self.intensity_ratio) and\
               ((self.from_charge == OUT_OF_RANGE_INT and self.to_charge == OUT_OF_RANGE_INT) or
                (self.from_charge == peak1.charge and self.to_charge == peak2.charge)):

                return abs(ppm_error(peak1.neutral_mass + self.offset, peak2.neutral_mass)) <= self.tolerance
            return False

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

    class MatchedSpectrum(object):

        def __init__(self, gsm):
            self.peak_match_map = gsm.peak_match_map
            self.peak_list = list(gsm)
            self.scan_time = gsm.scan_time
            self.peaks_explained = gsm.peaks_explained
            self.peaks_unexplained = gsm.peaks_unexplained
            self.id = gsm.id
            self.glycopeptide_sequence = gsm.glycopeptide_sequence

        def __iter__(self):
            for o in self.peak_list:
                yield o

        def peak_explained_by(self, peak_id):
            explained = set()
            try:
                matches = self.peak_match_map[peak_id]
                for match in matches:
                    explained.add(match["key"])
            except KeyError:
                pass
            return explained

        def __repr__(self):
            temp = "<MatchedSpectrum %s @ %d %d|%d>" % (
                self.glycopeptide_sequence, self.scan_time, self.peaks_explained,
                self.peaks_unexplained)
            return temp

        def reindex_peak_matches(self):
            mass_map = dict()
            for peak_idx, matches in self.peak_match_map.items():
                mass_map["%0.4f" % matches[0]['observed_mass']] = peak_idx
            for peak in self.peak_list:
                idx = mass_map.get("%0.4f" % peak.neutral_mass)
                if idx is None:
                    continue
                peak.id = idx
