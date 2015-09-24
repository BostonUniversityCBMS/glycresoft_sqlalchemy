from cpython.ref cimport PyObject
from cpython.string cimport PyString_AsString
from libc.stdlib cimport abort, malloc, free, realloc
from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_Append

import operator

cdef:
    object get_intensity = operator.attrgetter("intensity")


cpdef float ppm_error(float x, float y):
    return _ppm_error(x, y)

cdef inline float _ppm_error(float x, float y):
    return (x - y) / y

cpdef object tol_ppm_error(float x, float y, float tolerance):
    cdef float err
    err = (x - y) / y
    if abs(err) <= tolerance:
        return err
    else:
        return None


def test_compiled(training_spectrum, features):
    cdef:
        PeakStructArray* peaks
        PeakStructArray* matches
        MSFeatureStructArray* ms_features
        MSFeatureStruct cfeature
        MassOffsetFeature pfeature
        PeakStruct cpeak
        DPeak dpeak
        size_t i, j
        list py_peaks

    peaks = <PeakStructArray*>malloc(sizeof(PeakStructArray))
    py_peaks = list(training_spectrum)
    intensity_rank(py_peaks)
    peaks.size = len(py_peaks)
    peaks.peaks = <PeakStruct*>malloc(sizeof(PeakStruct) * peaks.size)
    for i in range(peaks.size):
        cpeak = peaks.peaks[i]
        dpeak = <DPeak>py_peaks[i]

        cpeak.neutral_mass = dpeak.neutral_mass
        cpeak.id = dpeak.id
        cpeak.charge = dpeak.charge
        cpeak.intensity = dpeak.intensity
        cpeak.rank = dpeak.rank
        cpeak.mass_charge_ratio = dpeak.mass_charge_ratio
        peaks.peaks[i] = cpeak

    ms_features = <MSFeatureStructArray*>malloc(sizeof(MSFeatureStructArray))
    ms_features.size = len(features)
    ms_features.features = <MSFeatureStruct*>malloc(sizeof(MSFeatureStruct) * ms_features.size)
    for i in range(ms_features.size):
        cfeature = ms_features.features[i]
        pfeature = <MassOffsetFeature>features[i]
        cfeature.offset = pfeature.offset
        cfeature.intensity_ratio = pfeature.intensity_ratio
        cfeature.from_charge = pfeature.from_charge
        cfeature.to_charge = pfeature.to_charge
        cfeature.tolerance = pfeature.tolerance
        cfeature.name = PyString_AsString(pfeature.name)
        ms_features.features[i] = cfeature

    for i in range(peaks.size):
        cpeak = peaks.peaks[i]
        for j in range(ms_features.size):
            cfeature = ms_features.features[j]
            matches = _search_spectrum(&cpeak, peaks, &cfeature)
            if matches.size > 0:
                print cfeature.name, matches.size




cdef inline int _intensity_ratio_function(PeakStruct* peak1, PeakStruct* peak2):
    cdef float ratio
    ratio = peak1.intensity / (peak2.intensity)
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

cdef inline bint feature_match(MSFeatureStruct* feature, PeakStruct* peak1, PeakStruct* peak2):
    if (feature.intensity_ratio == OUT_OF_RANGE_INT or _intensity_ratio_function(peak1, peak2) == feature.intensity_ratio) and\
       ((feature.from_charge == OUT_OF_RANGE_INT and feature.to_charge == OUT_OF_RANGE_INT) or
        (feature.from_charge == peak1.charge and feature.to_charge == peak2.charge)):
        return abs(_ppm_error(peak1.neutral_mass + feature.offset, peak2.neutral_mass)) <= feature.tolerance
    else:
        return False


cdef PeakStructArray* _search_spectrum(PeakStruct* peak, PeakStructArray* peak_list, MSFeatureStruct* feature):
    cdef:
        PeakStructArray* matches
        size_t i, j, n
        PeakStruct query_peak
        PeakStruct* temp
    matches = <PeakStructArray*>malloc(sizeof(PeakStructArray))
    matches.peaks = <PeakStruct*>malloc(sizeof(PeakStruct) * peak_list.size)
    matches.size = peak_list.size
    n = 0
    for i in range(peak_list.size):
        query_peak = peak_list.peaks[i]
        if feature_match(feature, peak, &query_peak):
            matches.peaks[n] = query_peak
            n += 1
    # temp = <PeakStruct*>malloc(sizeof(PeakStruct) * n)
    matches.peaks = <PeakStruct*>realloc(&matches.peaks, sizeof(PeakStruct) * n)
    #for i in range(n):
    #    temp[i] = matches.peaks[i]
    #free(matches.peaks)
    #matches.peaks = temp
    matches.size = n
    return matches


cdef int intensity_ratio_function(DPeak peak1, DPeak peak2):
    cdef float ratio
    ratio = peak1.intensity / (peak2.intensity)
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


cdef void intensity_rank(list peak_list, float minimum_intensity=100.):
    cdef:
        Py_ssize_t i = 0
        int step, rank, tailing
        DPeak p
    peak_list = sorted(peak_list, key=get_intensity, reverse=True)
    step = 0
    rank = 10
    tailing = 6
    for i in range(len(peak_list)):
        p = <DPeak>peak_list[i]
        if p.intensity < minimum_intensity:
            p.rank = 0
            continue
        step += 1
        if step == 10 and rank != 0:
            if rank == 1:
                if tailing != 0:
                    step = 0
                    tailing -= 1
                else:
                    step = 0
                    rank -= 1
            else:
                step = 0
                rank -= 1
        if rank == 0:
            break
        p.rank = rank


cdef int OUT_OF_RANGE_INT = -999

cdef class MassOffsetFeature(object):

    def __init__(self, offset, tolerance, name=None, intensity_ratio=OUT_OF_RANGE_INT, from_charge=OUT_OF_RANGE_INT, to_charge=OUT_OF_RANGE_INT):
        if name is None:
            name = "F:" + str(offset)
            if intensity_ratio is not OUT_OF_RANGE_INT:
                name += ", %r" % (intensity_ratio if intensity_ratio > OUT_OF_RANGE_INT else '')

        self.offset = offset
        self.tolerance = tolerance
        self.name = name
        self.intensity_ratio = intensity_ratio
        self.from_charge = from_charge
        self.to_charge = to_charge

    def __getstate__(self):
        return {
            "offset": self.offset,
            "tolerance": self.tolerance,
            "name": self.name,
            "intensity_ratio": self.intensity_ratio,
            "from_charge": self.from_charge,
            "to_charge": self.to_charge
        }

    def __setstate__(self, d):
        self.name = d['name']
        self.offset = d['offset']
        self.tolerance = d['tolerance']
        self.intensity_ratio = d['intensity_ratio']
        self.from_charge = d['from_charge']
        self.to_charge = d['to_charge']

    def __reduce__(self):
        return MassOffsetFeature, (0, 0), self.__getstate__()

    def __call__(self, DPeak query, DPeak peak):
        return self.test(query, peak)        

    cdef bint test(self, DPeak peak1, DPeak peak2):
        if (self.intensity_ratio == OUT_OF_RANGE_INT or intensity_ratio_function(peak1, peak2) == self.intensity_ratio) and\
           ((self.from_charge == OUT_OF_RANGE_INT and self.to_charge == OUT_OF_RANGE_INT) or
            (self.from_charge == peak1.charge and self.to_charge == peak2.charge)):
            return abs(ppm_error(peak1.neutral_mass + self.offset, peak2.neutral_mass)) <= self.tolerance
        else:
            return False

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash((self.name, self.offset, self.intensity_ratio, self.from_charge, self.to_charge))


cdef class DPeak(object):
    '''
    Defines a type for holding the same relevant information that the database model Peak does
    without the non-trivial overhead of descriptor access on the mapped object to check for
    updated data.
    '''


    def __init__(self, peak=None):
        if peak is not None:
            self.neutral_mass = PyFloat_AsDouble(peak.neutral_mass)
            self.id = PyInt_AsLong(peak.id)
            self.charge = PyInt_AsLong(peak.charge)
            self.intensity = PyFloat_AsDouble(peak.intensity)
        self.peak_relations = []

    def __repr__(self):
        return "<DPeak {} {} {}>".format(self.neutral_mass, self.charge, self.intensity)

    cdef PeakStruct* as_struct(self):
        cdef PeakStruct* result
        result = <PeakStruct*>malloc(sizeof(PeakStruct))
        result.neutral_mass = self.neutral_mass
        result.id = self.id
        result.charge = self.charge
        result.intensity = self.intensity
        result.rank = self.rank
        result.mass_charge_ratio = self.mass_charge_ratio
        return result

    def __getstate__(self):
        cdef dict d
        d = dict()
        d['neutral_mass'] = self.neutral_mass
        d['id'] = self.id
        d['charge'] = self.charge
        d['intensity'] = self.intensity
        d['rank'] = self.rank
        d['mass_charge_ratio'] = self.mass_charge_ratio
        d['peak_relations'] = self.peak_relations
        return d

    def __setstate__(self, dict d):
        self.neutral_mass = d['neutral_mass']
        self.id = d['id']
        self.charge = d['charge']
        self.intensity = d['intensity']
        self.rank = d['rank']
        self.mass_charge_ratio = d['mass_charge_ratio']
        self.peak_relations = d['peak_relations']

    def __reduce__(self):
        return DPeak, (None,), self.__getstate__()

    def __hash__(self):
        return hash((self.neutral_mass, self.intensity, self.charge))

    def __richmp__(self, other, int op):
        cdef bint res = True
        if op == 2:
            res &= self.neutral_mass == other.neutral_mass
            res &= self.intensity == other.intensity
            res &= self.charge == other.charge
            return res
        elif op == 3:
            return not self == other


cpdef DPeak DPeak_from_values(cls, float neutral_mass):
    cdef DPeak peak
    peak = DPeak()
    peak.neutral_mass = neutral_mass
    return peak


cpdef list search_spectrum(DPeak peak, list peak_list, MassOffsetFeature feature):
    cdef:
        list matches = []
        DPeak other_peak
        Py_ssize_t i
        object adder = matches.append
    for i in range(PyList_GET_SIZE(peak_list)):
        other_peak = <DPeak>PyList_GET_ITEM(peak_list, i)
        if feature.test(peak, other_peak):
            adder(other_peak)
    return matches

def pintensity_rank(list peak_list, float minimum_intensity=100.):
    intensity_rank(peak_list, minimum_intensity)

def pintensity_ratio_function(DPeak peak1, DPeak peak2):
    return intensity_ratio_function(peak1, peak2)


cdef class MatchedSpectrum(object):

    def __init__(self, gsm=None):
        if gsm is not None:
            self.peak_match_map = dict(gsm.peak_match_map)
            self.peak_list = list(gsm)
            self.scan_time = gsm.scan_time
            self.peaks_explained = gsm.peaks_explained
            self.peaks_unexplained = gsm.peaks_unexplained
            self.id = gsm.id
            self.glycopeptide_sequence = str(gsm.glycopeptide_sequence)

    def reindex_peak_matches(self):
        mass_map = dict()
        for peak_idx, matches in self.peak_match_map.items():
            mass_map["%0.4f" % matches[0]['observed_mass']] = peak_idx
        for peak in self.peak_list:
            idx = mass_map.get("%0.4f" % peak.neutral_mass)
            if idx is None:
                continue
            peak.id = idx


    def __iter__(self):
        cdef:
            Py_ssize_t i
            DPeak o
        for i in range(PyList_GET_SIZE(self.peak_list)):
            o = <DPeak>PyList_GET_ITEM(self.peak_list, i)
            yield o

    cpdef set peak_explained_by(self, object peak_id):
        cdef:
            set explained
            list matches
            dict match
            PyObject* temp
            Py_ssize_t i

        explained = set()
        temp = PyDict_GetItem(self.peak_match_map, peak_id)
        if temp == NULL:
            return explained
        matches = <list>temp
        for i in range(PyList_GET_SIZE(matches)):
            match = <dict>PyList_GET_ITEM(matches, i)
            explained.add(<str>PyDict_GetItem(match, "key"))

        return explained

    def __getstate__(self):
        d = {}
        d['peak_match_map'] = self.peak_match_map
        d['peak_list'] = self.peak_list
        d['scan_time'] = self.scan_time
        d['peaks_explained'] = self.peaks_explained
        d['peaks_unexplained'] = self.peaks_unexplained
        d['id'] = self.id
        d['glycopeptide_sequence'] = self.glycopeptide_sequence
        return d

    def __setstate__(self, d):
        self.peak_match_map = d['peak_match_map']
        self.peak_list = d['peak_list']
        self.scan_time = d['scan_time']
        self.peaks_explained = d['peaks_explained']
        self.peaks_unexplained = d['peaks_unexplained']
        self.id = d['id']
        self.glycopeptide_sequence = d['glycopeptide_sequence']

    def __reduce__(self):
        return MatchedSpectrum, (None,), self.__getstate__()

    def __repr__(self):
        temp = "<MatchedSpectrum %s @ %d %d|%d>" % (
            self.glycopeptide_sequence, self.scan_time, self.peaks_explained,
            self.peaks_unexplained)
        return temp
