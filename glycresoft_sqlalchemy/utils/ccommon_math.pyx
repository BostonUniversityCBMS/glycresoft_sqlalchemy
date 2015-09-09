from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE

import operator

cdef:
    object get_intensity = operator.attrgetter("intensity")


cpdef float ppm_error(float x, float y):
    return (x - y) / y


cpdef object tol_ppm_error(float x, float y, float tolerance):
    cdef float err
    err = (x - y) / y
    if abs(err) <= tolerance:
        return err
    else:
        return None


cpdef int mass_offset_match(float mass, DPeak peak, float offset=0., float tolerance=2e-5):
    return abs(ppm_error(mass + offset, peak.neutral_mass)) <= tolerance


cdef inline int intensity_ratio_function(DPeak peak1, DPeak peak2):
    cdef float ratio
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


cdef class MassOffsetFeature(object):
    cdef:
        public float offset
        public float tolerance
        public str name
        public int intensity_ratio
        public int intensity_rank

    def __init__(self, offset, tolerance, name=None, intensity_ratio=-1, intensity_rank=-1):
        if name is None:
            name = "F:" + str(offset)
            if intensity_ratio is not None:
                name += ", %r" % (intensity_ratio if intensity_ratio > 0 else '')
        if intensity_ratio is None:
            intensity_ratio = -1
        if intensity_rank is None:
            intensity_ratio = -1
        self.offset = offset
        self.tolerance = tolerance
        self.name = name
        self.intensity_ratio = intensity_ratio
        self.intensity_rank = intensity_rank

    def __getstate__(self):
        return {
            "offset": self.offset,
            "tolerance": self.tolerance,
            "name": self.name,
            "intensity_ratio": self.intensity_ratio,
            "intensity_rank": self.intensity_rank
        }

    def __setstate__(self, d):
        self.name = d['name']
        self.offset = d['offset']
        self.tolerance = d['tolerance']
        self.intensity_ratio = d['intensity_ratio']
        self.intensity_rank = d['intensity_rank']

    def __call__(self, DPeak query, DPeak peak):
        return self.test(query, peak)        

    cdef bint test(self, DPeak peak1, DPeak peak2):
        if self.intensity_ratio > -1 or intensity_ratio_function(peak1, peak2) == self.intensity_ratio:
            return abs(ppm_error(peak1.neutral_mass + self.offset, peak2.neutral_mass)) <= self.tolerance
        else:
            return False

    def __repr__(self):
        return self.name


cdef class PooledOffsetFeature(MassOffsetFeature):
    cdef:
        public list offsets

    def __init__(self, offsets, tolerance, name=None, intensity_ratio=-1, intensity_rank=-1):
        if name is None:
            name = "F:" + ', '.join(map(str, offsets))
            if intensity_ratio is not None:
                name += ", %r" % (intensity_ratio if intensity_ratio > 0 else '')
        if intensity_ratio is None:
            intensity_ratio = -1
        if intensity_rank is None:
            intensity_ratio = -1
        self.offsets = list(offsets)
        self.tolerance = tolerance
        self.name = name
        self.intensity_ratio = intensity_ratio
        self.intensity_rank = intensity_rank

    cdef bint test(self, DPeak peak1, DPeak peak2):
        cdef Py_ssize_t i = 0
        cdef float offset
        if self.intensity_ratio > 0 or intensity_ratio_function(peak1, peak2) == self.intensity_ratio:
            for i in range(len(self.offsets)):
                offset = self.offsets[i]
                if abs(ppm_error(peak1.neutral_mass + offset, peak2.neutral_mass)) <= self.tolerance:
                    return True
            return False
        else:
            return False

cdef class DPeak(object):
    '''
    Defines a type for holding the same relevant information that the database model Peak does
    without the non-trivial overhead of descriptor access on the mapped object to check for
    updated data.
    '''
    cdef:
        public float neutral_mass
        public long id
        public int charge
        public float intensity
        public int rank
        public float mass_charge_ratio
        public list peak_relations

    def __init__(self, peak=None):
        if peak is not None:
            self.neutral_mass = PyFloat_AsDouble(peak.neutral_mass)
            self.id = PyInt_AsLong(peak.id)
            self.charge = PyInt_AsLong(peak.charge)
            self.intensity = PyFloat_AsDouble(peak.intensity)
        self.peak_relations = []


cpdef DPeak DPeak_from_values(cls, float neutral_mass):
    cdef DPeak peak
    peak = DPeak()
    peak.neutral_mass = neutral_mass
    return peak


# DPeak.from_values = classmethod(from_values)


cdef class TheoreticalDPeak(DPeak):
    cdef:
        public str name


cpdef list search_spectrum(DPeak peak, list peak_list, MassOffsetFeature feature, float tolerance=2e-5):
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
