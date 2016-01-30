from cpython.ref cimport PyObject
from cpython.string cimport PyString_AsString, PyString_FromString
from libc.stdlib cimport abort, malloc, free, realloc, calloc
from libc.math cimport fabs
from libc.string cimport strcmp
from libc cimport *

cdef extern from * nogil:
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))
    int printf   (const char *template, ...)

from cython.parallel cimport parallel, prange # openmp must be enabled at compile time
cimport cython

from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble, PyFloat_FromDouble
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Values
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_Append


cpdef double ppm_error(double x, double y):
    return _ppm_error(x, y)


cdef double _ppm_error(double x, double y) nogil:
    return (x - y) / y



cdef int intensity_ratio_function(DPeak peak1, DPeak peak2):
    cdef double ratio
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


cdef int _intensity_ratio_function(PeakStruct* peak1, PeakStruct* peak2) nogil:
    cdef double ratio
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


cdef void _intensity_rank(PeakStructArray* peak_list, double minimum_intensity=100.) nogil:
    cdef:
        size_t i
        int step, rank, tailing
        PeakStruct p
    sort_by_intensity(peak_list)
    step = 0
    rank = 10
    tailing = 6
    for i in range(peak_list.size):
        p = peak_list.peaks[i]
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


cdef void intensity_rank(list peak_list, double minimum_intensity=100.):
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
            self.scan_peak_index = PyInt_AsLong(peak.scan_peak_index) 
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
        result.scan_peak_index = self.scan_peak_index
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


cdef PeakStructArray* unwrap_peak_list(list py_peaks):
    cdef:
        PeakStructArray* peaks
        PeakStruct cpeak
        DPeak dpeak
        size_t i, j

    j = PyList_GET_SIZE(py_peaks)
    peaks = <PeakStructArray*>malloc(sizeof(PeakStructArray))
    intensity_rank(py_peaks)
    peaks.peaks = <PeakStruct*>malloc(sizeof(PeakStruct) * j)
    peaks.size = j
    for i in range(j):
        dpeak = <DPeak>py_peaks[i]
        cpeak.neutral_mass = dpeak.neutral_mass
        cpeak.id = dpeak.id
        cpeak.charge = dpeak.charge
        cpeak.intensity = dpeak.intensity
        cpeak.rank = dpeak.rank
        cpeak.mass_charge_ratio = dpeak.mass_charge_ratio
        cpeak.scan_peak_index = dpeak.scan_peak_index
        peaks.peaks[i] = cpeak
    return peaks


# struct to cdef class wrapper
cdef DPeak wrap_peak(PeakStruct* peak):
    cdef DPeak dpeak = DPeak()
    dpeak.neutral_mass = peak.neutral_mass
    dpeak.charge = peak.charge
    dpeak.intensity = peak.intensity
    dpeak.rank = peak.rank
    dpeak.id = peak.id
    dpeak.scan_peak_index = peak.scan_peak_index
    return dpeak


# Sort PeakStructArray


cdef void sort_by_intensity(PeakStructArray* peak_list) nogil:
    qsort(peak_list.peaks, peak_list.size, sizeof(PeakStruct), compare_by_intensity)

cdef void sort_by_neutral_mass(PeakStructArray* peak_list) nogil:
    qsort(peak_list, peak_list.size, sizeof(PeakStruct), compare_by_neutral_mass)

cdef int compare_by_neutral_mass(const void * a, const void * b) nogil:
    if (<PeakStruct*>a).neutral_mass < (<PeakStruct*>b).neutral_mass:
        return -1
    elif (<PeakStruct*>a).neutral_mass == (<PeakStruct*>b).neutral_mass:
        return 0
    elif (<PeakStruct*>a).neutral_mass > (<PeakStruct*>b).neutral_mass:
        return 1

cdef int compare_by_intensity(const void * a, const void * b) nogil:
    if (<PeakStruct*>a).intensity < (<PeakStruct*>b).intensity:
        return -1
    elif (<PeakStruct*>a).intensity == (<PeakStruct*>b).intensity:
        return 0
    elif (<PeakStruct*>a).intensity > (<PeakStruct*>b).intensity:
        return 1


# Free functions

cdef void free_peak_struct_array(PeakStructArray* peaks) nogil:
    free(peaks.peaks)
    free(peaks)
