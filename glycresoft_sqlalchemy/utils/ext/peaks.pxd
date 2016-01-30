cpdef double ppm_error(double x, double y)

cdef double _ppm_error(double x, double y) nogil


cdef class DPeak(object):
    '''
    Defines a type for holding the same relevant information that the database model Peak does
    without the non-trivial overhead of descriptor access on the mapped object to check for
    updated data.
    '''
    cdef:
        public double neutral_mass
        public long id
        public long scan_peak_index
        public int charge
        public double intensity
        public int rank
        public double mass_charge_ratio
        public list peak_relations

        cdef PeakStruct* as_struct(self)


# Scalar Structs

cdef struct PeakStruct:
    double neutral_mass
    long id
    long scan_peak_index
    int charge
    double intensity
    int rank
    double mass_charge_ratio



# Array Structs

cdef struct PeakStructArray:
    PeakStruct* peaks
    Py_ssize_t size



# unwrap class

cdef PeakStructArray* unwrap_peak_list(list)


# wrap struct

cdef DPeak wrap_peak(PeakStruct* peak)


# Sort PeakStructArray
cdef void sort_by_intensity(PeakStructArray* peak_list) nogil
cdef void sort_by_neutral_mass(PeakStructArray* peak_list) nogil


# Free Functions
cdef void free_peak_struct_array(PeakStructArray* peaks) nogil
