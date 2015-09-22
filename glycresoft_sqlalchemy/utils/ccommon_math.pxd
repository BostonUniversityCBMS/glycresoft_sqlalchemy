
cpdef float ppm_error(float x, float y)

cpdef object tol_ppm_error(float x, float y, float tolerance)

cpdef int mass_offset_match(float mass, DPeak peak, float offset=*, float tolerance=*)


cdef int intensity_ratio_function(DPeak peak1, DPeak peak2)

cdef void intensity_rank(list peak_list, float minimum_intensity=*)


cdef class MassOffsetFeature(object):
    cdef:
        public float offset
        public float tolerance
        public str name
        public int intensity_ratio
        public int intensity_rank

    cdef bint test(self, DPeak peak1, DPeak peak2)

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

        cdef PeakStruct* as_struct(self)

cpdef DPeak DPeak_from_values(cls, float neutral_mass)


cdef struct PeakStruct:
    float neutral_mass
    long id
    int charge
    float intensity
    int rank
    float mass_charge_ratio

cpdef list search_spectrum(DPeak peak, list peak_list, MassOffsetFeature feature, float tolerance=?)
