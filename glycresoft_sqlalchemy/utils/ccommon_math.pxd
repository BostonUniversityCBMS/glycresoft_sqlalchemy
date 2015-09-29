
cpdef float ppm_error(float x, float y)

cpdef object tol_ppm_error(float x, float y, float tolerance)


cdef inline bint feature_match(MSFeatureStruct* feature, PeakStruct* peak1, PeakStruct* peak2) nogil

cdef int intensity_ratio_function(DPeak peak1, DPeak peak2)

cdef void intensity_rank(list peak_list, float minimum_intensity=*)
cdef void _intensity_rank(PeakStructArray* peak_list, float minimum_intensity=*) nogil


cdef class MassOffsetFeature(object):
    cdef:
        public float offset
        public float tolerance
        public str name
        public int intensity_ratio
        public int from_charge
        public int to_charge
        public str feature_type

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


cdef public struct PeakStruct:
    float neutral_mass
    long id
    int charge
    float intensity
    int rank
    float mass_charge_ratio


cdef public struct MSFeatureStruct:
    float offset
    float tolerance
    char* name
    int intensity_ratio
    int from_charge
    int to_charge
    char* feature_type


cdef public struct FragmentMatchStruct:
    float observed_mass
    float intensity
    char* key
    long peak_id


cdef public struct PeakStructArray:
    PeakStruct* peaks
    Py_ssize_t size


cdef public struct MSFeatureStructArray:
    MSFeatureStruct* features
    Py_ssize_t size


cdef public struct FragmentMatchStructArray:
    FragmentMatchStruct* matches
    size_t size

cpdef list search_spectrum(DPeak peak, list peak_list, MassOffsetFeature feature)


cdef class MatchedSpectrum(object):
    cdef:
        public dict peak_match_map
        #: Should be a list of DPeak instances
        public list peak_list
        public str glycopeptide_sequence
        public int scan_time
        public int peaks_explained
        public int peaks_unexplained
        public int id

    cpdef set peak_explained_by(self, object peak_id)


cdef public struct MatchedSpectrumStruct:
    FragmentMatchStructArray* peak_match_list
    PeakStructArray* peak_list
    char* glycopeptide_sequence
    int scan_time
    int peaks_explained
    int peaks_unexplained
    int id


cdef public struct MatchedSpectrumStructArray:
    MatchedSpectrumStruct* matches
    size_t size


cdef MatchedSpectrumStruct* unwrap_matched_spectrum(MatchedSpectrum ms)
cdef MSFeatureStructArray* unwrap_feature_functions(list features)


cdef FragmentMatchStructArray* matched_spectrum_struct_peak_explained_by(MatchedSpectrumStruct* ms, long peak_id) nogil

cdef void sort_by_intensity(PeakStructArray* peak_list) nogil

cdef int _intensity_ratio_function(PeakStruct* peak1, PeakStruct* peak2) nogil

cdef PeakStructArray* _search_spectrum(PeakStruct* peak, PeakStructArray* peak_list, MSFeatureStruct* feature) nogil
