from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem
from cpython.list cimport PyList_GET_ITEM
from glycresoft_sqlalchemy.utils.ccommon_math cimport (
    DPeak, MassOffsetFeature, ppm_error, intensity_rank,
    intensity_ratio_function, search_spectrum)


cdef class PeakRelation(object):
    cdef:
        public DPeak from_peak
        public DPeak to_peak
        public int from_charge
        public int to_charge
        public float intensity_ratio
        public str kind
        public MassOffsetFeature feature
        public bint same_terminal

    def __init__(self, DPeak from_peak, DPeak to_peak, MassOffsetFeature feature, float intensity_ratio, str kind=None):
        self.from_peak = from_peak
        self.to_peak = to_peak
        self.from_charge = from_peak.charge
        self.to_charge = to_peak.charge
        self.feature = feature
        self.intensity_ratio = intensity_ratio
        self.kind = kind

    def __repr__(self):
        template = "<PeakRelation {s.from_peak.neutral_mass}({s.from_charge}) ->" +\
            " {s.to_peak.neutral_mass}({s.to_charge}) by {s.feature.name} on {s.kind}>"
        return template.format(s=self)

cpdef tuple feature_function_estimator(list gsms, MassOffsetFeature feature_function, str kind='b'):
    cdef:
        DPeak peak
        DPeak match
        float total_on_kind_satisfied
        float total_off_kind_satisfied
        float total_on_kind_satisfied_normalized
        float total_off_kind_satisfied_normalized
        long total_on_kind
        long total_off_kind
        int rank
        list peak_relations
        list related
        list peaks, _peaks
        list matches
        Py_ssize_t i, j, k
        object spectrum
        object gsm
        PeakRelation pr
        bint is_on_kind
        set peak_explained_by
        list peak_explained_by_list
        dict db_search_match
        str match_kind

    total_on_kind_satisfied = 0
    total_off_kind_satisfied = 0
    total_on_kind = 0
    total_off_kind = 0
    peak_relations = []
    for i in range(len(gsms)):
        gsm = gsms[i]
        spectrum = gsm.spectrum
        peaks = list(spectrum)
        intensity_rank(peaks)
        _peaks = []
        # peaks = [p for p in peaks if p.rank > 0]
        for j in range(len(peaks)):
            peak = <DPeak>peaks[j]
            rank = peak.rank
            if rank > 0:
                _peaks.append(peak)
        peaks = _peaks
        related = []
        for j in range(len(peaks)):
            peak = peaks[j]
            peak_explained_by = gsm.peak_explained_by(peak.id)
            peak_explained_by_list = list(peak_explained_by)
            is_on_kind = False
            # is_on_kind = any([pg[0] == kind for pg in gsm.peak_explained_by(peak.id)])
            for k in range(len(peak_explained_by_list)):
                match_kind = <str>PyList_GET_ITEM(peak_explained_by_list, k)
                match_kind = match_kind[0]
                if match_kind == kind:
                    is_on_kind = True
                    break

            matches = search_spectrum(peak, peaks, feature_function)
            for k in range(len(matches)):
                match = <DPeak>PyList_GET_ITEM(matches, k)
                pr = PeakRelation(peak, match, feature_function, intensity_ratio_function(peak, match))
                related.append(pr)
                if is_on_kind:
                    total_on_kind_satisfied += 1
                    pr.kind = kind
                else:
                    total_off_kind_satisfied += 1
                    pr.kind = "Noise"
            if is_on_kind:
                total_on_kind += 1
            else:
                total_off_kind += 1
        if len(related) > 0:
            peak_relations.append((gsm, related))

    total_on_kind_satisfied_normalized = total_on_kind_satisfied / <float>(max(total_on_kind, 1))
    total_off_kind_satisfied_normalized = total_off_kind_satisfied / <float>(max(total_off_kind, 1))

    return total_on_kind_satisfied_normalized, total_off_kind_satisfied_normalized, peak_relations
