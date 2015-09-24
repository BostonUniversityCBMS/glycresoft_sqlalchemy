from cpython.sequence cimport PySequence_ITEM
from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.list cimport PyList_GET_ITEM, PyList_Append, PyList_GET_SIZE
from cpython.object cimport PyObject_CallFunctionObjArgs

from glycresoft_sqlalchemy.utils.ccommon_math cimport (
    DPeak, MassOffsetFeature, ppm_error, intensity_rank,
    intensity_ratio_function, search_spectrum, MatchedSpectrum)


cdef class PeakRelation(object):

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

    def __getstate__(self):
        d = {}
        d["from_peak"] = self.from_peak
        d["to_peak"] = self.to_peak
        d["from_charge"] = self.from_charge
        d["to_charge"] = self.to_charge
        d["feature"] = self.feature
        d["intensity_ratio"] = self.intensity_ratio
        d["kind"] = self.kind
        return d

    def __setstate__(self, d):
        self.from_peak = d["from_peak"]
        self.to_peak = d["to_peak"]
        self.from_charge = d["from_charge"]
        self.to_charge = d["to_charge"]
        self.feature = d["feature"]
        self.intensity_ratio = d["intensity_ratio"]
        self.kind = d["kind"]

    def __reduce__(self):
        return PeakRelation__new__, (PeakRelation,), self.__getstate__()


cdef object PeakRelation__new__ = PeakRelation.__new__


cdef inline PeakRelation make_peak_relation(DPeak from_peak, DPeak to_peak, MassOffsetFeature feature, float intensity_ratio, str kind=None):
    cdef PeakRelation pr
    pr = <PeakRelation>PyObject_CallFunctionObjArgs(PeakRelation__new__, <PyObject*>PeakRelation, NULL)
    pr.from_peak = from_peak
    pr.to_peak = to_peak
    pr.from_charge = from_peak.charge
    pr.to_charge = to_peak.charge
    pr.feature = feature
    pr.intensity_ratio = intensity_ratio
    pr.kind = kind
    return pr



cpdef dict search_features_on_spectrum(DPeak peak, list peak_list, list features, str kind=None):
    cdef:
        dict peak_match_relations = {}
        list match_list
        PyObject* temp
        Py_ssize_t peak_ind, feature_ind
        DPeak query_peak
        MassOffsetFeature feature
        PeakRelation pr

    for i in range(PyList_GET_SIZE(peak_list)):
        query_peak = <DPeak>PyList_GET_ITEM(peak_list, i)
        for j in range(PyList_GET_SIZE(features)):
            feature = <MassOffsetFeature>PyList_GET_ITEM(features, j)
            if feature.test(peak, query_peak):
                temp = PyDict_GetItem(peak_match_relations, feature)
                if temp == NULL:
                    match_list = []
                else:
                    match_list = <list>temp
                pr = make_peak_relation(peak, query_peak, feature, intensity_ratio_function(peak, query_peak), kind)
                match_list.append(pr)
                if temp == NULL:
                    PyDict_SetItem(peak_match_relations, feature, match_list)
    return peak_match_relations


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
        PeakRelation pr
        bint is_on_kind
        set peak_explained_by
        list peak_explained_by_list
        dict db_search_match
        str match_kind
        MatchedSpectrum gsm

    total_on_kind_satisfied = 0
    total_off_kind_satisfied = 0
    total_on_kind = 0
    total_off_kind = 0
    peak_relations = []
    for i in range(PyList_GET_SIZE(gsms)):
        gsm = <MatchedSpectrum>PyList_GET_ITEM(gsms, i)
        peaks = gsm.peak_list
        intensity_rank(peaks)
        _peaks = []
        # peaks = [p for p in peaks if p.rank > 0]
        for j in range(PyList_GET_SIZE(peaks)):
            peak = <DPeak>PyList_GET_ITEM(peaks, j)
            rank = peak.rank
            if rank > 0:
                PyList_Append(_peaks, peak)
        peaks = _peaks
        related = []
        for j in range(PyList_GET_SIZE(peaks)):
            peak = <DPeak>PyList_GET_ITEM(peaks, j)
            peak_explained_by = gsm.peak_explained_by(peak.id)
            peak_explained_by_list = list(peak_explained_by)
            is_on_kind = False
            # is_on_kind = any([pg[0] == kind for pg in gsm.peak_explained_by(peak.id)])
            for k in range(PyList_GET_SIZE(peak_explained_by_list)):
                match_kind = <str>PyList_GET_ITEM(peak_explained_by_list, k)
                match_kind = <str>PySequence_ITEM(match_kind, 0)
                if match_kind == kind:
                    is_on_kind = True
                    break

            matches = search_spectrum(peak, peaks, feature_function)
            for k in range(PyList_GET_SIZE(matches)):
                match = <DPeak>PyList_GET_ITEM(matches, k)
                pr = make_peak_relation(peak, match, feature_function, intensity_ratio_function(peak, match))
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
        if PyList_GET_SIZE(related) > 0:
            peak_relations.append((gsm, related))

    total_on_kind_satisfied_normalized = total_on_kind_satisfied / <float>(max(total_on_kind, 1))
    total_off_kind_satisfied_normalized = total_off_kind_satisfied / <float>(max(total_off_kind, 1))

    return total_on_kind_satisfied_normalized, total_off_kind_satisfied_normalized, peak_relations
