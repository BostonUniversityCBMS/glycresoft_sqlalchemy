from cpython.sequence cimport PySequence_ITEM
from cpython.string cimport PyString_AsString
from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.list cimport PyList_GET_ITEM, PyList_Append, PyList_GET_SIZE
from cpython.object cimport PyObject_CallFunctionObjArgs

from libc.stdlib cimport malloc, free, realloc


from glycresoft_sqlalchemy.utils.ccommon_math cimport (
    DPeak, MassOffsetFeature, ppm_error, intensity_rank,
    intensity_ratio_function, search_spectrum, MatchedSpectrum,
    MatchedSpectrumStruct, matched_spectrum_struct_peak_explained_by,
    _intensity_rank, _search_spectrum, MatchedSpectrumStructArray,
    MatchedSpectrumStruct, _intensity_ratio_function, unwrap_matched_spectrum,
    FragmentMatchStructArray, FragmentMatchStruct, PeakStructArray, unwrap_feature_functions,
    MSFeatureStructArray
    )


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
        return PeakRelation, (self.from_peak, self.to_peak, self.feature, self.intensity_ratio, self.kind)


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


cdef inline void make_peak_relation_struct(PeakStruct* from_peak, PeakStruct* to_peak, 
                                           MSFeatureStruct* feature, float intensity_ratio,
                                           char* kind, PeakRelationStruct* out) nogil:
    out.from_peak = from_peak
    out.to_peak = to_peak
    out.feature = feature
    out.intensity_ratio = intensity_ratio
    out.kind = kind


cdef inline void free_peak_relations_array(PeakRelationStructArray* relations) nogil:
    free(relations.relations)
    free(relations)


cdef FittedFeatureStruct* _feature_function_estimator(MatchedSpectrumStructArray* gsms, MSFeatureStruct* feature, char* kind):
    cdef:
        double total_on_kind_satisfied
        double total_off_kind_satisfied
        double total_on_kind_satisfied_normalized
        double total_off_kind_satisfied_normalized
        double temp
        double total_on_kind
        double total_off_kind
        int relation_count, peak_count
        bint is_on_kind
        size_t i, j, k, relations_counter, features_found_counter
        PeakStruct peak
        PeakStruct match
        PeakStruct stack_peak
        PeakRelationStructArray* related
        RelationSpectrumPairArray* peak_relations
        PeakStructArray* peaks
        PeakStructArray* matches
        PeakStruct* _peaks
        PeakRelationStruct* pr
        MatchedSpectrumStruct* gsm
        FragmentMatchStructArray* peak_explained_by_list
        FittedFeatureStruct* result

    total_on_kind_satisfied = 0
    total_off_kind_satisfied = 0
    total_on_kind = 0
    total_off_kind = 0
    features_found_counter = 0

    peak_relations = <RelationSpectrumPairArray*>malloc(sizeof(RelationSpectrumPairArray))
    peak_relations.pairs = <RelationSpectrumPair*>malloc(sizeof(RelationSpectrumPair) * gsms.size)
    peak_relations.size = gsms.size

    for i in range(gsms.size):
        print "On GSM", i
        gsm = &gsms.matches[i]
        peaks = gsm.peak_list
        _intensity_rank(peaks)
        peak_count = 0
        print "Filtering Peaks"
        _peaks = <PeakStruct*>malloc(sizeof(PeakStruct) * peaks.size)
        for j in range(peaks.size):
            stack_peak = peaks.peaks[j]
            if stack_peak.rank > 0:
                _peaks[peak_count] = stack_peak
                peak_count += 1
        print "Freed leftover peaks"
        free(peaks.peaks)
        peaks.peaks = _peaks
        peaks.size = peak_count

        # Pre-allocate space for a handful of relations assuming that there won't
        # be many per spectrum
        print "Allocating PeakRelationStructs"
        j = peaks.size / 5
        related = <PeakRelationStructArray*>malloc(sizeof(PeakRelationStructArray))
        related.relations = <PeakRelationStruct*>malloc(sizeof(PeakRelationStruct) * j)
        related.size = j

        # Search all peaks against all other peaks using the given feature
        for j in range(peaks.size):
            print "On Peak", j
            peak = peaks.peaks[j]

            # Determine if peak has been identified as on-kind or off-kind
            peak_explained_by_list = matched_spectrum_struct_peak_explained_by(gsm, peak.id)
            is_on_kind = False
            for k in range(peak_explained_by_list.size):
                print "Explained Ion", k
                if peak_explained_by_list.matches[k].key[0] == kind[0]:
                    is_on_kind = True
                    break
            print "Freeing list"
            free(peak_explained_by_list)

            # Perform search
            print "Performing matching"
            matches = _search_spectrum(&peak, peaks, feature)
            relations_counter = 0
            # For each match, construct a PeakRelationStruct and save it to `related`
            for k in range(matches.size):
                match = matches.peaks[k]
                print "On match", k, related.size
                # Make sure not to overflow `related.relations`. Should be rare because of pre-allocation
                if relations_counter == related.size:
                    related.relations = <PeakRelationStruct*>realloc(related.relations, related.size * 2)
                    related.size = related.size * 2
                make_peak_relation_struct(
                    &peak, &match, feature, _intensity_ratio_function(&peak, &match),
                    NULL, &related.relations[relations_counter])
                if is_on_kind:
                    total_on_kind_satisfied += 1
                    related.relations[relations_counter].kind = kind
                else:
                    total_off_kind_satisfied += 1
                    related.relations[relations_counter].kind = "Noise"

                # Increment counter
                relations_counter += 1

            related.relations = <PeakRelationStruct*>realloc(related.relations, relations_counter)
            related.size = relations_counter
            if is_on_kind:
                total_on_kind += 1
            else:
                total_off_kind += 1
            # Save results or free un-needed memory
            if relations_counter > 0:
                peak_relations.pairs[features_found_counter].matched_spectrum = gsm
                peak_relations.pairs[features_found_counter].relations = related
            else:
                free_peak_relations_array(related)

    # Ensure no division by zero
    total_on_kind = total_on_kind if total_on_kind > 0 else 1.0
    total_off_kind = total_off_kind if total_off_kind > 0 else 1.0

    # Compute feature parameters
    total_on_kind_satisfied_normalized = total_on_kind_satisfied / total_on_kind
    total_off_kind_satisfied_normalized = total_off_kind_satisfied / total_off_kind


    result = <FittedFeatureStruct*>malloc(sizeof(FittedFeatureStruct))

    result.on_kind = total_on_kind_satisfied_normalized
    result.off_kind = total_off_kind_satisfied_normalized
    result.relation_pairs = peak_relations

    return result

def test_compiled(list gsms, MassOffsetFeature feature, str kind):
    cdef:
        size_t i
        MatchedSpectrumStructArray* spectra
        MSFeatureStruct cfeature
        MSFeatureStructArray* cfeatures
        char* ckind
        FittedFeatureStruct* out
    spectra = <MatchedSpectrumStructArray*>malloc(sizeof(MatchedSpectrumStructArray))
    spectra.matches = <MatchedSpectrumStruct*>malloc(sizeof(MatchedSpectrumStruct) * len(gsms))
    spectra.size = len(gsms)

    print "Unwrapping gsms"
    for i in range(spectra.size):
        spectra.matches[i] = unwrap_matched_spectrum(gsms[i])[0]

    print "Unwrapping feature"
    cfeatures = unwrap_feature_functions([feature])
    print "1"
    cfeature = cfeatures.features[0]
    print "Unwrapping kind"
    ckind = PyString_AsString(kind)
    print ckind
    out = _feature_function_estimator(spectra, &cfeature, ckind)

