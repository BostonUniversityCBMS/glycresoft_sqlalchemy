import numpy as np

from glycresoft_sqlalchemy.utils import collectiontools
from collections import Counter, defaultdict

from .utils import (
    Sequence, MatchedSpectrum,
    chain_iterable, ppm_error, MassOffsetFeature,
    DPeak, intensity_ratio_function,
    intensity_rank, search_spectrum)


from glycresoft_sqlalchemy.utils.memoize import memoize
from glycresoft_sqlalchemy.structure.sequence import (
    Sequence, fragment_shift, Fragment, Composition,
    structure_constants)


# Lacking a reasonable definition of the "space between fragmentation sites"
SMALLEST_UNIT = 1000 * 2e-5
NOISE = "noise"

def preprocess_glycopeptide_spectrum_match(gsms):
    return map(MatchedSpectrum, gsms)


def _delta_series(sequence, mass_shift=0.0, step=1):
    sequence = Sequence(sequence)

    # The set of modification names.
    mod_dict = {}

    # The key is the position, the value is an array of fragments.
    # And the first element is always bare fragment.
    # The total number of HexNAc on the fragment should be recorded.

    sequence = sequence[::step]

    current_mass = 0
    for idx in range(len(sequence) - 1):
        for mod in sequence[idx][1]:
            mod_serial = mod.serialize()
            if mod_serial in mod_dict:
                mod_dict[mod_serial] += 1
            else:
                mod_dict[mod_serial] = 1

        if idx == 0:
            current_mass = sequence[0][0].mass + mass_shift
        else:
            current_mass = current_mass + sequence[idx][0].mass

        frag_dri = []
        # If incremental loss of HexNAc is not allowed, only one fragment of a given type is generated
        if not structure_constants.PARTIAL_HEXNAC_LOSS:
            frag = Fragment(str(mass_shift), idx + structure_constants.FRAG_OFFSET, (mod_dict).copy(), current_mass)
            frag_dri.append(frag)
            bare_dict = (mod_dict).copy()
            bare_dict["HexNAc"] = 0
            frag = Fragment(str(mass_shift), idx + structure_constants.FRAG_OFFSET, (bare_dict).copy(), current_mass)
            frag_dri.append(frag)
        # Else a fragment for each incremental loss of HexNAc must be generated
        else:
            frag = Fragment(str(mass_shift), idx + structure_constants.FRAG_OFFSET, (mod_dict).copy(), current_mass)
            frag_dri.extend(frag.partial_loss())
        yield frag_dri


def delta_finder(gsms, delta, step=1):
    total_sites = 0
    total_explained = 0

    total_sparsity = 0
    total_unexplained = 0
    proton = Composition("H+").mass
    for gsm in gsms:
        fragments = list(collectiontools.flatten(_delta_series(gsm.glycopeptide_sequence, delta, step)))
        n_frag_sites = len(fragments)
        peak_list = list(gsm)
        match_count = 0
        for fragment in fragments:
            for peak in gsm:
                if abs(ppm_error(fragment.mass - proton, peak.neutral_mass)) <= 2e-5:
                    match_count += 1

        total_explained += match_count
        total_sites += n_frag_sites
        total_unexplained += len(peak_list) - match_count
        total_sparsity += mass_accuracy_sparsity(gsm.glycopeptide_sequence, 'b') - n_frag_sites

    return total_explained / float(total_sites), total_unexplained / float(total_sparsity)


# BEGIN REFERENCE Estimation
# Simple one parameter estimator functions for learning the basic alpha, beta and p parameters


def offset_frequency(gsms, kind='b'):
    total_sites = 0
    total_explained = 0
    for gsm in gsms:
        n_frag_sites = count_fragmentation_sites(gsm.glycopeptide_sequence, kind)
        kind_explained = sum([1 for i in chain_iterable(gsm.peak_match_map.values()) if i['key'][0] == kind])
        total_sites += n_frag_sites
        total_explained += kind_explained
    return (total_explained)/float(total_sites)


def unknown_peak_rate(gsms, kind='b'):
    total_sparsity = 0
    total_unexplained = 0

    for gsm in gsms:
        sequence = gsm.glycopeptide_sequence
        n_frag_sites = count_fragmentation_sites(sequence, kind)
        kind_explained = sum([1 for i in chain_iterable(gsm.peak_match_map.values()) if i['key'][0] == kind])
        peaks_unexplained = gsm.peaks_unexplained + (gsm.peaks_explained - kind_explained)
        total_unexplained += peaks_unexplained
        total_sparsity += mass_accuracy_sparsity(sequence, kind) - n_frag_sites

    return total_unexplained / float(total_sparsity)


def count_fragmentation_sites(sequence, kind='b'):
    sequence = Sequence(sequence)
    fragmentation_sites = len(collectiontools.flatten(sequence.get_fragments(kind)))
    return fragmentation_sites


def prior_fragment_probability(gsms, kind='b'):
    hits = 0
    for gsm in gsms:
        sequence = Sequence(gsm.glycopeptide_sequence)
        random_mass = np.random.uniform(56., sequence.mass)
        for fragment in collectiontools.flatten(sequence.get_fragments(kind)):
            if abs(ppm_error(fragment.mass, random_mass)) <= 2e-5:
                hits += 1
            # elif fragment.mass - (random_mass + 230.) > 0:
            #     break
    return hits / float(len(gsms))

# END REFERENCE Estimation


def estimate_fragment_sparsity(sequence, kind='b'):
    return Sequence(sequence).mass / SMALLEST_UNIT


def mass_accuracy_sparsity(sequence, kind='b', tolerance=2e-5):
    mass_space, fragment_space = mass_dimensions(sequence, kind, tolerance)
    return mass_space - fragment_space


@memoize(200)
def mass_dimensions(sequence, kind='b', tolerance=2e-5):
    sequence = Sequence(sequence)
    mass_space = sequence.mass
    fragment_space = 0.
    for frag in chain_iterable(sequence.get_fragments(kind)):
        radius = tolerance * frag.mass
        fragment_space += 2 * radius
    return mass_space, fragment_space


def estimate_offset_parameters(gsms, kind='b'):

    total_sites = 0
    total_explained = 0
    total_sparsity = 0
    total_unexplained = 0
    mass_space_area = 0
    fragment_space_area = 0

    i = 0.
    for gsm in gsms:
        sequence = Sequence(gsm.glycopeptide_sequence)
        fragments = collectiontools.flatten(sequence.get_fragments(kind))

        n_frag_sites = len(fragments)

        kind_explained = sum([1 for d in chain_iterable(gsm.peak_match_map.values()) if d['key'][0] == kind])
        total_sites += n_frag_sites
        total_explained += kind_explained

        peaks_unexplained = gsm.peaks_unexplained + (gsm.peaks_explained - kind_explained)

        total_unexplained += peaks_unexplained
        total_sparsity += estimate_fragment_sparsity(sequence, kind) - n_frag_sites

        s, f = mass_dimensions(gsm.glycopeptide_sequence, kind=kind)
        mass_space_area += s
        fragment_space_area += f

        i += 1

    alpha = total_explained / float(total_sites)
    beta = total_unexplained / float(total_sparsity)
    prior_fragment_probability = fragment_space_area / mass_space_area
    return alpha, beta, prior_fragment_probability


def probability_of_peak_explained(offset_frequency, unknown_peak_rate, prior_fragment_probability):
    a = (prior_fragment_probability * offset_frequency)
    b = (1 - prior_fragment_probability) * unknown_peak_rate
    return a / (a + b)


def fragmentation_probability(peak, probability_of_fragment, features):
    probability_of_noise = 1 - probability_of_fragment
    numer = 1
    denom = 1
    if len(features) == 0:
        return probability_of_fragment
    for feature in features:
        match = probability_of_fragment * feature.on_kind
        total = match + (probability_of_noise * feature.off_kind)
        numer *= match
        denom *= total
    return numer / denom


def make_offset_function(offset=.0, tolerance=2e-5, name=None, **kwargs):
    return MassOffsetFeature(offset=offset, tolerance=tolerance, name=name, **kwargs)


class PeakRelation(object):
    def __init__(self, from_peak, to_peak, feature, intensity_ratio=None, kind=None):
        if intensity_ratio is None:
            intensity_ratio = intensity_ratio_function(from_peak, to_peak)
        self.from_peak = from_peak
        self.to_peak = to_peak
        self.feature = feature
        self.intensity_ratio = intensity_ratio
        self.same_terminal = None
        self.from_charge = from_peak.charge
        self.to_charge = to_peak.charge
        self.kind = kind or NOISE

    def __repr__(self):
        template = "<PeakRelation {s.from_peak.neutral_mass}({s.from_charge}) ->" +\
            " {s.to_peak.neutral_mass}({s.to_charge}) by {s.feature.name} on {s.kind}>"
        return template.format(s=self)


def feature_function_estimator(gsms, feature_function, kind='b'):
    total_on_kind_satisfied = 0.
    total_off_kind_satisfied = 0.
    total_on_kind = 0.
    total_off_kind = 0.
    peak_relations = []
    for gsm in gsms:
        peaks = list(map(DPeak, gsm))
        intensity_rank(peaks)
        peaks = [p for p in peaks if p.rank > 0]
        related = []
        for peak in peaks:
            is_on_kind = any(k[0] == kind for k in gsm.peak_explained_by(peak.id))
            matches = search_spectrum(peak, peaks, feature_function)
            for match in matches:
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

    total_on_kind_satisfied_normalized = total_on_kind_satisfied / max(total_on_kind, 1)
    total_off_kind_satisfied_normalized = total_off_kind_satisfied / max(total_off_kind, 1)

    return total_on_kind_satisfied_normalized, total_off_kind_satisfied_normalized, peak_relations


def search_features_on_spectrum(peak, peak_list, features):
    peak_match_relations = defaultdict(list)
    for query_peak in peak_list:
        for feature in features:
            if feature(peak, query_peak):
                match_list = peak_match_relations[peak.id]
                pr = PeakRelation(peak, query_peak, feature, intensity_ratio_function(peak, query_peak), feature.kind)
                match_list.append(pr)
    return peak_match_relations


def search_features(peak_list, features):
    peak_match_map = {}
    for peak in peak_list:
        r = search_features_on_spectrum(peak, peak_list, features)
        peak_match_map.update(r)
    return peak_match_map


class FittedFeature(object):
    def __init__(self, feature, kind, on_kind, off_kind, relations=None):
        if relations is None:
            relations = []
        self.feature = feature
        self.kind = kind
        self.on_kind = on_kind
        self.off_kind = off_kind
        self.relations = relations

    @property
    def name(self):
        return self.feature.name

    def __hash__(self):
        return hash((self.feature, self.kind))

    def __repr__(self):
        temp = "<FittedFeature ({feature}) u:{on_kind} v:{off_kind} @ {kind} {count_relations}>"
        return temp.format(
            feature=self.feature, on_kind=self.on_kind, off_kind=self.off_kind,
            kind=self.kind, count_relations=len(self.relations))

    def charge_relations(self):
        counter = Counter()
        for rel in self:
            counter[rel.from_charge, rel.to_charge] += 1
        return counter

    def intensity_ratio(self):
        counter = Counter()
        for rel in self:
            counter[intensity_ratio_function(rel.from_peak, rel.to_peak)] += 1
        return counter

    def charge_intensity_ratio(self):
        counter = Counter()
        for rel in self:
            counter[(rel.from_charge, rel.to_charge), intensity_ratio_function(rel.from_peak, rel.to_peak)] += 1
        return counter

    def peak_relations(self, include_noise=True):
        for spectrum_match, peak_relations in self.relations:
            for pr in peak_relations:
                if not include_noise and pr.kind == NOISE:
                    continue
                yield pr

    def __iter__(self):
        return self.peak_relations(False)

    def __call__(self, *args, **kwargs):
        return self.feature(*args, **kwargs)


class RelatedSpectrum(object):
    def __init__(self, spectrum, peak_relations=None):
        if peak_relations is None:
            peak_relations = {}
        self.peak_list = list(spectrum)
        self.peak_match_map = spectrum.peak_match_map
        self.peak_relation_map = peak_relations

    def peak_relationships(self, peak_id):
        rels = self.peak_relation_map[peak_id]
        return rels
