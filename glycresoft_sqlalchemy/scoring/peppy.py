import numpy as np
from scipy.misc import comb


def binomial_tail_probability(n, k, p):
    total = 0
    for i in range(k, n):
        total += comb(n, i) * (p ** i) * (1 - p) ** (n - i)
    return total


def peppy_binomial_fragments_matched(total_product_ion_count, count_product_ion_matches, ion_tolerance,
                                     precursor_mass):
    p = np.exp((np.log(ion_tolerance) + np.log(2)) +
               np.log(count_product_ion_matches) - np.log(precursor_mass))
    return binomial_tail_probability(total_product_ion_count, count_product_ion_matches, p)


def count_theoretical_product_ions(theoretical):
    total = len(theoretical.bare_b_ions)
    total += len(theoretical.bare_y_ions)
    total += len(theoretical.glycosylated_b_ions)
    total += len(theoretical.glycosylated_y_ions)
    total += len(theoretical.stub_ions)
    total += len(theoretical.oxonium_ions)
    return total


def medians(array):
    m1 = np.median(array)
    m2 = np.median(array[array > m1])
    m3 = np.median(array[array > m2])
    m4 = np.median(array[array > m3])
    return m1, m2, m3, m4


def peppy_binomial_intensity(peak_list, matched_peaks, total_product_ion_count):
    intensity_list = np.array([p.intensity for p in peak_list])
    m1, m2, m3, m4 = medians(intensity_list)

    matched_intensities = np.array(
        [match['intensity'] for p, matches in matched_peaks.items() for match in matches])
    counts = dict()
    last_count = total_product_ion_count
    next_count = (matched_intensities > m1).sum()
    counts[1] = binomial_tail_probability(last_count, next_count, 0.5)
    last_count = next_count

    next_count = (matched_intensities > m2).sum()
    counts[2] = binomial_tail_probability(last_count, next_count, 0.5)
    last_count = next_count

    next_count = (matched_intensities > m3).sum()
    counts[3] = binomial_tail_probability(last_count, next_count, 0.5)
    last_count = next_count

    next_count = (matched_intensities > m4).sum()
    counts[4] = binomial_tail_probability(last_count, next_count, 0.5)

    prod = 0
    for v in counts.values():
        prod += np.log(v)
    return np.exp(prod)


def peppy_score(spectrum, spectrum_match, match_tolerance):
    theoretical = spectrum_match.theoretical_glycopeptide
    total_theoretical_product_ions = count_theoretical_product_ions(
        theoretical)
    precursor_mass = theoretical.calculated_mass - \
        theoretical.glycan_composition.mass()
    fragment_match_component = peppy_binomial_fragments_matched(
        total_theoretical_product_ions,
        len(spectrum_match.peak_match_map),
        match_tolerance,
        precursor_mass
    )

    intensity_component = peppy_binomial_intensity(
        spectrum,
        spectrum_match.peak_match_map,
        total_theoretical_product_ions)
    return -np.log10(intensity_component * fragment_match_component)
