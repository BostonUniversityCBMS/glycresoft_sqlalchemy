from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, GlycopeptideMatch, object_session,
    HypothesisSampleMatch, GlycopeptideSpectrumMatch, TandemScan)

import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict


class PairedMatch(object):
    def __init__(self, subject=None, reference=None, scan_time=None):
        self.subject = subject
        self.reference = reference
        self.scan_time = scan_time

    def __repr__(self):
        return "<PairedMatch>\n%r\n%r" % (self.subject, self.reference)


class MatchComparator(object):
    def __init__(self, subject, reference):
        self.subject_session = subject
        self.reference_session = reference

    def get_reference_spectra(self, match):
        subject_spectra = match.spectrum_matches.filter(GlycopeptideSpectrumMatch.best_match).all()
        reference_spectra = self.reference_session.query(
            GlycopeptideSpectrumMatch).join(GlycopeptideMatch).filter(
            GlycopeptideSpectrumMatch.scan_time.in_([t.scan_time for t in subject_spectra]),
            GlycopeptideSpectrumMatch.best_match).filter(GlycopeptideMatch.is_not_decoy()).all()

        time_pairs = defaultdict(PairedMatch)
        for sm in subject_spectra:
            sm.label = "subject"
            time_pairs[sm.scan_time].subject = sm

        for sm in reference_spectra:
            sm.label = "reference"
            time_pairs[sm.scan_time].reference = sm
        return time_pairs


def describe_match(case):
    print case
    print case.mean_coverage, case.mean_hexnac_coverage
    print '\n'.join(' '.join(sorted(i['key'] for i in series)) for series in (
        case.bare_b_ions, case.bare_y_ions,
        case.glycosylated_b_ions, case.glycosylated_y_ions,
        case.stub_ions, case.oxonium_ions))

    fig = (pd.DataFrame(
        [{"time": h.scan_time, "peaks_matched": h.peaks_explained} for h in sorted(case.spectrum_matches, key=lambda x: x.scan_time)]).plot(
        x='time', y='peaks_matched', kind='scatter').get_figure())
    return fig


def describe_match_ipython(case):
    from IPython.display import display
    print case
    print case.mean_coverage, case.mean_hexnac_coverage
    print '\n'.join(' '.join(sorted(i['key'] for i in series)) for series in (
        case.bare_b_ions, case.bare_y_ions,
        case.glycosylated_b_ions, case.glycosylated_y_ions,
        case.stub_ions, case.oxonium_ions))

    fig = (pd.DataFrame(
        [{"time": h.scan_time, "peaks_matched": h.peaks_explained} for h in sorted(case.spectrum_matches, key=lambda x: x.scan_time)]).plot(
        x='time', y='peaks_matched', kind='scatter').get_figure())
    display(fig)
    plt.close(fig)
