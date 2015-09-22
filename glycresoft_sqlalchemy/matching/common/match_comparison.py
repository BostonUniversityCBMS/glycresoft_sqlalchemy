from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, GlycopeptideMatch, object_session,
    HypothesisSampleMatch, GlycopeptideSpectrumMatch, TandemScan)

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
