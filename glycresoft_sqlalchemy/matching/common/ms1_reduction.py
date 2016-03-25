import logging
try:
    logger = logging.getLogger("ms1_reduction")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

from glycresoft_sqlalchemy.data_model import (
    PipelineModule,
    SampleRun, Decon2LSPeakGroup,
    TandemScan, make_transient
)


class MS1Reduction(PipelineModule):
    def __init__(self, ms1_database_path, ms2_database_path, output_path, tolerance=2e-5):
        self.ms1 = self.manager_type(ms1_database_path)
        self.ms2 = self.manager_type(ms2_database_path)
        self.output = self.manager_type(output_path)
        self.tolerance = tolerance

    def reduce_peak_groups(self):
        reduced = set()
        tolerance = self.tolerance
        ms1 = self.ms1
        ms2 = self.ms2
        for tandem_scan in ms2.query(TandemScan):
            target_mass = tandem_scan.precursor_neutral_mass
            window_width = target_mass * tolerance
            lo = target_mass - window_width
            hi = target_mass + window_width
            rel_1_2 = ms1.query(Decon2LSPeakGroup).filter(
                Decon2LSPeakGroup.weighted_monoisotopic_mass.between(lo, hi)).all()
            if rel_1_2:
                reduced.update(rel_1_2)
        return reduced

    def save_extracted(self, sample_run, reduced_peaks):
        self.output.initialize()
        out = self.output.session()
        make_transient(sample_run)
        out.add(sample_run)
        list(map(make_transient, reduced_peaks))
        out.add_all(reduced_peaks)
        out.commit()
        self.inform("Extracted %d peak groups", len(reduced_peaks))
        return sample_run.id

    def run(self):
        peaks = self.reduce_peak_groups()
        sample_run = self.ms1.query(SampleRun).first()
        return self.save_extracted(sample_run, peaks)


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser("ms1_reduction")
    app.add_argument("ms1_database_path")
    app.add_argument("ms2_database_path")
    app.add_argument("output_path")

    args = app.parse_args()
    job = MS1Reduction(**args.__dict__)
    job.start()
