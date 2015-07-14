from glycresoft_sqlalchemy.data_model import (DatabaseManager, Hypothesis, TheoreticalGlycopeptide,
                                              TandemScan, Peak, SampleRun, PipelineModule, Protein)


def theoretical_sequence_to_spectrum(theoretical):
    ts = TandemScan()
    ts.precursor_neutral_mass = theoretical.calculated_mass
    ts.peaks = [Peak(neutral_mass=f['mass'], charge=1, intensity=1) for f in theoretical.fragments()]
    return ts


class SequenceDatabaseToSpectrumDatabase(PipelineModule):
    def __init__(self, database_path, hypothesis_id, output_path):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.output_manager = self.manager_type(output_path)

    def run(self):
        in_session = self.manager.session()
        self.output_manager.initialize()
        out_session = self.output_manager.session()
        hypothesis = in_session.query(Hypothesis).get(self.hypothesis_id)
        sample_run = SampleRun(name='converted-{}'.format(hypothesis.name))
        out_session.add(sample_run)
        out_session.flush()

        sample_run_id = sample_run.id
        counter = 0
        for theoretical in in_session.query(TheoreticalGlycopeptide).filter(
                TheoreticalGlycopeptide.protein_id == Protein.id,
                Protein.hypothesis_id == hypothesis.id):
            ts = theoretical_sequence_to_spectrum(theoretical)
            ts.sample_run_id = sample_run_id
            counter += 1
            out_session.add(ts)
            if counter % 1000 == 0:
                out_session.commit()
        out_session.commit()

        in_session.close()
        out_session.close()
