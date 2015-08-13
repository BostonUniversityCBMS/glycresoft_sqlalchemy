from glycresoft_sqlalchemy.data_model import HypothesisSampleMatch, SampleRun, SpectrumMatch, DatabaseManager
from sqlalchemy.orm import make_transient


def load(source, target, hsm_id, sample_run_id):
    source_db = DatabaseManager(source)
    target_db = DatabaseManager(target)
    source_session = source_db.session()
