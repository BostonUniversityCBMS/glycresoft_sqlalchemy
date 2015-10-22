from glycresoft_sqlalchemy.data_model import session, HypothesisSampleMatch
from glycresoft_sqlalchemy.report import export_csv

from .task_process import NullPipe, Message, Task


def taskmain(database_path, hypothesis_sample_match_id, filterfunc=lambda q: q,
             tempdir=None, comm=NullPipe(), **kwargs):
    comm.send(Message("Begin CSV export", type='info'))
    job = export_csv.CSVExportDriver(
        database_path, [hypothesis_sample_match_id], output_path=tempdir, filterfunc=filterfunc)
    paths = job.start()
    paths = [p for path_group in list(paths) for p in path_group]
    comm.send(Message({
        "hypothesis_sample_match_id": hypothesis_sample_match_id,
        "files": paths
        }, "files-to-download"))
    return


class ExportCSVTask(Task):
    def __init__(self, database_path, hypothesis_sample_match_id, filterfunc=lambda q: q,
                 tempdir=None, **kwargs):
        args = (database_path, hypothesis_sample_match_id, filterfunc, tempdir)
        hsm_name = session(database_path).query(HypothesisSampleMatch).get(hypothesis_sample_match_id).name
        name = "Export CSV {}".format(hsm_name or hypothesis_sample_match_id)
        kwargs.setdefault('name', name)
        super(ExportCSVTask, self).__init__(taskmain, args, **kwargs)
