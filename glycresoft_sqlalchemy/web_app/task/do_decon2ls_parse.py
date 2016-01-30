import os
from glycresoft_sqlalchemy.spectra.decon2ls_sa import Decon2LSIsosParser
from .task_process import NullPipe, Message, Task


def taskmain(database_path, isos_path, results_path, comm=None):
    # manager = DatabaseManager(database_path)
    if comm is None:
        comm = NullPipe()
    comm.send(Message("Begin conversion for %s." % isos_path))

    try:
        job = Decon2LSIsosParser(isos_path, results_path)
    except Exception, e:
        comm.send(Message(e, 'error'))
    else:
        comm.send(Message("Finished indexing %s" % os.path.basename(results_path), 'update'))
        comm.send(Message(job.sample_run.to_json(), "new-sample"))
    try:
        os.remove(isos_path)
    except Exception, e:
        comm.send(Message(["Could not remove isos file"]))


class Decon2LSIsosParseTask(Task):
    def __init__(self, database_path, isos_path, results_path, callback, **kwargs):
        args = (database_path, isos_path, results_path)
        job_name = "Decon2LS " + os.path.basename(results_path)
        kwargs.setdefault('name', job_name)
        super(Decon2LSIsosParseTask, self).__init__(taskmain, args, callback, **kwargs)
