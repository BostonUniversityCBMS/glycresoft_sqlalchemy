import os
from glycresoft_sqlalchemy.data_model import DatabaseManager
from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser
from .task_process import NullPipe, Message, Task


def taskmain(database_path, yaml_path, results_path, comm=None):
    # manager = DatabaseManager(database_path)
    if comm is None:
        comm = NullPipe()
    comm.send(Message("Begin conversion for %s." % yaml_path))
    try:
        BUPIDMSMSYamlParser(yaml_path, results_path)
    except Exception, e:
        comm.send(Message(e, 'error'))
    else:
        comm.send(Message("Finished indexing %s" % os.path.basename(results_path), 'update'))
    try:
        os.remove(yaml_path)
    except:
        comm.send(Message("Could not remove yaml file"))


class BUPIDYamlParseTask(Task):
    def __init__(self, database_path, yaml_path, results_path, callback, **kwargs):
        args = (database_path, yaml_path, results_path)
        job_name = "BUPID " + os.path.basename(results_path)
        kwargs.setdefault('name', job_name)
        super(BUPIDYamlParseTask, self).__init__(taskmain, args, callback, **kwargs)
