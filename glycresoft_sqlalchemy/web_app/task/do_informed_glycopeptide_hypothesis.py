from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder import integrated_omics
from .task_process import NullPipe, Message, Task
from ..common import get_random_string


def taskmain(database_path, hypothesis_name, protein_file, site_list_file,
             glycan_file, glycan_file_type, constant_modifications,
             comm=None, **kwargs):
    if comm is None:
        comm = NullPipe()
    manager = DatabaseManager(database_path)
    try:
        task = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
            database_path, hypothesis_name=hypothesis_name, mzid_path=protein_file,
            glycomics_path=glycan_file, glycomics_format=glycan_file_type,
            n_processes=kwargs.get("n_processes", 4))
        hypothesis_id = task.start()
        if task.status != 0:
            raise task.error
        session = manager.session()
        hypothesis = session.query(Hypothesis).get(hypothesis_id)
        comm.send(Message(hypothesis.to_json(), "new-hypothesis"))
        return hypothesis_id
    except Exception:
        comm.send(Message.traceback())
        raise


class IntegratedOmicsGlycopeptideHypothesisBuilderTask(Task):
    def __init__(self, database_path, hypothesis_name, protein_file, site_list_file,
                 glycan_file, glycan_file_type, constant_modifications,
                 callback, **kwargs):
        args = (database_path, hypothesis_name, protein_file, site_list_file,
                glycan_file, glycan_file_type, constant_modifications)
        if hypothesis_name == "":
            hypothesis_name = get_random_string(10)
        job_name = "Integrated-Omics Glycopeptide Hypothesis Builder " + hypothesis_name
        kwargs.setdefault('name', job_name)
        super(IntegratedOmicsGlycopeptideHypothesisBuilderTask, self).__init__(taskmain, args, callback, **kwargs)
