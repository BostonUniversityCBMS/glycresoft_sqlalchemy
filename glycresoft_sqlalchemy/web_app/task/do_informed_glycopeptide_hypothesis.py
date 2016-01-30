from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder import integrated_omics
from .task_process import NullPipe, Message, Task
from ..common import get_random_string


def taskmain(database_path, hypothesis_name, protein_ids, protein_file, site_list_file,
             glycan_options, maximum_glycosylation_sites,
             comm=None, **kwargs):
    if comm is None:
        comm = NullPipe()
    manager = DatabaseManager(database_path)
    try:
        kwargs.setdefault("n_processes", 4)
        task = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
            database_path, hypothesis_name=hypothesis_name, mzid_path=protein_file, protein_ids=protein_ids,
            n_processes=kwargs["n_processes"], maximum_glycosylation_sites=maximum_glycosylation_sites,
            **glycan_options)
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
    def __init__(self, database_path, hypothesis_name, protein_ids, protein_file, site_list_file,
                 glycan_options, maximum_glycosylation_sites,
                 callback, **kwargs):
        args = (database_path, hypothesis_name, protein_ids, protein_file, site_list_file,
                glycan_options, maximum_glycosylation_sites)
        if hypothesis_name == "":
            hypothesis_name = get_random_string(10)
        job_name = "Integrated-Omics Glycopeptide Hypothesis Builder " + hypothesis_name
        kwargs.setdefault('name', job_name)
        super(IntegratedOmicsGlycopeptideHypothesisBuilderTask, self).__init__(taskmain, args, callback, **kwargs)
