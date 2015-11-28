from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics
from .task_process import NullPipe, Message, Task


def taskmain(database_path, rules_table, constraints, hypothesis_name,
             reduction_type, derivatization_type, comm=None, **kwargs):
    if comm is None:
        comm = NullPipe()
    manager = DatabaseManager(database_path)
    manager.initialize()
    try:
        hypothesis_id = constrained_combinatorics.ConstrainedCombinatoricsGlycanHypothesisBuilder(
            database_path, rules_table=rules_table, constraints_list=constraints,
            reduction=reduction_type, derivatization=derivatization_type).start()
        session = manager()
        hypothesis = session.query(Hypothesis).get(hypothesis_id)
        hypothesis.name = hypothesis_name
        session.add(hypothesis)
        session.commit()
        comm.send(Message(hypothesis.to_json(), "new-hypothesis"))
        return hypothesis_id
    except:
        comm.send(Message.traceback())
        raise


class CombinatorialGlycanHypothesisTask(Task):
    def __init__(self, database_path, rules_table, constraints, hypothesis_name,
                 reduction_type, derivatization_type, callback, **kwargs):
        args = (database_path, rules_table, constraints, hypothesis_name, reduction_type, derivatization_type)
        job_name = "Combinatorial Glycan Hypothesis Builder " + hypothesis_name
        kwargs.setdefault('name', job_name)
        super(CombinatorialGlycanHypothesisTask, self).__init__(taskmain, args, callback, **kwargs)
