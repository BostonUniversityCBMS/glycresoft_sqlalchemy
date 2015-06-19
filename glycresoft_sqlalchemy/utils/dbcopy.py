import os
from ..data_model import Experiment, GlycopeptideMatch, DatabaseManager, Protein, TheoreticalGlycopeptide


def dbcopy(source_database_path, out_database_path, experiment_id):
    if not os.path.exists(out_database_path):
        out_manager = DatabaseManager(out_database_path)
        o = out_manager.connect()
        out_manager.initialize(o)
        out_session = out_manager.session(o)
    else:
        out_session = DatabaseManager(out_database_path).session()
    in_session = DatabaseManager(source_database_path).session()
    experiment = in_session.query(Experiment).filter(Experiment.id == experiment_id).first()
    proteins = experiment.proteins.values()

    in_session.expunge_all()

    out_session.merge(experiment)
    out_session.commit()
    print experiment.id
    for protein in proteins:
        out_session.merge(protein)
        out_session.commit()
        print protein.name
        for theoretical in in_session.query(TheoreticalGlycopeptide.id).filter(
                TheoreticalGlycopeptide.protein_id == protein.id):
            theoretical = in_session.query(
                TheoreticalGlycopeptide).filter(
                TheoreticalGlycopeptide.id == theoretical[0]).first()
            in_session.expunge(theoretical)
            out_session.merge(theoretical)
            out_session.commit()
        out_session.commit()

    out_session.commit()
    print out_session.query(TheoreticalGlycopeptide).count()
    out_session.commit()
    return out_session
