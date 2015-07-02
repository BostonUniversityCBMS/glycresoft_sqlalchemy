from sqlalchemy import func
from glycresoft_sqlalchemy.scoring import target_decoy


def test_main():
    job = target_decoy.TargetDecoyAnalyzer("./datafiles/integrated_omics_simple.db", 1, 2)
    session = job.manager.session()
    session.query(target_decoy.GlycopeptideMatch).update({"q_value": None}, synchronize_session=False)
    session.commit()
    before =  session.query(
        target_decoy.GlycopeptideMatch.q_value,
        func.count(target_decoy.GlycopeptideMatch.q_value)).group_by(
        target_decoy.GlycopeptideMatch.q_value).all()

    assert before == [(None, 0)]

    job.start()
    after = session.query(
        target_decoy.GlycopeptideMatch.q_value,
        func.count(target_decoy.GlycopeptideMatch.q_value)).group_by(
        target_decoy.GlycopeptideMatch.q_value).all()

    assert after != [(None, 0)]

if __name__ == '__main__':
    test_main()
