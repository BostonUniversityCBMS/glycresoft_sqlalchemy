from glycresoft_sqlalchemy.data_model import DatabaseManager, GlycopeptideMatch, Hypothesis
from sqlalchemy.sql.functions import max as sql_max, count


def format_q_value(pairs):
    for q_value, count in pairs:
        try:
            print "%0.3f %d" % (q_value, count)
        except:
            pass
    return "\n"


def main(database_path):
    session = DatabaseManager(database_path).session()
    for hypothesis in session.query(Hypothesis):
        print hypothesis
        if hypothesis.parameters.get("is_decoy", False):
            continue
        for protein in hypothesis.proteins.values():
            if protein.theoretical_glycopeptides.count() == 0:
                continue
            print protein
            print "Top Score: ", session.query(
                sql_max(GlycopeptideMatch.ms2_score)).filter(GlycopeptideMatch.protein_id == protein.id).first()
            print "Number of matches: ", session.query(GlycopeptideMatch).filter(
                GlycopeptideMatch.protein_id == protein.id).count()
            print "q-values:"
            format_q_value(session.query(
                GlycopeptideMatch.q_value,
                count(GlycopeptideMatch.q_value)).filter(
                GlycopeptideMatch.protein_id == protein.id).group_by(
                GlycopeptideMatch.q_value).all())
            print "\n"
    session.close()

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
