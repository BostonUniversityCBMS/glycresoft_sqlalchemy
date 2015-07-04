import os
import csv
import argparse
from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, Protein, GlycopeptideMatch


def export_matches_as_csv(database_path, hypothesis_id=1, output_path=None, include_ions_matched=True):
    manager = DatabaseManager(database_path)
    session = manager.session()
    gpms = session.query(GlycopeptideMatch).filter(
        GlycopeptideMatch.protein_id == Protein.id, Protein.hypothesis_id == hypothesis_id)
    if output_path is None:
        hypothesis = session.query(Hypothesis).get(hypothesis_id)
        output_path = database_path[:-3] + "_" + hypothesis.name + ".glycopeptide_matches.csv"
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        header = [
            "id", "ms1_score", "ms2_score", "q_value", "observed_mass", "volume", "ppm_error",
            "scan_id_range", "glycopeptide_sequence", "mean_coverage", "mean_hexnac_coverage",
            "stub_ions", "bare_b_ion_coverage", "bare_y_ion_coverage", "glycosylated_b_ion_coverage",
            "glycosylated_y_ion_coverage", "protein_name"
        ]
        if include_ions_matched:
            header.extend(["b_ions", "y_ions", 'stub_ions', 'oxonium_ions'])
        writer.writerow(header)
        for gpm in gpms:
            row = [
                gpm.id, gpm.ms1_score, gpm.ms2_score, gpm.q_value, gpm.observed_mass, gpm.volume, gpm.ppm_error,
                ' '.join(map(str, gpm.scan_id_range)), gpm.glycopeptide_sequence,
                gpm.mean_coverage, gpm.mean_hexnac_coverage,
                len(gpm.stub_ions),
                len(gpm.bare_b_ions),
                len(gpm.bare_y_ions),
                len(gpm.glycosylated_b_ions),
                len(gpm.glycosylated_y_ions),
                gpm.protein.name
            ]
            if include_ions_matched:
                row.extend([
                        ";".join([ion['key'] for series in [gpm.bare_b_ions, gpm.glycosylated_b_ions]
                                  for ion in series]),
                        ";".join([ion['key'] for series in [gpm.bare_y_ions, gpm.glycosylated_y_ions]
                                  for ion in series]),
                        ";".join([ion['key'] for ion in gpm.stub_ions]),
                        ";".join([ion['key'] for ion in gpm.oxonium_ions])
                    ])
            writer.writerow(row)
    return output_path


app = argparse.ArgumentParser("export-csv")

app.add_argument("database_path", help="path to the database file to analyze")
app.add_argument("-e", "--hypothesis-id", default=None, help="The hypothesis to analyze.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")


def main(database_path, hypothesis_id, out, include_ions_matched):
    return export_matches_as_csv(database_path, hypothesis_id, out, include_ions_matched)


def taskmain():
    args = app.parse_args()
    dbm = DatabaseManager(args.database_path)
    session = dbm.session()
    if args.hypothesis_id is None:
        hypothesis_ids = [x[0] for x in session.query(Hypothesis.id).all()]
    else:
        hypothesis_ids = [args.hypothesis_id]
    for hypothesis_id in hypothesis_ids:
        if args.out is not None and len(hypothesis_ids) > 1:
            parts = os.path.splitext(args.out)
            name = session.query(Hypothesis.name).get(hypothesis_id)[0]
            out = parts[0] + "_" + name + parts[1]
        else:
            out = args.out
        main(args.database_path, hypothesis_id, out)


if __name__ == '__main__':
    taskmain()
