import csv
from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, Protein, GlycopeptideMatch


def export_matches_as_csv(database_path, output_path=None):
    manager = DatabaseManager(database_path)
    session = manager.session()
    gpms = session.query(GlycopeptideMatch)
    if output_path is None:
        output_path = database_path[-3:] + ".glycopeptide_matches.csv"
    with open(output_path) as fh:
        writer = csv.writer(fh)
        header = [
            "id", "ms1_score", "ms2_score", "q_value", "observed_mass", "volume", "ppm_error",
            "scan_id_range", "glycopeptide_sequence", "mean_coverage", "mean_hexnac_coverage",
            "percent_bare_b_ion_coverage", "percent_bare_y_ion_coverage", "percent_glycosylated_b_ion_coverage",
            "percent_glycosylated_y_ion_coverage", "protein_name"
        ]
        writer.writerow(header)
        for gpm in gpms:
            row  = [
            gpm.id, gpm.ms1_score, gpm.ms2_score, gpm.q_value, gpm.observed_mass, gpm.volume, gpm.ppm_error,
            ' '.join(gpm.scan_id_range), gpm.glycopeptide_sequence, gpm.mean_coverage, gpm.mean_hexnac_coverage,
            len(gpm.stub_ions),
            len(gpm.bare_b_ions) / float(gpm.total_bare_b_ions_possible), 
            len(gpm.bare_y_ions) / float(gpm.total_bare_y_ions_possible),
            len(gpm.glycosylated_b_ions) / float(gpm.total_glycosylated_b_ions_possible),
            len(gpm.glycosylated_y_ions) / float(gpm.total_glycosylated_y_ions_possible),
            gpm.protein.name
            ]
            writer.writerow(row)
    return output_path


def main():
    import sys
    args = sys.argv[1:]    
    export_matches_as_csv(*args)

if __name__ == '__main__':
    main()
