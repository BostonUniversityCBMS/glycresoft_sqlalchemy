import os
import csv
import argparse
from glycresoft_sqlalchemy.data_model import (DatabaseManager, Hypothesis, Protein, TheoreticalGlycanComposition,
                                              GlycopeptideMatch, PipelineModule, HypothesisSampleMatch,
                                              TheoreticalGlycopeptideComposition)


def export_glycopeptide_ms2_matches(glycopeptides, output_path):
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        header = [
            "id", "ms1_score", "ms2_score", "q_value", "observed_mass", "volume", "ppm_error",
            "scan_id_range", "glycopeptide_sequence", "sequence_length", "mean_coverage", "mean_hexnac_coverage",
            "stub_ion_count", "bare_b_ion_coverage", "bare_y_ion_coverage", "glycosylated_b_ion_coverage",
            "glycosylated_y_ion_coverage", "protein_name", "b_ions", "y_ions", 'stub_ions', 'oxonium_ions',
            "start_position", "end_position"
        ]
        writer.writerow(header)
        for gpm in glycopeptides.yield_per(1000):
            row = [
                gpm.id, gpm.ms1_score, gpm.ms2_score, gpm.q_value, gpm.observed_mass, gpm.volume, gpm.ppm_error,
                ';'.join(map(str, gpm.scan_id_range)), gpm.glycopeptide_sequence, gpm.sequence_length,
                gpm.mean_coverage, gpm.mean_hexnac_coverage,
                len(gpm.stub_ions),
                len(gpm.bare_b_ions),
                len(gpm.bare_y_ions),
                len(gpm.glycosylated_b_ions),
                len(gpm.glycosylated_y_ions),
                gpm.protein.name,
                ";".join([ion['key'] for series in [gpm.bare_b_ions, gpm.glycosylated_b_ions]
                          for ion in series]),
                ";".join([ion['key'] for series in [gpm.bare_y_ions, gpm.glycosylated_y_ions]
                          for ion in series]),
                ";".join([ion['key'] for ion in gpm.stub_ions]),
                ";".join([ion['key'] for ion in gpm.oxonium_ions]),
                gpm.start_position, gpm.end_position

            ]
            writer.writerow(row)
    return output_path


def export_glycopeptide_ms1_matches_legacy(peak_group_matches, monosaccharide_identities, output_path):
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        headers = [
            "Score", "MassSpec MW", "Compound Key", "PeptideSequence", "PPM Error",
            "#ofAdduct", "#ofCharges", "#ofScans", "ScanDensity", "Avg A:A+2 Error",
            "A:A+2 Ratio", "Total Volume", "Signal to Noise Ratio", "Centroid Scan Error",
            "Centroid Scan", "MinScanNumber", "MaxScanNumber", "Hypothesis MW",
        ] + monosaccharide_identities + [
            "Adduct/Replacement", "PeptideModification", "PeptideMissedCleavage#", "#ofGlycanAttachmentToPeptide",
            "StartAA", "EndAA", "ProteinID"
        ]
        writer.writerow(headers)
        for pgm in peak_group_matches.yield_per(1000):
            theoretical_match = pgm.theoretical_match
            if theoretical_match is not None:
                glycan_composition = theoretical_match.glycan_composition_str[1:-1].split(';')
            else:
                glycan_composition = ['' for g in monosaccharide_identities]
            row = [
                pgm.ms1_score, pgm.weighted_monoisotopic_mass,
                theoretical_match.glycan_composition_str if pgm.matched else "",
                theoretical_match.modified_peptide_sequence if pgm.matched else "",
                pgm.ppm_error if pgm.matched else "",
                "", pgm.charge_state_count, pgm.scan_count, pgm.scan_density,
                pgm.a_peak_intensity_error, pgm.average_a_to_a_plus_2_ratio,
                pgm.total_volume, pgm.average_signal_to_noise, pgm.centroid_scan_error,
                pgm.centroid_scan_error, pgm.first_scan_id, pgm.last_scan_id,
                theoretical_match.calculated_mass if pgm.matched else ""
            ] + glycan_composition + [
                "/", theoretical_match.peptide_modifications if pgm.matched else "",
                theoretical_match.count_missed_cleavages if pgm.matched else "",
                theoretical_match.count_glycosylation_sites if pgm.matched else "",
                theoretical_match.start_position if pgm.matched else "",
                theoretical_match.end_position if pgm.matched else "",
                theoretical_match.protein.name if pgm.matched else "",
            ]
            writer.writerow(row)
        return output_path


def export_glycan_ms1_matches_legacy(peak_group_matches, monosaccharide_identities, output_path):
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        headers = [
            "Score", "MassSpec MW", "Compound Key", "PPM Error",
            "#ofAdduct", "#ofCharges", "#ofScans", "ScanDensity", "Avg A:A+2 Error",
            "A:A+2 Ratio", "Total Volume", "Signal to Noise Ratio", "Centroid Scan Error",
            "Centroid Scan", "MinScanNumber", "MaxScanNumber", "Hypothesis MW"
        ] + monosaccharide_identities + ["Adduct/Replacement", "ID"]

        writer.writerow(headers)

        for pgm in peak_group_matches.yield_per(1000):
            theoretical_match = pgm.theoretical_match
            if theoretical_match is not None:
                glycan_composition = [theoretical_match.glycan_composition[g] for g in monosaccharide_identities]
            else:
                glycan_composition = ['' for g in monosaccharide_identities]
            row = [
                pgm.ms1_score, pgm.weighted_monoisotopic_mass,
                theoretical_match.composition if pgm.matched else "",
                pgm.ppm_error if pgm.matched else "",
                pgm.mass_shift_count if (pgm.mass_shift is not None and pgm.mass_shift.mass != 0) else 0, pgm.charge_state_count, pgm.scan_count, pgm.scan_density,
                pgm.a_peak_intensity_error, pgm.average_a_to_a_plus_2_ratio,
                pgm.total_volume, pgm.average_signal_to_noise, pgm.centroid_scan_error,
                pgm.centroid_scan_error, pgm.first_scan_id, pgm.last_scan_id,
                theoretical_match.calculated_mass if pgm.matched else ""
            ] + glycan_composition + ['/' if pgm.mass_shift_type == None else "%s/%s" % (pgm.mass_shift.name, pgm.mass_shift.mass), pgm.id]
            writer.writerow(row)
        return output_path


class CSVExportDriver(PipelineModule):
    def __init__(self, database_path, hypothesis_sample_match_ids=None, output_path=None, filterfunc=lambda q: q):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()

        if hypothesis_sample_match_ids is None or len(hypothesis_sample_match_ids) == 0:
            hypothesis_sample_match_ids = [hsm.id for hsm in self.session.query(
                HypothesisSampleMatch.id)]
            self.hypothesis_sample_match_ids = hypothesis_sample_match_ids
        else:
            try:
                iter(hypothesis_sample_match_ids)
                self.hypothesis_sample_match_ids = hypothesis_sample_match_ids
            except:
                self.hypothesis_sample_match_ids = [hypothesis_sample_match_ids]

        if output_path is None:
            output_path = os.path.splitext(database_path)[0]
        self.output_path = output_path
        self.filterfunc = filterfunc

    def dispatch_export(self, hypothesis_sample_match_id):
        session = self.session
        hsm = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

        filterfunc = self.filterfunc

        def getname(hsm):
            try:
                return os.path.splitext(hsm.name)[0]
            except:
                return str(hsm.target_hypothesis.name) + "_on_" + str(hsm.sample_run_name)
        for res_type, query in hsm.results():
            print res_type
            if res_type == GlycopeptideMatch:
                output_path = self.output_path + '.{}.glycopeptide_matches.csv'.format(getname(hsm))
                # Only export target hypothesis
                export_glycopeptide_ms2_matches(filterfunc(query.filter(
                    GlycopeptideMatch.protein_id == Protein.id,
                    Protein.hypothesis_id == hsm.target_hypothesis_id)), output_path)
            elif (res_type == TheoreticalGlycopeptideComposition):
                output_path = self.output_path + '.{}.glycopeptide_compositions.csv'.format(getname(hsm))
                export_glycopeptide_ms1_matches_legacy(
                    filterfunc(query),
                    hsm.target_hypothesis.parameters['monosaccharide_identities'],
                    output_path)
            elif (res_type == TheoreticalGlycanComposition):
                output_path = self.output_path + ".{}.glycan_compositions.csv".format(getname(hsm))
                base_types = session.query(TheoreticalGlycanComposition.GlycanCompositionAssociation.base_type.distinct()).join(TheoreticalGlycanComposition).filter(
                    TheoreticalGlycanComposition.hypothesis_id == hsm.target_hypothesis.id)
                monosaccharide_identities = [b for q in base_types for b in q]
                export_glycan_ms1_matches_legacy(
                    filterfunc(query),
                    monosaccharide_identities,
                    output_path
                    )


            else:
                pass

    def run(self):
        return map(self.dispatch_export, self.hypothesis_sample_match_ids)

app = argparse.ArgumentParser("export-csv")

app.add_argument("database_path", help="path to the database file to export")
app.add_argument("-e", "--hypothesis-sample-match-id",  action="append", type=int, help="The hypothesis sample match to export.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")


def main(database_path, hypothesis_id, out):
    return CSVExportDriver(database_path, hypothesis_id, out).start()


def taskmain():
    args = app.parse_args()
    main(args.database_path, args.hypothesis_sample_match_id, args.out)


if __name__ == '__main__':
    taskmain()
