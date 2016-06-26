import os
import gc
import csv
import argparse
from glycresoft_sqlalchemy.data_model import (DatabaseManager, Hypothesis, Protein, TheoreticalGlycanComposition,
                                              GlycopeptideMatch, PipelineModule, HypothesisSampleMatch,
                                              TheoreticalGlycopeptideComposition)

from glypy.composition.glycan_composition import FrozenGlycanComposition


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

        def adduct_label(x):
            if x.mass_shift_type:
                return "%s:%d" % (x.mass_shift.name, x.mass_shift_count)
            return 'No Shift'

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
                pgm.modification_state_count, pgm.charge_state_count, pgm.scan_count, pgm.scan_density,
                pgm.a_peak_intensity_error, pgm.average_a_to_a_plus_2_ratio,
                pgm.total_volume, pgm.average_signal_to_noise, pgm.centroid_scan_error,
                pgm.centroid_scan_error, pgm.first_scan_id, pgm.last_scan_id,
                theoretical_match.calculated_mass if pgm.matched else ""
            ] + glycan_composition + [",".join(adduct_label(g) for g in pgm.subgroups), pgm.id]
            writer.writerow(row)
        return output_path


def export_glycan_ms1_matches_legacy_ungrouping(peak_group_matches, monosaccharide_identities, output_path):
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        headers = [
            "Score", "MassSpec MW", "Compound Key", "PPM Error",
            "#ofAdduct", "#ofCharges", "#ofScans", "ScanDensity", "Avg A:A+2 Error",
            "A:A+2 Ratio", "Total Volume", "Signal to Noise Ratio", "Centroid Scan Error",
            "Centroid Scan", "MinScanNumber", "MaxScanNumber", "Hypothesis MW"
        ] + monosaccharide_identities + ["Adduct/Replacement", "ID"]
        writer.writerow(headers)

        def adduct_label(x):
            if x.mass_shift_type:
                return "%s:%d" % (x.mass_shift.name, x.mass_shift_count)
            return 'No Shift'

        for pgm in peak_group_matches.yield_per(1000):
            theoretical_match = pgm.theoretical_match
            if theoretical_match is not None:
                glycan_composition = [theoretical_match.glycan_composition[g] for g in monosaccharide_identities]
            else:
                glycan_composition = ['' for g in monosaccharide_identities]
            for subgroup in pgm.subgroups:
                row = [
                    pgm.ms1_score, subgroup.weighted_monoisotopic_mass,
                    theoretical_match.composition if pgm.matched else "",
                    pgm.ppm_error if subgroup.matched else "",
                    pgm.modification_state_count, subgroup.charge_state_count,
                    subgroup.scan_count, subgroup.scan_density,
                    subgroup.a_peak_intensity_error, subgroup.average_a_to_a_plus_2_ratio,
                    subgroup.total_volume, subgroup.average_signal_to_noise, subgroup.centroid_scan_error,
                    subgroup.centroid_scan_error, subgroup.first_scan_id, subgroup.last_scan_id,
                    theoretical_match.calculated_mass if pgm.matched else ""
                ] + glycan_composition + [adduct_label(subgroup), subgroup.id]
                writer.writerow(row)
        return output_path


def export_glycan_composition_hypothesis(glycan_compositions, monosaccharide_identities, output_path):
    headers = ["theoretical mass", "composition"] + monosaccharide_identities
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        for row in glycan_compositions:
            writer.writerow([str(row.calculated_mass), str(row.composition)] + [
                str(row.glycan_composition[k]) for k in monosaccharide_identities])
    return output_path


def export_theoretical_glycopeptide_hypothesis(glycopeptides, monosaccharide_identities, output_path, session):
    headers = ["theoretical mass", "peptide sequence", "modifications", "glycan composition"] +\
        monosaccharide_identities + ["protein name", "start", "end", "missed cleavages", "glycosylation sites"]
    with open(output_path, 'wb') as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        i = 0
        n = glycopeptides.count()
        while i < n:
            for row in (glycopeptides.slice(i, i + 15000)):
                glycan_comp = FrozenGlycanComposition.parse(row.glycan_composition_str)
                writer.writerow([
                    str(row.calculated_mass), row.modified_peptide_sequence, str(row.peptide_modifications),
                    row.glycan_composition_str] + [str(glycan_comp[k]) for k in monosaccharide_identities] + [
                    str(row.protein.name), str(row.start_position), str(row.end_position),
                    str(row.count_missed_cleavages), str(row.count_glycosylation_sites)])
                i += 1
            session.expunge_all()
            print i, '/', n
    return output_path


class CSVExportDriver(PipelineModule):
    def __init__(self, database_path, hypothesis_ids=None, hypothesis_sample_match_ids=None,
                 output_path=None, filterfunc=lambda q: q):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()

        if hypothesis_sample_match_ids is None:
            self.hypothesis_sample_match_ids = []
        else:
            try:
                iter(hypothesis_sample_match_ids)
                self.hypothesis_sample_match_ids = hypothesis_sample_match_ids
            except:
                self.hypothesis_sample_match_ids = [hypothesis_sample_match_ids]

        if hypothesis_ids is None:
            self.hypothesis_ids = []
        else:
            try:
                iter(hypothesis_ids)
                self.hypothesis_ids = hypothesis_ids
            except:
                self.hypothesis_ids = [hypothesis_ids]

        if output_path is None:
            output_path = os.path.splitext(database_path)[0]
        self.output_path = output_path
        self.filterfunc = filterfunc

    def dispatch_export_hypothesis_sample_match(self, hypothesis_sample_match_id):
        session = self.session
        hsm = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

        filterfunc = self.filterfunc
        outputs = []

        def getname(hsm):
            try:
                return os.path.splitext(hsm.name)[0]
            except:
                return str(hsm.target_hypothesis.name) + "_on_" + str(hsm.sample_run_name)
        for res_type, query in hsm.results():
            if res_type == GlycopeptideMatch:
                output_path = self.output_path + '.{}.glycopeptide_matches.csv'.format(getname(hsm))
                outputs.append(output_path)
                # Only export target hypothesis
                export_glycopeptide_ms2_matches(filterfunc(query.filter(
                    GlycopeptideMatch.protein_id == Protein.id,
                    Protein.hypothesis_id == hsm.target_hypothesis_id)), output_path)
            elif (res_type == TheoreticalGlycopeptideComposition):
                output_path = self.output_path + '.{}.glycopeptide_compositions.csv'.format(getname(hsm))
                outputs.append(output_path)
                export_glycopeptide_ms1_matches_legacy(
                    filterfunc(query),
                    hsm.target_hypothesis.parameters['monosaccharide_identities'],
                    output_path)
            elif (res_type == TheoreticalGlycanComposition):
                output_path = self.output_path + ".{}.glycan_compositions.csv".format(getname(hsm))
                outputs.append(output_path)
                base_types = session.query(
                    TheoreticalGlycanComposition.GlycanCompositionAssociation.base_type.distinct()).join(
                    TheoreticalGlycanComposition).filter(
                    TheoreticalGlycanComposition.hypothesis_id == hsm.target_hypothesis.id)
                monosaccharide_identities = [b for q in base_types for b in q]
                export_glycan_ms1_matches_legacy(
                    filterfunc(query),
                    monosaccharide_identities,
                    output_path
                    )
            else:
                pass
        return outputs

    def dispatch_export_hypothesis(self, hypothesis_id):
        session = self.session
        hypothesis = session.query(Hypothesis).get(hypothesis_id)
        outputs = []

        def getname(hypothesis):
            try:
                return os.path.splitext(hypothesis.name)[0]
            except:
                return str(hypothesis.name)

        if hypothesis.theoretical_structure_type == TheoreticalGlycanComposition:
            output_path = self.output_path + ".{}.glycan_compositions.csv".format(getname(hypothesis))
            outputs.append(output_path)
            base_types = session.query(
                TheoreticalGlycanComposition.GlycanCompositionAssociation.base_type.distinct()).join(
                TheoreticalGlycanComposition).filter(
                TheoreticalGlycanComposition.hypothesis_id == hypothesis.id)
            monosaccharide_identities = [b for q in base_types for b in q]
            export_glycan_composition_hypothesis(
                session.query(TheoreticalGlycanComposition).filter(
                    TheoreticalGlycanComposition.hypothesis_id == hypothesis.id).yield_per(1000),
                monosaccharide_identities,
                output_path)
        elif hypothesis.theoretical_structure_type == TheoreticalGlycopeptideComposition:
            output_path = self.output_path + ".{}.glycopeptide_compositions.csv".format(getname(hypothesis))
            outputs.append(output_path)
            monosaccharide_identities = hypothesis.parameters['monosaccharide_identities']
            export_theoretical_glycopeptide_hypothesis(
                session.query(TheoreticalGlycopeptideComposition).join(Protein).filter(
                    Protein.hypothesis_id == hypothesis.id).yield_per(1000),
                monosaccharide_identities,
                output_path, session=session)
        else:
            print("Unsupported Hypothesis Type: %s" % hypothesis.theoretical_structure_type.__name__)
        return outputs

    def run(self):
        hsm_exports = map(self.dispatch_export_hypothesis_sample_match, self.hypothesis_sample_match_ids)
        hypothesis_exports = map(self.dispatch_export_hypothesis, self.hypothesis_ids)
        return hsm_exports, hypothesis_exports

app = argparse.ArgumentParser("export-csv")

app.add_argument("database_path", help="path to the database file to export")
app.add_argument("-d", "--hypothesis-id", action="append", type=int, help="The hypothesis to export.")
app.add_argument("-e", "--hypothesis-sample-match-id", action="append", type=int,
                 help="The hypothesis sample match to export.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")


def main(database_path, hypothesis_ids, hypothesis_sample_match_ids, out):
    return CSVExportDriver(
        database_path, hypothesis_ids=hypothesis_ids,
        hypothesis_sample_match_ids=hypothesis_sample_match_ids, output_path=out).start()


def taskmain():
    args = app.parse_args()
    main(args.database_path, args.hypothesis_id, args.hypothesis_sample_match_id, args.out)


if __name__ == '__main__':
    taskmain()
