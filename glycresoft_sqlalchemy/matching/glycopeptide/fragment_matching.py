import math
import multiprocessing
import functools
import itertools
import logging
import operator
try:
    logger = logging.getLogger("fragment_matching")
except:
    pass

from collections import defaultdict

from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser

from glycresoft_sqlalchemy import data_model as model
from glycresoft_sqlalchemy.data_model import PipelineModule, MSMSSqlDB, TandemScan
from glycresoft_sqlalchemy.utils.common_math import ppm_error


HypothesisSampleMatch = model.HypothesisSampleMatch
Hypothesis = model.Hypothesis
SampleRun = model.SampleRun
Protein = model.Protein
TheoreticalGlycopeptide = model.TheoreticalGlycopeptide
GlycopeptideSpectrumMatch = model.GlycopeptideSpectrumMatch
GlycopeptideMatch = model.GlycopeptideMatch


neutral_mass_getter = operator.attrgetter("neutral_mass")
key_getter = operator.itemgetter('key')

decon_format_lookup = {
    "bupid_yaml": BUPIDMSMSYamlParser,
    "db": MSMSSqlDB,
}

PROTON = 1.007276035

ms1_tolerance_default = 1e-5
ms2_tolerance_default = 2e-5


def yield_ids(session, hypothesis_id, chunk_size=1000, filter=lambda q: q):
    base_query = filter(session.query(TheoreticalGlycopeptide.id).filter(
        TheoreticalGlycopeptide.protein_id == Protein.id,
        Protein.hypothesis_id == hypothesis_id))
    chunk = []

    for item in base_query:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    yield chunk


def yield_scan_ids(session, sample_run_id, chunk_size=100, filter=lambda q: q):
    base_query = filter(session.query(TandemScan.id).filter(TandemScan.sample_run_id == sample_run_id))
    chunk = []

    for item in base_query:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    yield chunk


def batch_match_fragments(theoretical_ids, msmsdb_path, ms1_tolerance, ms2_tolerance,
                          database_manager, hypothesis_sample_match_id, sample_run_id,
                          hypothesis_id, intensity_threshold=0.0):
    try:
        session = database_manager.session()
        msmsdb = MSMSSqlDB(msmsdb_path)
        # Localized global references
        proton = PROTON
        lppm_error = ppm_error
        lfabs = math.fabs

        glycopeptide_matches_spectrum_matches = []

        for theoretical_id in theoretical_ids:
            theoretical = session.query(TheoreticalGlycopeptide).get(theoretical_id)
            if theoretical is None:
                raise ValueError("No theoretical glycopeptide with id %r" % theoretical_id)

            # Containers for global theoretical peak matches
            oxonium_ions = []
            bare_b_ions = []
            bare_y_ions = []
            glycosylated_b_ions = []
            glycosylated_y_ions = []
            stub_ions = []

            spectrum_matches = []

            query = msmsdb.ppm_match_tolerance_search(theoretical.calculated_mass, ms1_tolerance, sample_run_id)

            for spectrum in query:
                peak_list = spectrum.tandem_data
                peak_list = [p for p in peak_list if p.intensity >= intensity_threshold]
                peak_list = sorted(peak_list, key=neutral_mass_getter)
                peak_match_map = defaultdict(list)

                c = 0
                collect = oxonium_ions.append
                for theoretical_ion in theoretical.oxonium_ions:
                    query_mass = theoretical_ion['mass']
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        protonated_mass = observed_mass  # + proton
                        match_error = lppm_error(protonated_mass, query_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": protonated_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                            c += 1
                        elif protonated_mass > query_mass + 10:
                            break

                # If no oxonium ions were found, skip this spectrum
                if c < 1:
                    continue

                collect = bare_b_ions.append
                for theoretical_ion in theoretical.bare_b_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "intensity": peak.intensity,
                                      "observed_mass": observed_mass,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = bare_y_ions.append
                for theoretical_ion in theoretical.bare_y_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = glycosylated_b_ions.append
                for theoretical_ion in theoretical.glycosylated_b_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = glycosylated_y_ions.append
                for theoretical_ion in theoretical.glycosylated_y_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = stub_ions.append
                for theoretical_ion in theoretical.stub_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                spectrum_matches.append((spectrum, peak_match_map))

            if len(spectrum_matches) > 0:
                scan_ids = []
                # session = database_manager.session()
                for spectrum, peak_match_map in spectrum_matches:
                    scan_ids.append(spectrum.time)

                first_scan = min(scan_ids)
                last_scan = max(scan_ids)

                gpm = GlycopeptideMatch(
                    protein_id=theoretical.protein_id,
                    theoretical_glycopeptide_id=theoretical.id,
                    ms1_score=theoretical.ms1_score,
                    observed_mass=theoretical.observed_mass,
                    calculated_mass=theoretical.calculated_mass,
                    volume=theoretical.volume,
                    ppm_error=theoretical.ppm_error,

                    glycan_composition_str=theoretical.glycan_composition_str,
                    glycan_combination_id=theoretical.glycan_combination_id,

                    base_peptide_sequence=theoretical.base_peptide_sequence,
                    modified_peptide_sequence=theoretical.modified_peptide_sequence,
                    glycopeptide_sequence=theoretical.glycopeptide_sequence,
                    sequence_length=theoretical.sequence_length,
                    peptide_modifications=theoretical.peptide_modifications,
                    count_glycosylation_sites=theoretical.count_glycosylation_sites,
                    count_missed_cleavages=theoretical.count_missed_cleavages,
                    glycosylation_sites=theoretical.glycosylation_sites,
                    start_position=theoretical.start_position,
                    end_position=theoretical.end_position,
                    oxonium_ions=oxonium_ions,
                    stub_ions=stub_ions,
                    bare_b_ions=bare_b_ions,
                    bare_y_ions=bare_y_ions,
                    glycosylated_b_ions=glycosylated_b_ions,
                    glycosylated_y_ions=glycosylated_y_ions,
                    scan_id_range=scan_ids,
                    first_scan=first_scan,
                    last_scan=last_scan,
                    hypothesis_sample_match_id=hypothesis_sample_match_id
                )

                spectrum_matches_objects = []

                for spectrum, peak_match_map in spectrum_matches:
                    spectrum_match_inst = GlycopeptideSpectrumMatch(
                        scan_time=spectrum.time, peak_match_map=peak_match_map,
                        precursor_charge_state=spectrum.precursor_charge_state,
                        peaks_explained=len(peak_match_map),
                        peaks_unexplained=len(spectrum.tandem_data) - len(peak_match_map),
                        hypothesis_sample_match_id=hypothesis_sample_match_id,
                        hypothesis_id=hypothesis_id)
                    spectrum_matches_objects.append(spectrum_match_inst)
                glycopeptide_matches_spectrum_matches.append((gpm, spectrum_matches_objects))
        if len(glycopeptide_matches_spectrum_matches) > 0:
            session.add_all(g for g, spectrum_match_list in glycopeptide_matches_spectrum_matches)
            session.flush()
            for g, spectrum_match_list in glycopeptide_matches_spectrum_matches:
                for spectrum_match_inst in spectrum_match_list:
                    spectrum_match_inst.glycopeptide_match_id = g.id
                session.add_all(spectrum_match_list)
            session.commit()
        return len(theoretical_ids)
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise e
    finally:
        try:
            session.close()
        except:
            pass


def batch_match_theoretical_ions(scan_ids, msmsdb_path, ms1_tolerance, ms2_tolerance,
                                 database_manager, hypothesis_sample_match_id, sample_run_id,
                                 hypothesis_id, intensity_threshold=0.0):
    try:
        session = database_manager()
        msmsdb = MSMSSqlDB(msmsdb_path)()
        # Localized global references
        proton = PROTON
        lppm_error = ppm_error
        lfabs = math.fabs

        glycopeptide_matches_spectrum_matches = []

        for scan_id in scan_ids:
            scan_id = scan_id[0]
            spectrum = msmsdb.query(TandemScan).get(scan_id)

            # Containers for global theoretical peak matches
            oxonium_ions = []
            bare_b_ions = []
            bare_y_ions = []
            glycosylated_b_ions = []
            glycosylated_y_ions = []
            stub_ions = []

            spectrum_matches = []

            query = TheoreticalGlycopeptide.ppm_error_tolerance_search(
                session, spectrum.precursor_neutral_mass,
                ms1_tolerance, hypothesis_id)

            peak_list = spectrum.tandem_data
            peak_list = [p for p in peak_list if p.intensity >= intensity_threshold]
            peak_list = sorted(peak_list, key=neutral_mass_getter)

            for theoretical in query:
                peak_match_map = defaultdict(list)

                c = 0
                collect = oxonium_ions.append
                for theoretical_ion in theoretical.oxonium_ions:
                    query_mass = theoretical_ion['mass']
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        protonated_mass = observed_mass  # + proton
                        match_error = lppm_error(protonated_mass, query_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": protonated_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                            c += 1
                        elif protonated_mass > query_mass + 10:
                            break

                # If no oxonium ions were found, skip this spectrum
                if c < 1:
                    continue

                collect = bare_b_ions.append
                for theoretical_ion in theoretical.bare_b_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "intensity": peak.intensity,
                                      "observed_mass": observed_mass,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = bare_y_ions.append
                for theoretical_ion in theoretical.bare_y_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = glycosylated_b_ions.append
                for theoretical_ion in theoretical.glycosylated_b_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = glycosylated_y_ions.append
                for theoretical_ion in theoretical.glycosylated_y_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                collect = stub_ions.append
                for theoretical_ion in theoretical.stub_ions:
                    query_mass = theoretical_ion['mass']
                    deprotonated_mass = query_mass  # - proton
                    for peak in peak_list:
                        observed_mass = peak.neutral_mass
                        match_error = lppm_error(observed_mass, deprotonated_mass)
                        if lfabs(match_error) <= ms2_tolerance:
                            match = ({'key': theoretical_ion['key'],
                                      "observed_mass": observed_mass,
                                      "intensity": peak.intensity,
                                      'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                            collect(match)
                            peak_match_map[peak.scan_peak_index].append(match)
                        elif observed_mass > query_mass + 10:
                            break

                spectrum_matches.append((theoretical, peak_match_map))

            if len(spectrum_matches) > 0:
                for theoretical, peak_match_map in spectrum_matches:
                    spectrum_match_inst = GlycopeptideSpectrumMatch(
                        scan_time=spectrum.time, peak_match_map=peak_match_map,
                        precursor_charge_state=spectrum.precursor_charge_state,
                        peaks_explained=len(peak_match_map),
                        peaks_unexplained=len(spectrum.tandem_data) - len(peak_match_map),
                        hypothesis_sample_match_id=hypothesis_sample_match_id,
                        theoretical_glycopeptide_id=theoretical.id,
                        hypothesis_id=hypothesis_id)
                    glycopeptide_matches_spectrum_matches.append(spectrum_match_inst)

        if len(glycopeptide_matches_spectrum_matches) > 0:
            session.bulk_save_objects(glycopeptide_matches_spectrum_matches)
            session.commit()
        return len(scan_ids)
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise e
    finally:
        try:
            session.close()
        except:
            pass


def bundle_spectrum_matches(theoretical_glycopeptide_ids, database_manager, hypothesis_sample_match_id):
    session = database_manager()
    matches = []
    for theoretical_id in theoretical_glycopeptide_ids:
        theoretical_id = theoretical_id[0]

        spectrum_matches = session.query(GlycopeptideSpectrumMatch).filter(
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.theoretical_glycopeptide_id == theoretical_id).order_by(
            GlycopeptideSpectrumMatch.scan_time).all()

        if len(spectrum_matches) == 0:
            continue
        scan_ids = [scan.scan_time for scan in spectrum_matches]
        first_scan = scan_ids[0]
        last_scan = scan_ids[-1]

        theoretical = session.query(TheoreticalGlycopeptide).get(theoretical_id)

        gpm = GlycopeptideMatch(
                protein_id=theoretical.protein_id,
                theoretical_glycopeptide_id=theoretical.id,
                ms1_score=theoretical.ms1_score,
                observed_mass=theoretical.observed_mass,
                calculated_mass=theoretical.calculated_mass,
                volume=theoretical.volume,
                ppm_error=theoretical.ppm_error,

                glycan_composition_str=theoretical.glycan_composition_str,
                glycan_combination_id=theoretical.glycan_combination_id,

                base_peptide_sequence=theoretical.base_peptide_sequence,
                modified_peptide_sequence=theoretical.modified_peptide_sequence,
                glycopeptide_sequence=theoretical.glycopeptide_sequence,
                sequence_length=theoretical.sequence_length,
                peptide_modifications=theoretical.peptide_modifications,
                count_glycosylation_sites=theoretical.count_glycosylation_sites,
                count_missed_cleavages=theoretical.count_missed_cleavages,
                glycosylation_sites=theoretical.glycosylation_sites,
                start_position=theoretical.start_position,
                end_position=theoretical.end_position,
                # oxonium_ions=oxonium_ions,
                # stub_ions=stub_ions,
                # bare_b_ions=bare_b_ions,
                # bare_y_ions=bare_y_ions,
                # glycosylated_b_ions=glycosylated_b_ions,
                # glycosylated_y_ions=glycosylated_y_ions,
                scan_id_range=scan_ids,
                first_scan=first_scan,
                last_scan=last_scan,
                hypothesis_sample_match_id=hypothesis_sample_match_id)
        matches.append(gpm)
    session.bulk_save_objects(gpm)
    session.commit()
    return len(theoretical_glycopeptide_ids)


class IonMatching(PipelineModule):
    def __init__(self, database_path, hypothesis_id,
                 observed_ions_path,
                 observed_ions_type='bupid_yaml',
                 sample_run_id=None,
                 hypothesis_sample_match_id=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
                 intensity_threshold=150.0,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()
        self.hypothesis_id = hypothesis_id
        self.n_processes = n_processes

        self.ms1_tolerance = ms1_tolerance
        self.ms2_tolerance = ms2_tolerance
        self.observed_ions_path = observed_ions_path
        self.observed_ions_type = observed_ions_type
        self.intensity_threshold = intensity_threshold

        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.sample_run_id = sample_run_id

        self.msmsdb = MSMSSqlDB(observed_ions_path)

    def prepare_task_fn(self):
        task_fn = functools.partial(batch_match_fragments,
                                    msmsdb_path=self.msmsdb.path,
                                    ms1_tolerance=self.ms1_tolerance,
                                    ms2_tolerance=self.ms2_tolerance,
                                    database_manager=self.manager,
                                    hypothesis_sample_match_id=self.hypothesis_sample_match_id,
                                    sample_run_id=self.sample_run_id,
                                    hypothesis_id=self.hypothesis_id,
                                    intensity_threshold=self.intensity_threshold)
        return task_fn

    def stream_theoretical_glycopeptides(self, chunksize=500):
        session = self.manager.session()
        try:
            for name, protein_id in session.query(
                    Protein.name, Protein.id).filter(Protein.hypothesis_id == self.hypothesis_id):
                logger.info("Streaming %s (%d)", name, protein_id)
                theoretical_glycopeptide_ids = (session.query(
                       TheoreticalGlycopeptide.id).filter(TheoreticalGlycopeptide.protein_id == protein_id))
                chunk = []
                i = 0
                for theoretical_id in itertools.chain.from_iterable(theoretical_glycopeptide_ids):
                    chunk.append(theoretical_id)
                    i += 1
                    if i == chunksize:
                        yield chunk
                        chunk = []
                        i = 0
                yield chunk
        except Exception, e:
            logger.exception("An error occurred in stream_theoretical_glycopeptides, %r", locals(), exc_info=e)
            raise
        finally:
            session.close()

    def run(self):
        if self.hypothesis_sample_match_id is not None:
            for hsm in self.session.query(
                    HypothesisSampleMatch).filter(HypothesisSampleMatch.id == self.hypothesis_sample_match_id):
                hsm.parameters.update({
                        "ms1_ppm_tolerance": self.ms1_tolerance,
                        "ms2_ppm_tolerance": self.ms2_tolerance,
                        "intensity_threshold": self.intensity_threshold,
                        "observed_ions_path": self.observed_ions_path
                })
                self.session.add(hsm)
        self.session.commit()
        session = self.session
        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap(task_fn, self.stream_theoretical_glycopeptides(500)):
                cntr += res
                if (cntr - last) > 1000:
                    logger.info("%d Searches Complete." % cntr)
                    last = cntr
            pool.close()
            pool.join()
        else:
            for theoretical in self.stream_theoretical_glycopeptides():
                cntr += task_fn(theoretical)
                if (cntr - last) > 1000:
                    logger.info("%d Searches Complete." % cntr)
                    last = cntr
        session.commit()
        session.close()
