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

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, MSMSSqlDB, TandemScan, slurp, HypothesisSampleMatch,
    Protein, TheoreticalGlycopeptide, GlycopeptideSpectrumMatch,
    GlycopeptideMatch,
    )

from glycresoft_sqlalchemy.utils.common_math import ppm_error, median


neutral_mass_getter = operator.attrgetter("neutral_mass")
key_getter = operator.itemgetter('key')

decon_format_lookup = {
    "bupid_yaml": BUPIDMSMSYamlParser,
    "db": MSMSSqlDB,
}

PROTON = 1.007276035

ms1_tolerance_default = 1e-5
ms2_tolerance_default = 2e-5


def _sweep_solution(array, value, lo, hi, tolerance, verbose=False):
    best_index = -1
    best_error = float('inf')
    for i in range(hi - lo):
        target = array[lo + i]
        error = ppm_error(value, target.neutral_mass)
        abs_error = abs(error)
        if abs_error < tolerance and abs_error < best_error:
            best_index = lo + i
            best_error = abs_error
    if best_index == -1:
        return None
    else:
        return array[best_index]


def _binary_search(array, value, lo, hi, tolerance, verbose=False):
    if (hi - lo) < 5:
        return _sweep_solution(array, value, lo, hi, tolerance, verbose)
    else:
        mid = (hi + lo) / 2
        target = array[mid]
        target_value = target.neutral_mass
        error = ppm_error(value, target_value)

        if abs(error) <= tolerance:
            return _sweep_solution(array, value, max(mid - 5, lo), min(mid + 5, hi), tolerance, verbose)
        elif target_value > value:
            return _binary_search(array, value, lo, mid, tolerance, verbose)
        elif target_value < value:
            return _binary_search(array, value, mid, hi, tolerance, verbose)
    raise Exception("No recursion found!")


def binary_search(array, value, tolerance=2e-5, verbose=False):
    size = len(array)
    if size == 0:
        return None
    return _binary_search(array, value, 0, size, tolerance, verbose)


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

        theoreticals = slurp(session, TheoreticalGlycopeptide, theoretical_ids, flatten=False)

        for theoretical in theoreticals:

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
                if len(peak_list) == 0:
                    continue
                intensity_threshold = median(p.intensity for p in peak_list)
                peak_list = [p for p in peak_list if p.intensity >= intensity_threshold]
                peak_list = sorted(peak_list, key=neutral_mass_getter)
                peak_match_map = defaultdict(list)
                precursor_ppm_error = lppm_error(theoretical.calculated_mass, spectrum.precursor_neutral_mass)

                oxonium_ion_count = 0
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
                            oxonium_ion_count += 1
                        elif protonated_mass > query_mass + 10:
                            break

                # If no oxonium ions were found, skip this spectrum
                if oxonium_ion_count < 1:
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

                spectrum_matches.append((spectrum, peak_match_map, precursor_ppm_error, oxonium_ion_count))

            if len(spectrum_matches) > 0:
                scan_ids = []
                # session = database_manager.session()
                for spectrum, peak_match_map, precursor_ppm_error, oxcount in spectrum_matches:
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

                for spectrum, peak_match_map, precursor_ppm_error, oxcount in spectrum_matches:
                    spectrum_match_inst = GlycopeptideSpectrumMatch(
                        scan_time=spectrum.time, peak_match_map=peak_match_map,
                        precursor_charge_state=spectrum.precursor_charge_state,
                        precursor_ppm_error=precursor_ppm_error,
                        peaks_explained=len(peak_match_map) - oxcount,
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
                    spectrum_match_inst.theoretical_glycopeptide_id = g.theoretical_glycopeptide_id
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
            if len(peak_list) == 0:
                continue
            intensity_threshold = median(p.intensity for p in peak_list)
            peak_list = [p for p in peak_list if p.intensity >= intensity_threshold]
            peak_list = sorted(peak_list, key=neutral_mass_getter)

            for theoretical in query:
                precursor_ppm_error = lppm_error(theoretical.calculated_mass, spectrum.precursor_neutral_mass)
                peak_match_map = defaultdict(list)

                oxonium_ion_count = 0
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
                            oxonium_ion_count += 1
                        elif protonated_mass > query_mass + 10:
                            break

                # If no oxonium ions were found, skip this spectrum
                if oxonium_ion_count < 1:
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

                spectrum_matches.append((theoretical, peak_match_map, precursor_ppm_error, oxonium_ion_count))

            if len(spectrum_matches) > 0:
                for theoretical, peak_match_map, precursor_ppm_error, oxcount in spectrum_matches:
                    spectrum_match_inst = GlycopeptideSpectrumMatch(
                        scan_time=spectrum.time, peak_match_map=peak_match_map,
                        precursor_charge_state=spectrum.precursor_charge_state,
                        precursor_ppm_error=precursor_ppm_error,
                        peaks_explained=len(peak_match_map) - oxcount,
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


def search_spectrum(theoretical, spectrum):
    peak_list = spectrum.tandem_data
    peak_list = [p for p in peak_list if p.intensity >= 150.]
    peak_list = sorted(peak_list, key=neutral_mass_getter)
    peak_match_map = defaultdict(list)

    lppm_error = ppm_error
    lfabs = abs
    ms2_tolerance = 2e-5

    solution = defaultdict(list)
    for series_name, ions in theoretical.items():
        collect = solution[series_name].append

        for theoretical_ion in ions:
            query_mass = theoretical_ion.mass
            deprotonated_mass = query_mass
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion.name,
                              "intensity": peak.intensity,
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.scan_peak_index})
                    collect(match)
                    peak_match_map[peak.scan_peak_index].append(match)
                elif observed_mass > query_mass + 10:
                    break

    return solution, peak_match_map


class IonMatching(PipelineModule):
    def __init__(self, database_path, hypothesis_id,
                 observed_ions_path,
                 observed_ions_type='bupid_yaml',
                 sample_run_id=None,
                 hypothesis_sample_match_id=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
                 intensity_threshold=0.0,
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


class SpectrumMatching(PipelineModule):
    def __init__(self, database_path, hypothesis_id,
                 observed_ions_path,
                 observed_ions_type='bupid_yaml',
                 sample_run_id=None,
                 hypothesis_sample_match_id=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
                 intensity_threshold=0.0,
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
        task_fn = functools.partial(batch_match_theoretical_ions,
                                    msmsdb_path=self.msmsdb.path,
                                    ms1_tolerance=self.ms1_tolerance,
                                    ms2_tolerance=self.ms2_tolerance,
                                    database_manager=self.manager,
                                    hypothesis_sample_match_id=self.hypothesis_sample_match_id,
                                    sample_run_id=self.sample_run_id,
                                    hypothesis_id=self.hypothesis_id,
                                    intensity_threshold=self.intensity_threshold)
        return task_fn

    def stream_tandem_spectra(self, chunksize=100):
        session = self.msmsdb.session()
        try:
            scan_ids = session.query(TandemScan.id).filter(TandemScan.sample_run_id == self.sample_run_id).all()
            last = 0
            total = len(scan_ids)
            while last <= total:
                yield scan_ids[last:(last + chunksize)]
                last += chunksize

        except Exception, e:
            logger.exception("An error occurred in stream_tandem_spectra, %r", locals(), exc_info=e)
            raise
        finally:
            session.close()

    def run(self):
        if self.hypothesis_sample_match_id is not None:
            hsm = self.session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
            hsm.parameters.update({
                    "ms1_ppm_tolerance": self.ms1_tolerance,
                    "ms2_ppm_tolerance": self.ms2_tolerance,
                    "intensity_threshold": self.intensity_threshold,
                    "observed_ions_path": self.observed_ions_path
            })
            self.session.add(hsm)
            self.session.commit()

        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_tandem_spectra()):
                cntr += res
                if (cntr - last) > 100:
                    logger.info("%d Searches Complete." % cntr)
                    last = cntr
            pool.close()
            pool.join()
            pool.terminate()
        else:
            for theoretical in self.stream_tandem_spectra():
                cntr += task_fn(theoretical)
                if (cntr - last) > 100:
                    logger.info("%d Searches Complete." % cntr)
                    last = cntr
