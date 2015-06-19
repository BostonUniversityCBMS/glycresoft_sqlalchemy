import math
import multiprocessing
import functools
import itertools
import logging
import random
import operator
try:
    logger = logging.getLogger(__name__)
except:
    pass
from os.path import splitext
from collections import Counter, defaultdict

from glycresoft_ms2_classification.error_code_interface import NoIonsMatchedException
from glycresoft_ms2_classification.utils import collectiontools
from glycresoft_ms2_classification.ms import default_loader, default_serializer
from ..spectra import DeconIOBase, MSMSSqlDB, ObservedPrecursorSpectrum
from ..spectra.bupid_topdown_deconvoluter import BUPIDYamlParser
from glycresoft_ms2_classification.utils.parallel_opener import ParallelParser

from ..scoring import score_matches
from .. import data_model as model

Experiment = model.Experiment
ExperimentParameter = model.ExperimentParameter
Protein = model.Protein
TheoreticalGlycopeptide = model.TheoreticalGlycopeptide
SpectrumMatch = model.SpectrumMatch
GlycopeptideMatch = model.GlycopeptideMatch

neutral_mass_getter = operator.attrgetter("neutral_mass")
key_getter = operator.itemgetter('key')

decon_format_lookup = {
    "bupid_yaml": BUPIDYamlParser,
    "db": MSMSSqlDB,
    "default": default_loader
}

PROTON = 1.007276035

ms1_tolerance_default = 1e-5
ms2_tolerance_default = 2e-5


def ppm_error(x, y):
    return (x - y) / y


def match_fragments(theoretical, msmsdb_path, ms1_tolerance, ms2_tolerance, database_manager):
    msmsdb = MSMSSqlDB(msmsdb_path)

    session = database_manager.session()

    theoretical = session.query(TheoreticalGlycopeptide).filter(
        TheoreticalGlycopeptide.id == theoretical).first()

    # Containers for global theoretical peak matches
    oxonium_ions = []
    bare_b_ions = []
    bare_y_ions = []
    glycosylated_b_ions = []
    glycosylated_y_ions = []
    stub_ions = []

    # Localized global references
    proton = PROTON
    lppm_error = ppm_error
    lfabs = math.fabs

    spectrum_matches = []

    for spectrum in msmsdb.ppm_match_tolerance_search(theoretical.calculated_mass, ms1_tolerance):
        peak_list = spectrum.tandem_data
        peak_list = sorted(peak_list, key=neutral_mass_getter)
        peak_match_map = defaultdict(list)

        c = 0
        collect = oxonium_ions.append
        for theoretical_ion in theoretical.oxonium_ions:
            query_mass = theoretical_ion['mass']
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                protonated_mass = observed_mass + proton
                match_error = lppm_error(protonated_mass, query_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": protonated_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                    c += 1
                elif protonated_mass > query_mass + 10:
                    break

        # If no oxonium ions were found, skip this spectrum
        if c < 1:
            continue

        collect = bare_b_ions.append
        for theoretical_ion in theoretical.bare_b_ions:
            query_mass = theoretical_ion['mass']
            deprotonated_mass = query_mass - proton
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                elif observed_mass > query_mass + 10:
                    break

        collect = bare_y_ions.append
        for theoretical_ion in theoretical.bare_y_ions:
            query_mass = theoretical_ion['mass']
            deprotonated_mass = query_mass - proton
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                elif observed_mass > query_mass + 10:
                    break

        collect = glycosylated_b_ions.append
        for theoretical_ion in theoretical.glycosylated_b_ions:
            query_mass = theoretical_ion['mass']
            deprotonated_mass = query_mass - proton
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                elif observed_mass > query_mass + 10:
                    break

        collect = glycosylated_y_ions.append
        for theoretical_ion in theoretical.glycosylated_y_ions:
            query_mass = theoretical_ion['mass']
            deprotonated_mass = query_mass - proton
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                elif observed_mass > query_mass + 10:
                    break

        collect = stub_ions.append
        for theoretical_ion in theoretical.stub_ions:
            query_mass = theoretical_ion['mass']
            deprotonated_mass = query_mass - proton
            for peak in peak_list:
                observed_mass = peak.neutral_mass
                match_error = lppm_error(observed_mass, deprotonated_mass)
                if lfabs(match_error) <= ms2_tolerance:
                    match = ({'key': theoretical_ion['key'],
                              "observed_mass": observed_mass,
                              'ppm_error': match_error, "peak_id": peak.id})
                    collect(match)
                    peak_match_map[peak.id].append(match)
                elif observed_mass > query_mass + 10:
                    break

        spectrum_matches.append((spectrum, peak_match_map))

    if len(spectrum_matches) > 0:
        scan_ids = []
        # session = database_manager.session()
        for spectrum, peak_match_map in spectrum_matches:
            sm = SpectrumMatch(theoretical_glycopeptide_id=theoretical.id,
                               spectrum_id=spectrum.id, peak_match_map=peak_match_map)
            scan_ids.append(spectrum.id)
            session.add(sm)
        first_scan = min(scan_ids)
        last_scan = max(scan_ids)

        oxonium_ions = merge_ion_matches(oxonium_ions)
        bare_b_ions = merge_ion_matches(bare_b_ions)
        bare_y_ions = merge_ion_matches(bare_y_ions)
        glycosylated_b_ions = merge_ion_matches(glycosylated_b_ions)
        glycosylated_y_ions = merge_ion_matches(glycosylated_y_ions)
        stub_ions = merge_ion_matches(stub_ions)

        gpm = GlycopeptideMatch(
            id=theoretical.id,
            protein_id=theoretical.protein_id,
            ms1_score=theoretical.ms1_score,
            observed_mass=theoretical.observed_mass,
            calculated_mass=theoretical.calculated_mass,
            volume=theoretical.volume,
            ppm_error=theoretical.ppm_error,
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
            last_scan=last_scan
        )
        score_matches.apply(gpm, theoretical)
        session.add(gpm)
        session.commit()

    return 1


def merge_ion_matches(matches):
    groups = collectiontools.groupby(matches,
                                     key_getter)
    best_matches = []
    fabs = math.fabs
    for key, matched_key in groups.items():
        best_match = matched_key[0]
        best_ppm = fabs(best_match["ppm_error"])
        ppm_to_peak_id = {best_match["peak_id"]: best_match["ppm_error"]}
        for match in matched_key[1:]:
            ppm_to_peak_id[match["peak_id"]] = match["ppm_error"]
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
                best_ppm = fabs(match["ppm_error"])
        best_match["scan_map"] = ppm_to_peak_id
        best_matches.append(best_match)
    return best_matches


def _chunk_iter(iterable, size=50):
    results = []
    for entry in iterable:
        results.append(entry)
        if len(results) == size:
            yield results
            results = []
    yield results


class IonMatching(object):

    manager_type = model.DatabaseManager

    def __init__(self, database_path, observed_ions_path, observed_ions_type='bupid_yaml',
                 experiment_ids=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()
        if experiment_ids is None:
            experiment_ids = [exp_id for query in self.session.query(Experiment.id) for exp_id in query]
        self.experiment_ids = experiment_ids
        self.n_processes = n_processes

        self.ms1_tolerance = ms1_tolerance
        self.ms2_tolerance = ms2_tolerance
        self.observed_ions_path = observed_ions_path
        self.observed_ions_type = observed_ions_type

        if isinstance(observed_ions_path, str):
            if observed_ions_type != "db" and splitext(observed_ions_path)[1] != '.db':
                incoming_deconio = ParallelParser(decon_format_lookup[observed_ions_type], (observed_ions_path,))
                msmsdb = None
            else:
                incoming_deconio = None
                msmsdb = MSMSSqlDB(observed_ions_path)
        elif isinstance(observed_ions_path, MSMSSqlDB):
            incoming_deconio = None
            msmsdb = observed_ions_path

        self.incoming_deconio = incoming_deconio
        self.msmsdb = msmsdb

        for exp_id in self.experiment_ids:
            for exp_param in self.session.query(ExperimentParameter).filter(Experiment.id == exp_id):
                exp_param.value.update({
                        "ms1_ppm_tolerance": ms1_tolerance,
                        "ms2_ppm_tolerance": ms2_tolerance,
                        "observed_ions_path": observed_ions_path
                })
                self.session.add(exp_param)
        self.session.commit()

    def prepare_task_fn(self):
        task_fn = functools.partial(match_fragments,
                                    msmsdb_path=self.msmsdb.connection_string,
                                    ms1_tolerance=self.ms1_tolerance,
                                    ms2_tolerance=self.ms2_tolerance,
                                    database_manager=self.manager)
        return task_fn

    def stream_theoretical_glycopeptides(self, chunksize=50):
        session = self.manager.session()
        for exp_id in self.experiment_ids:
            for name, protein_id in session.query(
                    Protein.name, Protein.id).filter(Protein.experiment_id == exp_id):
                print(name)
                theoretical_glycopeptide_ids = (session.query(
                       TheoreticalGlycopeptide.id).filter(TheoreticalGlycopeptide.protein_id == protein_id))
                ## When the number of processes used is small, the more work done in the main process.
                # for theoretical_ids in _chunk_iter(
                #         itertools.chain.from_iterable(theoretical_glycopeptide_ids), chunksize):
                #     for theoretical in self.session.query(
                #             TheoreticalGlycopeptide).filter(TheoreticalGlycopeptide.id.in_(theoretical_ids)):
                #         yield theoretical
                ## When the number of processes is large, do more work in the child processes
                for theoretical_id in itertools.chain.from_iterable(theoretical_glycopeptide_ids):
                    yield theoretical_id
        session.close()

    def run(self):
        if self.msmsdb is None:
            decon_io = self.incoming_deconio.await()
            database = decon_io.to_db()
            self.msmsdb = database
        task_fn = self.prepare_task_fn()
        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap(task_fn, self.stream_theoretical_glycopeptides(500), 500):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Searches Complete." % cntr)
            pool.close()
            pool.join()
        else:
            for theoretical in self.stream_theoretical_glycopeptides():
                cntr += task_fn(theoretical)
                if cntr % 1000 == 0:
                    print "%d Searches Complete." % cntr
