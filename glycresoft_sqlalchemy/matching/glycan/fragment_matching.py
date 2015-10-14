import math
import multiprocessing
import functools
import itertools
import logging
import operator
try:
    logger = logging.getLogger(__name__)
except:
    pass
from os.path import splitext
from collections import defaultdict

from ..spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser

from ..data_model import (
    MS2GlycanHypothesis, MS2GlycanHypothesisSampleMatch, PipelineModule, MSMSSqlDB,
    TheoreticalGlycanStructure, GlycanStructureMatch, GlycanSpectrumMatch
    )

from ..utils.common_math import ppm_error

neutral_mass_getter = operator.attrgetter("neutral_mass")
key_getter = operator.itemgetter('key')

decon_format_lookup = {
    "bupid_yaml": BUPIDMSMSYamlParser,
    "db": MSMSSqlDB,
}

PROTON = 1.007276035

ms1_tolerance_default = 1e-5
ms2_tolerance_default = 2e-5


def match_fragments(theoretical, msmsdb_path, ms1_tolerance, ms2_tolerance,
                    database_manager, kinds, hypothesis_sample_match_id, sample_run_id,
                    hypothesis_id):
    '''
    *Task Function*

    Arguments
    ---------
    theoretical: int
        ID value for the theoretical sequence to search
    msmsdb_path: str
        Path to the ion database
    ms1_tolerance: float
        ppm mass error tolerance for MS1 matching
    ms2_tolerance: float
        ppm mass error tolerance for MS2 matching
    hypothesis_sample_match_id: int
        If set, associate all GlycopeptideMatches with the indicated HypothesisSampleMatch
    sample_run_id: int
        If set, only select ions from the indicated SampleRun

    '''
    try:
        msmsdb = MSMSSqlDB(msmsdb_path)
        session = database_manager.session()
        theoretical = session.query(TheoreticalGlycanStructure).get(theoretical)

        # Localized global references
        lppm_error = ppm_error
        lfabs = math.fabs

        spectrum_matches = []
        fragments = None
        ppm_errors = []
        observed_masses = []

        query = msmsdb.ppm_match_tolerance_search(theoretical.calculated_mass, ms1_tolerance, sample_run_id)

        matches = []
        collect = matches.append

        for spectrum in query:
            peak_list = spectrum.tandem_data
            peak_list = sorted(peak_list, key=neutral_mass_getter)
            peak_match_map = defaultdict(list)

            ppm_errors.append(ppm_error(theoretical.calculated_mass, spectrum.precursor_neutral_mass))
            observed_masses.append(spectrum.precursor_neutral_mass)

            c = 0
            if fragments is None:
                fragments = list(theoretical.fragments(kinds))
            for theoretical_ion in fragments:
                query_mass = theoretical_ion.mass
                for peak in peak_list:
                    observed_mass = peak.neutral_mass
                    match_error = lppm_error(observed_mass, query_mass)
                    if lfabs(match_error) <= ms2_tolerance:
                        c += 1
                        match = ({'key': theoretical_ion.name,
                                  "observed_mass": observed_mass,
                                  "intensity": peak.intensity,
                                  'ppm_error': match_error, "peak_id": peak.id})
                        collect(match)
                        peak_match_map[peak.id].append(match)

            if c > 0:
                spectrum_matches.append((spectrum, peak_match_map))

        if len(spectrum_matches) > 0:
            scan_ids = []
            # session = database_manager.session()
            for spectrum, peak_match_map in spectrum_matches:
                scan_ids.append(spectrum.time)

            first_scan = min(scan_ids)
            last_scan = max(scan_ids)

            structure_match = GlycanStructureMatch(
                theoretical_reference_id=theoretical.id,
                hypothesis_sample_match_id=hypothesis_sample_match_id,
                ms1_score=None,
                ms2_score=None,
                observed_mass=min(observed_masses),
                ppm_error=min(ppm_errors),
                volume=None,
                fragment_matches=matches,
                scan_id_range=scan_ids,
                first_scan=first_scan,
                last_scan=last_scan
            )

            session.add(structure_match)
            session.commit()

            for match in spectrum_matches:
                sm = GlycanSpectrumMatch(
                    glycopeptide_match_id=structure_match.id,
                    scan_time=spectrum.time, peak_match_map=peak_match_map,
                    peaks_explained=len(peak_match_map),
                    peaks_unexplained=len(spectrum.tandem_data) - len(peak_match_map),
                    hypothesis_sample_match_id=hypothesis_sample_match_id,
                    hypothesis_id=hypothesis_id)
                session.add(sm)
            session.commit()
            session.close()

    except Exception, e:
        logger.exception("An error occurred", exc_info=e)
        raise

    finally:
        session.close()
        msmsdb.close()


class IonMatching(PipelineModule):
    def __init__(self, database_path, hypothesis_id,
                 observed_ions_path,
                 observed_ions_type='bupid_yaml',
                 sample_run_id=None,
                 hypothesis_sample_match_id=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
                 kinds="BY",
                 n_processes=4):
        self.manager = self.manager_type(database_path)
        self.session = self.manager.session()
        self.hypothesis_id = hypothesis_id
        self.n_processes = n_processes

        self.ms1_tolerance = ms1_tolerance
        self.ms2_tolerance = ms2_tolerance
        self.observed_ions_path = observed_ions_path
        self.observed_ions_type = observed_ions_type

        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.sample_run_id = sample_run_id
        self.kinds = kinds

        if isinstance(observed_ions_path, str):
            if observed_ions_type != "db" and splitext(observed_ions_path)[1] != '.db':
                msmsdb = decon_format_lookup[observed_ions_type](observed_ions_path).to_db()
            else:
                msmsdb = MSMSSqlDB(observed_ions_path)
        elif isinstance(observed_ions_path, MSMSSqlDB):
            msmsdb = observed_ions_path

        self.msmsdb = msmsdb

        for hypothesis in self.session.query(MS2GlycanHypothesisSampleMatch).filter(
                MS2GlycanHypothesisSampleMatch.id == MS2GlycanHypothesisSampleMatch):
            hypothesis.parameters.update({
                    "ms1_ppm_tolerance": ms1_tolerance,
                    "ms2_ppm_tolerance": ms2_tolerance,
                    "observed_ions_path": observed_ions_path
            })
            self.session.add(hypothesis)
        self.session.commit()

    def prepare_task_fn(self):
        task_fn = functools.partial(match_fragments,
                                    msmsdb_path=self.msmsdb.path,
                                    ms1_tolerance=self.ms1_tolerance,
                                    ms2_tolerance=self.ms2_tolerance,
                                    database_manager=self.manager,
                                    hypothesis_sample_match_id=self.hypothesis_sample_match_id,
                                    sample_run_id=self.sample_run_id,
                                    hypothesis_id=self.hypothesis_id,
                                    kinds=self.kinds)
        return task_fn

    def stream_theoretical_glycan_structures(self):
        session = self.manager.session()
        try:
            for theoretical_id, in session.query(
                    TheoreticalGlycanStructure.id).filter(
                    TheoreticalGlycanStructure.hypothesis_id == self.hypothesis_id):
                yield theoretical_id
        except Exception, e:
            logger.exception("An error occurred in stream_theoretical_glycopeptides, %r", locals(), exc_info=e)
            raise
        finally:
            session.close()

    def run(self):
        session = self.session
        task_fn = self.prepare_task_fn()
        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap(task_fn, self.stream_theoretical_glycan_structures(), 500):
                cntr += res
                if cntr % 1000 == 0:
                    logger.info("%d Searches Complete." % cntr)
            pool.close()
            pool.join()
        else:
            for theoretical in self.stream_theoretical_glycan_structures():
                cntr += task_fn(theoretical)
                if cntr % 1000 == 0:
                    logger.info("%d Searches Complete." % cntr)
        session.commit()
        session.close()
