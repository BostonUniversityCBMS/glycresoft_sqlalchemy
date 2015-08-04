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

from ..utils import collectiontools
from ..spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser

from ..scoring import score_matches
from .. import data_model as model
from ..data_model import PipelineModule, MSMSSqlDB, TandemScan


HypothesisSampleMatch = model.HypothesisSampleMatch
Hypothesis = model.Hypothesis
SampleRun = model.SampleRun
Protein = model.Protein
TheoreticalGlycopeptide = model.TheoreticalGlycopeptide
SpectrumMatch = model.SpectrumMatch
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


def ppm_error(x, y):
    return (x - y) / y


def match_fragments(theoretical, msmsdb_path, ms1_tolerance, ms2_tolerance,
                    database_manager, hypothesis_sample_match_id, sample_run_id):
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

        query = msmsdb.ppm_match_tolerance_search(theoretical.calculated_mass, ms1_tolerance, sample_run_id)

        for spectrum in query:
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
                                  "intensity": peak.intensity,
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
                                  "intensity": peak.intensity,
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
                                  "intensity": peak.intensity,
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
                                  "intensity": peak.intensity,
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
                                  "intensity": peak.intensity,
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
                                  "intensity": peak.intensity,
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
                scan_ids.append(spectrum.id)

            first_scan = min(scan_ids)
            last_scan = max(scan_ids)

            oxonium_ions = merge_ion_matches(oxonium_ions)
            bare_b_ions = merge_ion_matches(bare_b_ions)
            bare_y_ions = merge_ion_matches(bare_y_ions)
            glycosylated_b_ions = merge_ion_matches(glycosylated_b_ions)
            glycosylated_y_ions = merge_ion_matches(glycosylated_y_ions)
            stub_ions = merge_ion_matches(stub_ions)

            gpm = GlycopeptideMatch(
                protein_id=theoretical.protein_id,
                theoretical_glycopeptide_id=theoretical.id,
                ms1_score=theoretical.ms1_score,
                observed_mass=theoretical.observed_mass,
                calculated_mass=theoretical.calculated_mass,
                volume=theoretical.volume,
                ppm_error=theoretical.ppm_error,
                glycan_composition_str=theoretical.glycan_composition_str,
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
            score_matches.evaluate(gpm, theoretical)
            session.add(gpm)
            session.commit()
            for spectrum, peak_match_map in spectrum_matches:
                sm = SpectrumMatch(glycopeptide_match_id=gpm.id,
                                   spectrum_id=spectrum.id, peak_match_map=peak_match_map)
                session.add(sm)
            session.commit()
            session.close()
        return 1

    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise e
    finally:
        try:
            session.close()
        except:
            pass


def merge_ion_matches(matches):
    """Group multiple matches to the same fragment

    Parameters
    ----------
    matches : list of dict
        The list of ion matches to group

    Returns
    -------
    list of dict: Merged ion matches
    """
    groups = collectiontools.groupby(matches,
                                     key_getter)
    best_matches = []
    fabs = math.fabs
    for key, matched_key in groups.items():
        best_match = matched_key[0]
        best_ppm = fabs(best_match["ppm_error"])
        peak_map = {best_match["peak_id"]: best_match}
        for match in matched_key[1:]:
            peak_map[match["peak_id"]] = match
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
                best_ppm = fabs(match["ppm_error"])
        best_match = best_match.copy()
        best_match["peak_map"] = peak_map
        best_matches.append(best_match)
    return best_matches


class IonMatching(PipelineModule):
    def __init__(self, database_path, hypothesis_id,
                 observed_ions_path,
                 observed_ions_type='bupid_yaml',
                 sample_run_id=None,
                 hypothesis_sample_match_id=None,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default,
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

        if isinstance(observed_ions_path, str):
            if observed_ions_type != "db" and splitext(observed_ions_path)[1] != '.db':
                msmsdb = decon_format_lookup[observed_ions_type](observed_ions_path).to_db()
            else:
                msmsdb = MSMSSqlDB(observed_ions_path)
        elif isinstance(observed_ions_path, MSMSSqlDB):
            msmsdb = observed_ions_path

        self.msmsdb = msmsdb

        for hypothesis in self.session.query(Hypothesis).filter(Hypothesis.id == hypothesis_id):
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
                                    sample_run_id=self.sample_run_id)
        return task_fn

    def stream_theoretical_glycopeptides(self, chunksize=50):
        session = self.manager.session()
        try:
            for name, protein_id in session.query(
                    Protein.name, Protein.id).filter(Protein.hypothesis_id == self.hypothesis_id):
                logger.info("Streaming %s (%d)", name, protein_id)
                theoretical_glycopeptide_ids = (session.query(
                       TheoreticalGlycopeptide.id).filter(TheoreticalGlycopeptide.protein_id == protein_id))
                for theoretical_id in itertools.chain.from_iterable(theoretical_glycopeptide_ids):
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
                    logger.info("%d Searches Complete." % cntr)
        session.commit()
        session.close()
