import csv
import os
import re
import datetime
import multiprocessing
import logging
import functools

from glycresoft_sqlalchemy.structure.modification import RestrictedModificationTable
from glycresoft_sqlalchemy.structure.modification import ModificationTable
from glycresoft_sqlalchemy.structure.sequence_space import SequenceSpace
from glycresoft_sqlalchemy.structure.stub_glycopeptides import StubGlycopeptide
from glycresoft_sqlalchemy.structure import constants
from glycresoft_sqlalchemy.proteomics import get_enzyme, msdigest_xml_parser

from glycresoft_sqlalchemy import data_model as model
from ..peptide_utilities import SiteListFastaFileParser
from glycresoft_sqlalchemy.data_model import (
    PipelineModule, Hypothesis, MS2GlycopeptideHypothesis,
    HypothesisSampleMatch, PeakGroupMatch, Protein,
    TheoreticalGlycopeptideGlycanAssociation,
    MS1GlycopeptideHypothesis,
    TheoreticalGlycopeptide, Hierarchy)

logger = logging.getLogger("search_space_builder")
mod_pattern = re.compile(r'(\d+)(\w+)')
g_colon_prefix = "G:"


def parse_site_file(site_lists):
    site_list_map = {}
    for entry in SiteListFastaFileParser(site_lists):
        site_list_map[entry['name']] = list(entry["glycosylation_sites"])
    return site_list_map


def _get_glycan_composition(mapping, glycan_identities):
    """Given a mapping and list of monosaccharide names, get the mapped integer counts
    for each monosaccharide
    """
    try:
        return {g: int(mapping[g_colon_prefix + g]) for g in glycan_identities}
    except:
        global g_colon_prefix
        g_colon_prefix = ""
        return {g: int(mapping[g_colon_prefix + g]) for g in glycan_identities}


class MS1GlycopeptideResult(object):
    '''
    Describes a row ofthe MS1 Results with format-agnostic basic
    construction and a CSV specific mapping
    '''

    @classmethod
    def interpolate(cls, source_type, monosaccharide_identities, *args, **kwargs):
        if source_type is MS1ResultsFile:
            return cls.from_csvdict(monosaccharide_identities, *args, **kwargs)
        elif source_type is MS1ResultsFacade:
            return cls.from_peak_group_match(*args, **kwargs)

    @classmethod
    def from_csvdict(cls, glycan_identities=None, **kwargs):
        if glycan_identities is None:
            glycan_identities = ()
        score = float(kwargs.get("Score", 0.))
        calculated_mass = float(kwargs.get("Hypothesis MW"))
        observed_mass = float(kwargs.get("MassSpec MW"))
        glycan_composition_str = kwargs.get("Compound Key")
        base_peptide_sequence = kwargs.get("PeptideSequence")
        peptide_modifications = kwargs.get("PeptideModification")
        count_missed_cleavages = int(kwargs.get("PeptideMissedCleavage#"))
        count_glycosylation_sites = int(kwargs.get("#ofGlycanAttachmentToPeptide"))
        ppm_error = float(kwargs.get("PPM Error"))
        volume = float(kwargs.get("Total Volume"))
        start_pos = int(kwargs.get("StartAA"))
        end_pos = int(kwargs.get("EndAA"))
        protein_name = kwargs.get("ProteinID", None)
        glycan_composition_map = _get_glycan_composition(kwargs, glycan_identities)
        return cls(
            score=score, calculated_mass=calculated_mass,
            observed_mass=observed_mass, glycan_composition_str=glycan_composition_str,
            base_peptide_sequence=base_peptide_sequence, modified_peptide_sequence=base_peptide_sequence,
            peptide_modifications=peptide_modifications,
            count_missed_cleavages=count_missed_cleavages, count_glycosylation_sites=count_glycosylation_sites,
            ppm_error=ppm_error, volume=volume, start_pos=start_pos, end_pos=end_pos,
            glycan_composition_map=glycan_composition_map, protein_name=protein_name,
            composition_id=None)

    @classmethod
    def from_peak_group_match(cls, peak_group_match=None, theoretical=None):
        pgm = peak_group_match
        theoretical = pgm.theoretical_match if theoretical is None else theoretical
        if theoretical is None:
            raise ValueError("PeakGroupMatch {} did not map to a real theoretical composition".format(pgm))
        return cls(
            score=pgm.ms1_score, calculated_mass=theoretical.calculated_mass,
            observed_mass=pgm.weighted_monoisotopic_mass,
            glycan_composition_str=theoretical.glycan_composition_str,
            base_peptide_sequence=theoretical.base_peptide_sequence,
            modified_peptide_sequence=theoretical.modified_peptide_sequence,
            glycopeptide_sequence=theoretical.glycopeptide_sequence,
            peptide_modifications=theoretical.peptide_modifications,
            count_missed_cleavages=theoretical.count_missed_cleavages,
            count_glycosylation_sites=theoretical.count_glycosylation_sites,
            ppm_error=pgm.ppm_error, volume=pgm.total_volume,
            start_pos=theoretical.start_position, end_pos=theoretical.end_position,
            protein_name=theoretical.protein.name, composition_id=theoretical.id)

    def __init__(self, score=None,
                 calculated_mass=None,
                 observed_mass=None,
                 glycan_composition_str=None,
                 base_peptide_sequence=None,
                 peptide_modifications=None,
                 count_missed_cleavages=None,
                 count_glycosylation_sites=None,
                 ppm_error=None,
                 volume=None,
                 start_pos=None,
                 end_pos=None,
                 glycan_composition_map=None,
                 protein_name=None,
                 fragments=None,
                 oxonium_ions=None,
                 id=None,
                 glycopeptide_sequence=None,
                 modified_peptide_sequence=None,
                 composition_id=None):
        self.ms1_score = score
        self.calculated_mass = calculated_mass
        self.observed_mass = observed_mass
        self.glycan_composition_str = glycan_composition_str
        self.base_peptide_sequence = base_peptide_sequence
        self.peptide_modifications = peptide_modifications
        self.count_missed_cleavages = count_missed_cleavages
        self.count_glycosylation_sites = count_glycosylation_sites
        self.ppm_error = ppm_error
        self.volume = volume
        self.start_position = start_pos
        self.end_position = end_pos
        self.glycan_composition_map = glycan_composition_map
        self.protein_name = protein_name
        self.id = id
        self.glycopeptide_sequence = glycopeptide_sequence
        self.modified_peptide_sequence = modified_peptide_sequence
        self.composition_id = composition_id

    @property
    def most_detailed_sequence(self):
        try:
            s = self.glycopeptide_sequence
            if s is not None and s != "":
                return s
        except:
            try:
                s = self.modified_peptide_sequence
                if s is not None and s != "":
                    return s
            except:
                return self.base_peptide_sequence

    def __repr__(self):
        return "MS1GlycopeptideResult(%s)" % str((self.most_detailed_sequence, self.peptide_modifications, self.glycan_composition_str))


def get_monosaccharide_identities(csv_columns):
    """Extract monosaccharide names from the column headers of MS1 Results CSV files

    Parameters
    ----------
    csv_columns: list
        A list of strings containing column names from the MS1 Results.

    Returns
    -------
    monosaccharide_identities : list
        The list of glycan identities
    """
    monosaccharide_identities = []
    extract_state = False
    for col in csv_columns:
        if col == "Hypothesis MW":
            extract_state = True
            continue
        elif col == "Adduct/Replacement":
            extract_state = False
            break
        elif extract_state:
            monosaccharide_identities.append(col.replace("G:", ""))
    logger.info("Glycan identities found: %s", monosaccharide_identities)
    return monosaccharide_identities


def get_peptide_modifications(peptide_mod_str, modification_table):
    """Maps the Protein Prospector modifications to :class:`.modifications.Modification`
    objects through `modification_table`

    Parameters
    ----------
    peptide_mod_str : str
        A string containing peptide modifications from ProteinProspector in the MS1 CSV
    modification_table : :class:`.modifications.ModificationTable`
        A mapping of possible :class:`.modifications.ModificationRule` functions provided
        by the user.

    Returns
    -------
    :class:`list` of :class:`.modifications.Modification` objects
    """
    items = mod_pattern.findall(peptide_mod_str)
    mod_list = []
    for i in items:
        if i[1] == '':
            continue
        for part in i[1].split(":"):
            mod = modification_table.get_modification(part, -1, int(i[0]))
            mod_list.append(mod)
    return mod_list


def get_search_space(ms1_result, glycan_sites, mod_list):
    """Create a :class:`.sequence_space.SequenceSpace` object from the collection
    of information interpolated from the MS1 Results

    Parameters
    ----------
    ms1_result : :class:`MS1GlycopeptideResult`
        The row of MS1 Results being operated on
    glycan_sites : list
        A list of sites along the parent sequence that are glycosylated
    mod_list : list
        The list of present modifications

    Returns
    -------
    :class:`.sequence_space.SequenceSpace`
        The generator of possible :class:`.sequence.Sequence` objects
    """

    seq_space = SequenceSpace(
        ms1_result.modified_peptide_sequence,
        ms1_result.glycan_composition_map,
        glycan_sites, mod_list)
    return seq_space


def generate_fragments(seq, ms1_result):
    """Consumes a :class:`.sequence.Sequence` object, and the contents of an MS1 Result row to
    generate the set of all theoretically observed fragments

    Parameters
    ----------
    seq: :class:`.sequence.Sequence`
        The binding of modifications to particular sites on a peptide sequence
    ms1_result: :class:`MS1GlycopeptideResult`
        Description of the precursor match

    Returns
    -------
    dict:
        Collection of theoretical ions from the given sequence,
        as well as the precursor information.
    """
    seq_mod = seq.get_sequence()
    fragments = zip(*map(seq.break_at, range(1, len(seq))))
    b_type = fragments[0]
    b_ions = []
    b_ions_hexnac = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if key == ("b1" or re.search(r'b1\+', key)) and constants.EXCLUDE_B1:
                # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                b_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]  # seq.get_fragments('Y')
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    pep_stubs = StubGlycopeptide(
        ms1_result.base_peptide_sequence,
        ms1_result.peptide_modifications,
        ms1_result.count_glycosylation_sites,
        seq.glycan)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()

    theoretical_glycopeptide = model.TheoreticalGlycopeptide(
        ms1_score=ms1_result.ms1_score,
        observed_mass=ms1_result.observed_mass,
        calculated_mass=ms1_result.calculated_mass,
        ppm_error=ms1_result.ppm_error,
        volume=ms1_result.volume,
        count_glycosylation_sites=ms1_result.count_glycosylation_sites,
        count_missed_cleavages=ms1_result.count_missed_cleavages,
        start_position=ms1_result.start_position,
        end_position=ms1_result.end_position,
        base_peptide_sequence=ms1_result.base_peptide_sequence,
        modified_peptide_sequence=seq_mod,
        peptide_modifications=ms1_result.peptide_modifications,
        glycopeptide_sequence=seq_mod + ms1_result.glycan_composition_str,
        sequence_length=len(seq),
        glycan_composition_str=ms1_result.glycan_composition_str,
        bare_b_ions=b_ions,
        bare_y_ions=y_ions,
        oxonium_ions=oxonium_ions,
        stub_ions=stub_ions,
        glycosylated_b_ions=b_ions_hexnac,
        glycosylated_y_ions=y_ions_hexnac,
        base_composition_id=ms1_result.composition_id
        )
    return theoretical_glycopeptide


def process_predicted_ms1_ion(row, modification_table, site_list_map,
                              renderer, database_manager, proteins):
    """Multiprocessing dispatch function to generate all theoretical sequences and their
    respective fragments from a given MS1 result

    Parameters
    ----------
    row: dict
        Line mapping from an MS1 results csv file
    modification_table: :class:`.modifications.RestrictedModificationTable`
        Modification table limited to only the rules specified by the user
    site_list: list of int
        List of putative glycosylation sites from the parent protein
    monosaccharide_identities: list
        List of glycan or monosaccaride names

    Returns
    -------
    list of dicts:
        List of each theoretical sequence and its fragment ions
    """
    try:
        session = database_manager.session()
        ms1_result = renderer.render(session, row)
        if (ms1_result.base_peptide_sequence == '') or (ms1_result.count_glycosylation_sites == 0):
            return 0

        # Compute the set of modifications that can occur.
        mod_list = get_peptide_modifications(
            ms1_result.peptide_modifications, modification_table)

        # Get the start and end positions of fragment relative to the
        glycan_sites = set(site_list_map.get(ms1_result.protein_name, [])).intersection(
            range(ms1_result.start_position, ms1_result.end_position))

        # No recorded sites, skip this component.
        if len(glycan_sites) == 0:
            return 0

        # Adjust the glycan_sites to relative position
        glycan_sites = [x - ms1_result.start_position for x in glycan_sites]
        ss = get_search_space(
            ms1_result, glycan_sites, mod_list)
        seq_list = ss.get_theoretical_sequence(ms1_result.count_glycosylation_sites)
        fragments = [generate_fragments(seq, ms1_result)
                     for seq in seq_list]
        i = 0
        for sequence in fragments:
            sequence.protein_id = proteins[ms1_result.protein_name].id
            session.add(sequence)
            i += 1
        session.commit()
        return i
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise


constructs = Hierarchy()


class TheoreticalSearchSpaceBuilder(PipelineModule):
    '''
    Describe the process of generating all theoretical sequences and their fragments
    from an MS1 Results CSV, a collection of constant and variable peptide modifications,
    and a list of putative glycosylation sites.

    Uses an :class:`sqlitedict.SqliteDict` instance as a storage backend.

    Includes a single- and multi-process compatible implementation. The more processes used,
    the more memory must be allocated to buffer results.
    '''

    HypothesisType = MS2GlycopeptideHypothesis

    @classmethod
    def from_hypothesis(cls, database_path, hypothesis_sample_match_id, n_processes=4):
        dbm = cls.manager_type(database_path)
        session = dbm.session()
        source_hypothesis_sample_match = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
        source_hypothesis = source_hypothesis_sample_match.target_hypothesis
        hypothesis_id = source_hypothesis.id
        variable_modifications = source_hypothesis.parameters.get("variable_modifications", [])
        constant_modifications = source_hypothesis.parameters.get("constant_modifications", [])
        enzyme = source_hypothesis.parameters.get("enzyme", [])
        inst = cls(
            database_path, database_path, enzyme=[enzyme], site_list=None,
            constant_modifications=constant_modifications,
            variable_modifications=variable_modifications,
            n_processes=n_processes, ms1_format=MS1ResultsFacade,
            source_hypothesis_id=hypothesis_id,
            source_hypothesis_sample_match_id=hypothesis_sample_match_id)
        session.close()
        return inst

    def __init__(self, ms1_results_file, db_file_name=None,
                 enzyme=None,
                 site_list=None,
                 constant_modifications=None,
                 variable_modifications=None, n_processes=4,
                 **kwargs):
        self.ms1_format = kwargs.get("ms1_format", MS1ResultsFile)

        if db_file_name is None:
            db_file_name = os.path.splitext(ms1_results_file)[0] + '.db'
        self.db_file_name = db_file_name

        self.manager = self.manager_type(db_file_name)
        self.manager.initialize()
        self.session = self.manager.session()

        self.ms1_results_file = ms1_results_file

        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")

        self.hypothesis = self.HypothesisType(name=u"target-{}".format(tag))
        self.session.add(self.hypothesis)
        self.session.commit()
        self.hypothesis_id = self.hypothesis.id
        logger.info("Building %r", self.hypothesis)

        # We are reading data from separate files
        if self.ms1_format is MS1ResultsFile:
            try:
                site_list_map = parse_site_file(site_list)
            except IOError, e:
                if isinstance(site_list, dict):
                    site_list_map = site_list
                else:
                    raise e
            self.glycosylation_site_map = site_list_map
            self.load_protein_from_sitelist(site_list_map)
            self.ms1_results_reader = self.ms1_format(self.ms1_results_file)
        # We are reading data from another hypothesis in the database
        else:
            source_hypothesis_id = kwargs['source_hypothesis_id']
            source_hypothesis_sample_match_id = kwargs['source_hypothesis_sample_match_id']
            self.load_protein_from_hypothesis(source_hypothesis_id)
            self.ms1_results_reader = self.ms1_format(
                self.manager,
                source_hypothesis_id,
                source_hypothesis_sample_match_id)
            site_list_map = self.glycosylation_site_map

        self.create_restricted_modification_table(constant_modifications, variable_modifications)

        self.n_processes = n_processes

        self.monosaccharide_identities = self.ms1_results_reader.monosaccharide_identities
        enzyme = map(get_enzyme, enzyme)

        self.hypothesis.parameters = {
            "monosaccharide_identities": self.monosaccharide_identities,
            "enzyme": enzyme,
            "site_list_map": site_list_map,
            "constant_modification_list": constant_modifications,
            "variable_modification_list": variable_modifications,
            "ms1_output_file": ms1_results_file,
            "enzyme": enzyme,
            "tag": tag,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS,
            "source_hypothesis_id": kwargs.get("source_hypothesis_id"),
            "source_hypothesis_sample_match_id": kwargs.get("source_hypothesis_sample_match_id")
        }
        self.session.add(self.hypothesis)

    def load_protein_from_sitelist(self, glycosylation_site_map):
        for name, sites in self.glycosylation_site_map.items():
            self.session.add(model.Protein(name=name, glycosylation_sites=sites, hypothesis_id=self.hypothesis.id))

        self.session.commit()

    def load_protein_from_hypothesis(self, hypothesis_id):
        session = self.session
        site_list_map = {}
        hid = self.hypothesis.id
        for protein in session.query(Protein).filter(Protein.hypothesis_id == hypothesis_id):
            self.session.add(model.Protein(
                name=protein.name,
                protein_sequence=protein.protein_sequence,
                glycosylation_sites=protein.glycosylation_sites,
                hypothesis_id=hid))
            site_list_map[protein.name] = protein.glycosylation_sites
        session.add(self.hypothesis)
        session.commit()
        self.glycosylation_site_map = site_list_map

    def create_restricted_modification_table(self, constant_modifications, variable_modifications):
        modification_table = RestrictedModificationTable.bootstrap(constant_modifications, variable_modifications)
        if constant_modifications is None and variable_modifications is None:
            modification_table = ModificationTable.bootstrap()
        self.modification_table = modification_table

    def prepare_task_fn(self):
        task_fn = functools.partial(process_predicted_ms1_ion, modification_table=self.modification_table,
                                    site_list_map=self.glycosylation_site_map,
                                    renderer=self.ms1_format,
                                    database_manager=self.manager, proteins=self.hypothesis.proteins)
        return task_fn

    def stream_results(self):
        for obj in self.ms1_results_reader:
            yield obj

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes
        '''
        task_fn = self.prepare_task_fn()
        cntr = 0
        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in worker_pool.imap_unordered(task_fn, self.stream_results(), chunksize=25):
                cntr += res
                if (cntr % 1000) == 0:
                    logger.info("Committing, %d records made", cntr)

            worker_pool.terminate()
        else:
            logger.debug("Building theoretical sequences sequentially")
            for row in self.ms1_results_reader:
                res = task_fn(row)
                cntr += res
                if (cntr % 1000) == 0:
                    logger.info("Committing, %d records made", cntr)
        return self.hypothesis_id


class MS1ResultsFile(object):
    renderer_id = "MS1ResultsFile"

    def __init__(self, file_path):
        self.file_path = file_path
        self.file_handle = open(file_path, 'rb')
        self.reader = csv.DictReader(self.file_handle)
        self.monosaccharide_identities = get_monosaccharide_identities(self.reader.fieldnames)

    @staticmethod
    def render(session, obj):
        return obj

    @classmethod
    def get_renderer(self):
        return lambda session, obj: obj

    def __iter__(self):
        for row in self.reader:
            try:
                if row["Compound Key"] != "":
                    yield MS1GlycopeptideResult.from_csvdict(self.monosaccharide_identities, **row)
            except:
                print row
                raise


class MS1ResultsFacade(object):
    renderer_id = "MS1ResultsFacade"
    '''
    A facade that behaves the same way MS1ResultsFile does, but uses a
    HypothesisSampleMatch and its target Hypothesis to construct input to produce
    new theoretical MS2 glycopeptides.
    '''
    def __init__(self, conn, hypothesis_id, hypothesis_sample_match_id):
        self.manager = conn
        session = self.manager.session()
        self.hypothesis_id = hypothesis_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.hypothesis_parameters = session.query(Hypothesis).get(hypothesis_id).parameters
        self.monosaccharide_identities = self.hypothesis_parameters["monosaccharide_identities"]

    @staticmethod
    def render(session, obj):
        return MS1GlycopeptideResult.from_peak_group_match(
            session.query(PeakGroupMatch).get(obj))

    @classmethod
    def get_renderer(self):
        def render(session, obj):
            return MS1GlycopeptideResult.from_peak_group_match(
                session.query(PeakGroupMatch).get(obj))
        return render

    def __iter__(self):
        session = self.manager.session()
        with session.no_autoflush:
            for pgm_id in session.query(PeakGroupMatch.id).filter(
                    PeakGroupMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                    PeakGroupMatch.matched):
                yield pgm_id[0]
        session.close()


parse_digest = msdigest_xml_parser.MSDigestParameters.parse
