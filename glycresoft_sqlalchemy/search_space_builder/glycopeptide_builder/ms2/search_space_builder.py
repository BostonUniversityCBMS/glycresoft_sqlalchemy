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
from glycresoft_sqlalchemy.structure import constants
from glycresoft_sqlalchemy.proteomics import get_enzyme, msdigest_xml_parser

from ..peptide_utilities import SiteListFastaFileParser
from ..utils import fragments

from glycresoft_sqlalchemy.utils import collectiontools

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, Hypothesis, MS2GlycopeptideHypothesis,
    HypothesisSampleMatch, PeakGroupMatchType, Protein,
    TheoreticalGlycopeptide, Hierarchy, MS1GlycopeptideHypothesisSampleMatch)


from glypy.composition.glycan_composition import FrozenGlycanComposition


logger = logging.getLogger("search_space_builder")
mod_pattern = re.compile(r'(\d+)([^\|]+)')
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
            protein_name=theoretical.protein.name, composition_id=theoretical.id,
            glycan_combination_id=theoretical.glycan_combination_id)

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
                 composition_id=None,
                 glycan_combination_id=None):
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
        self.glycan_combination_id = glycan_combination_id

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
        return "MS1GlycopeptideResult(%s)" % str((
            self.most_detailed_sequence, self.peptide_modifications, self.glycan_composition_str))


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
    if peptide_mod_str is None:
        return []
    items = mod_pattern.findall(peptide_mod_str)
    mod_list = []
    for i in items:
        if i[1] == '':
            continue
        for part in i[1].split(":"):
            mod = modification_table.get_modification(part, -1, int(i[0]))
            mod_list.append(mod)
    return mod_list


def extract_glycosites(ms1_result, site_list_map):
    glycan_sites = set(site_list_map.get(ms1_result.protein_name, [])).intersection(
            range(ms1_result.start_position, ms1_result.end_position))

    # No recorded sites, skip this component.
    if len(glycan_sites) == 0:
        return None

    # Adjust the glycan_sites to relative position
    glycan_sites = [x - ms1_result.start_position for x in glycan_sites]
    return glycan_sites


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
        ms1_result.glycopeptide_sequence,
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
    seq.glycan = FrozenGlycanComposition.parse(ms1_result.glycan_composition_str)
    seq_mod = seq.get_sequence(include_glycan=False)
    (oxonium_ions, b_ions, y_ions,
     b_ions_hexnac, y_ions_hexnac,
     stub_ions) = fragments(seq)

    # theoretical_glycopeptide = TheoreticalGlycopeptide(
    theoretical_glycopeptide = dict(
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
        base_composition_id=ms1_result.composition_id,
        glycan_combination_id=ms1_result.glycan_combination_id
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

        glycan_sites = extract_glycosites(ms1_result, site_list_map)
        if glycan_sites is None:
            return 0

        ss = get_search_space(
            ms1_result, glycan_sites, mod_list)
        seq_list = ss.get_theoretical_sequence(ms1_result.count_glycosylation_sites)
        fragments = [generate_fragments(seq, ms1_result)
                     for seq in seq_list]
        i = 0
        for sequence in fragments:
            sequence = TheoreticalGlycopeptide(**sequence)
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
    _do_remove_duplicates = True

    @classmethod
    def from_hypothesis_sample_match(cls, database_path, hypothesis_sample_match_id, n_processes=4, **kwargs):
        dbm = cls.manager_type(database_path)
        session = dbm.session()
        source_hypothesis_sample_match = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
        source_hypothesis = source_hypothesis_sample_match.target_hypothesis
        hypothesis_id = source_hypothesis.id
        variable_modifications = source_hypothesis.parameters.get("variable_modifications", [])
        constant_modifications = source_hypothesis.parameters.get("constant_modifications", [])
        enzyme = source_hypothesis.parameters.get("enzyme", [])

        builder = constructs[source_hypothesis_sample_match.__class__]

        name = Hypothesis.make_unique_name(
            session,
            kwargs.get("hypothesis_name", source_hypothesis.name + " (MS2)"))

        inst = builder(
            database_path, database_path, enzyme=[enzyme], site_list=None,
            constant_modifications=constant_modifications,
            variable_modifications=variable_modifications,
            n_processes=n_processes, ms1_format=MS1ResultsFacade,
            source_hypothesis_id=hypothesis_id,
            source_hypothesis_sample_match_id=hypothesis_sample_match_id,
            hypothesis_name=name)
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

        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications

        self._site_list_input = site_list

        self.create_restricted_modification_table(constant_modifications, variable_modifications)
        self.n_processes = n_processes

        self.enzyme = map(get_enzyme, enzyme)

        self.options = kwargs

    def bootstrap_hypothesis(self):
        kwargs = self.options
        name = kwargs.get("hypothesis_name")
        if name is None:
            tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
            name = u"target-{}".format(tag)
        self.hypothesis = self.HypothesisType(name=name)
        self.session.add(self.hypothesis)
        self.session.commit()
        self.hypothesis_id = self.hypothesis.id

    def bootstrap(self):
        self.bootstrap_hypothesis()
        self.configure_input_source()
        self.save_parameters()

    def save_parameters(self):
        kwargs = self.options
        self.hypothesis.parameters = {
            "monosaccharide_identities": self.monosaccharide_identities,
            "site_list_map": self.glycosylation_site_map,
            "constant_modification_list": self.constant_modifications,
            "variable_modification_list": self.variable_modifications,
            "ms1_output_file": self.ms1_results_file,
            "enzyme": self.enzyme,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS,
            "source_hypothesis_id": kwargs.get("source_hypothesis_id"),
            "source_hypothesis_sample_match_id": kwargs.get("source_hypothesis_sample_match_id")
        }
        self.session.add(self.hypothesis)
        self.session.commit()

    def configure_input_source(self):
        kwargs = self.options
        site_list = self._site_list_input
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

        self.monosaccharide_identities = self.ms1_results_reader.monosaccharide_identities

    def load_protein_from_sitelist(self, glycosylation_site_map):
        for name, sites in self.glycosylation_site_map.items():
            self.session.add(Protein(name=name, glycosylation_sites=sites, hypothesis_id=self.hypothesis.id))
        self.session.commit()

    def load_protein_from_hypothesis(self, hypothesis_id):
        session = self.session
        site_list_map = {}
        hid = self.hypothesis.id
        source_count = 0
        for protein in session.query(Protein).filter(Protein.hypothesis_id == hypothesis_id).group_by(
                Protein.id).all():
            source_count += 1
            prot = (Protein(
                name=protein.name,
                protein_sequence=protein.protein_sequence,
                glycosylation_sites=protein.glycosylation_sites,
                hypothesis_id=hid))
            session.add(prot)
            session.flush()
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
        assert self.hypothesis_id is not None
        task_fn = functools.partial(process_predicted_ms1_ion, modification_table=self.modification_table,
                                    site_list_map=self.glycosylation_site_map,
                                    renderer=self.ms1_format,
                                    database_manager=self.manager, proteins=self.hypothesis.proteins)
        return task_fn

    def stream_results(self):
        for obj in self.ms1_results_reader:
            yield obj

    def _remove_duplicates(self):
        if not self._do_remove_duplicates:
            return
        session = self.manager()
        hypothesis_id = self.hypothesis.id

        ids = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.from_hypothesis(hypothesis_id)).group_by(
            TheoreticalGlycopeptide.glycopeptide_sequence, TheoreticalGlycopeptide.protein_id)

        q = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == hypothesis_id,
            ~TheoreticalGlycopeptide.id.in_(ids.correlate(None)))

        conn = session.connection()
        conn.execute(TheoreticalGlycopeptide.__table__.delete(
            TheoreticalGlycopeptide.__table__.c.id.in_(q.selectable)))
        session.commit()

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes
        '''
        self.bootstrap()
        print self.glycosylation_site_map

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

        logger.info("Checking Integrity")
        self._remove_duplicates()

        final_count = self.session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == id).count()

        logger.info("%d Theoretical Glycopeptides created", final_count)
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
            session.query(PeakGroupMatchType).get(obj))

    @classmethod
    def get_renderer(self):
        def render(session, obj):
            return MS1GlycopeptideResult.from_peak_group_match(
                session.query(PeakGroupMatchType).get(obj))
        return render

    def __iter__(self):
        session = self.manager.session()
        with session.no_autoflush:
            for pgm_id in session.query(PeakGroupMatchType.id).filter(
                    PeakGroupMatchType.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                    PeakGroupMatchType.matched):
                yield pgm_id[0]
        session.close()


parse_digest = msdigest_xml_parser.MSDigestParameters.parse


@constructs.references(MS1GlycopeptideHypothesisSampleMatch)
class BatchingTheoreticalSearchSpaceBuilder(TheoreticalSearchSpaceBuilder):

    def stream_results(self, batch_size=1000):
        for chunk in collectiontools.chunk_iterator2(
                super(BatchingTheoreticalSearchSpaceBuilder, self).stream_results(), batch_size):
            yield chunk

    def prepare_task_fn(self):
        task_fn = functools.partial(batch_process_predicted_ms1_ion, modification_table=self.modification_table,
                                    site_list_map=self.glycosylation_site_map,
                                    renderer=self.ms1_format,
                                    database_manager=self.manager, proteins=self.hypothesis.proteins)
        return task_fn

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes
        '''
        self.bootstrap()
        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0
        step = 1000
        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in worker_pool.imap_unordered(task_fn, self.stream_results(4000), chunksize=1):
                cntr += res
                if (cntr > last + step):
                    last = cntr
                    logger.info("Committing, %d records made", cntr)

            worker_pool.terminate()
        else:
            logger.debug("Building theoretical sequences sequentially")
            for row in self.ms1_results_reader:
                res = task_fn(row)
                cntr += res
                if (cntr > last + step):
                    last = cntr
                    logger.info("Committing, %d records made", cntr)

        logger.info("Checking Integrity")
        self._remove_duplicates()

        final_count = self.session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id).count()

        logger.info("%d Theoretical Glycopeptides created", final_count)

        return self.hypothesis_id


def batch_process_predicted_ms1_ion(rows, modification_table, site_list_map,
                                    renderer, database_manager, proteins):
    """Multiprocessing dispatch function to generate all theoretical sequences and their
    respective fragments from a given MS1 result

    Parameters
    ----------
    rows: iterable
        Iterable of mappings from an MS1 results
    modification_table: :class:`.modifications.RestrictedModificationTable`
        Modification table limited to only the rules specified by the user
    site_list: list of int
        List of putative glycosylation sites from the parent protein
    monosaccharide_identities: list
        List of glycan or monosaccaride names

    Returns
    -------
    int:
        Count of new TheoreticalGlycopeptide items
    """
    session = database_manager.session()
    glycopeptide_acc = []
    try:
        i = 0

        for ms1_result in [renderer.render(session, row) for row in rows]:

            if (ms1_result.base_peptide_sequence == '') or (ms1_result.count_glycosylation_sites == 0):
                continue

            # Compute the set of modifications that can occur.
            mod_list = get_peptide_modifications(
                ms1_result.peptide_modifications, modification_table)

            glycan_sites = extract_glycosites(ms1_result, site_list_map)
            if glycan_sites is None:
                continue

            ss = get_search_space(
                ms1_result, glycan_sites, mod_list)
            seq_list = ss.get_theoretical_sequence(ms1_result.count_glycosylation_sites)
            fragments = [generate_fragments(seq, ms1_result)
                         for seq in seq_list]
            for sequence in fragments:
                sequence["protein_id"] = proteins[ms1_result.protein_name].id
                glycopeptide_acc.append(sequence)
                i += 1
        session.bulk_insert_mappings(TheoreticalGlycopeptide, glycopeptide_acc)
        session.commit()
        return i
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise
