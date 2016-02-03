import os
import csv
import logging
import functools
import itertools
import multiprocessing

from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, MS1GlycopeptideHypothesis, Protein, NaivePeptide, make_transient,
    TheoreticalGlycopeptideComposition, PipelineModule, PeptideBase, func,
    TheoreticalGlycanCombination, slurp)

from .include_glycomics import MS1GlycanImporter, MS1GlycanImportManager

from ..peptide_utilities import generate_peptidoforms, ProteinFastaFileParser, SiteListFastaFileParser
from ..glycan_utilities import get_glycan_combinations, GlycanCombinationProvider
from ..utils import flatten

from glycresoft_sqlalchemy.structure import composition
from glycresoft_sqlalchemy.utils.worker_utils import async_worker_pool
from glycresoft_sqlalchemy.utils.database_utils import toggle_indices

Composition = composition.Composition

logger = logging.getLogger("naive_glycopeptide_hypothesis")


def generate_glycopeptide_compositions(peptide, database_manager, hypothesis_id, max_sites=2):
    try:
        session = database_manager.session()
        peptide = session.query(NaivePeptide).get(peptide[0])
        i = 0

        glycopeptide_acc = []

        for glycan_set in get_glycan_combinations(
                database_manager.session(), min(peptide.count_glycosylation_sites, max_sites), hypothesis_id):

            glycan_composition_str = glycan_set.composition
            glycan_mass = glycan_set.dehydrated_mass()
            glycoform = TheoreticalGlycopeptideComposition(
                base_peptide_sequence=peptide.base_peptide_sequence,
                modified_peptide_sequence=peptide.modified_peptide_sequence,
                glycopeptide_sequence=peptide.modified_peptide_sequence + glycan_composition_str,
                protein_id=peptide.protein_id,
                start_position=peptide.start_position,
                end_position=peptide.end_position,
                peptide_modifications=peptide.peptide_modifications,
                count_missed_cleavages=peptide.count_missed_cleavages,
                count_glycosylation_sites=(glycan_set.count),
                glycosylation_sites=peptide.glycosylation_sites,
                sequence_length=peptide.sequence_length,
                calculated_mass=peptide.calculated_mass + glycan_mass,
                glycan_composition_str=glycan_composition_str,
                glycan_mass=glycan_mass,
                glycan_combination_id=glycan_set.id
            )

            glycopeptide_acc.append(glycoform)

            i += 1
            if i % 5000 == 0:
                logger.info("Flushing %d for %r-%r", i, peptide, peptide.protein)
                session.bulk_save_objects(glycopeptide_acc)
                session.commit()
                glycopeptide_acc = []

        session.bulk_save_objects(glycopeptide_acc)
        session.commit()
        session.close()
        return i
    except Exception, e:
        logger.exception("%r", locals(), exc_info=e)
        raise e
    finally:
        session.close()


def batch_generate_glycopeptide_compositions(peptide_ids, database_manager, hypothesis_id,
                                             glycan_combinator, max_sites=2):
    try:
        session = database_manager.session()
        i = 0
        glycopeptide_acc = []
        for peptide in slurp(session, NaivePeptide, peptide_ids, flatten=True):
            base_peptide_sequence = peptide.base_peptide_sequence
            modified_peptide_sequence = peptide.modified_peptide_sequence
            protein_id = peptide.protein_id
            start_position = peptide.start_position
            end_position = peptide.end_position
            peptide_modifications = peptide.peptide_modifications
            count_missed_cleavages = peptide.count_missed_cleavages
            glycosylation_sites = peptide.glycosylation_sites
            sequence_length = peptide.sequence_length
            calculated_mass = peptide.calculated_mass

            for glycan_set in glycan_combinator(min(peptide.count_glycosylation_sites, max_sites)):

                glycan_composition_str = glycan_set.composition
                glycan_mass = glycan_set.dehydrated_mass()
                glycoform = dict(
                    base_peptide_sequence=base_peptide_sequence,
                    modified_peptide_sequence=modified_peptide_sequence,
                    glycopeptide_sequence=modified_peptide_sequence + glycan_composition_str,
                    protein_id=protein_id,
                    start_position=start_position,
                    end_position=end_position,
                    peptide_modifications=peptide_modifications,
                    count_missed_cleavages=count_missed_cleavages,
                    count_glycosylation_sites=glycan_set.count,
                    glycosylation_sites=glycosylation_sites,
                    sequence_length=sequence_length,
                    calculated_mass=calculated_mass + glycan_mass,
                    glycan_composition_str=glycan_composition_str,
                    glycan_mass=glycan_mass,
                    glycan_combination_id=glycan_set.id
                )

                glycopeptide_acc.append(glycoform)

                i += 1
                if (i % 5000) == 0:
                    logging.info("Thunk, %d", i)
                    session.bulk_insert_mappings(TheoreticalGlycopeptideComposition, glycopeptide_acc)
                    # session.bulk_save_objects(glycopeptide_acc)
                    session.commit()
                    glycopeptide_acc = []

        session.bulk_insert_mappings(TheoreticalGlycopeptideComposition, glycopeptide_acc)
        # session.bulk_save_objects(glycopeptide_acc)

        session.commit()
        session.close()
        return i
    except Exception, e:
        logger.exception("%r", locals(), exc_info=e)
        raise e
    finally:
        session.close()


def digest_protein(protein, manager, constant_modifications, variable_modifications, enzyme, max_missed_cleavages):
    session = manager.session()
    try:
        peptidoforms = []
        i = 0
        logger.info("Digesting %r, %d glycosites", protein, len(protein.glycosylation_sites))
        for peptidoform in generate_peptidoforms(protein, constant_modifications,
                                                 variable_modifications, enzyme,
                                                 max_missed_cleavages):
            peptidoforms.append(peptidoform)
            i += 1
            if len(peptidoforms) > 10000:
                logger.info("Flushing %d peptidoforms for %r", i, protein)
                session.bulk_save_objects(peptidoforms)
                session.commit()
                session.expire_all()
                peptidoforms = []
        session.bulk_save_objects(peptidoforms)
        session.commit()
        logger.info("Digested %d peptidoforms for %s", i, protein)

        return 1
    except Exception, e:
        logger.exception("An exception occurred in digest_protein", exc_info=e)
        raise
    finally:
        session.close()


class NaiveGlycopeptideHypothesisBuilder(PipelineModule, MS1GlycanImportManager):
    HypothesisType = MS1GlycopeptideHypothesis

    def __init__(self, database_path, hypothesis_name, protein_file, site_list_file,
                 glycomics_path, glycomics_format, constant_modifications, variable_modifications,
                 enzyme, max_missed_cleavages=1, maximum_glycosylation_sites=2, n_processes=4, **kwargs):
        self.manager = self.manager_type(database_path)
        self.manager.initialize()

        MS1GlycanImportManager.__init__(
            self, glycomics_path, glycomics_format, None, maximum_glycosylation_sites,
            **kwargs)
        self.hypothesis_name = hypothesis_name
        self.protein_file = protein_file
        self.site_list_file = site_list_file
        self.glycomics_path = glycomics_path
        self.glycomics_format = glycomics_format
        self.session = self.manager.session()
        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications
        self.enzyme = enzyme
        self.max_missed_cleavages = max_missed_cleavages
        self.maximum_glycosylation_sites = maximum_glycosylation_sites
        self.options = kwargs
        self.n_processes = n_processes

    def bootstrap_hypothesis(self, session=None):
        if session is None:
            session = self.manager.session()

        hypothesis = None
        hypothesis = self.hypothesis = self.HypothesisType(name=self.hypothesis_name)
        hypothesis.parameters = ({
            "protein_file": self.protein_file,
            "site_list_file": self.site_list_file,
            "glycomics_path": self.glycomics_path,
            "glycomics_format": self.glycomics_format,
            "constant_modifications": self.constant_modifications,
            "variable_modifications": self.variable_modifications,
            "enzyme": self.enzyme,
            "max_missed_cleavages": self.max_missed_cleavages,
            "maximum_glycosylation_sites": self.maximum_glycosylation_sites
        })

        session.add(self.hypothesis)
        session.commit()
        self.hypothesis_id = hypothesis.id

    def prepare_task_fn(self):
        task_fn = functools.partial(
            batch_generate_glycopeptide_compositions,
            hypothesis_id=self.hypothesis_id,
            database_manager=self.manager,
            max_sites=self.maximum_glycosylation_sites,
            glycan_combinator=self.glycan_combinator)
        return task_fn

    def stream_proteins(self):
        return self.manager.session().query(Protein).filter(Protein.hypothesis_id == self.hypothesis.id)

    def batch_stream_peptides(self, chunk_size=20):
        ids = self.manager.session().query(NaivePeptide.id).filter(
                NaivePeptide.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis.id).all()
        total = len(ids)
        last = 0
        while last <= total:
            yield ids[last:(last + chunk_size)]
            last += chunk_size

    def digest_proteins(self, session):
        protein_digest_task = functools.partial(
            digest_protein,
            manager=self.manager,
            constant_modifications=self.constant_modifications,
            variable_modifications=self.variable_modifications,
            enzyme=self.enzyme,
            max_missed_cleavages=self.max_missed_cleavages)
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            async_worker_pool(pool, self.stream_proteins(), protein_digest_task)
            pool.terminate()
        else:
            for protein in self.stream_proteins():
                protein_digest_task(protein)
        session.commit()
        session.expire_all()
        session.query(NaivePeptide).filter(
            NaivePeptide.count_glycosylation_sites == None).delete("fetch")
        session.commit()
        session.expire_all()

        assert session.query(NaivePeptide).filter(
            NaivePeptide.count_glycosylation_sites == None).count() == 0

    def _remove_duplicates(self, session):
        logger.info("Checking integrity")
        ids = session.query(func.min(TheoreticalGlycopeptideComposition.id)).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis.id).group_by(
            TheoreticalGlycopeptideComposition.glycopeptide_sequence,
            TheoreticalGlycopeptideComposition.peptide_modifications,
            TheoreticalGlycopeptideComposition.glycosylation_sites,
            TheoreticalGlycopeptideComposition.protein_id)

        q = session.query(TheoreticalGlycopeptideComposition.id).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis.id,
            ~TheoreticalGlycopeptideComposition.id.in_(ids.correlate(None)))
        conn = session.connection()
        conn.execute(TheoreticalGlycopeptideComposition.__table__.delete(
            TheoreticalGlycopeptideComposition.__table__.c.id.in_(q.selectable)))
        session.commit()

    def run(self):
        session = self.session
        self.bootstrap_hypothesis(session)
        self.import_glycans()

        logger.info("Loaded glycans")

        site_list_map = {}
        if self.site_list_file is not None:
            for protein_entry in SiteListFastaFileParser(self.site_list_file):
                site_list_map[protein_entry['name']] = protein_entry['glycosylation_sites']
        for protein in ProteinFastaFileParser(self.protein_file):
            if protein.name in site_list_map:
                protein.glycosylation_sites = list(set(protein.glycosylation_sites) | site_list_map[protein.name])
            protein.hypothesis_id = self.hypothesis.id
            session.add(protein)
            logger.info("Loaded %s", protein)

        session.commit()

        self.digest_proteins(session)

        temporary_glycan_manager = self.mirror_glycans_to_temporary_storage()
        self.glycan_combinator = GlycanCombinationProvider(temporary_glycan_manager, self.hypothesis_id)

        index_controller = toggle_indices(self.manager.session(), TheoreticalGlycopeptideComposition)
        index_controller.drop()

        work_stream = self.batch_stream_peptides(2)
        task_fn = self.prepare_task_fn()

        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            async_worker_pool(pool, work_stream, task_fn)
            pool.terminate()
        else:
            for item in work_stream:
                cntr += task_fn(item)
                logger.info("%d done", cntr)

        session.commit()

        logger.info("Checking integrity")
        self._remove_duplicates(session)
        logger.info("Naive Hypothesis Complete. %d theoretical glycopeptide compositions genreated.",
                    session.query(TheoreticalGlycopeptideComposition).filter(
                        TheoreticalGlycopeptideComposition.protein_id == Protein.id,
                        Protein.hypothesis_id == self.hypothesis.id).count())

        index_controller.create()
        temporary_glycan_manager.clear()

        return self.hypothesis.id


class NaiveGlycopeptideHypothesisMS1LegacyCSV(PipelineModule):
    HypothesisType = MS1GlycopeptideHypothesis

    def __init__(self, database_path, hypothesis_id, protein_ids=None, output_path=None):
        self.manager = self.manager_type(database_path)

        self.hypothesis_id = hypothesis_id
        session = self.manager.session()
        if protein_ids is None:
            protein_ids = flatten(session.query(Protein.id).filter(Protein.hypothesis_id == hypothesis_id))
        elif isinstance(protein_ids[0], basestring):
            protein_ids = flatten(session.query(Protein.id).filter(Protein.name == name).first()
                                  for name in protein_ids)

        self.protein_ids = protein_ids
        if output_path is None:
            output_path = os.path.splitext(database_path)[0] + '.glycopeptides_compositions.csv'
        self.output_path = output_path
        hypothesis = session.query(self.HypothesisType).filter(self.HypothesisType.id == self.hypothesis_id).first()
        self.monosaccharide_identities = hypothesis.parameters['monosaccharide_identities']

        session.close()

    def stream_glycopeptides(self):
        session = self.manager.session()
        protein_ids = self.protein_ids
        try:
            for protein_id in protein_ids:
                logger.info("Streaming Protein ID %d", protein_id)
                for glycopeptide_composition in session.query(TheoreticalGlycopeptideComposition).filter(
                        TheoreticalGlycopeptideComposition.protein_id == protein_id).yield_per(1000):
                    session.expunge(glycopeptide_composition)
                    yield glycopeptide_composition
        finally:
            session.commit()
            session.close()

    def run(self):
        column_headers = self.compute_column_headers()
        with open(self.output_path, 'wb') as handle:
            writer = csv.writer(handle)
            writer.writerow(column_headers)
            task_fn = functools.partial(glycopeptide_to_columns,
                                        monosaccharide_identities=self.monosaccharide_identities)
            g = itertools.imap(task_fn, self.stream_glycopeptides())
            writer.writerows(g)

    def compute_column_headers(self):
        columns = ["Molecular Weight", "C", "Compositions"] + self.monosaccharide_identities +\
              ["Adduct/Replacement", "Adduct Amount", "Peptide Sequence", "Peptide Modification",
               "Peptide Missed Cleavage Number", "Number of Glycan Attachment to Peptide", "Start AA",
               "End AA", "ProteinID"]
        return columns


def glycopeptide_to_columns(glycopeptide, monosaccharide_identities):
    r = [glycopeptide.calculated_mass, 0, glycopeptide.glycan_composition_str] +\
        [glycopeptide.glycan_composition[m] for m in monosaccharide_identities] +\
        ["/0", 0,
         glycopeptide.modified_peptide_sequence, glycopeptide.peptide_modifications,
         glycopeptide.count_missed_cleavages,
         glycopeptide.count_glycosylation_sites,
         glycopeptide.start_position, glycopeptide.end_position, glycopeptide.protein_id]
    return r
