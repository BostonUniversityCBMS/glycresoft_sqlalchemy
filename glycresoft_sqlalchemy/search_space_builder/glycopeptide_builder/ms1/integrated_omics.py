import csv
import datetime
import logging
import itertools
import functools
import os
from multiprocessing import Pool

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, ExactMS1GlycopeptideHypothesis, Protein,
    func, TheoreticalGlycopeptideComposition)

from glycresoft_sqlalchemy.data_model import (
    InformedPeptide, InformedTheoreticalGlycopeptideComposition,
    )
from glycresoft_sqlalchemy.data_model import PipelineModule

from glycresoft_sqlalchemy.proteomics.enrich_peptides import EnrichDistinctPeptides
from glycresoft_sqlalchemy.utils.worker_utils import async_worker_pool
from glycresoft_sqlalchemy.utils.database_utils import toggle_indices

from .include_glycomics import MS1GlycanImporter, MS1GlycanImportManager
from .include_proteomics import ProteomeImporter

from ..glycan_utilities import GlycanCombinationProvider
from ..utils import flatten


from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.structure import modification
from glycresoft_sqlalchemy.structure import composition


logger = logging.getLogger("integrated_omics-ms1-search-space")

Sequence = sequence.Sequence
Modification = modification.Modification
NGlycanCoreGlycosylation = modification.NGlycanCoreGlycosylation
ModificationTable = modification.ModificationTable
Composition = composition.Composition

water = Composition("H2O").mass


def make_name(mzid_path, glycan_path):
    """Construct an hypothesis name tag from mzid and glycan file names

    Parameters
    ----------
    mzid_path : str
    glycan_path : str

    Returns
    -------
    str
    """
    mzid_part = os.path.splitext(os.path.basename(mzid_path))[0]
    glycan_part = os.path.splitext(os.path.basename(glycan_path))[0]
    tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
    return "integrated_omics-{}-{}-{}".format(mzid_part, glycan_part, tag)


def glycosylate_callback(peptide, glycan_combinator, position_selector, max_sites=2):
    """
    Generate all glycoform combinations of `peptide` with up to `max_sites` glycosylations

    Parameters
    ----------
    peptide : :class:`InformedPeptide`
        The peptide sequence to be glycosylated
    session : :class:`sqlalchemy.Session`
        An active database connection session to read :class:`TheoreticalGlycanCombination`
        with.
    hypothesis_id : int
        The hypothesis to read :class:`TheoreticalGlycanCombination` from.
    position_selector : function
        Function to call to use select glycosylation site combinations.
    max_sites : int, optional
        The maximum number of glycosylation sites to occupy.

    Returns
    -------
    result: list of :class:`PeptideSequence`
        Glycosylated product sequences
    """
    glycosylation_sites = peptide.glycosylation_sites
    n_sites = len(glycosylation_sites)
    modified_peptide_sequence = peptide.modified_peptide_sequence
    result = []
    for glycan_count in range(1, min(n_sites + 1, max_sites + 1)):
        for glycans in glycan_combinator(glycan_count):
            if len(result) > 10000:
                for item in result:
                    yield item
                result = []
            for sites in position_selector(glycosylation_sites, glycan_count):
                target = Sequence(modified_peptide_sequence)
                for site in sites:
                    for mod in target[site][1]:
                        target.drop_modification(site, mod)
                mass = target.mass
                for site in sites:
                    hexnac = Modification("HexNAc")
                    target.add_modification(site, hexnac)
                mass += glycans.dehydrated_mass()
                result.append((target, mass, glycans))

                # glycan_iter = iter(glycans)
                # for site in sites:
                #     glycan = glycan_iter.next()
                #     hexnac = Modification("HexNAc")
                #     hexnac.mass = glycan.calculated_mass - water
                #     for mod in target[site][1]:
                #         target.drop_modification(site, mod)
                #     target.add_modification(site, hexnac)
                # glycan_composition_string = glycans.composition
                # target.glycan = glycan_composition_string
                # result.append((target, glycans))

    # return result
    for item in result:
        yield item


def extract_peptides(session, ids):
    ids = [i[0] for i in ids]
    total = len(ids)
    last = 0
    step = 20
    results = []
    while last < total:
        results.extend(session.query(InformedPeptide).filter(
            InformedPeptide.id.in_(ids[last:last + step])))
        last += step
    return results


def backprop_informed_attributes(session, peptide, hypothesis_id):
    try:
        ids = session.query(TheoreticalGlycopeptideComposition.id).join(
            Protein, TheoreticalGlycopeptideComposition.protein_id == Protein.id).filter(
            TheoreticalGlycopeptideComposition.modified_peptide_sequence == peptide.modified_peptide_sequence,
            Protein.id == peptide.protein_id,
            TheoreticalGlycopeptideComposition.protein_id == peptide.protein_id,
            Hypothesis.id == hypothesis_id
            ).all()
        peptide_score = peptide.peptide_score
        peptide_id = peptide.id
        other = peptide.other
        args = [{"id": i[0], "other": other, "peptide_score": peptide_score,
                 "base_peptide_id": peptide_id} for i in ids]
        session.execute(InformedTheoreticalGlycopeptideComposition.__table__.insert(), args)
    except Exception, e:
        logging.exception("An error ocurred in backprop_informed_attributes for %r, %r", peptide, hypothesis_id,
                          exc_info=e)
        raise e


def batch_make_theoretical_glycopeptides(peptide_ids, position_selector, database_manager,
                                         hypothesis_id, glycan_combinator, max_sites=2):
    try:
        session = database_manager.session()
        glycopeptide_acc = []
        i = 0

        peptides = extract_peptides(session, peptide_ids)

        for peptide in peptides:

            protein_id = peptide.protein_id
            base_peptide_sequence = peptide.base_peptide_sequence
            modified_peptide_sequence = peptide.modified_peptide_sequence
            start_position = peptide.start_position
            end_position = peptide.end_position
            count_glycosylation_sites = peptide.count_glycosylation_sites
            count_missed_cleavages = peptide.count_missed_cleavages

            glycoforms = glycosylate_callback(peptide, glycan_combinator, position_selector, max_sites=max_sites)
            for glycoform, total_mass, glycans in glycoforms:
                glycan_ids = glycans.id
                informed_glycopeptide = dict(
                    protein_id=protein_id,
                    base_peptide_sequence=base_peptide_sequence,
                    modified_peptide_sequence=modified_peptide_sequence,
                    glycopeptide_sequence=str(glycoform) + (glycans.composition),
                    calculated_mass=total_mass,
                    glycosylation_sites=list(glycoform.n_glycan_sequon_sites),
                    start_position=start_position,
                    end_position=end_position,
                    count_glycosylation_sites=count_glycosylation_sites,
                    count_missed_cleavages=count_missed_cleavages,
                    glycan_mass=glycans.calculated_mass,
                    glycan_composition_str=glycans.composition,
                    glycan_combination_id=glycan_ids,
                    sequence_type="InformedTheoreticalGlycopeptideComposition"
                    )

                glycopeptide_acc.append(informed_glycopeptide)
                i += 1
                if i % 50000 == 0:
                    session.bulk_insert_mappings(TheoreticalGlycopeptideComposition, glycopeptide_acc)
                    session.commit()
                    glycopeptide_acc = []
                    i = 0

            session.bulk_insert_mappings(TheoreticalGlycopeptideComposition, glycopeptide_acc)
            session.commit()
            glycopeptide_acc = []
            backprop_informed_attributes(session, peptide, hypothesis_id)
            session.commit()

        i = 0
        return len(peptide_ids)
    except Exception, e:
        logging.exception("An error occurred in make_theoretical_glycopeptide", exc_info=e)

        raise


def load_proteomics(database_path, mzid_path, hypothesis_id=None, hypothesis_type=ExactMS1GlycopeptideHypothesis,
                    include_all_baseline=True, target_proteins=None):
    if target_proteins is None:
        target_proteins = []
    task = ProteomeImporter(
        database_path, mzid_path, hypothesis_id=hypothesis_id, hypothesis_type=hypothesis_type,
        include_all_baseline=include_all_baseline, target_proteins=target_proteins)
    task.start()
    return task.hypothesis_id


def load_glycomics_naive(database_path, glycan_path, hypothesis_id, format='txt', **kwargs):
    task = MS1GlycanImporter(database_path, glycan_path, hypothesis_id, format, **kwargs)
    task.start()
    return task.hypothesis_id


class IntegratedOmicsMS1SearchSpaceBuilder(PipelineModule, MS1GlycanImportManager):
    hypothesis_type = ExactMS1GlycopeptideHypothesis

    def __init__(self, database_path, hypothesis_id=None, protein_ids=None, mzid_path=None,
                 glycomics_path=None, glycomics_format='txt', hypothesis_name=None, maximum_glycosylation_sites=2,
                 n_processes=4, **kwargs):
        self.manager = self.manager_type(database_path)
        self.manager.initialize()

        MS1GlycanImportManager.__init__(
            self, glycomics_path, glycomics_format, None, maximum_glycosylation_sites,
            **kwargs)

        self.hypothesis_id = hypothesis_id
        self.mzid_path = mzid_path

        self.protein_ids = protein_ids
        self.n_processes = n_processes
        self.hypothesis_name = hypothesis_name
        self.include_all_baseline = kwargs.get("include_all_baseline", True)
        self.options = kwargs

    def bootstrap_hypothesis(self):
        """
        Perform any required initialization steps to prepare the :class:`Hypothesis`
        instance this `PipelineModule` will populate.

        This step may create :class:`Protein`, :class:`TheoreticalGlycanComposition` and
        :class:`TheoreticalGlycanCombination` entities in addition to a :class:`Hypothesis`
        if one does not already exist for the given :attr:`hypothesis_id`.

        Subtasks
            1. :class:`ProteomeImporter`
            2. :class:`MS1GlycanImporter`

        This step may also update :attr:`protein_ids` if they were not
        provided, or provided as a list of names instead of id numbers.
        """
        session = self.manager.session()

        hypothesis = None
        if self.hypothesis_id is not None:
            hypothesis = session.query(Hypothesis).get(self.hypothesis_id)

        if hypothesis is None:
            if self.hypothesis_name is None:
                self.hypothesis_name = make_name(self.mzid_path, self.glycomics_path)
            hypothesis = self.hypothesis_type(name=self.hypothesis_name)
            session.add(hypothesis)
            session.commit()
            self.hypothesis_id = hypothesis.id
            load_proteomics(
                self.database_path, self.mzid_path,
                hypothesis_id=hypothesis.id,
                hypothesis_type=self.hypothesis_type,
                include_all_baseline=self.include_all_baseline,
                target_proteins=self.protein_ids)
            self.import_glycans()

        if self.protein_ids is None:
            protein_ids = flatten(session.query(Protein.id).filter(Protein.hypothesis_id == self.hypothesis_id))
        elif isinstance(self.protein_ids[0], basestring):
            def tryfind(name):
                result = session.query(Protein.id).filter(
                    Protein.name == name,
                    Protein.hypothesis_id == self.hypothesis_id).first()
                if result is None:
                    logger.exception("No protein could be found for name %s", name)
                    return None
                else:
                    return result[0]
            protein_ids = [id for id in [tryfind(name) for name in self.protein_ids] if id is not None]
        else:
            protein_ids = self.protein_ids
        session.close()

        self.protein_ids = protein_ids

    def stream_peptides(self, chunk_size=50):
        session = self.manager.session()
        protein_ids = self.protein_ids
        try:

            for protein_id in protein_ids:
                protein = session.query(Protein).get(protein_id)
                logger.info("Streaming Protein %s", protein.name)
                all_ids = session.query(InformedPeptide.id).filter(
                        InformedPeptide.protein_id == protein_id,
                        InformedPeptide.count_glycosylation_sites > 0).all()
                total = len(all_ids)
                last = 0
                while last <= total:
                    yield all_ids[last:(last + chunk_size)]
                    last += chunk_size

        finally:
            session.commit()
            session.close()

    def prepare_task_fn(self):
        task_fn = functools.partial(
            batch_make_theoretical_glycopeptides,
            hypothesis_id=self.hypothesis_id,
            position_selector=itertools.combinations,
            database_manager=self.manager,
            max_sites=self.maximum_glycosylation_sites,
            glycan_combinator=self.glycan_combinator)
        return task_fn

    def _remove_duplicates(self, session):
        ids = session.query(
            func.min(
                TheoreticalGlycopeptideComposition.id)).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id).group_by(
            TheoreticalGlycopeptideComposition.glycopeptide_sequence,
            TheoreticalGlycopeptideComposition.start_position,
            TheoreticalGlycopeptideComposition.protein_id
            )

        q = session.query(TheoreticalGlycopeptideComposition.id).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id,
            ~TheoreticalGlycopeptideComposition.id.in_(ids.correlate(None)))

        # Delete the base class, and let the cascade through the foreign key delete
        # the subclass.
        logger.info("Removing Duplicates")
        conn = session.connection()
        conn.execute(TheoreticalGlycopeptideComposition.__table__.delete(
            TheoreticalGlycopeptideComposition.__table__.c.id.in_(q.selectable)))
        session.commit()

        logger.info("%d TheoreticalGlycopeptideCompositions created", session.query(
            TheoreticalGlycopeptideComposition).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id).count())

        session.close()

    def run(self):
        self.bootstrap_hypothesis()
        subjob = EnrichDistinctPeptides(self.database_path, self.hypothesis_id, self.protein_ids, max_distance=0)
        subjob.start()

        temporary_glycan_manager = self.mirror_glycans_to_temporary_storage()
        self.glycan_combinator = GlycanCombinationProvider(temporary_glycan_manager, self.hypothesis_id)

        task_fn = self.prepare_task_fn()
        index_controller = toggle_indices(self.manager.session(), TheoreticalGlycopeptideComposition)
        index_controller.drop()
        cntr = 0
        if self.n_processes > 1:
            pool = Pool(self.n_processes)
            async_worker_pool(pool, self.stream_peptides(), task_fn)
            pool.terminate()
        else:
            for peptide in self.stream_peptides():
                cntr += task_fn(peptide)
                logger.info("Completed %d sequences", cntr)
        logger.info("Checking Integrity")
        session = self.manager.session()

        self._remove_duplicates(session)

        index_controller.create()
        temporary_glycan_manager.clear()

        return self.hypothesis_id


class IntegratedOmicsMS1LegacyCSV(PipelineModule):

    def __init__(self, database_path, hypothesis_id, protein_ids=None, output_path=None):
        self.manager = self.manager_type(database_path)

        self.hypothesis_id = hypothesis_id
        session = self.manager()
        if protein_ids is None:
            protein_ids = flatten(session.query(Protein.id).filter(Protein.hypothesis_id == hypothesis_id))
        elif isinstance(protein_ids[0], basestring):
            protein_ids = flatten(session.query(Protein.id).filter(Protein.name == name).first()
                                  for name in protein_ids)
        logger.info("Protein IDs %r", protein_ids)

        self.protein_ids = protein_ids
        if output_path is None:
            output_path = os.path.splitext(database_path)[0] + '.informed_glycopeptides.csv'
        self.output_path = output_path
        hypothesis = session.query(Hypothesis).filter(Hypothesis.id == self.hypothesis_id).first()
        self.monosaccharide_identities = hypothesis.parameters['monosaccharide_identities']

        session.close()

    def stream_glycopeptides(self):
        session = self.manager.session()
        protein_ids = self.protein_ids
        try:
            for protein_id in protein_ids:
                logger.info("Streaming Protein ID %d", protein_id)
                for informed_glycopeptide in session.query(InformedTheoreticalGlycopeptideComposition).filter(
                        InformedTheoreticalGlycopeptideComposition.protein_id == protein_id).group_by(
                        InformedTheoreticalGlycopeptideComposition.modified_peptide_sequence,
                        InformedTheoreticalGlycopeptideComposition.glycan_composition_str,
                        InformedTheoreticalGlycopeptideComposition.start_position):
                    session.expunge(informed_glycopeptide)
                    yield informed_glycopeptide
        finally:
            session.commit()
            session.close()

    def run(self):
        column_headers = self.compute_column_headers()
        with open(self.output_path, 'wb') as handle:
            writer = csv.writer(handle)
            writer.writerow(column_headers)
            g = itertools.imap(glycopeptide_to_columns, self.stream_glycopeptides())
            writer.writerows(g)
        return self.output_path

    def compute_column_headers(self):
        columns = ["Molecular Weight", "C", "Compositions"] + self.monosaccharide_identities +\
              ["Adduct/Replacement", "Adduct Amount", "Peptide Sequence", "Peptide Modification",
               "Peptide Missed Cleavage Number", "Number of Glycan Attachment to Peptide", "Start AA",
               "End AA", "ProteinID"]
        return columns


def glycopeptide_to_columns(glycopeptide):
    r = [glycopeptide.calculated_mass, 0, glycopeptide.glycan_composition_str] +\
         glycopeptide.glycan_composition_str[1:-1].split(";") +\
        ["/0", 0,
         glycopeptide.modified_peptide_sequence, "", 0, glycopeptide.glycopeptide_sequence.count("HexNAc"),
         glycopeptide.start_position, glycopeptide.end_position, glycopeptide.protein_id]
    return r
