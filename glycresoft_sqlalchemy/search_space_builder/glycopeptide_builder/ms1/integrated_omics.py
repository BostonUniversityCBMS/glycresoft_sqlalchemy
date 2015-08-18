import csv
import datetime
import logging
import itertools
import functools
import os
from multiprocessing import Pool

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, MS1GlycopeptideHypothesis, Protein,
    Glycan, DatabaseManager, TheoreticalGlycopeptideCompositionGlycanAssociation,
    func, TheoreticalGlycopeptideComposition)

from glycresoft_sqlalchemy.data_model.informed_proteomics import (InformedPeptide, InformedTheoreticalGlycopeptideComposition)
from glycresoft_sqlalchemy.data_model import PipelineModule

from glycresoft_sqlalchemy.proteomics.enrich_peptides import EnrichDistinctPeptides
from glycresoft_sqlalchemy.utils.worker_utils import async_worker_pool
from glycresoft_sqlalchemy.utils.database_utils import toggle_indices

from .include_glycomics import MS1GlycanImporter
from .include_proteomics import ProteomeImporter

from ..glycan_utilities import get_glycan_combinations, merge_compositions
from ..utils import flatten


from glycresoft_sqlalchemy.structure import sequence
from glycresoft_sqlalchemy.structure import modification
from glycresoft_sqlalchemy.structure import composition
from glycresoft_sqlalchemy.structure import constants


logger = logging.getLogger(__name__)

Sequence = sequence.Sequence
Modification = modification.Modification
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


def glycosylate_callback(peptide, session, hypothesis_id, position_selector, max_sites=2):
    n_sites = len(peptide.glycosylation_sites)
    result = []
    for glycan_count in range(1, min(n_sites + 1, max_sites)):
        for glycans in get_glycan_combinations(session, glycan_count, hypothesis_id):
            for sites in position_selector(peptide.glycosylation_sites, glycan_count):
                target = Sequence(peptide.modified_peptide_sequence)
                glycan_composition = []
                glycan_iter = iter(glycans)
                for site in sites:
                    glycan = glycan_iter.next()
                    hexnac = Modification("HexNAc")
                    hexnac.mass = glycan.theoretical_mass - water
                    for mod in target[site][1]:
                        target.drop_modification(site, mod)
                    target.add_modification(site, hexnac)
                    glycan_composition.append(glycan.glycan_composition)
                glycan_composition_string = merge_compositions(glycan_composition)
                target.glycan = glycan_composition_string
                result.append((target, [g.id for g in glycans]))
    return result


def make_theoretical_glycopeptide(peptide, position_selector, database_manager, hypothesis_id):
    """Multiprocessing task to create exact glycopeptides from `peptide`. An exact glycopeptide
    has all of its non-glycan modifications positioned and unambiguous.

    Parameters
    ----------
    peptide : InformedPeptide
        The template peptide onto which glycans are added
    position_selector : function
        A function that controls the manner in which combinations or
        permutations of glycans are added to the template
    database_manager : DatabaseManager
        The manager of database connections
    hypothesis_id : int
        The hypothesis to look up information in.

    Returns
    -------
    int: The number  of glycopeptides produced.
    """
    try:
        session = database_manager.session()
        conn = session.connection()
        glycoforms = glycosylate_callback(peptide, session, hypothesis_id, position_selector)
        for glycoform, glycan_ids in glycoforms:
            informed_glycopeptide = InformedTheoreticalGlycopeptideComposition(
                protein_id=peptide.protein_id,
                base_peptide_sequence=peptide.base_peptide_sequence,
                modified_peptide_sequence=peptide.modified_peptide_sequence,
                glycopeptide_sequence=str(glycoform),
                calculated_mass=glycoform.mass,
                glycosylation_sites=list(glycoform.n_glycan_sequon_sites),
                peptide_score=peptide.peptide_score,
                start_position=peptide.start_position,
                end_position=peptide.end_position,
                count_glycosylation_sites=peptide.count_glycosylation_sites,
                other=peptide.other,
                glycan_composition_str=glycoform.glycan,
                base_peptide_id=peptide.id)
            session.add(informed_glycopeptide)
            session.flush()
            conn.execute(
                TheoreticalGlycopeptideCompositionGlycanAssociation.insert(),
                [{'peptide_id': informed_glycopeptide.id, 'glycan_id': gid} for gid in glycan_ids])
        if len(glycoforms) > 0:
            session.commit()
            session.close()
        return len(glycoforms)
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise


def load_proteomics(database_path, mzid_path):
    task = ProteomeImporter(database_path, mzid_path)
    task.start()
    return task.hypothesis_id


def load_glycomics_naive(database_path, glycan_path, hypothesis_id):
    task = MS1GlycanImporter(database_path, glycan_path, hypothesis_id, "hypothesis")
    task.start()
    return task.hypothesis_id


class IntegratedOmicsMS1SearchSpaceBuilder(PipelineModule):
    HypothesisType = MS1GlycopeptideHypothesis

    def __init__(self, database_path, hypothesis_id, protein_ids=None, n_processes=4):
        self.manager = self.manager_type(database_path)
        if not os.path.exists(database_path):
            self.manager.initialize()
            self.manager = self.manager_type(database_path)

        session = self.manager.session()
        hypothesis = session.query(Hypothesis).filter(Hypothesis.id == hypothesis_id).first()
        hypothesis.hypothesis_type = self.HypothesisType.__name__

        self.hypothesis_id = hypothesis.id
        if protein_ids is None:
            protein_ids = flatten(session.query(Protein.id).filter(Protein.hypothesis_id == hypothesis_id))
        elif isinstance(protein_ids[0], basestring):
            protein_ids = flatten(session.query(Protein.id).filter(Protein.name == name).first()
                                  for name in protein_ids)
        session.close()

        self.protein_ids = protein_ids

        self.n_processes = n_processes

    def stream_peptides(self):
        session = self.manager.session()
        protein_ids = self.protein_ids
        try:
            for protein_id in protein_ids:
                logger.info("Streaming Protein ID %d", protein_id)
                for informed_peptide in session.query(InformedPeptide).filter(
                        InformedPeptide.protein_id == protein_id,
                        InformedPeptide.count_glycosylation_sites > 0):
                    session.expunge(informed_peptide)
                    yield informed_peptide
        finally:
            session.commit()
            session.close()

    def prepare_task_fn(self):
        task_fn = functools.partial(
            make_theoretical_glycopeptide,
            hypothesis_id=self.hypothesis_id,
            position_selector=itertools.combinations,
            database_manager=self.manager)
        return task_fn

    def run(self):
        subjob = EnrichDistinctPeptides(self.database_path, self.hypothesis_id, self.protein_ids, max_distance=0)
        subjob.start()

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
        ids = session.query(
            func.min(
                TheoreticalGlycopeptideComposition.id)).filter(
            TheoreticalGlycopeptideComposition.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id).group_by(
            TheoreticalGlycopeptideComposition.glycopeptide_sequence,
            TheoreticalGlycopeptideComposition.start_position,
            TheoreticalGlycopeptideComposition.protein_id)

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
        index_controller.create()
        return self.hypothesis_id


class IntegratedOmicsMS1LegacyCSV(PipelineModule):

    def __init__(self, database_path, hypothesis_id, protein_ids=None, output_path=None):
        self.manager = self.manager_type(database_path)

        self.hypothesis_id = hypothesis_id
        session = self.manager.session()
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
