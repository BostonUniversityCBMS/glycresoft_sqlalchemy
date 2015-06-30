import os
import csv
import logging
import functools
import itertools
import multiprocessing

from ..data_model import (DatabaseManager, Hypothesis, Protein, Glycan, NaivePeptide,
                          TheoreticalGlycopeptideComposition, PipelineModule, PeptideBase)
from .include_glycomics import MS1GlycanImporter
from .peptide_utilities import generate_peptidoforms, ProteinFastaFileParser, SiteListFastaFileParser
from .glycan_utilities import get_glycan_combinations, merge_compositions
from .utils import flatten

from glycresoft_ms2_classification.structure import composition

Composition = composition.Composition

format_mapping = {
    "csv": "hypothesis",
    "hypothesis": None
}

logger = logging.getLogger("naive_glycopeptide_hypothesis")


def generate_glycopeptide_compositions(peptide, database_manager, hypothesis_id):
    try:
        session = database_manager.session()
        peptide = session.query(NaivePeptide).get(peptide[0])
        # logger.info("Working on %r", peptide)
        i = 0
        for glycan_set in get_glycan_combinations(
                database_manager.session(), peptide.count_glycosylation_sites, hypothesis_id):
            glycan_composition_str = merge_compositions([(g.composition) for g in glycan_set])
            glycan_mass = sum([(g.mass) for g in glycan_set]) - (Composition("H2O").mass * len(glycan_set))
            glycoform = TheoreticalGlycopeptideComposition(
                base_peptide_sequence=peptide.base_peptide_sequence,
                protein_id=peptide.protein_id,
                start_position=peptide.start_position,
                end_position=peptide.end_position,
                peptide_modifications=peptide.peptide_modifications,
                count_missed_cleavages=peptide.count_missed_cleavages,
                count_glycosylation_sites=peptide.count_glycosylation_sites,
                glycosylation_sites=peptide.glycosylation_sites,
                sequence_length=peptide.sequence_length,
                calculated_mass=peptide.calculated_mass + glycan_mass,
                glycan_composition_str=glycan_composition_str,
                glycan_mass=glycan_mass,
            )
            session.add(glycoform)
            i += 1
            if i % 100000 == 0:
                session.commit()
        session.commit()
        session.close()
        return i
    except Exception, e:
        logger.exception("%r", locals(), exc_info=e)
        raise e
    finally:
        session.close()


class NaiveGlycopeptideHypothesisBuilder(PipelineModule):
    def __init__(self, database_path, hypothesis_name, protein_file, site_list_file,
                 glycan_file, glycan_file_type, constant_modifications, variable_modifications,
                 enzyme, max_missed_cleavages=1, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.manager.initialize()
        self.protein_file = protein_file
        self.site_list_file = site_list_file
        self.glycan_file = glycan_file
        self.glycan_file_type = glycan_file_type
        self.session = session = self.manager.session()
        self.hypothesis = Hypothesis(name=hypothesis_name)
        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications
        self.enzyme = enzyme
        self.max_missed_cleavages = max_missed_cleavages
        self.hypothesis.parameters = ({
            "protein_file": protein_file,
            "site_list_file": site_list_file,
            "glycan_file": glycan_file,
            "glycan_file_type": glycan_file_type,
            "constant_modifications": constant_modifications,
            "variable_modifications": variable_modifications,
            "enzyme": enzyme,
            "max_missed_cleavages": max_missed_cleavages
        })
        session.add(self.hypothesis)
        session.commit()
        self.n_processes = n_processes

    def prepare_task_fn(self):
        return functools.partial(generate_glycopeptide_compositions,
                                 database_manager=self.manager,
                                 hypothesis_id=self.hypothesis.id)

    def stream_peptides(self):
        for item in self.manager.session().query(NaivePeptide.id).filter(
                NaivePeptide.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis.id):
            yield item

    def run(self):
        session = self.session
        if self.glycan_file_type == "csv":
            importer = MS1GlycanImporter(self.database_path, self.glycan_file, self.hypothesis.id, 'hypothesis')
            importer.run()
        else:
            raise NotImplementedError()
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
        for protein in session.query(Protein).filter(Protein.hypothesis_id == self.hypothesis.id):
            peptidoforms = list(generate_peptidoforms(protein, self.constant_modifications,
                                                      self.variable_modifications, self.enzyme,
                                                      self.max_missed_cleavages))
            session.add_all(peptidoforms)
            logger.info("Digested %s", protein)
            session.commit()

        session.query(PeptideBase).filter(PeptideBase.count_glycosylation_sites == None).delete('fetch')
        session.commit()

        work_stream = self.stream_peptides()
        task_fn = self.prepare_task_fn()

        cntr = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            # apply_async gets around the block-wise batch wait that all of the *map* variants have
            work_queue = []
            for item in work_stream:
                work_queue.append(pool.apply_async(task_fn, [item]))
            last_length = len(work_queue)
            while len(work_queue) > 0:
                next_round = []
                for i in range(len(work_queue)):
                    task = work_queue[i]
                    task.wait(0.1)
                    if task.ready():
                        cntr += task.get()
                    else:
                        next_round.append(task)
                work_queue = next_round
                if last_length != len(work_queue):
                    last_length = len(work_queue)
                    logger.info("%d tasks remaining", last_length)

        else:
            for item in work_stream:
                cntr += task_fn(item)
                logger.info("%d done", cntr)
        logger.info("Naive Hypothesis Complete. %d theoretical glycopeptide compositions genreated.",
                    session.query(TheoreticalGlycopeptideComposition).filter(
                        TheoreticalGlycopeptideComposition.protein_id == Protein.id,
                        Protein.hypothesis_id == self.hypothesis.id).count())
        session.commit()
        return self.hypothesis.id


class NaiveGlycopeptideHypothesiMS1LegacyCSV(PipelineModule):
    manager_type = DatabaseManager

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
            output_path = os.path.splitext(database_path)[0] + '.glycopeptides.csv'
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
                for glycopeptide_composition in session.query(TheoreticalGlycopeptideComposition).filter(
                        TheoreticalGlycopeptideComposition.protein_id == protein_id):
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
            g = itertools.imap(glycopeptide_to_columns, self.stream_glycopeptides())
            writer.writerows(g)

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
