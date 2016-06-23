import os
import datetime
import logging
try:
    logger = logging.getLogger("include_proteomics")
except:
    pass


from ..peptide_utilities import SiteListFastaFileParser, generate_peptidoforms
from glycresoft_sqlalchemy.structure import sequence, modification

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, Protein, make_transient, InformedPeptide,
    ExactMS1GlycopeptideHypothesis, func, _TemplateNumberStore)

from glycresoft_sqlalchemy.proteomics.mzid_sa import Proteome as MzIdentMLProteome

from glycresoft_sqlalchemy.structure.sequence import find_n_glycosylation_sequons

from glycresoft_sqlalchemy.utils.database_utils import temp_table


def make_base_sequence(peptide, constant_modifications, modification_table):
    seq = peptide.base_peptide_sequence
    seq = sequence.parse(seq)
    for mod in {modification_table[const_mod]
                for const_mod in constant_modifications}:
        for site in mod.find_valid_sites(seq):
            seq.add_modification(site, mod.name)
    make_transient(peptide)
    peptide.id = None
    peptide.modified_peptide_sequence = str(seq)
    peptide.calculated_mass = seq.mass
    return peptide


def find_best_peptides(session, hypothesis_id):
    q = session.query(
        InformedPeptide.id, InformedPeptide.peptide_score,
        InformedPeptide.modified_peptide_sequence, InformedPeptide.protein_id).join(
        Protein).filter(Protein.hypothesis_id == hypothesis_id).yield_per(10000)
    keepers = dict()
    for id, score, modified_peptide_sequence, protein_id in q:
        try:
            old_id, old_score = keepers[modified_peptide_sequence, protein_id]
            if score > old_score:
                keepers[modified_peptide_sequence, protein_id] = id, score
        except KeyError:
            keepers[modified_peptide_sequence, protein_id] = id, score
    return keepers


def store_best_peptides(session, keepers):
    table = temp_table(_TemplateNumberStore)
    conn = session.connection()
    table.create(conn)
    payload = [{"value": x[0]} for x in keepers.values()]
    conn.execute(table.insert(), payload)
    session.commit()
    return table


def remove_duplicates(session, hypothesis_id):
    keepers = find_best_peptides(session, hypothesis_id)
    table = store_best_peptides(session, keepers)
    ids = session.query(table.c.value)
    q = session.query(InformedPeptide.id).filter(
        InformedPeptide.protein_id == Protein.id,
        Protein.hypothesis_id == hypothesis_id,
        ~InformedPeptide.id.in_(ids.correlate(None)))

    session.execute(InformedPeptide.__table__.delete(
        InformedPeptide.__table__.c.id.in_(q.selectable)))
    conn = session.connection()
    table.drop(conn)
    session.commit()


class ProteomeImporter(PipelineModule):

    def __init__(self, database_path, mzid_path, glycosylation_sites_file=None,
                 hypothesis_id=None, constant_modifications=("Carbamidomethyl (C)",),
                 target_proteins=None,
                 hypothesis_type=ExactMS1GlycopeptideHypothesis, peptide_type=InformedPeptide,
                 baseline_missed_cleavages=1, include_all_baseline=False):
        self.manager = self.manager_type(database_path)
        self.mzid_path = mzid_path
        self.hypothesis_id = hypothesis_id
        self.glycosylation_sites_file = glycosylation_sites_file
        self.constant_modifications = constant_modifications
        self.hypothesis_type = hypothesis_type
        self.peptide_type = peptide_type
        self.baseline_missed_cleavages = baseline_missed_cleavages
        self.include_all_baseline = include_all_baseline
        self.target_proteins = target_proteins or {}

    def build_baseline_peptides(self, session, protein):
        if len(self.enzymes):
            if (len(self.target_proteins) == 0) or (protein.id in self.target_proteins) or\
                    (protein.name in self.target_proteins):
                logger.info("Building All Baseline Peptides For %r", protein)
                for peptidoform in generate_peptidoforms(
                        protein, self.constant_modifications, [],
                        self.enzymes[0], missed_cleavages=self.baseline_missed_cleavages,
                        peptide_class=self.peptide_type,
                        **{"count_variable_modifications": 0, "peptide_score": 0}):
                    peptidoform.count_variable_modifications = 0
                    if peptidoform.modified_peptide_sequence is None:
                        peptidoform.modified_peptide_sequence = peptidoform.base_peptide_sequence
                    session.add(peptidoform)
                protein.informed_peptides.filter(
                    self.peptide_type.count_glycosylation_sites == None).delete("fetch")
                session.commit()

    def build_unmodified_peptides(self, session):
        modification_table = modification.RestrictedModificationTable.bootstrap(
            list(self.constant_modifications), [])
        basic_peptides = []
        for peptide in session.query(self.peptide_type).filter(
                self.peptide_type.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis_id).group_by(
                self.peptide_type.base_peptide_sequence):
            peptide = make_base_sequence(peptide, self.constant_modifications, modification_table)
            basic_peptides.append(peptide)
            if len(basic_peptides) > 1000:
                session.add_all(basic_peptides)
                session.commit()
                basic_peptides = []

        session.add_all(basic_peptides)
        session.commit()

    def _display_protein_peptide_counts(self, session):
        peptide_counts = session.query(Protein.name, func.count(InformedPeptide.id)).filter(
            Protein.hypothesis_id == self.hypothesis_id).join(InformedPeptide).group_by(
            InformedPeptide.protein_id, Protein.name).all()

        logger.info("Peptide Counts: %r", peptide_counts)

    def _remove_duplicates(self, session):
        remove_duplicates(session, self.hypothesis_id)
        logger.info("Saved %d InformedPeptides", session.query(InformedPeptide).filter(
            InformedPeptide.from_hypothesis(self.hypothesis_id)).count())

    def _bootstrap_hypothesis(self, session):
        self.manager.initialize()
        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
        name = 'target-{}-{}'.format(
                os.path.splitext(
                    os.path.basename(self.mzid_path)
                    )[0],
                tag)
        hypothesis = self.hypothesis_type(name=name, parameters={"protein_file": self.mzid_path})
        session.add(hypothesis)
        session.commit()
        self.hypothesis_id = hypothesis.id
        return hypothesis

    def load_proteome(self, session):
        mzident_parser = MzIdentMLProteome(self.manager.path, self.mzid_path, self.hypothesis_id)
        self._display_protein_peptide_counts(session)
        self.constant_modifications = mzident_parser.constant_modifications
        self.enzymes = mzident_parser.enzymes

        logger.info("Constant Modifications: %r", self.constant_modifications)
        logger.info("Enzyme: %r", self.enzymes)
        if True:  # self.glycosylation_sites_file is None:
            for protein in session.query(Protein).filter(Protein.hypothesis_id == self.hypothesis_id):
                session.add(protein)
                session.flush()
                if self.include_all_baseline:
                    self.build_baseline_peptides(session, protein)

        # else:
        #     site_list_gen = SiteListFastaFileParser(self.glycosylation_sites_file)
        #     site_list_map = {d['name']: d["glycosylation_sites"] for d in site_list_gen}
        #     for protein in session.query(Protein):
        #         try:
        #             protein.glycosylation_sites = site_list_map[protein.name]
        #         except KeyError:
        #             protein.glycosylation_sites = find_n_glycosylation_sequons(protein.protein_sequence)
        #         session.add(protein)
        #         session.flush()
        #         if self.include_all_baseline:
        #             self.build_baseline_peptides(session, protein)

        if not self.include_all_baseline:
            logger.info("Building Unmodified Reference Peptides")
            self.build_unmodified_peptides(session)

        assert session.query(InformedPeptide).filter(
            InformedPeptide.protein_id == None).count() == 0

    def run(self):
        try:
            session = self.manager.session()
            if self.hypothesis_id is None:
                hypothesis = self._bootstrap_hypothesis(session)
            else:
                hypothesis = session.query(self.hypothesis_type).get(self.hypothesis_id)

            self.load_proteome(session)

            session.add(hypothesis)
            session.commit()

            self._remove_duplicates(session)
            session.commit()
            self._display_protein_peptide_counts(session)

            hypothesis.parameters.update({
                "protein_file": self.mzid_path,
                "enzyme": self.enzymes,
                "constant_modifications": self.constant_modifications,
                "include_all_baseline": self.include_all_baseline,
                "variable_modifications": []
            })

            session.add(hypothesis)
            session.close()
        except Exception, e:
            logger.exception("%r", locals(), exc_info=e)
            raise e
