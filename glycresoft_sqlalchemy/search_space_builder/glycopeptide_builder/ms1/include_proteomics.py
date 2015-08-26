import os
import datetime
import logging
try:
    logger = logging.getLogger("include_proteomics")
except:
    pass


from ..peptide_utilities import SiteListFastaFileParser
from glycresoft_sqlalchemy.structure import sequence, modification

from glycresoft_sqlalchemy.data_model import MS1GlycopeptideHypothesis, Protein
from glycresoft_sqlalchemy.data_model import PipelineModule, make_transient, InformedPeptide

from glycresoft_sqlalchemy.proteomics.mzid_sa import Proteome as MzIdentMLProteome

from glycresoft_sqlalchemy.structure.sequence import find_n_glycosylation_sequons


def make_base_sequence(peptide, constant_modifications, modification_table):
    seq = peptide.base_peptide_sequence
    seq = sequence.Sequence(seq)
    for mod in {modification_table[const_mod]
                for const_mod in constant_modifications}:
        for site in mod.find_valid_sites(seq):
            seq.add_modification(site, mod.name)
    make_transient(peptide)
    peptide.id = None
    peptide.modified_peptide_sequence = str(seq)
    peptide.calculated_mass = seq.mass
    return peptide


class ProteomeImporter(PipelineModule):

    def __init__(self, database_path, mzid_path, glycosylation_sites_file=None,
                 hypothesis_id=None, constant_modifications=("Carbamidomethyl (C)",)):
        self.manager = self.manager_type(database_path)
        self.mzid_path = mzid_path
        self.hypothesis_id = hypothesis_id
        self.glycosylation_sites_file = glycosylation_sites_file
        self.constant_modifications = constant_modifications

    def run(self):
        try:
            session = self.manager.session()
            if self.hypothesis_id is None:
                self.manager.initialize()
                tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
                name = 'target-{}-{}'.format(
                        os.path.splitext(
                            os.path.basename(self.mzid_path)
                            )[0],
                        tag)
                hypothesis = MS1GlycopeptideHypothesis(name=name, parameters={"mzid_path": self.mzid_path})
                session.add(hypothesis)
                session.commit()
                self.hypothesis_id = hypothesis.id
            session.close()
            MzIdentMLProteome(self.manager.path, self.mzid_path, self.hypothesis_id)

            session = self.manager.session()
            if self.glycosylation_sites_file is None:
                for protein in session.query(Protein):
                    protein.glycosylation_sites = find_n_glycosylation_sequons(protein.protein_sequence)
                    session.add(protein)
            else:
                site_list_gen = SiteListFastaFileParser(self.glycosylation_sites_file)
                site_list_map = {d['name']: d["glycosylation_sites"] for d in site_list_gen}
                for protein in session.query(Protein):
                    try:
                        protein.glycosylation_sites = site_list_map[protein.name]
                    except KeyError:
                        protein.glycosylation_sites = find_n_glycosylation_sequons(protein.protein_sequence)
                    session.add(protein)

            session.add(hypothesis)
            session.commit()

            modification_table = modification.RestrictedModificationTable.bootstrap(list(self.constant_modifications), [])
            basic_peptides = []
            for peptide in session.query(InformedPeptide).filter(
                    InformedPeptide.protein_id == Protein.id, Protein.hypothesis_id == self.hypothesis_id).group_by(
                    InformedPeptide.base_peptide_sequence):
                peptide = make_base_sequence(peptide, self.constant_modifications, modification_table)
                basic_peptides.append(peptide)
                if len(basic_peptides) > 1000:
                    session.add_all(basic_peptides)
                    session.commit()
                    basic_peptides = []

            session.add_all(basic_peptides)
            session.commit()

            session.close()
        except Exception, e:
            logger.exception("%r", locals(), exc_info=e)
            raise e
