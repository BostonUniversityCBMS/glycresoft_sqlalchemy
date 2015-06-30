import os
import datetime
import logging
try:
    logger = logging.getLogger("include_proteomics")
except:
    pass
from glycresoft_ms2_classification.proteomics import mzid
from glycresoft_ms2_classification.structure.parser import strip_modifications

from ..data_model import DatabaseManager, Hypothesis, Protein
from ..data_model.informed_proteomics import InformedPeptide
from ..data_model import PipelineModule


class ProteomeImporter(PipelineModule):
    manager_type = DatabaseManager

    def __init__(self, database_path, mzid_path, hypothesis_id=None):
        self.manager = self.manager_type(database_path)
        self.mzid_path = mzid_path
        self.hypothesis_id = hypothesis_id

    def run(self):
        session = self.manager.session()
        try:
            if self.hypothesis_id is None:
                self.manager.initialize()
                tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
                name = 'target-{}-{}'.format(
                        os.path.splitext(
                            os.path.basename(self.mzid_path)
                            )[0],
                        tag)
                hypothesis = Hypothesis(name=name)
                session.add(hypothesis)
                session.commit()
                self.hypothesis_id = hypothesis.id

            proteome = mzid.extract_proteins(self.mzid_path)

            for reference_protein in proteome:
                protein = Protein(name=reference_protein.accession,
                                  protein_sequence=reference_protein.sequence,
                                  other=reference_protein.metadata,
                                  glycosylation_sites=list(reference_protein.n_glycan_sequon_sites),
                                  hypothesis_id=self.hypothesis_id)
                session.add(protein)
                session.commit()
                for reference_peptide in reference_protein:
                    if reference_peptide.mod_index.get("HexNAc", 0) > 0:
                        reference_peptide.deglycosylate()
                    n_glycosites = reference_peptide.n_glycan_sequon_sites
                    peptide = InformedPeptide(
                        protein_id=protein.id,
                        base_peptide_sequence=strip_modifications(str(reference_peptide)),
                        modified_peptide_sequence=str(reference_peptide),
                        calculated_mass=reference_peptide.mass,
                        glycosylation_sites=list(n_glycosites),
                        count_glycosylation_sites=len(n_glycosites),
                        peptide_score=reference_peptide.peptide_score,
                        start_position=reference_peptide.start,
                        end_position=reference_peptide.end,
                        other=reference_peptide.metadata)
                    session.add(peptide)

                session.commit()
            session.commit()
        except Exception, e:
            logger.exception("%r", locals(), exc_info=e)
            raise e
        session.close()
