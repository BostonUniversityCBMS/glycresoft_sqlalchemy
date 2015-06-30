import logging
try:
    logger = logging.getLogger("include_glycomics")
except:
    pass
from ..data_model import DatabaseManager, Glycan, TheoreticalGlycopeptide, Protein, GlycopeptideMatch, Hypothesis
from ..data_model import PipelineModule

from glycresoft_ms2_classification.structure import glycans as ms1_glycomics

glycan_file_type_map = {
    "hypothesis": ms1_glycomics.GlycanHypothesis,
    "ms1_matches": ms1_glycomics.GlycanMS1Results
}

glycan_composition_string = ms1_glycomics.Glycan.glycan_composition_string


class MS1GlycanImporter(PipelineModule):
    manager_type = DatabaseManager

    def __init__(self, database_path, glycan_path, hypothesis_id, glycan_file_type="hypothesis"):
        self.manager = self.manager_type(database_path)
        self.glycan_file_type = glycan_file_type
        self.glycan_reader = glycan_file_type_map[glycan_file_type](glycan_path)
        self.hypothesis_id = hypothesis_id

        session = self.manager.session()
        hypothesis = session.query(Hypothesis).filter(Hypothesis.id == hypothesis_id).first()
        hypothesis.parameters['monosaccharide_identities'] = self.glycan_reader.glycan_identities
        session.add(hypothesis)
        session.commit()

    def run(self):
        session = self.manager.session()
        try:
            for glycan in self.glycan_reader.glycans:
                g = Glycan(composition=glycan_composition_string(glycan.composition),
                           mass=glycan.mass, other=glycan.data,
                           hypothesis_id=self.hypothesis_id)
                session.add(g)
            session.commit()
        except Exception, e:
            logger.exception("%r", locals(), exc_info=e)
            raise
        finally:
            session.close()


def map_glycans_by_composition(database_path, hypothesis_id,
                               target_tables=(TheoreticalGlycopeptide, GlycopeptideMatch)):
    session = DatabaseManager(database_path).session()
    for composition, glycan_id in (session.query(Glycan.composition, Glycan.id).filter(
            Glycan.hypothesis_id == hypothesis_id)):
        for table in target_tables:
            session.query(table).filter(
                table.glycan_composition_str == composition,
                table.protein_id == Protein.id,
                Protein.hypothesis_id == hypothesis_id).update(
                {"glycan_id": glycan_id}, synchronize_session=False)
            session.commit()
    session.commit()
    session.close()
