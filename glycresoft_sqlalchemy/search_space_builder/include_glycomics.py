from ..data_model import DatabaseManager, Glycan, TheoreticalGlycopeptide, Protein, GlycopeptideMatch

from glycresoft_ms2_classification.structure import glycans as ms1_glycomics

glycan_file_type_map = {
    "hypothesis": ms1_glycomics.GlycanHypothesis,
    "ms1_matches": ms1_glycomics.GlycanMS1Results
}


class MS1GlycanImporter(object):
    manager_type = DatabaseManager

    def __init__(self, database_path, glycan_path, experiment_id, glycan_file_type="hypothesis"):
        self.manager = self.manager_type(database_path)
        self.glycan_file_type = glycan_file_type
        self.glycan_reader = glycan_file_type_map[glycan_file_type]
        self.experiment_id = experiment_id

    def run(self):
        session = self.manager.session()
        for glycan in self.glycan_reader.glycans:
            g = Glycan(composition=glycan.composition,
                       mass=glycan.mass, other=glycan.data,
                       experiment_id=self.experiment_id)
            session.add(g)
        session.commit()
        session.close()


def map_glycans_by_composition(database_path, experiment_id,
                               target_tables=(TheoreticalGlycopeptide, GlycopeptideMatch)):
    session = DatabaseManager(database_path).session()
    for composition, glycan_id in (session.query(Glycan.composition, Glycan.id).filter(
            Glycan.experiment_id == experiment_id)):
        for table in target_tables:
            session.query(table).filter(
                table.glycan_composition_str == composition,
                table.protein_id == Protein.id,
                Protein.experiment_id == experiment_id).update(
                {"glycan_id": glycan_id}, synchronize_session=False)
            session.commit()
    session.commit()
    session.close()
