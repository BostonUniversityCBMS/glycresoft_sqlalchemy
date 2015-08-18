import logging
try:
    logger = logging.getLogger("include_glycomics")
except:
    pass
from glycresoft_sqlalchemy.data_model import DatabaseManager, Glycan, TheoreticalGlycopeptide, Protein, GlycopeptideMatch, Hypothesis
from glycresoft_sqlalchemy.data_model import PipelineModule

from glycresoft_sqlalchemy.structure import glycans as ms1_glycomics
from ..glycan_builder.composition_source import GlycanCompositionHypothesisBuilder


glycan_file_type_map = {
   "txt": GlycanCompositionHypothesisBuilder
}


class MS1GlycanImporter(PipelineModule):
    manager_type = DatabaseManager

    def __init__(self, database_path, glycan_path, hypothesis_id, glycan_file_type="hypothesis"):
        self.manager = self.manager_type(database_path)
        self.glycan_file_type = glycan_file_type
        self.hypothesis_id = hypothesis_id
        self.glycan_path = glycan_path
