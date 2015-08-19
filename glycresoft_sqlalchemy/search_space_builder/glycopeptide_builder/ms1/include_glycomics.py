import logging
try:
    logger = logging.getLogger("include_glycomics")
except:
    pass

from glycresoft_sqlalchemy.data_model import PipelineModule

from ...glycan_builder.composition_source import GlycanCompositionHypothesisBuilder


glycan_file_type_map = {
   "txt": lambda database_path, glycan_path, id, options: GlycanCompositionHypothesisBuilder(
    database_path, glycan_path, "txt", hypothesis_id=id, **options)
}


class MS1GlycanImporter(PipelineModule):

    def __init__(self, database_path, glycan_path, hypothesis_id, glycan_file_type="txt", **kwargs):
        self.manager = self.manager_type(database_path)
        self.glycan_file_type = glycan_file_type
        self.hypothesis_id = hypothesis_id
        self.glycan_path = glycan_path
        self.options = kwargs

    def run(self):
        loader = glycan_file_type_map[self.glycan_file_type](
            self.database_path, self.glycan_path, self.hypothesis_id,
            self.options)
        loader.start()
