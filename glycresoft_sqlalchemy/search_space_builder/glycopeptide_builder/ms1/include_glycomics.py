import logging
import os

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, DatabaseManager, Hypothesis)
from ...glycan_builder.composition_source import (
    GlycanCompositionHypothesisBuilder, OtherGlycanHypothesisGlycanHypothesisBuilder)
from ..glycan_utilities import create_combinations

from glycresoft_sqlalchemy.utils.collectiontools import decoratordict
from glycresoft_sqlalchemy.utils.database_utils import database_exists as db_exists


class UnknownGlycomicsFormatError(Exception):
    pass


class InvalidGlycomicsConfigurationError(Exception):
    pass


class MissingGlycomicsDataError(Exception):
    pass


logger = logging.getLogger("include_glycomics")

glycan_file_type_map = {
   "txt": lambda database_path, glycan_path, id, options: GlycanCompositionHypothesisBuilder(
    database_path, glycan_path, "txt", hypothesis_id=id, **options),
   "hypothesis": lambda database_path, glycan_path, id, options: OtherGlycanHypothesisGlycanHypothesisBuilder(
    database_path, hypothesis_id=id, source_path=glycan_path, **options)
}


class MS1GlycanImporter(PipelineModule):

    def __init__(self, database_path, glycan_path, hypothesis_id, glycan_file_type="txt",
                 combination_size=2, **kwargs):
        self.manager = self.manager_type(database_path)
        self.glycan_file_type = glycan_file_type
        self.hypothesis_id = hypothesis_id
        self.glycan_path = glycan_path
        self.combination_size = combination_size
        self.options = kwargs

    def run(self):
        session = self.manager()
        loader = glycan_file_type_map[self.glycan_file_type](
            self.database_path, self.glycan_path, self.hypothesis_id,
            self.options)
        loader.start(verbose=False)
        self.inform("Building glycan combinations")
        count = create_combinations(session, self.combination_size, self.hypothesis_id)
        logger.info("Produced %d glycan composition combinations", count)
        if count == 0:
            raise MissingGlycomicsDataError("No glycan combinations were produced")


keyword_validator = decoratordict()


@keyword_validator("txt")
def validate_txt(manager):
    exists = os.path.exists(manager.glycomics_path)
    return exists, "Invalid glycan text file path"


@keyword_validator("hypothesis")
def validate_hypothesis(manager):
    id = manager.options.get("source_hypothesis_id")
    if id is None:
        return False, "Missing glycan source hypothesis id"
    db_manager = DatabaseManager(manager.glycomics_path)
    connection_manager = db_manager.connection_manager
    url = connection_manager.construct_url()
    exists = db_exists(url)
    if exists:
        session = db_manager()
        hypothesis = session.query(Hypothesis).get(id)
        return hypothesis is not None, "Source hypothesis does not exist"
    else:
        return False, "Source database does not exist"


class MS1GlycanImportManager(object):
    def __init__(self, glycomics_path, glycomics_format, hypothesis_id, maximum_glycosylation_sites, **kwargs):
        self.glycomics_path = glycomics_path
        self.glycomics_format = glycomics_format
        self.hypothesis_id = hypothesis_id
        self.maximum_glycosylation_sites = maximum_glycosylation_sites
        self.options = kwargs

        self.validate()

    def validate(self):
        valid, reason = keyword_validator[self.glycomics_format](self)
        if not valid:
            raise InvalidGlycomicsConfigurationError(reason)

    def import_glycans(self):
        session = self.manager()
        hypothesis = session.query(Hypothesis).get(self.hypothesis_id)
        hypothesis.parameters.update({
            "glycomics_path": self.glycomics_path,
            "glycomics_format": self.glycomics_format,
            "maximum_glycosylation_sites": self.maximum_glycosylation_sites
            })
        session.add(hypothesis)
        session.commit()
        session.close()
        if self.glycomics_format == "txt":
            importer = MS1GlycanImporter(
                self.database_path, self.glycomics_path, self.hypothesis_id, 'txt',
                combination_size=self.maximum_glycosylation_sites)
            importer.start()
        elif self.glycomics_format == "hypothesis":
            importer = MS1GlycanImporter(
                self.database_path, self.glycomics_path, self.hypothesis_id, "hypothesis",
                source_hypothesis_id=self.options.get("source_hypothesis_id"),
                combination_size=self.maximum_glycosylation_sites)
            importer.start()
        elif self.glycomics_format == "hypothesis-sample-match":
            importer = MS1GlycanImporter(
                self.database_path, self.glycomics_path, self.hypothesis_id, "hypothesis-sample-match",
                source_hypothesis_sample_match_id=self.options.get("source_hypothesis_sample_match_id"),
                combination_size=self.maximum_glycosylation_sites)
            importer.start()
        else:
            raise UnknownGlycomicsFormatError(self.glycomics_format)
