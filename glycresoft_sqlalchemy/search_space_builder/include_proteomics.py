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
from ..proteomics.mzid_sa import Proteome as MzIdentMLProteome

class ProteomeImporter(PipelineModule):
    manager_type = DatabaseManager

    def __init__(self, database_path, mzid_path, hypothesis_id=None):
        self.manager = self.manager_type(database_path)
        self.mzid_path = mzid_path
        self.hypothesis_id = hypothesis_id

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
                hypothesis = Hypothesis(name=name, parameters={"mzid_path": self.mzid_path})
                session.add(hypothesis)
                session.commit()
                self.hypothesis_id = hypothesis.id
            session.close()
            MzIdentMLProteome(self.manager.path, self.mzid_path, self.hypothesis_id)
        except Exception, e:
            logger.exception("%r", locals(), exc_info=e)
            raise e
