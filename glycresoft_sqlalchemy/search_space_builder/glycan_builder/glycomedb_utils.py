import logging
import os
import datetime
import gzip
import logging
import requests
from lxml import etree

from glypy.utils import StringIO
from glypy.io import glycoct
from glypy import Glycan, GlycanComposition

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, Taxon, ReferenceDatabase, ReferenceAccessionNumber, TheoreticalGlycanComposition,
    TheoreticalGlycanStructure, MS1GlycanHypothesis, MS2GlycanHypothesis, func)


logger = logging.getLogger("glycomedb_utils")

def timestamp():
    return datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")


cache_name = "glycome-db-cache"


class GlycomeDBDownloader(PipelineModule):
    def __init__(self, database_path, drop_stems=True, drop_positions=True):
        self.manager = self.manager_type(database_path)
        self.drop_stems = drop_stems
        self.drop_positions = drop_positions

    def run(self):
        self.manager.initialize()
        if os.path.exists(cache_name):
            data_source = open(cache_name, "rb")
        else:
            response = requests.get(u'http://www.glycome-db.org/http-services/getStructureDump.action?user=eurocarbdb')
            response.raise_for_status()
            open(cache_name, "wb").write(response.content)
            data_source = open(cache_name, "rb")

        handle = gzip.GzipFile(fileobj=data_source)
        xml = etree.parse(handle)
        session = self.manager.session()
        hypothesis = MS2GlycanHypothesis(name="GlycomeDB-" + timestamp())
        session.add(hypothesis)
        glycomedb = ReferenceDatabase.get(session, name="Glycome-DB")
        session.add(glycomedb)
        session.commit()
        hypothesis_id = hypothesis.id
        i = 0

        logger.info("Parsing database structures")
        drop_stems = self.drop_stems
        drop_positions = self.drop_positions
        CompositionTaxonomyAssociationTable = TheoreticalGlycanComposition.TaxonomyAssociationTable
        for structure in xml.iterfind(".//structure"):
            try:
                accession = structure.attrib['id']
                i += 1
                glycoct_str = structure.find("sequence").text
                taxa = [Taxon.get(session, t.attrib['ncbi']) for t in structure.iterfind(".//taxon")]
                glycan = glycoct.loads(glycoct_str)
                if (glycoct.loads(str(glycan)).mass() - glycan.mass()) > 0.00001:
                    raise Exception("Mass did not match on reparse, %f" % (glycoct.loads(str(glycan)).mass() - glycan.mass()))

                composition = GlycanComposition.from_glycan(glycan)
                if drop_stems:
                    composition.drop_stems()
                if drop_positions:
                    composition.drop_positions()
                composition.drop_configurations()
                composition.collapse()

                reduction = "ReducedEnd" if glycan.reducing_end else None

                record = TheoreticalGlycanStructure(
                    glycoct=glycoct_str,
                    composition=composition.serialize(),
                    reduction=reduction,
                    theoretical_mass=glycan.mass(),
                    hypothesis_id=hypothesis_id)
                record.taxa = taxa
                record.references = [ReferenceAccessionNumber.get(session, accession, glycomedb.id)]
                session.add(record)
                if i % 1000 == 0:
                    session.commit()
                    self.inform("Commit %s", i)
            except glycoct.GlycoCTSectionUnsupported:
                pass
            except Exception, e:
                logger.exception("%s", accession, exc_info=e)
        session.commit()

        self.inform("Linking Compositions")
        i = 0
        for composition, in session.query(
                func.distinct(TheoreticalGlycanStructure.composition)).filter(
                TheoreticalGlycanStructure.hypothesis_id == hypothesis_id):
            i += 1
            logger.info("Creating Composition %s", composition)
            tgc = TheoreticalGlycanComposition(
                composition=composition,
                hypothesis_id=hypothesis_id,
                theoretical_mass=GlycanComposition.parse(composition).mass())
            session.add(tgc)
            session.flush()
            logger.info("Relating Composition to Structures %s", tgc)
            session.query(TheoreticalGlycanStructure).filter(
                TheoreticalGlycanStructure.composition == composition).update({
                "composition_reference_id": tgc.id
                }, synchronize_session=False)
            logger.info("Assigning Taxa")
            taxa_assoc = [{"taxon_id": id[0], "entity_id": tgc.id} for id in session.query(func.distinct(Taxon.id)).join(
                TheoreticalGlycanStructure.taxa).filter(
                TheoreticalGlycanStructure.hypothesis_id == hypothesis_id)]
            session.execute(CompositionTaxonomyAssociationTable.insert(), taxa_assoc)
            if i % 1000 == 0:
                session.commit()
        session.commit()
