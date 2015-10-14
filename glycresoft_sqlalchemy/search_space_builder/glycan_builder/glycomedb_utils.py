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
from glypy.algorithms import subtree_search

from taxonomylite import Taxonomy

from sqlalchemy import select

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, make_transient, StructureMotif, TheoreticalGlycanStructureToMotifTable,
    PipelineModule, Taxon, ReferenceDatabase, ReferenceAccessionNumber, TheoreticalGlycanComposition,
    TheoreticalGlycanStructure, MS1GlycanHypothesis, MS2GlycanHypothesis, func)

from glycresoft_sqlalchemy.utils.database_utils import temp_table
from glycresoft_sqlalchemy.utils import appdir

from glycresoft_sqlalchemy.search_space_builder.glycan_builder import registry

logger = logging.getLogger("glycomedb_utils")


motif_families = {
    "N-Linked Glycans": "%N-Glycan%",
    "O-Linked Glycans": "%O-Glycan%",
}


def timestamp():
    return datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")

set_reducing_end = registry.set_reducing_end
derivatize = registry.derivatize

cache_name = "glycome-db-cache"
reference_hypothesis_name_prefix = "GlycomeDB-Reference-"


class GlycomeDBDownloader(PipelineModule):
    def __init__(self, database_path, drop_stems=True, drop_positions=True, verbose=False):
        self.manager = self.manager_type(database_path)
        self.drop_stems = drop_stems
        self.drop_positions = drop_positions
        self.verbose = verbose

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
        hypothesis = MS2GlycanHypothesis(name=reference_hypothesis_name_prefix + timestamp())
        session.add(hypothesis)
        glycomedb = ReferenceDatabase.get(session, name="Glycome-DB")
        session.add(glycomedb)
        session.commit()
        hypothesis_id = hypothesis.id
        i = 0

        motifs = session.query(StructureMotif).all()

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
                    calculated_mass=glycan.mass(),
                    hypothesis_id=hypothesis_id)
                record.taxa = taxa
                for motif in motifs:
                    if motif.matches(record):
                        record.motifs.append(motif)
                record.references = [ReferenceAccessionNumber.get(session, accession, glycomedb.id)]
                session.add(record)
                if i % 1000 == 0:
                    session.commit()
                    self.inform("Commit %s", i)
            except glycoct.GlycoCTSectionUnsupported:
                pass
            except Exception, e:
                if self.verbose:
                    logger.exception("%s", accession, exc_info=e)
        session.commit()
        self.hypothesis_id = hypothesis_id

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
                calculated_mass=GlycanComposition.parse(composition).mass())
            session.add(tgc)
            session.flush()
            logger.info("Relating Composition to Structures %s", tgc)
            session.query(TheoreticalGlycanStructure).filter(
                TheoreticalGlycanStructure.composition == composition).update({
                "composition_reference_id": tgc.id
                }, synchronize_session=False)
            logger.info("Assigning Taxa")
            taxa_assoc = [{"taxon_id": id[0], "entity_id": tgc.id} for id in session.query(
                func.distinct(Taxon.id)).join(
                TheoreticalGlycanStructure.taxa).filter(
                TheoreticalGlycanStructure.hypothesis_id == hypothesis_id)]
            session.execute(CompositionTaxonomyAssociationTable.insert(), taxa_assoc)
            if i % 1000 == 0:
                session.commit()
        session.commit()
        return hypothesis_id


class TaxonomyFilter(PipelineModule):
    def __init__(
            self, database_path, hypothesis_id, taxonomy_path=None,
            taxa_ids=None, include_descendent_taxa=True, motif_family=None):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.taxa_ids = taxa_ids
        self.include_descendent_taxa = include_descendent_taxa
        self.motif_family = motif_family

        self.taxonomy_path = taxonomy_path

    def run(self):
        return self.stream_chosen_glycans()

    def stream_chosen_glycans(self):
        taxonomy = None
        if self.include_descendent_taxa:
            if self.taxonomy_path is None or not os.path.exists(self.taxonomy_path):
                if self.taxonomy_path is None:
                    self.taxonomy_path = "taxonomy.db"
                self.inform("Downloading NCBI Taxonomic Hierarchy From Source")
                taxonomy = Taxonomy.from_source()
            else:
                taxonomy = Taxonomy(self.taxonomy_path)

        session = self.manager.session()

        TempTaxonTable = temp_table(Taxon)
        conn = session.connection()
        TempTaxonTable.create(conn)

        for taxon_id in self.taxa_ids:
            if self.include_descendent_taxa:
                all_ids = taxonomy.children(taxon_id, deep=True)
            else:
                all_ids = [taxon_id]

            conn.execute(TempTaxonTable.insert(), {"id": i for i in all_ids})
        session.commit()
        conn = session.connection()
        chosen_taxa = select([TempTaxonTable.c.id])

        last_composition_reference_id = None
        composition_reference = None

        query = session.query(func.distinct(TheoreticalGlycanStructure.id)).join(
            TheoreticalGlycanStructure.TaxonomyAssociationTable).filter(
            TheoreticalGlycanStructure.TaxonomyAssociationTable.c.taxon_id.in_(chosen_taxa)).order_by(
            TheoreticalGlycanStructure.composition_reference_id)

        if self.motif_family is not None:
            query = query.join(TheoreticalGlycanStructure.motifs).filter(StructureMotif.name.like(
                motif_families[self.motif_family]))
        i = 0
        for id, in query:
            i += 1
            if i % 1000 == 0:
                self.inform("%d Glycans processed" % i)
            structure = session.query(TheoreticalGlycanStructure).get(id)
            if structure.composition_reference_id != last_composition_reference_id:
                last_composition_reference_id = structure.composition_reference_id
                composition_reference = session.query(TheoreticalGlycanComposition).get(last_composition_reference_id)
                yield (TheoreticalGlycanComposition, composition_reference)
            yield TheoreticalGlycanStructure, structure
        TempTaxonTable.drop(conn)
        session.commit()


@registry.composition_source_type.register("glycome-db")
class GlycomeDBHypothesis(PipelineModule):
    batch_size = 5000

    def __init__(
            self, database_path, hypothesis_id=None, glycomedb_path=None,
            taxonomy_path=None, taxa_ids=None, include_descendent_taxa=False, include_structures=True,
            motif_family=None, reduction=None, derivatization=None):
        self.manager = self.manager_type(database_path)
        self.glycomedb_path = glycomedb_path
        self.taxonomy_path = taxonomy_path
        self.taxa_ids = taxa_ids
        self.include_descendent_taxa = include_descendent_taxa
        self.hypothesis_id = hypothesis_id
        self.include_structures = include_structures
        self.motif_family = motif_family
        self.reduction = reduction
        self.derivatization = derivatization

    def run(self):
        self.manager.initialize()
        session = self.manager.session()
        hypothesis = None
        if self.hypothesis_id is None:
            if self.include_structures:
                hypothesis = MS2GlycanHypothesis()
            else:
                hypothesis = MS1GlycanHypothesis()
            session.add(hypothesis)
            session.commit()
            self.hypothesis_id = hypothesis.id
        else:
            hypothesis = session.query(Hypothesis).get(self.hypothesis_id)

        hypothesis.parameters = hypothesis.parameters or {}
        hypothesis.parameters['taxa_ids'] = self.taxa_ids
        hypothesis.parameters['taxa_include_descendent_taxa'] = self.include_descendent_taxa
        hypothesis.parameters['include_structures'] = self.include_structures
        hypothesis.parameters['motif_family'] = self.motif_family

        self.fetch_relevant_glycans()

    def fetch_relevant_glycans(self):
        session = self.manager.session()
        glycomedb_reference = self.resolve_glycomedb()
        taxonomy_filter = TaxonomyFilter(
            self.glycomedb_path, glycomedb_reference.id, self.taxonomy_path,
            self.taxa_ids, self.include_descendent_taxa, self.motif_family)

        i = 0

        reduction = self.reduction
        derivatization = self.derivatization

        last_composition_reference_id = None
        for record_type, record in taxonomy_filter.start():
            taxa = list(record.taxa)
            if record_type is TheoreticalGlycanComposition:
                if reduction is not None or derivatization is not None:
                    composition = record.glycan_composition.copy()
                    if reduction is not None:
                        set_reducing_end(composition, reduction)
                    if derivatization is not None:
                        derivatize(composition, derivatization)
                    new_record = TheoreticalGlycanComposition(
                        calculated_mass=composition.mass(),
                        composition=record.composition,
                        derivatization=derivatization,
                        reduction=reduction,
                        hypothesis_id=self.hypothesis_id
                        )
                else:
                    new_record = TheoreticalGlycanComposition(
                        calculated_mass=record.calculated_mass,
                        composition=record.composition,
                        derivatization=record.derivatization,
                        reduction=record.reduction,
                        hypothesis_id=self.hypothesis_id)
            elif record_type is TheoreticalGlycanStructure:
                if not self.include_structures:
                    continue
                if reduction is not None or derivatization is not None:
                    structure = record.structure().clone()
                    if reduction is not None:
                        set_reducing_end(structure, reduction)
                    if derivatization is not None:
                        derivatize(structure, derivatization)
                    new_record = TheoreticalGlycanStructure(
                        calculated_mass=structure.mass(),
                        composition=record.composition,
                        derivatization=derivatization,
                        reduction=reduction,
                        glycoct=str(structure),
                        hypothesis_id=self.hypothesis_id
                        )
                else:
                    new_record = TheoreticalGlycanStructure(
                        calculated_mass=record.calculated_mass,
                        composition=record.composition,
                        derivatization=record.derivatization,
                        reduction=record.reduction,
                        glycoct=record.glycoct,
                        hypothesis_id=self.hypothesis_id)
            else:
                raise Exception("Unknown Record Type: %s" % record_type)
            session.add(new_record)
            session.flush()

            if record_type is TheoreticalGlycanComposition:
                last_composition_reference_id = new_record.id
            elif record_type is TheoreticalGlycanStructure:
                new_record.composition_reference_id = last_composition_reference_id
                session.execute(
                    TheoreticalGlycanStructureToMotifTable.insert(),
                    [{"glycan_id": new_record.id, "motif_id": motif.id} for motif in record.motifs])

            session.execute(record_type.TaxonomyAssociationTable.insert(), [
                {"taxon_id": t.id, "entity_id": new_record.id} for t in taxa])
            if len(record.references) != 0:
                session.execute(record_type.ReferenceAccessionAssocationTable.insert(), [
                    {"entity_id": new_record.id, "accession_code": r.id, "database_id": r.database_id} for r in record.references
                    ])
            i += 1

            if i % self.batch_size == 0:
                # self.inform("%d records processed" % i)
                print ("%d records processed" % i)
                session.commit()
        session.commit()

    def resolve_glycomedb(self):
        try:
            glycomedb_manager = self.manager_type(self.glycomedb_path)
            try:
                glycomedb_session = glycomedb_manager.session()
                glycomedb_reference = glycomedb_session.query(
                    Hypothesis).filter(Hypothesis.name.like(reference_hypothesis_name_prefix + "%")).first()
                if glycomedb_reference is None:
                    raise Exception("No Reference Found")
                return glycomedb_reference
            except Exception, e:
                logger.exception(
                    "An error occurred while resolving the GlycomeDB reference database",
                    exc_info=e)
                job = GlycomeDBDownloader(self.glycomedb_path)
                hypothesis_id = job.start()
                glycomedb_reference = glycomedb_manager.session().query(Hypothesis).get(hypothesis_id)
                if glycomedb_reference is None:
                    raise Exception("No Reference Found")
                return glycomedb_reference
        except Exception, e:
            logger.exception(
                ("An error occured while resolving the GlycomeDB reference database."
                 " The reference could not be created. Check to see if %s is a valid file path or"
                 " database URI.") % self.glycomedb_path, exc_info=e)
            raise e
