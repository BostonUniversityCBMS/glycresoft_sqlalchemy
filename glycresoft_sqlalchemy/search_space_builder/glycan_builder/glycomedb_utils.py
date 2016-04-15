import logging
import os
import datetime
import gzip
import logging
import requests
from lxml import etree

from glypy.io import glycoct
from glypy import GlycanComposition
from glypy.composition.glycan_composition import FrozenGlycanComposition

from taxonomylite import Taxonomy

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, make_transient, StructureMotif, TheoreticalGlycanStructureToMotifTable,
    TheoreticalGlycanCompositionToMotifTable, PipelineModule, Taxon, ReferenceDatabase,
    ReferenceAccessionNumber, TheoreticalGlycanComposition, TheoreticalGlycanStructure,
    MS1GlycanHypothesis, MS2GlycanHypothesis, func)

from glycresoft_sqlalchemy.utils.database_utils import temp_table
from glycresoft_sqlalchemy.utils.data_files import taxonomylite_store, glycomedb_store, glycomedb_download_cache

from glycresoft_sqlalchemy.search_space_builder.glycan_builder import registry

logger = logging.getLogger("glycomedb_utils")


motif_families = {
    "N-Linked": "%N-Glycan%",
    "O-Linked": "%O-Glycan%",
}


def timestamp():
    return datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")

set_reducing_end = registry.set_reducing_end
derivatize = registry.derivatize

cache_name = glycomedb_download_cache
reference_hypothesis_name_prefix = "GlycomeDB-Reference-"


class GlycomeDBDownloader(PipelineModule):
    def __init__(self, database_path, drop_stems=True, drop_positions=True, verbose=False):
        self.manager = self.manager_type(database_path)
        self.drop_stems = drop_stems
        self.drop_positions = drop_positions
        self.verbose = verbose

    def map_motifs_from_structure_to_composition(self, session, composition):
        motif_ids = {x[1] for x in session.query(
                TheoreticalGlycanStructure.id, TheoreticalGlycanStructureToMotifTable.c.motif_id).filter(
                TheoreticalGlycanStructure.composition_reference_id == composition.id,
                TheoreticalGlycanStructure.hypothesis_id == self.hypothesis_id,
                ).join(
                TheoreticalGlycanStructure.motifs).group_by(TheoreticalGlycanStructureToMotifTable.c.motif_id)}
        return tuple(motif_ids)

    def migrate_motifs_to_compositions(self, session, composition):
        motif_ids = self.map_motifs_from_structure_to_composition(session, composition)
        session.execute(
            TheoreticalGlycanCompositionToMotifTable.insert(),
            [{"glycan_id": composition.id, "motif_id": m_id} for m_id in motif_ids])

    def migrate_taxa_to_compositions(self, session, composition):
        taxa_assoc = [{"taxon_id": id[0], "entity_id": composition.id} for id in session.query(
                func.distinct(Taxon.id)).join(
                TheoreticalGlycanStructure.taxa).filter(
                TheoreticalGlycanStructure.hypothesis_id == self.hypothesis_id)]
        session.execute(TheoreticalGlycanComposition.TaxonomyAssociationTable.insert(), taxa_assoc)

    def run(self):
        self.manager.initialize()
        logger.debug("Checking %s for downloaded data", cache_name)
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
        self.hypothesis_id = hypothesis.id
        i = 0

        motifs = session.query(StructureMotif).all()

        drop_stems = self.drop_stems
        drop_positions = self.drop_positions

        taxa = {int(t.attrib['ncbi']) for t in xml.iterfind(".//taxon")}
        [Taxon.get(session, tid) for tid in taxa]
        session.flush()


        logger.info("Parsing database structures")
        taxon_acc = []
        motif_acc = []
        for structure in xml.iterfind(".//structure"):
            try:
                accession = structure.attrib['id']
                glycoct_str = structure.find("sequence").text
                taxa = [int(t.attrib['ncbi']) for t in structure.iterfind(".//taxon")]
                glycan = glycoct.loads(glycoct_str)
                if (glycoct.loads(str(glycan)).mass() - glycan.mass()) > 0.00001:
                    # Parity Error
                    continue

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
                    hypothesis_id=self.hypothesis_id)
                record.references = [ReferenceAccessionNumber.get(session, accession, glycomedb.id)]
                session.add(record)
                session.flush()
                for motif in motifs:
                    if motif.matches(record):
                        motif_acc.append({"motif_id": motif.id, "glycan_id": record.id})

                taxon_acc.extend({"taxon_id": tid, "entity_id": record.id} for tid in taxa)
                i += 1
                if (i % 100) == 0:
                    session.commit()
                    if taxon_acc:
                        session.execute(TheoreticalGlycanStructure.TaxonomyAssociationTable.insert(),
                                        taxon_acc)
                    if motif_acc:
                        session.execute(TheoreticalGlycanStructureToMotifTable.insert(), motif_acc)
                    taxon_acc = []
                    motif_acc = []
                    session.commit()
                    self.inform("Commit %r", i)
            except (KeyboardInterrupt, IntegrityError):
                raise
            except (glycoct.GlycoCTSectionUnsupported, IndexError):
                pass
            except Exception, e:
                logger.exception("%s", accession, exc_info=e)
                pass
                if isinstance(e, KeyboardInterrupt):
                    raise
        if taxon_acc:
            session.execute(TheoreticalGlycanStructure.TaxonomyAssociationTable.insert(),
                            taxon_acc)
        if motif_acc:
            session.execute(TheoreticalGlycanStructureToMotifTable.insert(), motif_acc)

        self.inform("Imported %d structures", i)
        session.commit()
        self.inform("Linking Compositions")
        i = 0
        for composition, in session.query(
                func.distinct(TheoreticalGlycanStructure.composition)).filter(
                TheoreticalGlycanStructure.hypothesis_id == self.hypothesis_id):
            i += 1
            # logger.info("Creating Composition %s", composition)
            tgc = TheoreticalGlycanComposition(
                composition=composition,
                hypothesis_id=self.hypothesis_id,
                calculated_mass=FrozenGlycanComposition.parse(composition).mass())
            session.add(tgc)
            session.flush()
            # logger.info("Relating Composition to Structures %s", tgc)
            session.query(TheoreticalGlycanStructure).filter(
                TheoreticalGlycanStructure.composition == composition).update({
                "composition_reference_id": tgc.id
                }, synchronize_session=False)
            # logger.info("Assigning Taxa")
            self.migrate_taxa_to_compositions(session, tgc)
            # logger.info("Assigning Motifs")
            # self.migrate_motifs_to_compositions(session, tgc)

            if i % 100 == 0:
                session.commit()
                self.inform("Commit %s", i)
        session.commit()
        return self.hypothesis_id


class TaxonomyFilter(PipelineModule):
    def __init__(
            self, database_path, hypothesis_id, taxonomy_path=None,
            taxa_ids=None, include_descendent_taxa=True, motif_family=None):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.taxa_ids = taxa_ids
        self.include_descendent_taxa = include_descendent_taxa
        self.motif_family = motif_family

        self.taxonomy_path = taxonomy_path or taxonomylite_store
        self.temp_taxonomy_table = None
        self._all_taxonomy_ids = None

    def run(self):
        return self.stream_chosen_glycans()

    def resolve_taxonomy_db(self):
        taxonomy = None
        if self.include_descendent_taxa:
            if self.taxonomy_path is None:
                self.taxonomy_path = taxonomylite_store
            if not os.path.exists(self.taxonomy_path):
                logger.info("Downloading NCBI Taxonomic Hierarchy From Source")
                if self.taxonomy_path is not None:
                    taxonomy = Taxonomy.from_source(store_path=self.taxonomy_path)
                else:
                    taxonomy = Taxonomy.from_source()
            else:
                taxonomy = Taxonomy(self.taxonomy_path)
        return taxonomy

    def get_all_taxonomies(self):
        if self._all_taxonomy_ids is not None:
            return self._all_taxonomy_ids
        taxonomy = self.resolve_taxonomy_db()
        results = []
        for taxon_id in self.taxa_ids:
            if self.include_descendent_taxa:
                all_ids = [taxon_id] + taxonomy.children(taxon_id, deep=True)
            else:
                all_ids = [taxon_id]
            results.extend(all_ids)
        results = list(set(results))
        self._all_taxonomy_ids = results
        return results

    def prefetch_linked_taxonomy_ids(self):
        session = self.manager()
        TempTaxonTable = self.construct_taxonomy_set(session)
        chosen_taxa = select([TempTaxonTable.c.id])
        query = session.query(TheoreticalGlycanStructure.id, Taxon.id).join(
            TheoreticalGlycanStructure.TaxonomyAssociationTable).filter(
            TheoreticalGlycanStructure.TaxonomyAssociationTable.c.taxon_id.in_(chosen_taxa)).order_by(
            TheoreticalGlycanStructure.composition_reference_id).group_by(Taxon.id)

        if self.motif_family is not None:
            query = query.join(TheoreticalGlycanStructure.motifs).filter(StructureMotif.name.like(
                motif_families[self.motif_family]))
        return tuple(set(int(x[1]) for x in query))

    def construct_taxonomy_set(self, session):
        if self.temp_taxonomy_table is not None:
            return self.temp_taxonomy_table

        conn = session.connection()
        self.temp_taxonomy_table = temp_table(Taxon)
        self.temp_taxonomy_table.create(conn)
        all_ids = self.get_all_taxonomies()
        conn.execute(self.temp_taxonomy_table.insert(), [{"id": i} for i in set(all_ids)])
        session.commit()

        return self.temp_taxonomy_table

    def stream_chosen_glycans(self):
        session = self.manager.session()

        TempTaxonTable = self.construct_taxonomy_set(session)

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

        def slurp(iterable, session, chunk_size=100):
            all_items = [x[0] for x in iterable]
            total = len(all_items)
            last = 0
            while last <= total:
                chunk = all_items[last:(last + chunk_size)]
                chunk_objs = session.query(TheoreticalGlycanStructure).filter(
                    TheoreticalGlycanStructure.id.in_(chunk)).all()
                for obj in chunk_objs:
                    yield obj
                last += chunk_size

        for structure in slurp(query, session):
            i += 1
            if i % 100 == 0:
                self.inform("%d Glycans processed" % i)
            # structure = session.query(TheoreticalGlycanStructure).get(id)
            if structure.composition_reference_id != last_composition_reference_id:
                last_composition_reference_id = structure.composition_reference_id
                composition_reference = session.query(TheoreticalGlycanComposition).get(last_composition_reference_id)
                yield (TheoreticalGlycanComposition, composition_reference)
            yield TheoreticalGlycanStructure, structure
        conn = session.connection()
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
                hypothesis = MS2GlycanHypothesis(name=timestamp())
            else:
                hypothesis = MS1GlycanHypothesis(name=timestamp())
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

    def make_new_record(self, record_type, record):
        if record_type is TheoreticalGlycanComposition:
            if self.reduction is not None or self.derivatization is not None:
                composition = record.glycan_composition.copy()
                if self.reduction is not None:
                    set_reducing_end(composition, self.reduction)
                if self.derivatization is not None:
                    derivatize(composition, self.derivatization)
                new_record = TheoreticalGlycanComposition(
                    calculated_mass=composition.mass(),
                    composition=record.composition,
                    derivatization=self.derivatization,
                    reduction=self.reduction,
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
                return
            if self.reduction is not None or self.derivatization is not None:
                structure = record.structure().clone()
                if self.reduction is not None:
                    set_reducing_end(structure, self.reduction)
                if self.derivatization is not None:
                    derivatize(structure, self.derivatization)
                new_record = TheoreticalGlycanStructure(
                    calculated_mass=structure.mass(),
                    composition=record.composition,
                    derivatization=self.derivatization,
                    reduction=self.reduction,
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
        return new_record

    def fetch_relevant_glycans(self):
        session = self.manager.session()
        glycomedb_reference = self.resolve_glycomedb()
        taxonomy_filter = TaxonomyFilter(
            self.glycomedb_path, glycomedb_reference.id, self.taxonomy_path,
            self.taxa_ids, self.include_descendent_taxa, self.motif_family)

        i = 0

        motif_acc = []
        reference_acc = {
            TheoreticalGlycanComposition: [],
            TheoreticalGlycanStructure: []
        }

        taxonomy_acc = {
            TheoreticalGlycanComposition: [],
            TheoreticalGlycanStructure: []
        }

        last_composition_reference_id = None
        taxa_ids = taxonomy_filter.prefetch_linked_taxonomy_ids()
        [Taxon.get(session, taxid) for taxid in taxa_ids]
        session.commit()
        for record_type, record in taxonomy_filter.start():
            new_record = self.make_new_record(record_type, record)
            if new_record is None:
                continue
            session.add(new_record)
            session.flush()

            if record_type is TheoreticalGlycanComposition:
                last_composition_reference_id = new_record.id
            elif record_type is TheoreticalGlycanStructure:
                new_record.composition_reference_id = last_composition_reference_id
                motif_acc.extend([{"glycan_id": new_record.id, "motif_id": motif.id} for motif in record.motifs])
                if len(motif_acc) > 1000:
                    session.execute(TheoreticalGlycanStructureToMotifTable.insert(), motif_acc)
                    motif_acc = []

            i += 1

            # taxonomy_acc[record_type].extend({"taxon_id": t.id, "entity_id": new_record.id} for t in record.taxa)
            if len(taxonomy_acc[record_type]) > 1000:
                session.execute(record_type.TaxonomyAssociationTable.insert(), taxonomy_acc[record_type])
                taxonomy_acc[record_type] = []
            if len(record.references) != 0:
                reference_acc[record_type].extend(
                    {"entity_id": new_record.id, "accession_code": r.id, "database_id": r.database_id}
                    for r in record.references)
                if len(reference_acc[record_type]) > 1000:
                    session.execute(
                        record_type.ReferenceAccessionAssocationTable.insert(), reference_acc[record_type])
                    reference_acc[record_type] = []

            if i % self.batch_size == 0:
                self.inform("%d records processed", i)
                session.commit()

        session.execute(TheoreticalGlycanStructureToMotifTable.insert(), motif_acc)
        for record_type, acc in reference_acc.items():
            if len(acc) == 0:
                continue
            session.execute(
                        record_type.ReferenceAccessionAssocationTable.insert(),
                        acc)
        for record_type, acc in taxonomy_acc.items():
            if len(acc) == 0:
                continue
            session.execute(
                record_type.TaxonomyAssociationTable.insert(),
                acc)
        self.inform("%d records processed", i)
        session.commit()
        return self.hypothesis_id

    def download_glycomedb(self):
        job = GlycomeDBDownloader(self.glycomedb_path)
        hypothesis_id = job.start()
        glycomedb_manager = self.manager_type(self.glycomedb_path)
        glycomedb_reference = glycomedb_manager.session().query(Hypothesis).get(hypothesis_id)
        if glycomedb_reference is None:
            raise Exception("No Reference Found")
        return glycomedb_reference

    def resolve_glycomedb(self):
        try:
            if self.glycomedb_path is None:
                self.glycomedb_path = glycomedb_store

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
                return self.download_glycomedb()
        except Exception, e:
            logger.exception(
                ("An error occured while resolving the GlycomeDB reference database."
                 " The reference could not be created. Check to see if %s is a valid file path or"
                 " database URI.") % self.glycomedb_path, exc_info=e)
            raise e
