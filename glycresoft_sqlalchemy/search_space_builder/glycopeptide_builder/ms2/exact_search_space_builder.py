import multiprocessing
import logging
import functools

from glycresoft_sqlalchemy.structure.sequence import Sequence, strip_modifications

from glypy import GlycanComposition

from .search_space_builder import TheoreticalSearchSpaceBuilder, constructs, BatchingTheoreticalSearchSpaceBuilder
from ..utils import WorkItemCollectionFlat as WorkItemCollection, fragments

from glycresoft_sqlalchemy.data_model import (
    TheoreticalGlycopeptide, Protein,
    ExactMS1GlycopeptideHypothesisSampleMatch,
    ExactMS2GlycopeptideHypothesis)

from sqlalchemy import func


logger = logging.getLogger("search_space_builder")


def generate_fragments(seq, ms1_result):
    """Construct an instance of :class:`TheoreticalGlycopeptide` and
    compute all backbone fragments of interest.

    Parameters
    ----------
    seq : Sequnce
        An instance of :class:`glycresoft_sqlalchemy.structure.sequence.Sequence` which
        will be fragmented.
    ms1_result : MS1GlycopeptideResult
        Provides context for the created TheoreticalGlycopeptide instance

    Returns
    -------
    TheoreticalGlycopeptide
    """
    seq.glycan = GlycanComposition.parse(ms1_result.glycan_composition_str)
    seq_mod = seq.get_sequence(include_glycan=False)
    (oxonium_ions, b_ions, y_ions,
     b_ions_hexnac, y_ions_hexnac,
     stub_ions) = fragments(seq)

    # theoretical_glycopeptide = TheoreticalGlycopeptide(
    theoretical_glycopeptide = dict(
        ms1_score=ms1_result.ms1_score,
        observed_mass=ms1_result.observed_mass,
        calculated_mass=ms1_result.calculated_mass,
        ppm_error=ms1_result.ppm_error,
        volume=ms1_result.volume,
        count_glycosylation_sites=ms1_result.count_glycosylation_sites,
        count_missed_cleavages=ms1_result.count_missed_cleavages,
        start_position=ms1_result.start_position,
        end_position=ms1_result.end_position,
        base_peptide_sequence=strip_modifications(ms1_result.base_peptide_sequence),
        modified_peptide_sequence=seq_mod,
        peptide_modifications=ms1_result.peptide_modifications,
        glycopeptide_sequence=seq_mod + ms1_result.glycan_composition_str,
        sequence_length=len(seq),
        glycan_composition_str=ms1_result.glycan_composition_str,
        bare_b_ions=b_ions,
        bare_y_ions=y_ions,
        oxonium_ions=oxonium_ions,
        stub_ions=stub_ions,
        glycosylated_b_ions=b_ions_hexnac,
        glycosylated_y_ions=y_ions_hexnac,
        protein_id=ms1_result.protein_name,
        base_composition_id=ms1_result.composition_id,
        glycan_combination_id=ms1_result.glycan_combination_id
        )

    return theoretical_glycopeptide


def from_sequence(ms1_result, database_manager, protein_map, source_type):
    """Produce an instance of :class:`.TheoreticalGlycopeptide` from an MS1
    match with an exact sequence definition. This is a **task function** intended
    for execution in the main loop of `run` or in one of the worker processes therein.

    This function assumes that all peptide backbone modifications have already
    been placed exactly, directly writing the sequence into the more detailed
    structure, unlike the equivalent function used by `TheoreticalSearchSpaceBuilder`
    and `PoolingTheoreticalSearchSpaceBuilder`.

    Parameters
    ----------
    ms1_result : MS1GlycopeptideResult or int
        A generic structure representing an MS1 match from LC-MS
        database search. The structure is expected to reference an
        exact glycopeptide as opposed to an ambiguous one requiring
        combinatorial expansion.

        Alternatively may be an int object corresponding to the primary key
        :attr:`id` of a :class:`PeakGroupMatch` object which will be used
        to construct an instance of :class:`MS1GlycopeptideResult`.
    database_manager : DatabaseManager
        Provides connection to database
    protein_map : dict
        Maps :attr:`.Protein.name` to :attr:`.Protein.id`
    source_type : MS1ResultsFile or MS1ResultsFacade
        The source of ms1_result, defining a :method:`render` method which
        converts the raw input into a database

    Returns
    -------
    TheoreticalGlycopeptide:
        The single glycopeptide with pre-computed fragment ion masses. The structure
        has not yet been persisted to the database and must be written there by the
        receiving process.
    """
    try:
        session = database_manager.session()
        ms1_result = source_type.render(session, ms1_result)
        if len(ms1_result.base_peptide_sequence) == 0:
            return None
        seq = Sequence(ms1_result.most_detailed_sequence)
        seq.glycan = ''
        product = generate_fragments(seq, ms1_result)
        if not isinstance(product.protein_id, int):
            product.protein_id = protein_map[product.protein_id]
        # print product.glycans.all(), product.glycans.count()
        session.close()
        return product
    except Exception, e:
        logging.exception("An error occurred, %r", ms1_result, exc_info=e)
        raise


#  @constructs.references(ExactMS1GlycopeptideHypothesisSampleMatch)
class ExactSearchSpaceBuilder(TheoreticalSearchSpaceBuilder):
    HypothesisType = ExactMS2GlycopeptideHypothesis

    """Summary"""
    def __init__(self, ms1_results_file, db_file_name, enzyme, site_list=None,
                 n_processes=4, **kwargs):
        try:
            kwargs.pop("variable_modifications")
            kwargs.pop("constant_modifications")
        except:
            pass

        super(ExactSearchSpaceBuilder, self).__init__(
            ms1_results_file, db_file_name,
            constant_modifications=[], variable_modifications=[], enzyme=enzyme,
            site_list=site_list, n_processes=n_processes, **kwargs)

    def prepare_task_fn(self):
        """Construct the partial function to be applied to each input in this step of
        of the pipeline. The base function used is :func:`from_sequence`

        Returns
        -------
        function

        See Also
        --------
        from_sequence
        """
        protein_map = dict(self.session.query(Protein.name, Protein.id).filter(
            Protein.hypothesis_id == self.hypothesis.id))
        return functools.partial(
            from_sequence, database_manager=self.manager, protein_map=protein_map, source_type=self.ms1_format)

    def run(self):
        """Main work method.

        Applies the function returned by :meth:`prepare_task_fn` to each item to be processed.

        If :attr:`n_processes` > 1, this will be done using multiple worker processes through a
        `multiprocessing.Pool` object with :attr:`n_processes` workers.

        Returns
        -------
        int : The :attr:`id` of the constructed :class:`.MS2GlycopeptideHypothesis`
        """
        self.bootstrap()

        session = self.session
        cntr = 0
        last = 0
        commit_interval = 1000
        task_fn = self.prepare_task_fn()
        accumulator = WorkItemCollection(session)

        session.add(self.hypothesis)
        session.commit()
        id = self.hypothesis.id

        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for theoretical in pool.imap_unordered(task_fn, self.ms1_results_reader, chunksize=100):
                if theoretical is None:
                    continue
                accumulator.add(theoretical)
                cntr += 1
                if cntr > (last + commit_interval):
                    accumulator.commit()
                    logger.info("%d records handled", cntr)
                    last = cntr
            pool.terminate()
        for ms1_result in self.ms1_results_reader:
            theoretical = task_fn(ms1_result)
            if theoretical is None:
                continue
            accumulator.add(theoretical)
            cntr += 1
            if cntr > last + commit_interval:
                accumulator.commit()
                logger.info("%d records handled", cntr)
                last = cntr

        accumulator.commit()

        logger.info("Checking integrity")

        # Remove duplicates
        ids = session.query(func.min(TheoreticalGlycopeptide.id)).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == id).group_by(
            TheoreticalGlycopeptide.glycopeptide_sequence,
            TheoreticalGlycopeptide.protein_id)

        q = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == id,
            ~TheoreticalGlycopeptide.id.in_(ids.correlate(None)))
        conn = session.connection()
        conn.execute(TheoreticalGlycopeptide.__table__.delete(
            TheoreticalGlycopeptide.__table__.c.id.in_(q.selectable)))
        session.commit()
        final_count = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == id).count()

        logger.info("%d Theoretical Glycopeptides created", final_count)

        session.close()
        return id


def batch_from_sequence(ms1_results, database_manager, protein_map, source_type):
    try:
        session = database_manager.session()
        glycopeptide_acc = []
        i = 0
        for ms1_result in [source_type.render(session, ms1_result) for ms1_result in ms1_results]:

            if len(ms1_result.base_peptide_sequence) == 0:
                return None
            seq = Sequence(ms1_result.most_detailed_sequence)
            seq.glycan = ''
            product = generate_fragments(seq, ms1_result)
            if not isinstance(product["protein_id"], int):
                product["protein_id"] = protein_map[product["protein_id"]]
            glycopeptide_acc.append(product)
            i += 1
        session.bulk_insert_mappings(TheoreticalGlycopeptide, glycopeptide_acc)
        session.commit()
        return i
    except Exception, e:
        logging.exception("An error occurred in batch_from_sequence", exc_info=e)
        raise e


@constructs.references(ExactMS1GlycopeptideHypothesisSampleMatch)
class BatchingExactSearchSpaceBuilder(BatchingTheoreticalSearchSpaceBuilder):
    HypothesisType = ExactMS2GlycopeptideHypothesis

    def __init__(self, ms1_results_file, db_file_name, enzyme, site_list=None,
                 n_processes=4, **kwargs):
        try:
            kwargs.pop("variable_modifications")
            kwargs.pop("constant_modifications")
        except:
            pass

        super(BatchingExactSearchSpaceBuilder, self).__init__(
            ms1_results_file, db_file_name,
            constant_modifications=[], variable_modifications=[], enzyme=enzyme,
            site_list=site_list, n_processes=n_processes, **kwargs)

    def prepare_task_fn(self):
        """Construct the partial function to be applied to each input in this step of
        of the pipeline. The base function used is :func:`from_sequence`

        Returns
        -------
        function

        See Also
        --------
        from_sequence
        """
        protein_map = dict(self.session.query(Protein.name, Protein.id).filter(
            Protein.hypothesis_id == self.hypothesis.id))
        return functools.partial(
            batch_from_sequence, database_manager=self.manager, protein_map=protein_map, source_type=self.ms1_format)

    def run(self):
        self.bootstrap()

        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0
        step = 1000

        session = self.session

        session.add(self.hypothesis)
        session.commit()
        id = self.hypothesis.id

        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in worker_pool.imap_unordered(task_fn, self.stream_results(), chunksize=1):
                cntr += res
                if (cntr > last + step):
                    last = cntr
                    logger.info("Committing, %d records made", cntr)

            worker_pool.terminate()
        else:
            logger.debug("Building theoretical sequences sequentially")
            for row in self.ms1_results_reader:
                res = task_fn(row)
                cntr += res
                if (cntr > last + step):
                    last = cntr
                    logger.info("Committing, %d records made", cntr)

        logger.info("Checking integrity")

        self._remove_duplicates()
        final_count = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.protein_id == Protein.id,
            Protein.hypothesis_id == id).count()

        logger.info("%d Theoretical Glycopeptides created", final_count)

        session.close()
        return id
