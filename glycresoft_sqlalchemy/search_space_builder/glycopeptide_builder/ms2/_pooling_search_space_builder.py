import multiprocessing
import logging
import functools

from .search_space_builder import (MS1GlycopeptideResult, constructs,
                                   get_peptide_modifications, get_search_space,
                                   generate_fragments, TheoreticalSearchSpaceBuilder)
from ..utils import WorkItemCollectionFlat as WorkItemCollection
from glycresoft_sqlalchemy.data_model import MS1GlycopeptideHypothesisSampleMatch

logger = logging.getLogger("search_space_builder")


def process_predicted_ms1_ion(row, database_manager, modification_table, site_list_map, proteins, source_type):
    """Multiprocessing dispatch function to generate all theoretical sequences and their
    respective fragments from a given MS1 result

    Parameters
    ----------
    row: dict
        Line mapping from an MS1 results csv file
    modification_table: :class:`.modifications.RestrictedModificationTable`
        Modification table limited to only the rules specified by the user
    site_list: list of int
        List of putative glycosylation sites from the parent protein
    monosaccharide_identities: list
        List of glycan or monosaccaride names

    Returns
    -------
    list of dicts:
        List of each theoretical sequence and its fragment ions
    """
    try:
        session = database_manager.session()
        ms1_result = source_type.render(session, row)

        if (ms1_result.base_peptide_sequence == '') or (ms1_result.count_glycosylation_sites == 0):
            return []

        # Compute the set of modifications that can occur.
        mod_list = get_peptide_modifications(
            ms1_result.peptide_modifications, modification_table)

        # Get the start and end positions of fragment relative to the
        glycan_sites = set(site_list_map.get(ms1_result.protein_name, [])).intersection(
            range(ms1_result.start_position, ms1_result.end_position))

        # No recorded sites, skip this component.
        if len(glycan_sites) == 0:
            return []

        # Adjust the glycan_sites to relative position
        glycan_sites = [x - ms1_result.start_position for x in glycan_sites]
        ss = get_search_space(
            ms1_result, glycan_sites, mod_list)
        seq_list = ss.get_theoretical_sequence(ms1_result.count_glycosylation_sites)
        fragments = [generate_fragments(seq, ms1_result)
                     for seq in seq_list]

        for sequence in fragments:
            sequence.protein_id = proteins[ms1_result.protein_name]

        return fragments
    except Exception, e:
        logger.exception("An error occurred, %r", locals(), exc_info=e)
        raise e


#@constructs.references(MS1GlycopeptideHypothesisSampleMatch)
class PoolingTheoreticalSearchSpaceBuilder(TheoreticalSearchSpaceBuilder):
    '''
    Describe the process of generating all theoretical sequences and their fragments
    from an MS1 Results CSV, a collection of constant and variable peptide modifications,
    and a list of putative glycosylation sites.

    Uses an :class:`sqlitedict.SqliteDict` instance as a storage backend.

    Includes a single- and multi-process compatible implementation. The more processes used,
    the more memory must be allocated to buffer results.
    '''

    def __init__(self, ms1_results_file, db_file_name=None,
                 enzyme=None,
                 site_list=None,
                 constant_modifications=None,
                 variable_modifications=None,
                 n_processes=4,
                 commit_checkpoint=1000, **kwargs):
        super(PoolingTheoreticalSearchSpaceBuilder, self).__init__(
            ms1_results_file, db_file_name,
            constant_modifications=constant_modifications,
            enzyme=enzyme, site_list=site_list,
            variable_modifications=variable_modifications, n_processes=n_processes, **kwargs)
        self.commit_checkpoint = commit_checkpoint

    def prepare_task_fn(self):
        task_fn = functools.partial(process_predicted_ms1_ion,
                                    database_manager=self.manager,
                                    modification_table=self.modification_table,
                                    site_list_map=self.glycosylation_site_map,
                                    source_type=self.ms1_format,
                                    proteins={k: v.id for k, v in self.hypothesis.proteins.items()})
        return task_fn

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes in a pool, or
        in the main process if :attr:`n_processes` == 1.

        All records produced in worker processes are pooled in the main process, and
        are commit to the database approximately every :attr:`commit_checkpoint` records
        or so to prevent bottlenecks.
        '''
        task_fn = self.prepare_task_fn()
        cntr = 0
        checkpoint = 0
        accumulator = []
        accumulator = WorkItemCollection(self.session)
        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in worker_pool.imap_unordered(task_fn, self.stream_results(), chunksize=500):
                [accumulator.add(i) for i in res]
                cntr += len(res)
                if cntr >= checkpoint:
                    logger.info("Committing, %d records made", cntr)
                    accumulator.commit()
                    checkpoint = cntr + self.commit_checkpoint
            worker_pool.terminate()
        else:
            logger.debug("Building theoretical sequences sequentially")
            for row in self.ms1_results_reader:
                res = task_fn(row)
                [accumulator.add(i) for i in res]
                cntr += len(res)
                if cntr >= checkpoint:
                    logger.info("Committing, %d records made", cntr)
                    accumulator.commit()
                    checkpoint = cntr + self.commit_checkpoint
        accumulator.commit()
        return self.hypothesis_id
