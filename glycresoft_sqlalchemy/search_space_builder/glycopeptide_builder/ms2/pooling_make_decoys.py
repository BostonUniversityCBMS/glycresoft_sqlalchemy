import multiprocessing
import functools
import itertools
import logging

from glycresoft_sqlalchemy.data_model import TheoreticalGlycopeptide

from .make_decoys import (DecoySearchSpaceBuilder, strip_modifications, decoy_type_map,
                          reverse_preserve_sequon, fragments)

logger = logging.getLogger(__name__)


def make_decoy(theoretical_sequence_id, prefix_len=0, suffix_len=1,
               protein_decoy_map=None, database_manager=None,
               permute_fn=reverse_preserve_sequon):
    try:
        session = database_manager.session()

        theoretical_sequence = session.query(TheoreticalGlycopeptide).get(theoretical_sequence_id)

        permuted_sequence = permute_fn(theoretical_sequence.glycopeptide_sequence,
                                       prefix_len=prefix_len, suffix_len=suffix_len)

        (oxonium_ions, bare_b_ions, bare_y_ions, glycosylated_b_ions,
            glycosylated_y_ions, stub_ions) = fragments(permuted_sequence)

        decoy = TheoreticalGlycopeptide(
            ms1_score=theoretical_sequence.ms1_score,
            observed_mass=theoretical_sequence.observed_mass,
            calculated_mass=theoretical_sequence.calculated_mass,
            ppm_error=theoretical_sequence.ppm_error,
            volume=theoretical_sequence.volume,
            count_glycosylation_sites=theoretical_sequence.count_glycosylation_sites,
            count_missed_cleavages=theoretical_sequence.count_missed_cleavages,
            start_position=theoretical_sequence.start_position,
            end_position=theoretical_sequence.end_position,
            base_peptide_sequence=strip_modifications(str(permuted_sequence)),
            modified_peptide_sequence=str(permuted_sequence),
            peptide_modifications=theoretical_sequence.peptide_modifications,
            glycopeptide_sequence=str(permuted_sequence) + theoretical_sequence.glycan_composition_str,
            sequence_length=len(permuted_sequence),
            glycan_composition_str=theoretical_sequence.glycan_composition_str,
            bare_b_ions=bare_b_ions,
            bare_y_ions=bare_y_ions,
            oxonium_ions=oxonium_ions,
            stub_ions=stub_ions,
            glycosylated_b_ions=glycosylated_b_ions,
            glycosylated_y_ions=glycosylated_y_ions,
            protein_id=protein_decoy_map[theoretical_sequence.protein_id]
        )

        return decoy
    except Exception, e:
        logger.exception("An error occurred in make_decoy\n%r", locals(), exc_info=e)
        raise e
    finally:
        session.close()


class PoolingDecoySearchSpaceBuilder(DecoySearchSpaceBuilder):
    def __init__(self, database_path, prefix_len=0, suffix_len=1, hypothesis_ids=None,
                 n_processes=4, decoy_type=0, commit_checkpoint=1000):
        super(PoolingDecoySearchSpaceBuilder, self).__init__(
            database_path=database_path, prefix_len=prefix_len,
            suffix_len=suffix_len, hypothesis_ids=hypothesis_ids,
            decoy_type=decoy_type, n_processes=n_processes)
        self.commit_checkpoint = commit_checkpoint

    def prepare_task_fn(self):
        return functools.partial(make_decoy, prefix_len=self.prefix_len, suffix_len=self.suffix_len,
                                 protein_decoy_map=self.protein_decoy_map, database_manager=self.manager,
                                 permute_fn=decoy_type_map[self.decoy_type])

    def run(self):
        task_fn = self.prepare_task_fn()
        cntr = 0
        session = self.manager.session()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap_unordered(task_fn, self.stream_theoretical_glycopeptides(), chunksize=500):
                cntr += 1
                session.add(res)
                if cntr % self.commit_checkpoint == 0:
                    session.commit()
                    logger.info("%d Decoys Complete." % cntr)
            pool.terminate()
        else:
            for res in itertools.imap(task_fn, self.stream_theoretical_glycopeptides()):
                cntr += 1
                session.add(res)
                if cntr % self.commit_checkpoint == 0:
                    session.commit()
                    logger.info("%d Decoys Complete." % cntr)
        session.commit()
        return self.decoy_hypothesis_ids
