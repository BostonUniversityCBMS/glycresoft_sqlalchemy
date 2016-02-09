import itertools
import functools
import multiprocessing
import logging

from sqlalchemy import distinct

from glycresoft_sqlalchemy.data_model import (
    TheoreticalGlycopeptide, GlycopeptideMatch, PipelineModule, GlycopeptideSpectrumMatch)

from glycresoft_sqlalchemy.scoring import simple_scoring_algorithm

logger = logging.getLogger('spectrum_assignment')


def score_spectrum_matches(times, database_manager, scorer, score_parameters):
    try:
        scores_collection = []
        spectrum_match_collection = []
        session = database_manager()
        for time in times:
            matches = session.query(GlycopeptideSpectrumMatch).filter(
                GlycopeptideSpectrumMatch.scan_time == time[0]).all()
            best_score = -(float('inf'))
            best_match = None
            for match in matches:
                score = scorer(match.as_match_like(), match.theoretical_glycopeptide, **score_parameters)
                scores_collection.append(score)
                spectrum_match_collection.append(match)
                if score.value > best_score:
                    best_score = score.value
                    best_match = match
            best_match.best_match = True
        session.bulk_save_objects(scores_collection)
        session.bulk_save_objects(spectrum_match_collection)
        session.commit()
        return len(times)
    except Exception, e:
        logger.exception("An exception occcured in score_spectrum_matches", exc_info=e)
    finally:
        session.close()


class SpectrumAssigner(PipelineModule):
    def __init__(self, database_path, hypothesis_id, hypothesis_sample_match_id, scorer=None,
                 score_parameters=None, n_processes=4):
        if scorer is None:
            scorer = simple_scoring_algorithm.SimpleSpectrumScorer()
        if score_parameters is None:
            score_parameters = {}

        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.scorer = scorer
        self.score_parameters = score_parameters
        self.n_processes = n_processes

    def prepare_task_fn(self):
        return functools.partial(
            score_spectrum_matches,
            database_manager=self.manager,
            scorer=self.scorer,
            score_parameters=self.score_parameters
            )

    def stream_spectrum_time_id(self, chunk_size=20):
        session = self.manager()
        all_times = session.query(distinct(GlycopeptideSpectrumMatch.scan_time)).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.hypothesis_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).all()
        final = len(all_times)
        last = 0
        while 1:
            next_chunk = all_times[last:(last + chunk_size)]
            if last <= final:
                yield next_chunk
                last += chunk_size
            else:
                break

    def run(self):
        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for res in pool.imap(task_fn, self.stream_spectrum_time_id()):
                cntr += res
                if (cntr - last) > 1000:
                    logger.info("%d Assignments Complete." % cntr)
                    last = cntr
            pool.terminate()

        else:
            for step in self.self.stream_spectrum_time_id():
                cntr += task_fn(step)
                if (cntr - last) > 1000:
                    logger.info("%d Assignments Complete." % cntr)
                    last = cntr


def theoretical_to_match(theoretical, hypothesis_sample_match_id):
    match = GlycopeptideMatch(
        protein_id=theoretical.protein_id,
        theoretical_glycopeptide_id=theoretical.id,
        ms1_score=theoretical.ms1_score,
        observed_mass=theoretical.observed_mass,
        calculated_mass=theoretical.calculated_mass,
        volume=theoretical.volume,
        ppm_error=theoretical.ppm_error,

        glycan_composition_str=theoretical.glycan_composition_str,
        glycan_combination_id=theoretical.glycan_combination_id,

        base_peptide_sequence=theoretical.base_peptide_sequence,
        modified_peptide_sequence=theoretical.modified_peptide_sequence,
        glycopeptide_sequence=theoretical.glycopeptide_sequence,
        sequence_length=theoretical.sequence_length,
        peptide_modifications=theoretical.peptide_modifications,
        count_glycosylation_sites=theoretical.count_glycosylation_sites,
        count_missed_cleavages=theoretical.count_missed_cleavages,
        glycosylation_sites=theoretical.glycosylation_sites,
        start_position=theoretical.start_position,
        end_position=theoretical.end_position,

        hypothesis_sample_match_id=hypothesis_sample_match_id)
    return match


def build_matches(theoretical_ids, database_manager, hypothesis_sample_match_id, score_name):
    session = database_manager()

    def best_scoring_key(match):
        return match.scores[score_name]

    summary_matches = []

    for theoretical_id in theoretical_ids:
        theoretical = session.query(TheoreticalGlycopeptide).get(theoretical_id[0])

        spectrum_matches = session.query(GlycopeptideSpectrumMatch).filter(
            GlycopeptideSpectrumMatch.theoretical_glycopeptide_id == theoretical_id[0],
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.best_match).order_by(GlycopeptideSpectrumMatch.scan_time.asc()).all()

        match_total = theoretical_to_match(theoretical, hypothesis_sample_match_id)
        match_total.scan_id_range = [s.scan_time for s in spectrum_matches]
        match_total.first_scan = match_total.scan_id_range[0]
        match_total.last_scan = match_total.scan_id_range[-1]

        best_scoring = max(spectrum_matches, key=best_scoring_key)
        ion_series_matches = simple_scoring_algorithm.split_ion_list(best_scoring.ion_matches())
        match_total.bare_b_ions = ion_series_matches.get("bare_b_ions", [])
        match_total.bare_y_ions = ion_series_matches.get("bare_y_ions", [])
        match_total.glycosylated_b_ions = ion_series_matches.get("glycosylated_b_ions", [])
        match_total.glycosylated_y_ions = ion_series_matches.get("glycosylated_y_ions", [])
        match_total.oxonium_ions = ion_series_matches.get("oxonium_ions", [])
        match_total.stub_ions = ion_series_matches.get("stub_ions", [])

        score = best_scoring_key(best_scoring)
        match_total.ms2_score = score
        summary_matches.append(match_total)

    session.add_all(summary_matches)
    session.flush()
    conn = session.connection()
    T_GlycopeptideSpectrumMatch = GlycopeptideSpectrumMatch.__table__
    for summary in summary_matches:
        conn.execute(T_GlycopeptideSpectrumMatch.update().where(
            (T_GlycopeptideSpectrumMatch.c.theoretical_glycopeptide_id == summary.theoretical_glycopeptide_id) and
            (T_GlycopeptideSpectrumMatch.c.hypothesis_sample_match_id == summary.hypothesis_sample_match_id)).values(
            glycopeptide_match_id=summary.id))
    session.commit()
    return len(theoretical_ids)


class SummaryMatchBuilder(PipelineModule):
    def __init__(self, database_path, hypothesis_id, hypothesis_sample_match_id, score_name, n_processes=4):
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.score_name = score_name
        self.n_processes = n_processes

    def prepare_task_fn(self):
        return functools.partial(
            build_matches,
            database_manager=self.manager,
            hypothesis_id=self.hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            score_name=self.score_name)

    def stream_theoretical_ids(self, chunk_size=500):
        session = self.manager()

        ids = session.query(TheoreticalGlycopeptide.id).filter(
            TheoreticalGlycopeptide.hypothesis_id == self.hypothesis_id).all()

        step = chunk_size
        last = 0
        final = len(ids)
        while last <= final:
            yield ids[last:(last + step)]
            last += step

    def run(self):
        task_fn = self.prepare_task_fn()
        cntr = 0
        last = 0

        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for increment in pool.imap_unordered(task_fn, self.stream_theoretical_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d records completed", cntr)
                    last = cntr
            pool.terminate()
        else:
            for increment in itertools.imap(task_fn, self.stream_theoretical_ids()):
                cntr += increment
                if cntr - last > 1000:
                    logger.info("%d records completed", cntr)
                    last = cntr
