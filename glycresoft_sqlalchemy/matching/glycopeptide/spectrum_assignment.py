import itertools
import functools
import operator
import multiprocessing
import logging

from sqlalchemy import distinct, bindparam

from glycresoft_sqlalchemy.data_model import (
    TheoreticalGlycopeptide, GlycopeptideMatch, PipelineModule,
    GlycopeptideSpectrumMatch, HypothesisSampleMatch, slurp)

from glycresoft_sqlalchemy.scoring import simple_scoring_algorithm, target_decoy

from glycresoft_sqlalchemy.utils import database_utils

from glycresoft_sqlalchemy.utils import collectiontools


logger = logging.getLogger('spectrum_assignment')


def score_spectrum_matches(times, database_manager, scorer, score_parameters, hypothesis_id,
                           hypothesis_sample_match_id):
    try:
        scores_collection = []
        spectrum_match_collection = []
        session = database_manager()

        def multislurp(session, times_ids):
            return [slurp(session, GlycopeptideSpectrumMatch, ids, flatten=False) for ids in times_ids]

        for matches in multislurp(session, times):
            # matches = session.query(GlycopeptideSpectrumMatch).filter(
            #     GlycopeptideSpectrumMatch.scan_time == time[0]).all()
            best_score = -(float('inf'))
            best_match = []
            for match in matches:
                match.best_match = False
                score = scorer(match.as_match_like(), match.theoretical_glycopeptide, **score_parameters)
                scores_collection.append(score)
                spectrum_match_collection.append(match)
                if score.value > best_score:
                    best_score = score.value
                    best_match = [match]
                elif abs(score.value - best_score) <= 1e-6:
                    best_match.append(match)

            for case in best_match:
                case.best_match = True
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
            score_parameters=self.score_parameters,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            hypothesis_id=self.hypothesis_id
            )

    def stream_spectrum_time_id(self, chunk_size=200):
        session = self.manager()

        chunker = session.query(GlycopeptideSpectrumMatch.id, GlycopeptideSpectrumMatch.scan_time).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.hypothesis_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).order_by(
            GlycopeptideSpectrumMatch.scan_time).all()

        get1 = operator.itemgetter(1)
        get0 = operator.itemgetter(0)

        groups = collectiontools.groupby(chunker, get1).items()

        # all_times = session.query(distinct(GlycopeptideSpectrumMatch.scan_time)).filter(
        #     GlycopeptideSpectrumMatch.hypothesis_id == self.hypothesis_id,
        #     GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id).all()
        final = len(groups)
        last = 0
        while 1:
            next_chunk = groups[last:(last + chunk_size)]
            next_chunk = [list(map(get0, group)) for group in map(get1, next_chunk)]
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
            for step in self.stream_spectrum_time_id():
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


def multislurp(session, theoretical_ids):
    theoreticals = slurp(
        session, TheoreticalGlycopeptide,
        [theoretical_id for theoretical_id, spectrum_match_ids in theoretical_ids],
        flatten=False)

    match_ids = [x for theoretical_id, group in theoretical_ids for x in group]

    matches = slurp(
        session, GlycopeptideSpectrumMatch,
        match_ids, flatten=False)

    theoreticals = collectiontools.groupby(theoreticals, operator.attrgetter("id"))
    matches = collectiontools.groupby(matches, lambda x: x.theoretical_glycopeptide_id)
    indices = theoreticals.keys()
    return zip([theoreticals[i][0] for i in indices], [matches[i] for i in indices])


def build_matches(theoretical_ids, database_manager, hypothesis_sample_match_id, score_name):
    session = database_manager()

    def best_scoring_key(match):
        return match.scores[score_name]

    summary_matches = []
    spectrum_matches_by_theoretical_id = {}

    def multislurp(session, theoretical_ids):
        theoreticals = slurp(
            session, TheoreticalGlycopeptide,
            [theoretical_id for theoretical_id, spectrum_match_ids in theoretical_ids],
            flatten=False)

        match_ids = [x for theoretical_id, group in theoretical_ids for x in group]

        matches = slurp(
            session, GlycopeptideSpectrumMatch,
            match_ids, flatten=False)

        theoreticals = collectiontools.groupby(theoreticals, operator.attrgetter("id"))
        matches = collectiontools.groupby(matches, lambda x: x.theoretical_glycopeptide_id)
        indices = theoreticals.keys()
        return zip([theoreticals[i][0] for i in indices], [matches[i] for i in indices])

    try:
        for theoretical, spectrum_matches in multislurp(session, theoretical_ids):
            theoretical_id = theoretical.id
            if len(spectrum_matches) == 0:
                continue

            match_total = theoretical_to_match(theoretical, hypothesis_sample_match_id)
            match_total.scan_id_range = [s.scan_time for s in spectrum_matches]
            match_total.first_scan = match_total.scan_id_range[0]
            match_total.last_scan = match_total.scan_id_range[-1]

            spectrum_matches_by_theoretical_id[theoretical_id] = spectrum_matches

            best_scoring = max(spectrum_matches, key=best_scoring_key)

            ion_series_matches = simple_scoring_algorithm.split_ion_list(best_scoring.ion_matches())
            match_total.bare_b_ions = ion_series_matches.get("bare_b_ions", [])
            match_total.bare_y_ions = ion_series_matches.get("bare_y_ions", [])
            match_total.glycosylated_b_ions = ion_series_matches.get("glycosylated_b_ions", [])
            match_total.glycosylated_y_ions = ion_series_matches.get("glycosylated_y_ions", [])
            match_total.oxonium_ions = ion_series_matches.get("oxonium_ions", [])
            match_total.stub_ions = ion_series_matches.get("stub_ions", [])

            match_total.mean_coverage = simple_scoring_algorithm.mean_coverage(match_total)
            match_total.mean_hexnac_coverage = simple_scoring_algorithm.mean_hexnac_coverage(match_total, theoretical)

            score = best_scoring_key(best_scoring)
            match_total.ms2_score = score
            try:
                match_total.q_value = best_scoring.scores['q_value']
            except KeyError:
                pass
            summary_matches.append(match_total)
        session.add_all(summary_matches)
        session.commit()

        spectrum_match_updates = []

        for summary in summary_matches:
            matches = spectrum_matches_by_theoretical_id[summary.theoretical_glycopeptide_id]
            for match in matches:
                if match.id is None:
                    continue
                spectrum_match_updates.append({"b_id": match.id, "glycopeptide_match_id": summary.id})

        if spectrum_match_updates:
            conn = session.connection()
            T_GlycopeptideSpectrumMatch = GlycopeptideSpectrumMatch.__table__
            stmt = T_GlycopeptideSpectrumMatch.update().where(
                T_GlycopeptideSpectrumMatch.c.id == bindparam("b_id")).values(
                glycopeptide_match_id=bindparam("glycopeptide_match_id"))
            conn.execute(stmt, spectrum_match_updates)
            session.commit()
        return len(theoretical_ids)
    except KeyboardInterrupt, e:
        logger.exception("An exception occcured in build_matches", exc_info=e)
        raise e
        return -1


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
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            score_name=self.score_name)

    def stream_theoretical_ids(self, chunk_size=100):
        session = self.manager()

        index = database_utils.get_index(GlycopeptideSpectrumMatch.__table__, "theoretical_glycopeptide_id")
        force_index = "indexed by " + index.name

        ids_query = session.query(
            GlycopeptideSpectrumMatch.theoretical_glycopeptide_id, GlycopeptideSpectrumMatch.id).filter(
            GlycopeptideSpectrumMatch.hypothesis_id == self.hypothesis_id,
            GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            GlycopeptideSpectrumMatch.best_match)
        ids_query = ids_query.with_hint(GlycopeptideSpectrumMatch, force_index, "sqlite")
        getter = operator.itemgetter(0)
        ids_grouper = collectiontools.groupby(ids_query, getter).items()
        result = []
        for key, grouped in ids_grouper:
            gsm_ids = [gsm_id for key, gsm_id in grouped]
            result.append((key, gsm_ids))

        ids = result

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
        session = self.manager()
        index_controller = database_utils.toggle_indices(session, GlycopeptideMatch)
        logger.info("Begin GlycopeptideMatch Creation")
        index_controller.drop()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for increment in pool.imap_unordered(task_fn, self.stream_theoretical_ids()):
                cntr += increment
                if cntr - last > 200:
                    logger.info("%d records completed", cntr)
                    last = cntr
            pool.terminate()
        else:
            for increment in itertools.imap(task_fn, self.stream_theoretical_ids()):
                cntr += increment
                if cntr - last > 200:
                    logger.info("%d records completed", cntr)
                    last = cntr
        logger.info("End GlycopeptideMatch Creation. Reconstructing Indices")
        index_controller.create()


class SpectrumMatchAnalyzer(PipelineModule):
    def __init__(self, database_path, hypothesis_sample_match_id, scorer=None, score_parameters=None, n_processes=4):
        if scorer is None:
            scorer = simple_scoring_algorithm.SimpleSpectrumScorer()
        if score_parameters is None:
            score_parameters = {}

        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.scorer = scorer
        self.score_parameters = score_parameters
        self.n_processes = n_processes

        session = self.manager()

        hsm = session.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
        self.target_hypothesis_id = hsm.target_hypothesis_id
        self.decoy_hypothesis_id = hsm.decoy_hypothesis_id

        session.close()

    def do_assignment(self):
        target = SpectrumAssigner(
            self.database_path, hypothesis_id=self.target_hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            scorer=self.scorer, score_parameters=self.score_parameters, n_processes=self.n_processes)
        target.start()

        decoy = SpectrumAssigner(
            self.database_path, hypothesis_id=self.decoy_hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            scorer=self.scorer, score_parameters=self.score_parameters,
            n_processes=self.n_processes)
        decoy.start()

        spectrum_tda = target_decoy.TargetDecoySpectrumMatchAnalyzer(
            self.database_path, target_hypothesis_id=self.target_hypothesis_id,
            decoy_hypothesis_id=self.decoy_hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id, score=self.scorer.score_name)
        spectrum_tda.start()

    def do_summarize(self):
        task = SummaryMatchBuilder(
            self.database_path, hypothesis_id=self.target_hypothesis_id,
            hypothesis_sample_match_id=self.hypothesis_sample_match_id,
            score_name=self.scorer.score_name, n_processes=self.n_processes)
        task.start()

        # task = SummaryMatchBuilder(
        #     self.database_path, hypothesis_id=self.decoy_hypothesis_id,
        #     hypothesis_sample_match_id=self.hypothesis_sample_match_id,
        #     score_name=self.scorer.score_name, n_processes=self.n_processes)
        # task.start()

    def start(self):
        self.do_assignment()
        self.do_summarize()
