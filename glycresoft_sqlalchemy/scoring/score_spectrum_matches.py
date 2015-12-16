import itertools
import functools
import multiprocessing
try:
    import logging
    logger = logging.getLogger('score_spectrum_matches')
except:
    pass
from ..data_model import (
    TheoreticalGlycopeptide, GlycopeptideMatch, PipelineModule, GlycopeptideSpectrumMatch, func, Protein)
from .simple_scoring_algorithm import evaluate, merge_ion_matches, split_ion_list


def score(matched, theoretical, **parameters):
    evaluate(matched, theoretical, **parameters)


class MatchScorer(PipelineModule):
    def __init__(self, database_path, scorer=None, score_parameters=None):
        if scorer is None:
            scorer = score
        self.manager = self.manager_type(database_path)
        self.score_parameters = score_parameters or {}
        self.scorer = scorer

    def run(self):
        session = self.manager.session()
        for i in session.query(GlycopeptideMatch.id):
            i = i[0]
            q = session.query(GlycopeptideMatch, TheoreticalGlycopeptide).filter(
                GlycopeptideMatch.id == i, GlycopeptideMatch.theoretical_glycopeptide_id == TheoreticalGlycopeptide.id)
            for match, theoretical in q:
                self.scorer(match, theoretical, **self.score_parameters)
                session.add(match)

        session.commit()


class SimpleSpectrumAssignment(PipelineModule):
    def __init__(self, database_path, hypothesis_id, hypothesis_sample_match_id,
                 scorer=None, score_parameters=None, n_processes=4):
        if scorer is None:
            scorer = score
        self.manager = self.manager_type(database_path)
        self.hypothesis_id = hypothesis_id
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.scorer = scorer
        self.score_parameters = score_parameters or {}
        self.n_processes = n_processes

    def prepare_task_fn(self):
        return functools.partial(
            score_on_limited_peak_matches,
            scorer=self.scorer,
            database_manager=self.manager, score_parameters=self.score_parameters)

    def stream_glycopeptide_match_ids(self):
        session = self.manager.session()
        for id, in session.query(GlycopeptideMatch.id).filter(
                GlycopeptideMatch.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis_id,
                GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id):
            yield id
        session.close()

    def run(self):
        session = self.manager.session()
        GSM = GlycopeptideSpectrumMatch
        logger.info("Marking all GlycopeptideSpectrumMatches")
        session.query(GSM).filter(
            GSM.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
            GSM.hypothesis_id == self.hypothesis_id).update(
            {"best_match": False}, synchronize_session=False)
        session.commit()
        logger.info("Calculating best matches.")
        accepted_peaks_map = session.query(
            GSM.scan_time, func.max(GSM.peaks_explained)).group_by(
            GSM.scan_time).filter(
            GSM.hypothesis_id == self.hypothesis_id,
            GSM.hypothesis_sample_match_id == self.hypothesis_sample_match_id)
        i = 0
        logger.info("Updating best_match flags")
        for scan_time, peaks_explained in accepted_peaks_map:
            i += 1
            session.query(GSM).filter(
                GSM.scan_time == scan_time,
                GSM.peaks_explained == peaks_explained,
                GSM.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                GSM.hypothesis_id == self.hypothesis_id).update(
                {"best_match": True}, synchronize_session=False)
            if i % 1000 == 0:
                logger.info("%d scans updated", i)
        session.commit()

        task_fn = self.prepare_task_fn()

        cntr = 0
        accum = []
        logger.info("Begin Scoring")
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for result in pool.imap_unordered(task_fn, self.stream_glycopeptide_match_ids()):
                accum.append(result)
                cntr += 1
                if cntr % 1000 == 0:
                    for a in accum:
                        session.merge(a)
                    logger.info("%d matches scored", cntr)
                    session.commit()
                    accum = []
            pool.terminate()
            del pool
        else:
            for id in self.stream_glycopeptide_match_ids():
                result = task_fn(id)
                accum.append(result)
                cntr += 1
                if cntr % 1000 == 0:
                    for a in accum:
                        session.merge(a)
                    logger.info("%d matches scored", cntr)
                    session.commit()
                    accum = []
        for a in accum:
            session.merge(a)
        session.commit()
        accum = []
        session.commit()
        session.close()


def score_on_limited_peak_matches(glycopeptide_match, database_manager, scorer, score_parameters):
    try:
        session = database_manager.session()
        with session.no_autoflush:
            glycopeptide_match = session.query(GlycopeptideMatch).get(glycopeptide_match)
            ion_matches = split_ion_list(merge_ion_matches(itertools.chain.from_iterable(itertools.chain.from_iterable(
                [gsm.peak_match_map.values() for gsm in glycopeptide_match.spectrum_matches
                 if gsm.best_match]))))

            for series, value in ion_matches.items():
                setattr(glycopeptide_match, series, value)
            theoretical_sequence = glycopeptide_match.theoretical_reference
            glycopeptide_match.ms2_score = 0.0
            scorer(glycopeptide_match, theoretical_sequence, **score_parameters)
            return glycopeptide_match
    except Exception, e:
        logger.exception("An error occurred processing %r", glycopeptide_match, exc_info=e)
        raise e
    finally:
        session.close()


def score_spectrum_match(glycopeptide_spectrum_match_id, database_manager, scorer, score_parameters):
    try:
        session = database_manager()
        glycopeptide_spectrum_match = session.query(glycopeptide_spectrum_match_id)
        match_like = glycopeptide_spectrum_match.as_match_like()
        theoretical_sequence = glycopeptide_spectrum_match.glycopeptide_match.theoretical_reference
        scorer(match_like, theoretical_sequence, **score_parameters)

    except Exception, e:
        logger.exception("An error occurred processing %r", glycopeptide_spectrum_match, exc_info=e)
        raise e
    finally:
        session.close()
