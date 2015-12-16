import logging
from functools import partial
from multiprocessing import Pool
import itertools

from glycresoft_sqlalchemy.data_model import (
    HypothesisSampleMatch, PipelineModule, GlycopeptideMatch,
    GlycopeptideSpectrumMatch)

from glycresoft_sqlalchemy.scoring import target_decoy
from glycresoft_sqlalchemy.scoring import score_spectrum_matches, simple_scoring_algorithm

logger = logging.getLogger("rescore")


def do_glycopeptide_match_scoring(glycopeptide_match_ids, database_manager, scorer, score_parameters):
    session = database_manager()
    collection = []
    try:
        for glycopeptide_match_id in glycopeptide_match_ids:
            match = session.query(GlycopeptideMatch).get(glycopeptide_match_id)
            theoretical = match.theoretical_reference
            scorer(match, theoretical, **score_parameters)
            collection.append(match)
        session.add_all(collection)
        session.commit()
        return len(collection)
    except Exception, e:
        logger.exception("An error occurred processing %r", match, exc_info=e)
        raise e
    finally:
        session.close()


def do_glycopeptide_spectrum_match_scoring(glycopeptide_spectrum_match_ids, database_manager, scorer, score_parameters):
    session = database_manager()
    collection = []
    try:
        for glycopeptide_spectrum_match_id in glycopeptide_spectrum_match_ids:
            spectrum_match = session.query(GlycopeptideSpectrumMatch).get(glycopeptide_spectrum_match_id)
            theoretical = spectrum_match.glycopeptide_match.theoretical_reference
            match = spectrum_match.as_match_like()
            score = scorer(match, theoretical, **score_parameters)
            collection.append(score)
        session.add_all(collection)
        session.commit()

        return len(collection)
    except Exception, e:
        logger.exception("An error occurred processing %r", match, exc_info=e)
        raise e
    finally:
        session.close()


def yield_ids(session, base_query, chunk_size=200):
    base_query = tuple(base_query)
    last = 0
    final = len(base_query)
    while 1:
        next_chunk = base_query[last:(last + chunk_size)]
        if last <= final:
            yield next_chunk
            last += chunk_size
        else:
            break


class RescoreHypothesisSampleMatch(PipelineModule):
    def __init__(self, database_path, hypothesis_sample_match_id, scorer=None,
                 score_parameters=None, n_processes=4, **kwargs):
        if scorer is None:
            scorer = score_spectrum_matches.score
        if score_parameters is None:
            score_parameters = {}
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.scorer = scorer
        self.score_parameters = score_parameters
        self.n_processes = n_processes
        self.options = kwargs

    def stream_ids(self, is_decoy=False):
        session = self.manager()
        if is_decoy:
            q = session.query(GlycopeptideMatch.id).filter(
                GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                ~GlycopeptideMatch.is_not_decoy())
        else:
            q = session.query(GlycopeptideMatch.id).filter(
                GlycopeptideMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                GlycopeptideMatch.is_not_decoy())
        for chunk in yield_ids(session, q):
            yield chunk
        session.close()

    def prepare_scoring_fn(self):
        return partial(
            do_glycopeptide_match_scoring, database_manager=self.manager,
            scorer=self.scorer, score_parameters=self.score_parameters)

    def do_score(self, session, is_decoy=False):
        fn = self.prepare_scoring_fn()
        cntr = 0
        last = 0
        if self.n_processes > 1:
            pool = Pool(self.n_processes)
            for res in pool.imap_unordered(fn, self.stream_ids(is_decoy=is_decoy), 1):
                cntr += res
                if cntr - last > 1000:
                    logger.info("Scored %d matches", cntr)
                    last = cntr
            pool.close()
            pool.terminate()
        else:
            for res in itertools.imap(fn, self.stream_ids(is_decoy=is_decoy)):
                cntr += res
                if cntr - last > 1000:
                    logger.info("Scored %d matches", cntr)
                    last = cntr

    def run(self):
        session = self.manager()
        hsm = session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
        self.do_score(session, is_decoy=False)
        self.do_score(session, is_decoy=True)

        job = target_decoy.TargetDecoyAnalyzer.from_hypothesis_sample_match(self.database_path, hsm)
        job.start()


class RescoreSpectrumHypothesisSampleMatch(PipelineModule):
    def __init__(self, database_path, hypothesis_sample_match_id, scorer=None,
                 score_parameters=None, n_processes=4, **kwargs):
        if scorer is None:
            scorer = simple_scoring_algorithm.SimpleSpectrumScorer()
        if score_parameters is None:
            score_parameters = {}
        self.manager = self.manager_type(database_path)
        self.hypothesis_sample_match_id = hypothesis_sample_match_id
        self.scorer = scorer
        self.score_parameters = score_parameters
        self.n_processes = n_processes
        self.options = kwargs

    def stream_ids(self, is_decoy=False):
        session = self.manager()
        if is_decoy:
            q = session.query(GlycopeptideSpectrumMatch.id).filter(
                (GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id) & (
                    GlycopeptideSpectrumMatch.is_decoy()))
        else:
            q = session.query(GlycopeptideSpectrumMatch.id).filter(
                GlycopeptideSpectrumMatch.hypothesis_sample_match_id == self.hypothesis_sample_match_id,
                GlycopeptideSpectrumMatch.is_not_decoy())
        for chunk in yield_ids(session, q):
            yield chunk
        session.close()

    def prepare_scoring_fn(self):
        return partial(
            do_glycopeptide_spectrum_match_scoring, database_manager=self.manager,
            scorer=self.scorer, score_parameters=self.score_parameters)

    def do_score(self, session, is_decoy=False):
        fn = self.prepare_scoring_fn()
        cntr = 0
        last = 0
        if self.n_processes > 1:
            pool = Pool(self.n_processes)
            for res in pool.imap_unordered(fn, self.stream_ids(is_decoy=is_decoy), 1):
                cntr += res
                if cntr - last > 1000:
                    logger.info("Scored %d matches", cntr)
                    last = cntr
            pool.close()
            pool.terminate()
        else:
            for res in itertools.imap(fn, self.stream_ids(is_decoy=is_decoy)):
                cntr += res
                if cntr - last > 1000:
                    logger.info("Scored %d matches", cntr)
                    last = cntr

    def run(self):
        session = self.manager()
        hsm = session.query(HypothesisSampleMatch).get(self.hypothesis_sample_match_id)
        self.do_score(session, is_decoy=False)
        self.do_score(session, is_decoy=True)
