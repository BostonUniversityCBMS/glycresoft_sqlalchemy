import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.scoring import pair_counting, simple_scoring_algorithm
from glycresoft_sqlalchemy.scoring import rescore


def testmain():
    db_file_name = "./datafiles/integrated_omics_simple.db"
    # frequency_counter = pair_counting.pickle.load(open('datafiles/Phil-82-Training-Data/pair_counts.pkl'))
    # scorer = pair_counting.FrequencyScorer(frequency_counter)

    # scorer = simple_scoring_algorithm.SimpleSpectrumScorer()
    scorer = simple_scoring_algorithm.SimpleSpectrumScorer2()

    job = rescore.RescoreSpectrumHypothesisSampleMatch(db_file_name, hypothesis_sample_match_id=2, scorer=scorer)
    job.start()


if __name__ == '__main__':
    testmain()
