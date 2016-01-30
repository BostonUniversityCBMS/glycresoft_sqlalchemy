import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass

from glycresoft_sqlalchemy.matching.glycopeptide.pipeline import GlycopeptideFragmentMatchingPipeline
from glycresoft_sqlalchemy.matching.glycopeptide.spectrum_assignment import SpectrumAssigner


def test_main():
    db_file_name = "./datafiles/scaling_complexity_test_data/large-dataset-comparison.db"  # "./datafiles/integrated_omics_simple.db"

    hypothesis_id = 5
    hypothesis_sample_match_id = 4

    job = SpectrumAssigner(db_file_name, hypothesis_id, hypothesis_sample_match_id)
    job.start()



if __name__ == '__main__':
    test_main()
