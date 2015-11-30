import os
import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy import data_model
from glycresoft_sqlalchemy.search_space_builder import exact_search_space_builder, pooling_make_decoys
from glycresoft_sqlalchemy.matching import matching
from glycresoft_sqlalchemy.scoring import target_decoy, score_spectrum_matches
from glycresoft_sqlalchemy.search_space_builder import integrated_omics


def test_main():
    db_file_name = "./datafiles/integrated_omics_simple.db"
    # try:
    #     os.remove(db_file_name)
    #     pass
    # except:
    #     pass

    hypothesis_sample_match_id = 1

    job = exact_search_space_builder.ExactSearchSpaceBuilder.from_hypothesis_sample_match(
        db_file_name, hypothesis_sample_match_id, 4)
    hypothesis_id = job.start()


if __name__ == '__main__':
    test_main()
