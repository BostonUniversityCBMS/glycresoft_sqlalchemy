from glycresoft_sqlalchemy.search_space_builder import make_decoys

import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")


def test_main():
    job = make_decoys.DecoySearchSpaceBuilder("./datafiles/ResultOf20140918_01_isos.db", hypothesis_ids=[1], n_processes=4)
    job.session.query(make_decoys.Hypothesis.id).filter(make_decoys.Hypothesis.name.like("decoy")).delete('fetch')
    job.session.commit()
    job.start()

if __name__ == '__main__':
    test_main()
