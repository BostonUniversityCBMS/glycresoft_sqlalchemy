from glycresoft_sqlalchemy.search_space_builder import make_decoys

import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")

if __name__ == '__main__':
    job = make_decoys.DecoySearchSpaceBuilder("./datafiles/ResultOf20140918_01_isos.db", experiment_ids=[1], n_processes=4)
    job.session.query(make_decoys.Experiment.id).filter(make_decoys.Experiment.name.like("decoy")).delete('fetch')
    job.session.commit()
    job.run()
