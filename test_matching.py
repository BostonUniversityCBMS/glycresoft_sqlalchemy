import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")
from glycresoft_sqlalchemy.matching import matching

if __name__ == '__main__':
    job = matching.IonMatching("./datafiles/ResultOf20140918_01_isos.db", "./datafiles/20140918_01.yaml.db",
                               experiment_ids=None, n_processes=6)
    job.session.query(matching.GlycopeptideMatch).delete()
    job.session.query(matching.SpectrumMatch).delete()
    job.session.commit()
    job.run()
