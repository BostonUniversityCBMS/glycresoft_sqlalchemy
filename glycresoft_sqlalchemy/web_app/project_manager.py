import os
import shutil
import logging
from glycresoft_sqlalchemy.data_model import DatabaseManager, MSMSSqlDB, Hypothesis, SampleRun, HypothesisSampleMatch

logger = logging.getLogger("project_manager")

path = os.path


def logmethod(method):
    name = method.__name__
    def wrapper(*args, **kwargs):
        logger.info("%s called with (%r, %r)", name, args, kwargs)
        result = method(*args, **kwargs)
        logger.info("%s returned %r", name, result)
        return result
    return wrapper


class ProjectManager(DatabaseManager):
    def __init__(self, database_path, sample_dir=None, results_dir=None):
        super(ProjectManager, self).__init__(database_path)
        if sample_dir is None:
            sample_dir = path.join(path.dirname(database_path), 'sample_dir')
        if results_dir is None:
            results_dir = path.join(path.dirname(database_path), 'results_dir')
        self.sample_dir = sample_dir
        self.results_dir = results_dir
        self._ensure_paths_exist()
        self.sample_submanagers = []
        for sample in os.listdir(self.sample_dir):
            manager = MSMSSqlDB(sample)
            try:
                manager.connect()
                self.sample_submanagers.append(manager)
            except:
                pass

    @logmethod
    def _ensure_paths_exist(self):
        try:
            os.makedirs(self.sample_dir)
        except Exception, e:
            # logger.exception("An error occurred", exc_info=e)
            pass
        try:
            os.makedirs(self.results_dir)
        except Exception, e:
            # logger.exception("An error occurred", exc_info=e)
            pass

    @logmethod
    def find_sample(self, name):
        for manager in self.sample_submanagers:
            match = manager.session().query(SampleRun).filter(SampleRun.name == name).first()
            if match is None:
                continue
            else:
                return match, manager
        raise KeyError(name)

    @logmethod
    def add_sample_path(self, fname):
        name = path.basename(fname)
        dest = path.join(self.sample_dir, name)
        shutil.copy(fname, dest)
        self.sample_submanagers.append(MSMSSqlDB(dest))

    def get_sample_path(self, name):
        return path.join(self.sample_dir, name)

    @logmethod
    def add_sample_stream(self, fname, stream):
        dest = path.join(self.sample_dir, fname)
        with open(dest, 'wb') as fh:
            for line in stream:
                fh.write(line)
        self.sample_submanagers.append(MSMSSqlDB(dest))

    @logmethod
    def hypotheses(self):
        session = self.session()
        return session.query(Hypothesis)

    @logmethod
    def samples(self):
        results = []
        for manager in self.sample_submanagers:
            results.extend(manager.session().query(SampleRun).all())
        return results
