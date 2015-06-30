import logging
import time
from .data_model import DatabaseManager

logger = logging.getLogger("pipeline_module")


class Pipeline(object):
    def __init__(self, steps=None, index=0):
        if steps is None:
            steps = []
        self.steps = []
        self.index = 0

    def start(self, at=None):
        if at is None:
            at = self.index
        for i in range(at, len(self.steps)):
            step = self.steps[i]
            step.start()


class PipelineModule(object):
    manager_type = DatabaseManager

    def start(self, *args, **kwargs):
        self._begin(*args, **kwargs)
        try:
            out = self.run()
        except Exception, e:
            logger.exception("An error occurred: %r", e, exc_info=e)
            out = self.error = e
            self.status = -1
        else:
            self.status = 0
        self._end(*args, **kwargs)
        return out

    def _begin(self, *args, **kwargs):
        self.start_time = time.time()
        logger.info("Begin %s", self.__class__.__name__)

    def set_runner(self, callable):
        self.run = callable

    def _end(self, *args, **kwargs):
        self.end_time = time.time()
        logger.info("End %s", self.__class__.__name__)
        logger.info(self.summarize())

    def summarize(self):
        chunks = [
            "Started at %s." % time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self.start_time)),
            "Ended at %s." % time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self.end_time)),
            "Total time elapsed: %s" % time.strftime("%j:%H:%M:%S", time.gmtime(self.end_time - self.start_time)),
            "%s completed successfully." % self.__class__.__name__ if self.status == 0 else
            "%s failed with exit code %d, error message %r" % (self.__class__.__name__, self.status, self.error)
        ]
        return '\n'.join(chunks)

    @property
    def database_path(self):
        return self.manager.path