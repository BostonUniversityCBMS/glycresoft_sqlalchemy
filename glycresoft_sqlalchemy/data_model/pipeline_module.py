import types
import logging
import time
import datetime
import pprint

from sqlalchemy import (PickleType, Numeric, Unicode, Table, DateTime, func,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)
from sqlalchemy.orm import relationship

from .connection import DatabaseManager
from .base import Base2 as Base

logger = logging.getLogger("pipeline_module")


debug = True


class User(Base):
    __tablename__ = "User"

    id = Column(Integer, primary_key=True)
    email_address = Column(Unicode(128))
    password = Column(Unicode(128))
    messages = relationship("Message", lazy='dynamic')


class Message(Base):
    __tablename__ = "Message"

    id = Column(Integer, primary_key=True)
    owner = Column(Integer, ForeignKey(User.id), index=True)
    message_type = Column(Unicode(64))
    contents = Column(PickleType)
    created = Column(DateTime, index=True, default=func.now())

    def __repr__(self):
        return "<{s.__class__.__name__}:{s.message_type} {s.contents} @ {s.created}>".format(s=self)


class JobState(Base):
    __tablename__ = "JobState"

    id = Column(Integer, primary_key=True)
    name = Column(Unicode(128))
    state = Column(Unicode(128))


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


class GeneratorWithCallback(object):
    def __init__(self, generator, callback):
        self.generator = generator
        self.callback = callback

    def __iter__(self):
        for item in self.generator:
            yield item
        self.callback()


class PipelineModule(object):
    '''
    Represent a single step in a program pipeline. This is a base class for all
    other pipeline steps. It provides basic state tracking information about the
    task run, and sets up some common logging events.

    Attributes
    ----------
    manager_type: DatabaseManager
        A subclass of :class:`DatabaseManager` which provides :class:`Session` objects
        from a bound database. Does not contain any open connections itself, so it may
        be safely shared with worker processes.
    start_time: float
        The instant the task began, as returned by :func:`time.time`
    end_time: float
        The instant the task ended, as returned by :func:`time.time`
    status: int
        An equivalent to a return code
    error: Exception or None
        Any Exception that was unhandled by the task and caused it to terminate.
        The exception does not propagate, but is logged, and a status code of -1
        is set.
    '''
    manager_type = DatabaseManager
    error = None
    raise_on_error = True

    def start(self, *args, **kwargs):
        self._begin(*args, **kwargs)
        try:
            out = self.run()
        except KeyboardInterrupt if debug else Exception, e:
            logger.exception("An error occurred: %r", e, exc_info=e)
            out = self.error = e
            self.status = -1
            if self.raise_on_error:
                raise e
        else:
            self.status = 0
        if isinstance(out, types.GeneratorType):
            return GeneratorWithCallback(out, lambda: self._end(*args, **kwargs))
        else:
            self._end(*args, **kwargs)
            return out

    def _begin(self, verbose=True, *args, **kwargs):
        self.start_time = datetime.datetime.now()
        if verbose:
            logger.info("Begin %s\n%s\n", self.__class__.__name__, pprint.pformat(self.__dict__))

    def set_runner(self, callable):
        self.run = callable

    def inform(self, *args, **kwargs):
        now = datetime.datetime.now() - self.start_time
        logger.info(
            "%s: %s time elapsed.", self.__class__.__name__, str(now))
        logger.info(*args, **kwargs)

    def _end(self, verbose=True, *args, **kwargs):
        self.end_time = datetime.datetime.now()
        if verbose:
            logger.info("End %s", self.__class__.__name__)
            logger.info(self.summarize())

    def summarize(self):
        chunks = [
            "Started at %s." % self.start_time,
            "Ended at %s." % self.end_time,
            "Total time elapsed: %s" % (self.end_time - self.start_time),
            "%s completed successfully." % self.__class__.__name__ if self.status == 0 else
            "%s failed with exit code %d, error message %r" % (self.__class__.__name__, self.status, self.error)
        ]
        return '\n'.join(chunks)

    @property
    def database_path(self):
        return self.manager.path
