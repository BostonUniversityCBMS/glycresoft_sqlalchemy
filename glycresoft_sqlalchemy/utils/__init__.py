import os
import logging
import urlparse
import urllib

try:
    import cPickle as pickle
except:
    import pickle


from glypy.utils import classproperty
from glypy.utils.enum import Enum

from .vendor import appdir, pager, sqlitedict


def get_scale():
    return int(os.getenv('GLYCRESOFT_SCALE', 1))


def simple_repr(self):  # pragma: no cover
    template = "{self.__class__.__name__}({d})"
    d = [
        "%s=%r" % (k, v) if v is not self else "(...)" for k, v in sorted(
            self.__dict__.items(), key=lambda x: x[0])
        if not k.startswith("_") and not callable(v)]
    return template.format(self=self, d=', '.join(d))


def path2url(path):
    return urlparse.urljoin(
      'file:', urllib.pathname2url(path))


def file_uri(path):
    return path2url(os.path.abspath(path))


class ItemsAsAttributes(object):

    __getitem__ = object.__getattribute__
    __setitem__ = object.__setattr__

    __repr__ = simple_repr


class Bundle(ItemsAsAttributes):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            self[k] = v


def setup_logging():
    try:
        logging.basicConfig(level=logging.DEBUG, filename='glycresoft-log', filemode='w',
                            format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                            datefmt="%H:%M:%S")
        logging.captureWarnings(True)

        fmt = logging.Formatter(
            "%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s", "%H:%M:%S")
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logging.getLogger().addHandler(handler)

        warner = logging.getLogger('py.warnings')
        warner.setLevel("CRITICAL")

    except Exception, e:
        logging.exception("Error, %r", e, exc_info=e)
        raise e


def setup_logging_child():
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                            datefmt="%H:%M:%S")
        logging.captureWarnings(True)

        fmt = logging.Formatter(
            "%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s", "%H:%M:%S")
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logging.getLogger().addHandler(handler)

        warner = logging.getLogger('py.warnings')
        warner.setLevel("CRITICAL")

    except Exception, e:
        logging.exception("Error, %r", e, exc_info=e)
        raise e
