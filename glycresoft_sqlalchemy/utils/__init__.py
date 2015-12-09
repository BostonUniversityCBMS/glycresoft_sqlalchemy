try:
    import cPickle as pickle
except:
    import pickle

from glypy.utils import classproperty


def simple_repr(self):  # pragma: no cover
    template = "{self.__class__.__name__}({d})"
    d = {"%s=%r" % (k, v) for k, v in sorted(self.__dict__.items(), key=lambda x: x[0]) if not k.startswith("_")}
    return template.format(self=self, d=', '.join(d))


class ItemsAsAttributes(object):

    __getitem__ = object.__getattribute__
    __setitem__ = object.__setattr__

    __repr__ = simple_repr


class Bundle(ItemsAsAttributes):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            self[k] = v
