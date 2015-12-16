from glycresoft_sqlalchemy.utils import Bundle
from sqlalchemy.ext.declarative import declarative_base

Base2 = declarative_base()


class _Namespace(Bundle):
    def __init__(self, *args, **kwargs):
        self.initialization_list = []
        Bundle.__init__(self, *args, **kwargs)

    def initializer(self, fn):
        self.initialization_list.append(fn)
        return fn

    def __str__(self):
        return "Namespace(%s)" % self.__dict__

Namespace = _Namespace()


class Hierarchy(dict):
    def __init__(self):
        dict.__init__(self)

    def references(self, referent):
        def wrapper(cls):
            self[referent] = cls
            return cls
        return wrapper
    register = references
