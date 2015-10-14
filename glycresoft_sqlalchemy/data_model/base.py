from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class _Namespace(object):
    def __str__(self):
        return "Namespace(%s)" % self.__dict__

Namespace = _Namespace()
Namespace.initialization_list = []


class Hierarchy(dict):
    def __init__(self):
        dict.__init__(self)

    def references(self, referent):
        def wrapper(cls):
            self[referent] = cls
            return cls
        return wrapper
    register = references
