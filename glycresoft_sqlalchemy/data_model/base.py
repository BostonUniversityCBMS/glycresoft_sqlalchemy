from glycresoft_sqlalchemy.utils import Bundle
from sqlalchemy.ext.declarative import declarative_base
from glycresoft_sqlalchemy.utils import database_utils

Base = declarative_base()
Base2 = Base


def _base_eq(self, other):
    try:
        if other is None:
            return False
        return self.id == other.id and self.__class__ == other.__class__
    except AttributeError:
        return False


def _base_hash(self):
    return hash((self.id, self.__class__))


def _base_get_uri(self):
    return database_utils.get_uri_for_instance(self)


Base.__eq__ = _base_eq
Base.__hash__ = _base_hash
Base.uri = _base_get_uri


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


def slurp(session, model, ids, flatten=True):
    if flatten:
        ids = [j for i in ids for j in i]
    total = len(ids)
    last = 0
    step = 100
    results = []
    while last < total:
        results.extend(session.query(model).filter(
            model.id.in_(ids[last:last + step])))
        last += step
    return results
