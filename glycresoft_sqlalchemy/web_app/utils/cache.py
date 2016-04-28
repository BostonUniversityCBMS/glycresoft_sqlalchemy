import operator

from contextlib import contextmanager
from functools import partial


class LRUDict(dict):
    def __init__(self, size=20, **kwargs):
        dict.__init__(self)
        self.counts = dict()
        self.size = size
        for key, value in kwargs.items():
            self[key] = value

    def __getitem__(self, key):
        value = dict.__getitem__(self, key)
        self.counts[key] += 1
        return value

    def __setitem__(self, key, value):
        if key in self:
            dict.__setitem__(self, key, value)
        else:
            if len(self) >= self.size:
                self.shrink()
            dict.__setitem__(self, key, value)
            self.counts[key] = 1

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        del self.counts[key]

    def shrink(self):
        size = self.size
        total_items = len(self)
        order = sorted(self.counts.items(), key=operator.itemgetter(1))
        i = 0
        while size <= total_items:
            self.popitem(order[i][0])
            self.counts.popitem(order[i][0])
            i += 1
            total_items -= 1


class ApplicationDataCache(object):
    def __init__(self, size=20):
        self.app_data = LRUDict(size)

    def __getitem__(self, key):
        return self.app_data[key]

    def __setitem__(self, key, value):
        self.app_data[key] = value

    def cache(self, obj, callpath, *args, **kwargs):
        if callable(callpath):
            fn = callpath
        else:
            getter = operator.attrgetter(callpath)
            fn = getter(obj)
        source_id, source_type = obj.id, obj.__class__.__name__

        try:
            value = self[source_id, source_type, callpath, tuple(args), frozenset(kwargs.items())]
            return value
        except KeyError:
            value = fn(*args, **kwargs)
            self[source_id, source_type, callpath, tuple(args), frozenset(kwargs.items())] = value
            return value

    @contextmanager
    def transaction(self):
        yield
        try:
            self.app_data.commit()
        except:
            pass

    def invalidate(self, key):
        del self.app_data[key]


class CachablePartialFunction(partial):
    def __init__(self, func, *args, **kwargs):
        super(CachablePartialFunction, self).__init__(func, *args, **kwargs)
        self._hash = hash((func, tuple(args), frozenset(kwargs.items())))

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (self.func == other.func) and (self.args == other.args) and (self.keywords == other.keywords)
