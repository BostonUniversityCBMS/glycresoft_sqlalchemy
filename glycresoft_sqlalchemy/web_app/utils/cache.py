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
    def __init__(self):
        self.app_data = LRUDict()

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

    def cache_querycount(self, source, query=None, filter=None):
        source_id, source_type = source.id, source.__class__.__name__

        try:
            return self[(source_id, source_type, query, filter)]
        except KeyError:
            if query is not None:
                if isinstance(query, basestring):
                    query_obj = getattr(source, query)
                else:
                    query_obj = query
            value = self[(source_id, source_type, query, filter)] = filter(query_obj).count()
            return value

    @contextmanager
    def transaction(self):
        yield
        try:
            self.app_data.commit()
        except:
            pass


class CachablePartialFunction(partial):
    def __init__(self, func, *args, **kwargs):
        super(CachablePartialFunction, self).__init__(func, *args, **kwargs)
        self._hash = hash((func, tuple(args), frozenset(kwargs.items())))

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (self.func == other.func) and (self.args == other.args) and (self.keywords == other.keywords)


class MonosaccharideFilterSet(object):
    def __init__(self, constraints=None):
        if constraints is None:
            constraints = tuple()
        self.constraints = constraints

    def __eq__(self, other):
        return self.constraints == other.constraints

    def __hash__(self):
        return hash(self.constraints)

    def items(self):
        for constraint in self.constraints:
            yield constraint.monosaccharide, constraint

    @classmethod
    def fromdict(cls, filters):
        filters = [MonosaccharideFilter(monosaccharide=monosaccharide, **constraint)
                   for monosaccharide, constraint in filters.items()]
        filters = tuple(filters)
        return cls(filters)

    def __repr__(self):
        return "MonosaccharideFilterSet(%r)" % (self.constraints,)


class MonosaccharideFilter(object):
    def __init__(self, monosaccharide, minimum, maximum, include):
        self.monosaccharide = monosaccharide
        self.minimum = minimum
        self.maximum = maximum
        self.include = include

    def __eq__(self, other):
        return self.monosaccharide == other.monosaccharide and\
            self.minimum == other.minimum and\
            self.maximum == other.maximum and\
            self.include == other.include

    def __hash__(self):
        return hash((self.monosaccharide, self.minimum, self.maximum,
                     self.include))

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return "MonosaccharideFilter(monosaccharide=%r, minimum=%d, maximum=%d, include=%r)" % (
            self.monosaccharide, self.minimum, self.maximum, self.include)
