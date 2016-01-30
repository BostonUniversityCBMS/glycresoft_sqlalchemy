from collections import defaultdict
import itertools

from glypy.utils.multimap import MultiMap

from .vendor import sqlitedict


try:
    range = xrange
except Exception, e:
    pass

try:
    import cPickle as pickle
except:
    import pickle as pickle


class DefaultSqliteDict(sqlitedict.SqliteDict):
    def __init__(self, *args, **kwargs):
        factory_fn = kwargs.pop("default", None)
        sqlitedict.SqliteDict.__init__(self, *args, **kwargs)
        self.default_factory = factory_fn

    def __getitem__(self, key):
        try:
            return sqlitedict.SqliteDict.__getitem__(self, key)
        except KeyError, e:
            if self.default_factory is not None:
                self.__setitem__(key, self.default_factory())
                return self[key]
            else:
                raise e


class SqliteSet(object):
    def __init__(self, iterable=()):
        self.store = sqlitedict.open()
        self.store.clear()
        for x in iterable:
            self.add(x)

    def add(self, obj):
        value = pickle.dumps(obj, -1)
        self.store[value] = 1

    def update(self, objs):
        for obj in objs:
            value = pickle.dumps(obj, -1)
            self.store[value] = 1

    def __contains__(self, obj):
        value = pickle.dumps(obj, -1)
        return value in self.store

    def __iter__(self):
        for x in self.store:
            yield pickle.loads(x)

    def __repr__(self):
        return '{%s}' % ', '.join(map(str, self.store.keys()))


def _identity(i):
    return i


def groupby(ungrouped_list, key_fn=_identity, transform_fn=_identity, kind=list):
    if kind == list:
        return groupby_list(ungrouped_list, key_fn, transform_fn)
    elif kind == set:
        return groupby_set(ungrouped_list, key_fn, transform_fn)
    else:
        raise TypeError("Collection {kind} is not defined".format(kind=kind))


def groupby_list(ungrouped_list, key_fn, transform_fn):
    groups = defaultdict(list)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].append(transform_fn(item))
    return groups


def groupby_set(ungrouped_list, key_fn, transform_fn):
    groups = defaultdict(set)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].add(transform_fn(item))
    return groups


def _repeat_label(group_size):
    label = 0
    while(True):
        for value in itertools.repeat(label, group_size):
            yield value
        label += 1


def chunk_iterator(iterator, chunk_size):
    labeler = _repeat_label(chunk_size)
    for label, chunk in itertools.groupby(iterator, lambda x: labeler.next()):
        yield list(chunk)


def chunk_iterator2(iterator, chunk_size):
    current_chunk = []
    i = 0
    for o in iterator:
        current_chunk.append(o)
        i += 1
        if i > chunk_size:
            yield current_chunk
            current_chunk = []
            i = 0
    yield current_chunk


def invert_index(list_of_lists, key_fn=_identity):
    groups = defaultdict(list)
    for chunk in list_of_lists:
        keys = set()
        for item in chunk:
            keys.add(key_fn(item))
        for key in keys:
            groups[key].append(chunk)

    return groups


def descending_combination_counter(counter):
    keys = counter.keys()
    count_ranges = map(lambda x: range(x + 1), counter.values())
    for combination in itertools.product(*count_ranges):
        yield dict(zip(keys, combination))


def flatten(iterable):
    return tuple(itertools.chain.from_iterable(iterable))


class decoratordict(dict):
    def __call__(self, key):
        def wrapper(f):
            self[key] = f
            return f
        return wrapper
