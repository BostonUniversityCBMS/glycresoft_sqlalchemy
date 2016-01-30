import operator
import logging
import random
import string
from hashlib import sha512


logger = logging.getLogger('wcom')


def logmethod(method):
    name = method.__name__

    def wrapper(*args, **kwargs):
        logger.info("%s called with (%r, %r)", name, args, kwargs)
        result = method(*args, **kwargs)
        logger.info("%s returned %r", name, result)
        return result
    return wrapper


SIMPLE_CHARS = string.ascii_letters + string.digits


def get_random_string(length=24):
    return ''.join(random.choice(SIMPLE_CHARS) for i in xrange(length))

random_string = get_random_string


def get_random_hash(length=24):
    hash = sha512()
    hash.update(get_random_string())
    return hash.hexdigest()[:length]


def intify(value, default=0):
    try:
        return int(value)
    except:
        return default


def logerror(f):
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception, e:
            logger.exception("An error occured in %s", f, exc_info=e)
            raise e
    return wrapper


def listify(value):
    if isinstance(value, (list, tuple)):
        return list(value)
    else:
        return [value]


def remove_empty_rows(*columns):
    lengths = map(len, columns)
    assert all(i == lengths[0] for i in lengths), "Not all columns are the same length"
    keep = []
    for i in range(lengths[0]):
        drop = False
        for column in columns:
            if column[i] in [None, "", " "]:
                drop = True
                break
        if not drop:
            keep.append(i)
    if len(keep) == 0:
        return [[] for i in lengths]
    op = operator.itemgetter(*keep)
    return [listify(op(column)) for column in columns]
