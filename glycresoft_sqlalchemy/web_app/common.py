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


def get_random_hash(length=24):
    hash = sha512()
    hash.update(get_random_string())
    return hash.hexdigest()[:length]
