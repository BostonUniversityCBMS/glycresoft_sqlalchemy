import logging

logger = logging.getLogger('wcom')


def logmethod(method):
    name = method.__name__

    def wrapper(*args, **kwargs):
        logger.info("%s called with (%r, %r)", name, args, kwargs)
        result = method(*args, **kwargs)
        logger.info("%s returned %r", name, result)
        return result
    return wrapper
