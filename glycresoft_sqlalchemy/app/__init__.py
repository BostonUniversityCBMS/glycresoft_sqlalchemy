import sys
import logging
from contextlib import contextmanager

try:
    import matplotlib
    matplotlib.use("agg")
    logging.basicConfig(level=logging.DEBUG, filename='glycresoft-log', filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
    logging.captureWarnings(True)

    fmt = logging.Formatter(
        "%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s", "%H:%M:%S")
    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    logging.getLogger().addHandler(handler)

    warner = logging.getLogger('py.warnings')
    warner.setLevel("CRITICAL")

except Exception, e:
    logging.exception("Error, %r", e, exc_info=e)
    raise e


@contextmanager
def let(obj):
    try:
        yield obj
    finally:
        pass


def fail(message):
    print(message)
    sys.exit(1)
