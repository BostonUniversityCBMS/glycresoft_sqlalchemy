import os
import json
import logging
from ConfigParser import ConfigParser

# from ..ms import constants as ms_constants
# from ..structure import constants as structure_constants
# from ..prediction_tools import constants as prediction_constants

logger = logging.getLogger("config_loader")

namespace_lookup = {
    # "ms": ms_constants,
    # "spectra": ms_constants,
    # "structure": structure_constants,
    # "prediction": prediction_constants
}


def inject_option(namespace, key, value):
    if hasattr(namespace, key):
        setattr(namespace, key, value)
    else:
        d = namespace.__dict__
        key_map = {k.lower(): k for k in d}
        if key in key_map:
            setattr(namespace, key_map[key], value)
        else:
            raise AttributeError("Could not find {key} in {namespace}".format(
                key=key, namespace=namespace))


def inject_options(namespace, options):
    for option, value in options.items():
        inject_option(namespace, option, value)


def read_json(path):
    conf = json.load(open(path))
    for section, options in conf.items():
        try:
            namespace = namespace_lookup[section.lower()]
            inject_options(namespace, options)
        except KeyError, e:
            logger.error("An error occurred during configuration: Unknown Section: %s", section, exc_info=e)
            raise


def convert(val):
    rval = val
    try:
        rval = int(val)
    except:
        try:
            rval = float(val)
        except:
            pass
    if val.lower() in {"yes", "true"}:
        rval = True
    elif val.lower() in {"no", "false"}:
        rval = False
    return rval


def read_ini(path):
    reader = ConfigParser()
    reader.read(path)
    for section in reader.sections():
        options = {k: convert(v) for k, v in reader.items(section)}
        try:
            namespace = namespace_lookup[section.lower()]
            inject_options(namespace, options)
        except KeyError, e:
            logger.error("An error occurred during configuration: Unknown Section: %s", section, exc_info=e)
            raise


def gather():
    cfg = {}
    for section in ["spectra", "structure", "prediction"]:
        ns = namespace_lookup[section]
        cfg[section] = ns.__dict__
    return cfg


def write_json(path):
    json.dump(gather(), open(path, 'w'))


def write_ini(path):
    options = gather()
    writer = ConfigParser()
    for section in options:
        writer.add_section(section.lower())
        for name, value in options[section].items():
            writer.set(section.lower(), name.lower(), str(value))
    writer.write(open(path, 'w'))


def load(path):
    try:
        read_json(path)
    except:
        read_ini(path)
