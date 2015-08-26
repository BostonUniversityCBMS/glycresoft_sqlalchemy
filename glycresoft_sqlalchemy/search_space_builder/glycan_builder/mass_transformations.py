from glypy.composition import composition_transform
from glypy.structure.monosaccharide import ReducedEnd
from glypy.structure import Substituent
import functools

# Mass Transform Functions


def derivatize(glycan, derivative):
    return composition_transform.derivatize(glycan, derivative)


def reduce(glycan, reducing_end_type):
    if reducing_end_type is not None:
        glycan.set_reducing_end(reducing_end_type())
    return glycan

derivatization_map = {
    None: None,
    True: None,
    False: None,
    "Permethylation": "methyl",
    "Acetylation": "acetyl",
}

reduction_map = {
    None: None,
    False: False,
    True: ReducedEnd,
    "Reduced": ReducedEnd,
    "Free": None,
    "Deutero-Reduced": lambda: (ReducedEnd("H[2]H"))
}


def get_derivatization(derivatization_method):
    try:
        method = derivatization_map[derivatization_method]
    except KeyError:
        try:
            method = Substituent(name=derivatization_method)
        except:
            method = Substituent(composition=derivatization_method)
    derivatizer = functools.partial(derivatize, derivatize=method)
    return derivatizer


def get_reduction(reduction_method):
    try:
        method = reduction_map[reduction_method]
    except KeyError:
        try:
            def method():
                return ReducedEnd(reduction_method)
            method()
        except:
            raise Exception("Could not resolve Reduction Method: %s" % reduction_method)
    reducer = functools.partial(reduce, reducing_end_type=method)
    return reducer
