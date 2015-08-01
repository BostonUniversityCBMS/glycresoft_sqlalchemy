from itertools import cycle
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib.path import Path
from matplotlib.colors import cnames, hex2color


def lighten(rgb, factor=0.25):
    '''Given a triplet of rgb values, lighten the color by `factor`%'''
    factor += 1
    return [min(c * factor, 1) for c in rgb]


def darken(rgb, factor=0.25):
    '''Given a triplet of rgb values, darken the color by `factor`%'''
    factor = 1 - factor
    return [(c * factor) for c in rgb]


colors = cycle([hex2color(cnames[name]) for name in ("red", "blue", "yellow", "purple", "navy", "grey")])


color_name_map = {
    "HexNAc": hex2color(cnames["steelblue"])
}


def get_color(name):
    """Given a name, find the color mapped to that name, or
    select the next color from the `colors` generator and assign
    it to the name and return the new color.

    Parameters
    ----------
    name : object
        Any hashable object, usually a string

    Returns
    -------
    tuple: RGB triplet
    """
    try:
        return color_name_map[name]
    except KeyError:
        color_name_map[name] = colors.next()
        return color_name_map[name]
