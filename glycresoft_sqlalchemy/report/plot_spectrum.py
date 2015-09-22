try:   # pragma: no cover
    from cStringIO import StringIO
except:  # pragma: no cover
    try:
        from StringIO import StringIO
    except:
        from io import StringIO
try:  # pragma: no cover
    from lxml import etree as ET
except ImportError:  # pragma: no cover
    try:
        from xml.etree import cElementTree as ET
    except:
        from xml.etree import ElementTree as ET

from itertools import cycle
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib.path import Path
from matplotlib.colors import cnames, hex2color, Normalize
from matplotlib import cm as colormap

import operator
import numpy as np

from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Decon2LSPeakGroup, PeakGroupMatch,
    Hypothesis, Decon2LSPeak)

from ..spectra.spectrum_model import mass_charge_ratio


mz_getter = operator.attrgetter("mass_charge_ratio")


def plot_scan(scan, ax=None, logscale=False, color='black', **kwargs):
    peaks = list(scan)
    if peaks[0].charge is not None:
        for p in peaks:
            p.mass_charge_ratio = mass_charge_ratio(p.neutral_mass, p.charge)
    else:
        for p in peaks:
            p.mass_charge_ratio = p.neutral_mass

    peaks.sort(key=mz_getter)

    if ax is None:
        figure, ax = plt.subplots(1)
        label_plot(ax)
    mzs = [p.mass_charge_ratio for p in peaks]
    intensities = [p.intensity for p in peaks]
    if logscale:
        intensities = np.log(intensities)

    kwargs.setdefault("width", 0.1)
    kwargs.setdefault("edgecolor", color)
    ax.bar(mzs, intensities, **kwargs)
    return ax


def label_plot(ax, **kwargs):
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative Intensity")
    ax.autoscale()
