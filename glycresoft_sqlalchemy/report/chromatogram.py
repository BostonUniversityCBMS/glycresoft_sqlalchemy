import argparse
import os
import re

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

from ..utils.collectiontools import groupby


intensity_getter = operator.attrgetter('monoisotopic_intensity')
scan_time_getter = operator.attrgetter('scan.time')


def draw_chromatogram(peak_list, ax=None, set_bounds=True, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)
    if isinstance(peak_list, (Decon2LSPeakGroup, PeakGroupMatch)):
        time, abundance_over_time = get_chromatogram_peak_data(peak_list.peak_data)
    elif isinstance(peak_list[0], Decon2LSPeak):
        time, abundance_over_time = get_chromatogram(peak_list)
    else:
        time, abundance_over_time = map(np.array, peak_list)

    color = kwargs.pop('color', "blue")
    alpha = kwargs.pop("alpha", 0.4)
    linewidth = kwargs.pop("linewidth", 0.1)

    ax.fill_between(time, abundance_over_time, alpha=alpha)
    ax.plot(time, abundance_over_time, color=color, alpha=alpha + 0.2, linewidth=linewidth)

    if set_bounds:
        ax.set_xlim(0, time.max() + 200)
        ax.set_ylim(0, abundance_over_time.max() + 200)
        ax.set_ylabel("Relative Intensity")
        ax.set_xlabel("Scan Number")
        ax.xaxis.set_tick_params(width=0)
        ax.set_title("Extracted Chromatogram")

    return ax


def get_chromatogram(peak_list):
    intensity = map(intensity_getter, peak_list)
    scans = map(scan_time_getter, peak_list)

    scan_groups = groupby(zip(scans, intensity), key_fn=operator.itemgetter(0), transform_fn=operator.itemgetter(1))
    scans = []
    intensity = []
    for scan, peak_group in scan_groups.items():
        scans.append(scan)
        intensity.append(sum(peak_group))

    intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
    time = np.arange(0, scans.max() + 1)
    abundance_over_time = np.zeros_like(time)
    abundance_over_time[scans] = intensity

    return time, abundance_over_time


def get_chromatogram_peak_data(peak_data):
    scans = peak_data['scan_times']
    intensity = peak_data['intensities']

    scan_groups = groupby(zip(scans, intensity), key_fn=operator.itemgetter(0), transform_fn=operator.itemgetter(1))
    scans = []
    intensity = []
    for scan, peak_group in scan_groups.items():
        scans.append(scan)
        intensity.append(sum(peak_group))

    intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
    time = np.arange(0, scans.max() + 1)
    abundance_over_time = np.zeros_like(time)
    abundance_over_time[scans] = intensity

    return time, abundance_over_time
