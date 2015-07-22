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

from glycresoft_sqlalchemy.data_model import DatabaseManager, Decon2LSPeakGroup, PeakGroupMatch, Hypothesis
from ..utils.collectiontools import groupby


def draw_chromatogram(peak_grouping, ax=None, set_bounds=True, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)
    intensity = peak_grouping.peak_data['intensities']
    scans = peak_grouping.peak_data['scan_ids']
    if len(intensity) != len(scans):
        # Don't try to collapse multiple peaks per scan yet
        print "Multiple peaks per scan detected"
        return ax
    color = kwargs.get('color')
    intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
    time = np.arange(0, scans.max() + 1)
    abundance_over_time = np.zeros_like(time)
    abundance_over_time[scans] = intensity
    ax.fill_between(time, abundance_over_time, alpha=0.5, color=color)
    ax.plot(time, abundance_over_time, color=color, alpha=0.7)
    if set_bounds:
        ax.set_xlim(0, time.max() + 200)
        ax.set_ylim(0, abundance_over_time.max() + 200)
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Relative Intensity")
    return ax


intensity_getter = operator.attrgetter('monoisotopic_intensity')
scan_id_getter = operator.attrgetter('scan_id')


def draw_chromatogram_peak_list(peak_list, ax=None, set_bounds=True, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)
    intensity = map(intensity_getter, peak_list)
    scans = map(scan_id_getter, peak_list)

    scan_groups = groupby(zip(scans, intensity), key_fn=operator.itemgetter(0))
    scans = []
    intensity = []
    for scan, peak_group in scan_groups.items():
        scans.append(scan)
        intensity.append(sum(*peak_group))

    intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
    time = np.arange(0, scans.max() + 1)
    abundance_over_time = np.zeros_like(time)
    abundance_over_time[scans] = intensity

    color = kwargs.get('color')

    ax.fill_between(time, abundance_over_time, alpha=0.5)
    ax.plot(time, abundance_over_time, color=color, alpha=0.7)

    if set_bounds:
        ax.set_xlim(0, time.max() + 200)
        ax.set_ylim(0, abundance_over_time.max() + 200)

    return ax


def get_chromatogram(peak_list):
    intensity = map(intensity_getter, peak_list)
    scans = map(scan_id_getter, peak_list)

    scan_groups = groupby(zip(scans, intensity), key_fn=operator.itemgetter(0))
    scans = []
    intensity = []
    for scan, peak_group in scan_groups.items():
        scans.append(scan)
        intensity.append(sum(*peak_group))

    intensity, scans = map(np.array, zip(*sorted(zip(intensity, scans), key=lambda x: x[1])))
    time = np.arange(0, scans.max() + 1)
    abundance_over_time = np.zeros_like(time)
    abundance_over_time[scans] = intensity

    return time, abundance_over_time
