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
from matplotlib.colors import cnames, hex2color

import numpy as np

from glycresoft_sqlalchemy.data_model import DatabaseManager, Decon2LSPeakGroup, PeakGroupMatch, Hypothesis


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
    points = np.zeros_like(time)
    points[scans] = intensity
    ax.fill_between(time, points, alpha=0.5, color=color)
    ax.plot(time, points, color=color, alpha=0.7)
    if set_bounds:
        ax.set_xlim(0, time.max() + 200)
        ax.set_ylim(0, points.max() + 200)
    return ax
