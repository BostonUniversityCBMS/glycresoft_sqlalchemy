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

from glycresoft_sqlalchemy.data_model import DatabaseManager, GlycopeptideMatch, Protein, Hypothesis
from glycresoft_ms2_classification.structure import sequence


def clean_file_name(file_name):
    name = re.sub(r"\|", "_", file_name)
    return name


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


def span_overlap(a, b):
    return a.spans(b.start_position) or a.spans(b.end_position) or\
           b.spans(a.start_position) or b.spans(a.end_position)


def layout_layers(gpms):
    layers = [[]]
    gpms.sort(key=lambda x: x.ms2_score, reverse=True)
    for gpm in gpms:
        placed = False
        for layer in layers:
            collision = False
            for member in layer:
                if span_overlap(gpm, member):
                    collision = True
                    break
            if not collision:
                layer.append(gpm)
                placed = True
                break
        if not placed:
            layers.append([gpm])
    return layers


class IDMapper(dict):
    def __init__(self):
        dict.__init__(self)
        self.counter = 0

    def add(self, key, value, meta):
        label = key % self.counter
        value.set_gid(label)
        self[label] = meta
        self.counter += 1
        return label


def draw_layers(layers, protein):
    figure, ax = plt.subplots(1, 1)
    i = 0
    id_mapper = IDMapper()
    for i, aa in enumerate(protein.protein_sequence):
        ax.text(-0.2 + i, 5.5, aa, family="monospace", fontsize=4)
    rect = mpatches.Rectangle((0, 0), i, 5, facecolor='red', alpha=0.5)
    id_mapper.add("main-sequence-%d", rect, {'main-sequence': True})
    ax.add_patch(rect)
    if i != 0:
        defer_x_bounds = False
        ax.set_xlim(-10, i + 10)
    else:
        defer_x_bounds = True
    layer_height = 6
    row_width = 80
    y_step = -8
    cur_y = -5
    cur_position = 0
    next_row = cur_position + row_width

    # Track the most extreme sequence for guessing the plot bounds when the
    # parent sequence is not available
    max_position = 0
    for layer in layers:
        for gpm in layer:
            max_position = max(max_position, gpm.end_position)
            rect = mpatches.Rectangle(
                (gpm.start_position, cur_y), width=gpm.sequence_length, height=layer_height,
                facecolor='lightseagreen', edgecolor='darkseagreen',
                alpha=min(max(gpm.ms2_score * 2, 0.2), 0.8))
            label = id_mapper.add("glycopeptide-%d", rect, {
                "sequence": gpm.glycopeptide_sequence,
                "start-position": gpm.start_position,
                "end-position": gpm.end_position,
                "ms2-score": gpm.ms2_score,
                "q-value": gpm.q_value
            })
            ax.add_patch(rect)
            seq = sequence.Sequence(gpm.glycopeptide_sequence)
            for i, pos in enumerate(seq):
                if len(pos[1]) > 0:
                    color = get_color(pos[1][0].name)
                    facecolor, edgecolor = lighten(color), darken(color)
                    mod_patch = mpatches.Rectangle(
                        (gpm.start_position + i, cur_y), width=1, height=layer_height,
                        facecolor=facecolor, edgecolor=edgecolor,
                        alpha=min(max(gpm.ms2_score * 2, 0.4), 0.8))
                    id_mapper.add("modification-%d", mod_patch, {"modification-type": pos[1][0].name, "parent": label})
                    ax.add_patch(mod_patch)
                    ax.text(gpm.start_position + i + .13, cur_y, str(pos[1][0])[0], fontsize=5., family='monospace')
        cur_y += y_step
    ax.set_ylim(cur_y - 100, 50)
    if defer_x_bounds:
        ax.set_xlim(-10, max_position + 10)
    ax.axis('off')

    return ax, id_mapper


def plot_glycoforms(protein, filterfunc=lambda x: x):
    gpms = filterfunc(protein.glycopeptide_matches).all()
    if len(gpms) == 0:
        return None, None
    layers = layout_layers(gpms)
    ax, id_mapper = draw_layers(layers, protein)
    return ax, id_mapper


def plot_glycoforms_svg(protein, filterfunc=lambda x: x):
    ax, id_mapper = plot_glycoforms(protein, filterfunc)

    fig = ax.get_figure()
    fig.tight_layout(pad=0.2)
    buff = StringIO()
    fig.savefig(buff, format='svg')
    root, ids = ET.XMLID(buff.getvalue())
    root.attrib['class'] = 'plot-glycoforms-svg'
    for id, attributes in id_mapper.items():
        element = ids[id]
        element.attrib.update({("data-" + k): str(v) for k, v in attributes.items()})
        element.attrib['class'] = id.rsplit('-')[0]
    svg = ET.tostring(root)
    plt.close()
    return svg


def main(database_path, protein_id, filterfunc=lambda x: x, save_path=None, **kwargs):
    dbm = DatabaseManager(database_path)
    s = dbm.session()
    protein = s.query(Protein).get(protein_id)

    if save_path is None:
        path_template = "{base}_{hypothesis}_{protein}_glycoforms.{format}"
        save_path = path_template.format(
            base=os.path.splitext(database_path)[0], hypothesis=protein.hypothesis.name,
            protein=clean_file_name(protein.name), format=kwargs.get("format", "png"))
        kwargs.setdefault("format", "png")
    # Special handling for SVG to assign element ids
    if kwargs.get('format') == 'svg' and kwargs.get("include_data", True):
        svg = plot_glycoforms_svg(protein, filterfunc)
        with open(save_path, 'wb') as fh:
            fh.write(svg)
        return
    else:
        ax, id_mapper = plot_glycoforms(protein, filterfunc=filterfunc)
        if ax is None:
            return
        fig = ax.get_figure()
        fig.set_size_inches(16, 12)
        fig.savefig(save_path, bbox_inches='tight', **kwargs)
        plt.close()

app = argparse.ArgumentParser("draw_glycoforms")
app.add_argument("database_path", help="path to the database file to analyze")
app.add_argument("-p", "--protein-id", default=None, help="The protein to analyze.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")
app.add_argument("-f", "--format", default='png', choices=["png", "jpeg", "svg", 'pdf'], help="Image file format")


def taskmain():
    args = app.parse_args()
    session = DatabaseManager(args.database_path).session()
    if args.protein_id is None:
        proteins = [j for i in session.query(Protein.id) for j in i]
    else:
        proteins = [args.protein_id]

    for protein_id in proteins:
        if args.out is not None and len(proteins) > 1:
            parts = os.path.splitext(args.out)
            name = session.query(Protein.name).get(protein_id)[0]
            out = parts[0] + "_" + clean_file_name(name) + parts[1]
        else:
            out = (args.out)
        main(args.database_path, protein_id, save_path=out, format=args.format)

if __name__ == '__main__':
    taskmain()
