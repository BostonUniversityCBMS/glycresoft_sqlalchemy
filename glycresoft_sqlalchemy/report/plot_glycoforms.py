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

import matplotlib
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib.textpath import TextPath

from ..data_model import DatabaseManager, GlycopeptideMatch, Protein
from ..structure import sequence
from .colors import lighten, darken, get_color


font_options = font_manager.FontProperties(family='monospace')


def clean_file_name(file_name):
    name = re.sub(r"\|", "_", file_name)
    return name


def span_overlap(a, b):
    return a.spans(b.start_position) or a.spans(b.end_position - 1) or\
           b.spans(a.start_position) or b.spans(a.end_position - 1)


def layout_layers(gpms):
    '''
    Produce a non-overlapping stacked layout of individual peptide-like
    identifications across a protein sequence.
    '''
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
    '''
    A dictionary-like container which uses a format-string
    key pattern to generate unique identifiers for each entry

    Key Pattern: '<type-name>-%d'

    Associates each generated id with a dictionary of metadata and
    sets the `gid` of the passed `matplotlib.Artist` to the generated
    id. Only the id and metadata are stored.

    Used to preserve a mapping of metadata to artists for later SVG
    serialization.
    '''
    def __init__(self):
        dict.__init__(self)
        self.counter = 0

    def add(self, key, value, meta):
        label = key % self.counter
        value.set_gid(label)
        self[label] = meta
        self.counter += 1
        return label


def draw_layers(layers, protein, scale_factor=1.0, **kwargs):
    '''
    Render fixed-width stacked peptide identifications across
    a protein. Each shape is rendered with a unique identifier.
    '''
    figure, ax = plt.subplots(1, 1)
    id_mapper = IDMapper()
    i = 0
    row_width = 70
    layer_height = 0.56 * scale_factor
    y_step = (layer_height + 0.15) * -scale_factor
    cur_y = -3

    cur_position = 0

    mod_text_x_offset = 0.15 * scale_factor
    sequence_font_size = 6. * scale_factor
    mod_font_size = 2.08 * scale_factor
    mod_text_y_offset = 0.1 * scale_factor
    mod_width = 0.5 * scale_factor
    mod_x_offset = 0.25 * scale_factor
    total_length = len(protein.protein_sequence or '')
    protein_pad = -0.365 * scale_factor
    peptide_pad = protein_pad * (2./3.) * scale_factor

    glycosites = set(protein.glycosylation_sites)
    for layer in layers:
        layer.sort(key=lambda x: x.start_position)

    while cur_position < total_length:
        next_row = cur_position + row_width
        i = -2
        text_path = TextPath(
            (protein_pad + i, layer_height + .2 + cur_y),
            str(cur_position + 1), size=sequence_font_size/7.5, prop=font_options, stretch=1000)
        patch = mpatches.PathPatch(text_path, facecolor='grey', lw=0.04)
        ax.add_patch(patch)

        i = row_width + 2
        text_path = TextPath(
            (protein_pad + i, layer_height + .2 + cur_y),
            str(next_row), size=sequence_font_size/7.5, prop=font_options, stretch=1000)
        patch = mpatches.PathPatch(text_path, facecolor='grey', lw=0.04)
        ax.add_patch(patch)

        for i, aa in enumerate(protein.protein_sequence[cur_position:next_row]):
            text_path = TextPath(
                (protein_pad + i, layer_height + .2 + cur_y),
                aa, size=sequence_font_size/7.5, prop=font_options, stretch=1000)
            color = 'red' if any(
                (((i + cur_position) in glycosites),
                 ((i + cur_position - 1) in glycosites),
                 ((i + cur_position - 2) in glycosites))
                ) else 'black'
            patch = mpatches.PathPatch(text_path, facecolor=color, lw=0.04)
            ax.add_patch(patch)

        for layer in layers:
            c = 0
            for gpm in layer:
                if gpm.start_position < cur_position:
                    continue
                elif gpm.start_position > next_row:
                    break
                c += 1
                rect = mpatches.Rectangle(
                    (peptide_pad + gpm.start_position - cur_position, cur_y),
                    width=gpm.sequence_length, height=layer_height,
                    facecolor='lightblue', edgecolor='black', linewidth=0.15,
                    alpha=min(max(gpm.ms2_score * 2, 0.2), 0.8))
                label = id_mapper.add("glycopeptide-%d", rect, {
                    "sequence": gpm.glycopeptide_sequence,
                    "start-position": gpm.start_position,
                    "end-position": gpm.end_position,
                    "ms2-score": gpm.ms2_score,
                    "q-value": gpm.q_value,
                    "record-id": gpm.id,
                    "calculated-mass": gpm.calculated_mass,
                    "spectra-count": gpm.spectrum_matches.filter_by(best_match=True).count()
                })
                ax.add_patch(rect)
                seq = sequence.Sequence(gpm.glycopeptide_sequence)
                for i, pos in enumerate(seq):
                    if len(pos[1]) > 0:
                        color = get_color(pos[1][0].name)
                        facecolor, edgecolor = lighten(color), darken(color, 0.6)

                        mod_patch = mpatches.Rectangle(
                            (gpm.start_position - cur_position + i - mod_x_offset, cur_y),
                            width=mod_width, height=layer_height, alpha=0.4,
                            facecolor=facecolor, edgecolor=edgecolor, linewidth=0.5,
                            #alpha=min(max(gpm.ms2_score * 2, 0.3), 0.6)
                            )
                        id_mapper.add("modification-%d", mod_patch,
                                      {"modification-type": pos[1][0].name, "parent": label})
                        ax.add_patch(mod_patch)
                        text_path = TextPath(
                            (gpm.start_position - cur_position + i - mod_text_x_offset, cur_y + mod_text_y_offset),
                            str(pos[1][0])[0], size=mod_font_size/4.5, prop=font_options)
                        patch = mpatches.PathPatch(text_path, facecolor='black', lw=0.04)
                        ax.add_patch(patch)
                        # ax.text(gpm.start_position - cur_position + i - text_shift, cur_y + 0.1,
                        #         str(pos[1][0])[0], fontsize=mod_font_size, family='monospace')
            if c > 0:
                cur_y += y_step
        cur_y += y_step * 6
        cur_position = next_row

    ax.set_ylim(cur_y - 50, 50)
    ax.set_xlim(-30, row_width + 30)
    ax.axis('off')
    return ax, id_mapper


def plot_glycoforms(protein, filterfunc=lambda x: x.filter(GlycopeptideMatch.ms2_score > 0.2), **kwargs):
    gpms = filterfunc(protein.glycopeptide_matches).all()
    layers = layout_layers(gpms)
    ax, id_mapper = draw_layers(layers, protein, **kwargs)
    return ax, id_mapper


def plot_glycoforms_svg(protein, filterfunc=lambda x: x.filter(GlycopeptideMatch.ms2_score > 0.2), scale=1.5, **kwargs):
    '''
    A specialization of :func:`plot_glycoforms` which adds additional features to SVG images, such
    adding shape metadata to XML tags and properly configuring the viewport and canvas for the figure's
    dimensions.
    '''
    ax, id_mapper = plot_glycoforms(protein, filterfunc)
    old_size = matplotlib.rcParams["figure.figsize"]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.autoscale()

    x_size = sum(map(abs, xlim))
    y_size = sum(map(abs, ylim))

    aspect_ratio = x_size / y_size
    canvas_x = 8.
    canvas_y = canvas_x / aspect_ratio

    matplotlib.rcParams["figure.figsize"] = canvas_x, canvas_y

    fig = ax.get_figure()
    fig.tight_layout(pad=0.2)
    fig.patch.set_visible(False)
    fig.set_figwidth(canvas_x)
    fig.set_figheight(canvas_y)

    ax.patch.set_visible(False)
    buff = StringIO()
    fig.savefig(buff, format='svg')
    root, ids = ET.XMLID(buff.getvalue())
    root.attrib['class'] = 'plot-glycoforms-svg'
    for id, attributes in id_mapper.items():
        element = ids[id]
        element.attrib.update({("data-" + k): str(v) for k, v in attributes.items()})
        element.attrib['class'] = id.rsplit('-')[0]
    min_x, min_y, max_x, max_y = map(int, root.attrib["viewBox"].split(" "))
    min_x += 100
    max_x += 200
    view_box = ' '.join(map(str, (min_x, min_y, max_x, max_y)))
    root.attrib["viewBox"] = view_box
    width = float(root.attrib["width"][:-2]) * 2
    root.attrib["width"] = "%dpt" % width

    height = width / (aspect_ratio)

    root.attrib["height"] = "%dpt" % (height * 1.2)
    root.attrib["preserveAspectRatio"] = "xMinYMin meet"
    root[1].attrib["transform"] = "scale(%f)" % scale
    svg = ET.tostring(root)
    plt.close()

    matplotlib.rcParams["figure.figsize"] = old_size
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
