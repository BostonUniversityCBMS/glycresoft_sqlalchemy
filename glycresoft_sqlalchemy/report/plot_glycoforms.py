import argparse
import os

from itertools import cycle
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib.path import Path
from matplotlib.colors import cnames, hex2color

from glycresoft_sqlalchemy.data_model import DatabaseManager, GlycopeptideMatch, Protein, Hypothesis
from glycresoft_ms2_classification.structure import sequence


def lighten(rgb, factor=0.25):
    factor += 1
    return [min(c * factor, 1) for c in rgb]

def darken(rgb, factor=0.25):
    factor = 1 - factor
    return [(c * factor) for c in rgb]


colors = cycle([hex2color(cnames[name]) for name in ("red", "blue", "yellow", "purple", "navy", "grey")])


color_name_map = {
    "HexNAc": hex2color(cnames["steelblue"])
}


def get_color(name):
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


def draw_layers(layers, protein):
    figure, ax = plt.subplots(1,1)
    for i, aa in enumerate(protein.protein_sequence):
        ax.text(-0.2 + i, 5.5, aa, family="monospace", fontsize=4)
    rect = mpatches.Rectangle((0, 0), i, 5, facecolor='red', alpha=0.5)
    ax.add_patch(rect)
    ax.set_xlim(-10, i + 10)
    layer_height = 6
    row_width = 80
    y_step = -8
    
    cur_y = -5
    cur_position = 0
    next_row = cur_position + row_width
    for layer in layers:
        for gpm in layer:
            rect = mpatches.Rectangle((gpm.start_position, cur_y), width=gpm.sequence_length, height=layer_height, 
                                      facecolor='lightseagreen', edgecolor='darkseagreen', alpha=min(max(gpm.ms2_score * 2, 0.2), 0.8))
            ax.add_patch(rect)
            seq = sequence.Sequence(gpm.glycopeptide_sequence)
            for i, pos in enumerate(seq):
                if len(pos[1]) > 0:
                    color = get_color(pos[1][0].name)
                    facecolor, edgecolor = lighten(color), darken(color)
                    mod_patch = mpatches.Rectangle((gpm.start_position + i, cur_y), width=1, height=layer_height,
                                                   facecolor=facecolor, edgecolor=edgecolor,
                                                   alpha=min(max(gpm.ms2_score * 2, 0.4), 0.8))
                    ax.add_patch(mod_patch)
                    ax.text(gpm.start_position + i + .02, cur_y, str(pos[1][0])[0], fontsize=5., family='monospace')
        cur_y += y_step
    ax.set_ylim(cur_y-100, 50)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return ax


def plot_glycoforms(protein, filterfunc=lambda x: x):
    gpms = filterfunc(protein.glycopeptide_matches).all()
    layers = layout_layers(gpms)
    ax = draw_layers(layers, protein)
    return ax


def main(database_path, protein_id, filterfunc=lambda x: x, save_path=None, **kwargs):
    dbm = DatabaseManager(database_path)
    s = dbm.session()
    protein = s.query(Protein).get(protein_id)

    if save_path is None:
        save_path = protein.hypothesis.name + "_" + protein.name + "_glycoforms.png"
        kwargs['format'] = 'png'

    ax = plot_glycoforms(protein, filterfunc=filterfunc)
    plt.savefig(save_path, bbox_inches='tight', pad_inches=0.2, **kwargs)
    plt.close()

app = argparse.ArgumentParser("draw_glycoforms")
app.add_argument("database_path", help="path to the database file to analyze")
app.add_argument("-p", "--protein-id", default=None, help="The protein to analyze.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")

def taskmain():
    args = app.parse_args()
    session = DatabaseManager(args.database_path).session()
    if args.protein_id is None:
        proteins = [j for i in session.query(Protein.id) for j in i]
    else:
        proteins = [args.protein_id]

    for protein_id in proteins
        if args.out is not None and len(proteins) > 1:
            parts = os.path.splitext(args.out)
            name = session.query(Protein.name).get(protein_id)[0]
            out = parts[0] + "_" + name + parts[1]
        else:
            out = args.out
        main(args.database_path, protein_id, out)

if __name__ == '__main__':
    taskmain()
