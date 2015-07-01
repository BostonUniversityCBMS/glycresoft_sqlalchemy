import itertools

from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, GlycopeptideMatch, Protein

import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

matplotlib.rcParams['figure.figsize'] = 10, 8


def fetch(database_path, hypothesis_id):
    session = DatabaseManager(database_path).session()
    qry = session.query(GlycopeptideMatch.ms2_score,
                        GlycopeptideMatch.ms1_score,
                        GlycopeptideMatch.volume,
                        GlycopeptideMatch.q_value,
                        Protein.name).filter(
                        GlycopeptideMatch.protein_id == Protein.id,
                        Protein.hypothesis_id == hypothesis_id)
    return qry


def plot(qry):
    colors = itertools.cycle(["red", "blue", "green", "grey", "gold", "purple"])
    frame = pd.DataFrame(iter(qry), columns=["MS2 Score", "MS1 Score", "Abundance", "q-value", "Protein Name"])
    fig, axes = plt.subplots(2)
    ax1, ax2 = axes
    color_map = {pn: colors.next() for pn in set(frame["Protein Name"])}
    frame.plot(x='MS2 Score', y='MS1 Score', alpha=0.5, kind='scatter',
               s=(frame.Abundance / frame.Abundance.max()) * 200, c=frame["Protein Name"].map(color_map),
               ax=ax1)
    frame.plot(y='q-value', x='MS2 Score', alpha=0.5, kind='scatter',
               s=(frame.Abundance / frame.Abundance.max()) * 200, c=frame["Protein Name"].map(color_map),
               ax=ax2)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-0.5, 1)

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    handles, labels = (ax1.get_legend_handles_labels())
    handles += tuple(
        mpatches.Patch(color=color_map[n], label=n) for n in color_map) + (mpatches.Patch(color="white", label='Size ~ Volume'), ) 
    legend = ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
                        ncol=3, fancybox=True, handles=handles)
    return axes, legend


def plot_scores(database_path, hypothesis_id):
    return plot(fetch(database_path, hypothesis_id))


def main(database_path, hypothesis_id, save_path=None, **kwargs):
    if save_path is None:
        session = DatabaseManager(database_path).session()
        exp = session.query(Hypothesis).get(hypothesis_id)
        save_path = exp.name + "_score_scatterplot.png"
        kwargs["format"] = 'png'
    ax, legend = plot_scores(database_path, hypothesis_id)
    plt.savefig(save_path, bbox_extra_artists=(legend,), bbox_inches='tight', pad_inches=0.2, **kwargs)
    plt.close()

if __name__ == '__main__':
    import argparse
    import os
    app = argparse.ArgumentParser("score_scatterplot")
    app.add_argument("database_path", help="path to the database file to analyze")
    app.add_argument("-e", "--hypothesis-id", default=None, help="The hypothesis to analyze.")
    app.add_argument("-o", "--out", default=None, help="Where to save the result")

    args = app.parse_args()
    session = DatabaseManager(args.database_path).session()
    if args.hypothesis_id is None:
        hypothesiss = [j for i in session.query(Hypothesis.id) for j in i]
    else:
        hypothesiss = [args.hypothesis_id]
    for hypothesis_id in hypothesiss:
        if args.out is not None and len(hypothesiss) > 1:
            parts = os.path.splitext(args.out)
            name = session.query(Hypothesis.name).get(hypothesis_id)[0]
            out = parts[0] + "_" + name + parts[1]
        else:
            out = args.out
        main(args.database_path, hypothesis_id, out)

