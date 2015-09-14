import argparse
import os
from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, GlycopeptideMatch, Protein

import pandas as pd
import matplotlib
from matplotlib import pyplot as plt

# matplotlib.rcParams['figure.figsize'] = 10, 8


def fetch(database_path, hypothesis_id):
    session = DatabaseManager(database_path).session()
    qry = session.query(GlycopeptideMatch.q_value,
                        GlycopeptideMatch.protein_id,
                        Protein.name).filter(
                        GlycopeptideMatch.protein_id == Protein.id,
                        Protein.hypothesis_id == hypothesis_id)
    return qry


def plot(qry):
    hypothesis = pd.DataFrame(iter(qry), columns=["mass", "protein_id", "protein_name"])
    groups = hypothesis.groupby("protein_name")
    for protein_id, group in groups:
        ax = group['mass'].hist(label=protein_id, binsize=0.01, alpha=0.5)

    ax.set_title("q-value frequency")
    ax.set_xlabel("q value")
    ax.set_ylabel("Counts")
    handles, labels = (ax.get_legend_handles_labels())
    legend = ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
                       ncol=3, fancybox=True, handles=handles)

    return ax, legend


def histogram_complexity(database_path, hypothesis_id):
    query = fetch(database_path, hypothesis_id)
    return plot(query)


def main(database_path, hypothesis_id, save_path=None, **kwargs):
    if save_path is None:
        session = DatabaseManager(database_path).session()
        exp = session.query(Hypothesis).get(hypothesis_id)
        save_path = exp.name + "_q_values.png"
        kwargs["format"] = 'png'
    ax, legend = histogram_complexity(database_path, hypothesis_id)
    plt.savefig(save_path, bbox_extra_artists=(legend,), bbox_inches='tight', pad_inches=0.2, **kwargs)
    plt.close()

app = argparse.ArgumentParser("q_value_vis")
app.add_argument("database_path", help="path to the database file to analyze")
app.add_argument("-e", "--hypothesis-id", default=None, help="The hypothesis to analyze.")
app.add_argument("-o", "--out", default=None, help="Where to save the result")


def taskmain():
    args = app.parse_args()
    session = DatabaseManager(args.database_path).session()
    if args.hypothesis_id is None:
        hypotheses = [j for i in session.query(Hypothesis.id) for j in i]
    else:
        hypotheses = [args.hypothesis_id]
    for hypothesis_id in hypotheses:
        if args.out is not None and len(hypotheses) > 1:
            parts = os.path.splitext(args.out)
            name = session.query(Hypothesis.name).get(hypothesis_id)[0]
            out = parts[0] + "_" + name + parts[1]
        else:
            out = args.out
        main(args.database_path, hypothesis_id, out)

if __name__ == '__main__':
    taskmain()
