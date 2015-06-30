from glycresoft_sqlalchemy.web_app import report
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, Response
import argparse

from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch
from glycresoft_sqlalchemy.data_model import informed_proteomics


app = Flask(__name__)
report.prepare_environment(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'


def connect_db():
    g.db = DatabaseManager(DATABASE).session()


@app.route("/")
def index():
    hypotheses = g.db.query(Hypothesis).all()
    return render_template("index.templ", hypotheses=hypotheses)


@app.route("/hypothesis")
def show_hypotheses():
    return render_template("show_hypotheses.templ", hypotheses=g.db.query(Hypothesis).all())


@app.route("/glycan_search_space")
def build_naive_glycan_search():
    return render_template("glycan_search_space.templ")


@app.route("/glycan_search_space", methods=["POST"])
def build_naive_glycan_search_process():
    print request.values
    return jsonify(**dict(request.values))


@app.route("/glycopeptide_search_space")
def build_naive_glycopeptide_search_space():
    return render_template("glycopeptide_search_space.templ")


@app.route("/match_samples")
def match_samples():
    return render_template("match_samples.templ")


@app.before_request
def before_request():
    print session
    connect_db()


@app.teardown_request
def teardown_request(exception):
    db = getattr(g, 'db', None)
    if db is not None:
        db.close()


@app.context_processor
def inject_model():
    return {
        "Hypothesis": Hypothesis,
        "Protein": Protein,
        "TheoreticalGlycopeptide": TheoreticalGlycopeptide,
        "GlycopeptideMatch": GlycopeptideMatch
    }


@app.context_processor
def inject_functions():
    def query(args):
        return g.db.query(args)
    return locals()

parser = argparse.ArgumentParser('view-results')
parser.add_argument("results_database")


def main(results_database):
    global DATABASE
    DATABASE = results_database
    app.debug = DEBUG
    app.secret_key = SECRETKEY
    app.run()

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.results_database)
