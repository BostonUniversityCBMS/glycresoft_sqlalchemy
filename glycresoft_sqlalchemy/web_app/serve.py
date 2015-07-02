import logging
try:
    logging.basicConfig(level='DEBUG')
except:
    pass

from glycresoft_sqlalchemy.web_app import report
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, \
     Response

from werkzeug import secure_filename
import argparse

from glycresoft_sqlalchemy.data_model import Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch
from project_manager import ProjectManager


app = Flask(__name__)
report.prepare_environment(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'

manager = None


def connect_db():
    g.db = manager.session()


@app.route("/")
def index():
    return render_template("index.templ")


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


@app.route("/add_sample", methods=["POST"])
def post_add_sample():
    print dir(request.files['observed-ions-file'])
    run_name = request.values['sample_name']
    secure_name = secure_filename(run_name) + '.unprocessed'
    path = manager.get_sample_path(secure_name)
    print path
    request.files['observed-ions-file'].save(path)
    return redirect("/")


@app.route("/add_sample")
def add_sample():
    return render_template("add_sample_form.templ")


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
        "GlycopeptideMatch": GlycopeptideMatch,
        "Manager": manager
    }


@app.context_processor
def inject_functions():
    def query(args):
        return g.db.query(args)
    return locals()

parser = argparse.ArgumentParser('view-results')
parser.add_argument("results_database")


def main(results_database):
    global DATABASE, manager
    DATABASE = results_database
    manager = ProjectManager(DATABASE)
    app.debug = DEBUG
    app.secret_key = SECRETKEY
    app.run()

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.results_database)
