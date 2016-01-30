import logging
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, \
     Response

import argparse

from werkzeug.contrib.profiler import ProfilerMiddleware
from werkzeug.wsgi import LimitedStream

from glycresoft_sqlalchemy.web_app.project_manager import ProjectManager
from glycresoft_sqlalchemy.web_app.utils.cache import ApplicationDataCache
from glycresoft_sqlalchemy.web_app import report

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, Protein, TheoreticalGlycopeptide, PeakGroupMatchType, TheoreticalGlycopeptideComposition,
    GlycopeptideMatch)


from glycresoft_sqlalchemy.web_app.services.server_sent_events import server_sent_events
from glycresoft_sqlalchemy.web_app.services.json_api import api as json_api
from glycresoft_sqlalchemy.web_app.services.preferences import app_config as preferences
from glycresoft_sqlalchemy.web_app.services.task_management import task_actions
from glycresoft_sqlalchemy.web_app.services.file_exports import file_exports
from glycresoft_sqlalchemy.web_app.services.sample_management import sample_management
from glycresoft_sqlalchemy.web_app.services.view_database_search_results import view_database_search_results
from glycresoft_sqlalchemy.web_app.services.make_glycopeptide_hypothesis import make_glycopeptide_hypothesis
from glycresoft_sqlalchemy.web_app.services.peak_group_matching import peak_group_matching
from glycresoft_sqlalchemy.web_app.services.tandem_glycopeptide_search import tandem_glycopeptide_search
from glycresoft_sqlalchemy.web_app.services.make_glycan_hypothesis import make_glycan_hypothesis
from glycresoft_sqlalchemy.web_app.services.view_hypothesis import view_hypothesis

from glycresoft_sqlalchemy.web_app.utils.pagination import paginate


class StreamConsumingMiddleware(object):

    def __init__(self, app):
        self.app = app

    def __call__(self, environ, start_response):
        stream = LimitedStream(environ['wsgi.input'],
                               int(environ['CONTENT_LENGTH'] or 0))
        environ['wsgi.input'] = stream
        app_iter = self.app(environ, start_response)
        try:
            stream.exhaust()
            for event in app_iter:
                yield event
        finally:
            if hasattr(app_iter, 'close'):
                app_iter.close()



app = Flask(__name__)
app.wsgi_app = StreamConsumingMiddleware(app.wsgi_app)
# app.wsgi_app = ProfilerMiddleware(app.wsgi_app, profile_dir='profiling')
app.config['PROPAGATE_EXCEPTIONS'] = True
report.prepare_environment(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'
SERVER = None


manager = None


app.register_blueprint(server_sent_events)
app.register_blueprint(json_api)
app.register_blueprint(preferences)
app.register_blueprint(task_actions)
app.register_blueprint(file_exports)
app.register_blueprint(sample_management)
app.register_blueprint(view_database_search_results)
app.register_blueprint(make_glycopeptide_hypothesis)
app.register_blueprint(make_glycan_hypothesis)
app.register_blueprint(peak_group_matching)
app.register_blueprint(tandem_glycopeptide_search)
app.register_blueprint(view_hypothesis)


@app.route("/ms1_or_ms2_choice")
def branch_ms1_ms2():
    # This would be better done completely client-side, but
    # requires some templating engine/cache on the client
    ms1_choice = request.values.get("ms1_choice")
    ms2_choice = request.values.get("ms2_choice")
    return render_template("components/ms1_or_ms2_choice.templ",
                           ms1_choice=ms1_choice, ms2_choice=ms2_choice)


# ----------------------------------------
# Server Shutdown
# ----------------------------------------

def shutdown_server():
    if no_gevent:
        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')

        func()
    else:
        SERVER.stop()


@app.route('/internal/shutdown', methods=['POST'])
def shutdown():
    g.manager.halting = True
    shutdown_server()
    g.manager.stoploop()
    return Response("Should be dead")

# ----------------------------------------
#
# ----------------------------------------


@app.route("/internal/show_cache")
def show_cache():
    print cache.app_data
    return Response("Printed")


def connect_db():
    g.manager = manager
    g.db = manager.session()


@app.route("/")
def index():
    return render_template("index.templ")


@app.before_request
def before_request():
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
        "Manager": manager,
        "PeakGroupMatchType": PeakGroupMatchType,
        "TheoreticalGlycopeptideComposition": TheoreticalGlycopeptideComposition,
    }


@app.context_processor
def inject_functions():
    def query(args):
        return g.db.query(args)
    return {
        "query": query,
        "paginate": paginate
    }

parser = argparse.ArgumentParser('view-results')
parser.add_argument("store_path")
parser.add_argument("-n", "--no-execute-tasks", action="store_true", required=False, default=False)
parser.add_argument("--external", action='store_true', required=False, default=False,
                    help='Let non-host machines connect to the server')
parser.add_argument("-p", "--port", required=False, type=int, default=5000, help="The port on which to run the server")


DEBUG = False


def setup_logging():
    try:
        logging.basicConfig(
            level=logging.INFO, filename='glycresoft-log', filemode='w',
            format="%(asctime)s - %(processName)s:%(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
            datefmt="%H:%M:%S")
        fmt = logging.Formatter(
            "%(asctime)s - %(processName)s:%(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s", "%H:%M:%S")
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logging.getLogger().addHandler(handler)
    except Exception, e:
        logging.exception("Error, %r", e, exc_info=e)
        raise e


def run(store_path, external, no_execute_tasks, port, **kwargs):
    global DATABASE, manager, CAN_EXECUTE, SERVER
    host = None
    if external:
        host = "0.0.0.0"
    DATABASE = store_path
    CAN_EXECUTE = not no_execute_tasks

    manager = ProjectManager(DATABASE)

    app.debug = DEBUG
    app.secret_key = SECRETKEY
    # setup_logging()
    app.run(host=host, use_reloader=False, threaded=True, debug=DEBUG, port=port, passthrough_errors=True)


def main():
    args = parser.parse_args()
    run(**args.__dict__)

if __name__ == "__main__":
    main()
