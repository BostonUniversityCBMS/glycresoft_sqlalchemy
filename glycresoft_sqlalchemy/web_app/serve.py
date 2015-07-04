import logging
try:
    logging.basicConfig(level='DEBUG')
except:
    pass
import json
import functools
from glycresoft_sqlalchemy.web_app import report
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, \
     Response

from werkzeug import secure_filename
import argparse

from glycresoft_sqlalchemy.data_model import Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch
from glycresoft_sqlalchemy.web_app.project_manager import ProjectManager

from glycresoft_sqlalchemy.web_app.task.do_bupid_yaml_parse import BUPIDYamlParseTask
from glycresoft_sqlalchemy.web_app.task.task_process import QueueEmptyException

app = Flask(__name__)
report.prepare_environment(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'

manager = None


def connect_db():
    g.db = manager.session()


def message_queue_stream():
    """Implement a simple Server Side Event (SSE) stream based on the
    stream of events emit from the :attr:`TaskManager.messages` queue of `manager`.

    These messages are handled on the client side.

    Yields
    ------
    str: Formatted Server Side Event Message

    References
    ----------
    [1] - http://stackoverflow.com/questions/12232304/how-to-implement-server-push-in-flask-framework
    """
    payload = 'id: {id}\nevent: {event_name}\ndata: {data}\n\n'
    i = 0
    yield payload.format(id=i, event_name='begin-stream', data=json.dumps('Starting Stream'))
    i += 1
    while True:
        try:
            message = manager.messages.get(True, 3)
            event = payload.format(
                id=i, event_name=message.type,
                data=json.dumps(message.message))
            i += 1
            print message, event
            yield event
        except KeyboardInterrupt:
            break
        except QueueEmptyException, e:
            # Send a comment to keep the connection alive
            yield ":no events\n\n"
        except Exception, e:
            logging.exception("An error occurred in message_queue_stream", exc_info=e)


@app.route("/")
def index():
    return render_template("index.templ")


@app.route("/tasks")
def get_tasks(self):
    return jsonify(**{t.id: t.to_json() for t in manager.tasks.values()})


@app.route('/stream')
def message_stream():

    return Response(message_queue_stream(),
                    mimetype="text/event-stream")


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
    """Handle an uploaded sample file

    Returns
    -------
    TYPE : Description
    """
    run_name = request.values['sample_name']
    secure_name = secure_filename(run_name)
    path = manager.get_temp_path(secure_name)
    request.files['observed-ions-file'].save(path)
    dest = manager.get_sample_path(run_name)
    # Construct the task with a callback to add the processed sample
    # to the set of project samples
    task = BUPIDYamlParseTask(
        manager.path,
        path,
        dest,
        functools.partial(manager.add_sample, path=dest))
    manager.add_task(task)
    return redirect("/")


@app.route("/add_sample")
def add_sample():
    return render_template("add_sample_form.templ")


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
    app.run(use_reloader=False, threaded=True)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.results_database)
