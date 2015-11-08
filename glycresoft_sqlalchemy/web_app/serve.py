import logging
import os
import base64
import json
import functools
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, \
     Response

from werkzeug import secure_filename
import argparse
import random

from glycresoft_sqlalchemy.web_app.project_manager import ProjectManager
from glycresoft_sqlalchemy.web_app import report

from glycresoft_sqlalchemy.data_model import (
    Hypothesis, Protein, TheoreticalGlycopeptide, PeakGroupMatch, TheoreticalGlycopeptideComposition,
    GlycopeptideMatch, HypothesisSampleMatch, MassShift, json_type, make_transient, MS1GlycanHypothesisSampleMatch,
    MS2GlycanHypothesisSampleMatch, MS1GlycopeptideHypothesisSampleMatch, MS2GlycopeptideHypothesisSampleMatch)

from glycresoft_sqlalchemy.report import microheterogeneity
from glycresoft_sqlalchemy.utils.database_utils import get_or_create
from glycresoft_sqlalchemy.web_app.task.do_bupid_yaml_parse import BUPIDYamlParseTask
from glycresoft_sqlalchemy.web_app.task.do_decon2ls_parse import Decon2LSIsosParseTask
from glycresoft_sqlalchemy.web_app.task.do_ms2_search import TandemMSGlycoproteomicsSearchTask
from glycresoft_sqlalchemy.web_app.task.do_ms1_peak_group_matching import LCMSSearchTask
from glycresoft_sqlalchemy.web_app.task.do_naive_glycopeptide_hypothesis import NaiveGlycopeptideHypothesisBuilderTask
from glycresoft_sqlalchemy.web_app.task.do_export_csv import ExportCSVTask
from glycresoft_sqlalchemy.web_app.task.task_process import QueueEmptyException
from glycresoft_sqlalchemy.web_app.task.dummy import DummyTask

from glycresoft_sqlalchemy.web_app.utils.pagination import paginate


app = Flask(__name__)
report.prepare_environment(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'

manager = None


JSONEncoderType = json_type.new_alchemy_encoder()


def connect_db():
    g.db = manager.session()


def message_queue_stream():
    """Implement a simple Server Side Event (SSE) stream based on the
    stream of events emit from the :attr:`TaskManager.messages` queue of `manager`.

    These messages are handled on the client side.

    At the moment, messages are not "addressed" to a particular recipient. If multiple users
    are connected at once, who receives which message is undefined. A solution to this would
    be to create labeled queues, but this requires a user identification system.

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
    yield payload.format(id=i - 1, event_name='update', data=json.dumps('Initialized'))
    i += 1
    while True:
        try:
            message = manager.messages.get(True, 1)
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
            if random.random() > 0.8:
                yield payload.format(id=i, event_name='tick', data=json.dumps('Tick'))
        except Exception, e:
            logging.exception("An error occurred in message_queue_stream", exc_info=e)


@app.route("/")
def index():
    return render_template("index.templ")


@app.route('/stream')
def message_stream():
    return Response(message_queue_stream(),
                    mimetype="text/event-stream")


# ----------------------------------------
#           Settings and Preferences
# ----------------------------------------


@app.route("/preferences")
def show_preferences():
    preferences = request.values
    return render_template("components/preferences.templ", **preferences)


@app.route("/preferences", methods=["POST"])
def update_preferences():
    preferences = request.values
    print "Minimum Score:", preferences["minimum-score"]
    return jsonify(**dict(preferences.items()))


@app.route("/internal/update_settings", methods=["POST"])
def update_settings():
    '''
    TODO
    ----
    Diff incoming settings with server-side settings and
    send back the union of the settings to the client.
    '''
    settings = request.values
    print settings
    return jsonify(**settings)


@app.route("/test/task-test")
def test_page():
    return render_template("test.templ")


@app.route("/internal/log/<task_id>")
def send_log(task_id):
    return Response("<pre>%s</pre>" % open(manager.get_task_path(task_id + '.log'), 'r').read(), mimetype='application/text')


# ---------------------------------------
#        File Download Operations
# ---------------------------------------


@app.route("/internal/file_download/<b64path>")
def download_file(b64path):
    path = base64.b64decode(b64path)
    name = os.path.basename(path)
    if manager.temp_dir in path:
        def yielder():
            for line in open(path):
                yield line
        return Response(yielder(), mimetype="application/octet-stream",
                        headers={"Content-Disposition": "attachment; filename=%s;" % name})


@app.route("/view_database_search_results/<int:id>", methods=["POST"])
def view_database_search_results(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    results_type = hsm.results_type
    if PeakGroupMatch == results_type:
        return view_composition_database_search_results(id)
    elif GlycopeptideMatch == results_type:
        return view_tandem_glycopeptide_database_search_results(id)
    return Response("No Match, " + results_type.__name__)

# -------------------------------------------------------------------------
#           View Tandem Glycopeptide Database Search Results
# -------------------------------------------------------------------------


# Dispatch through view_database_search_results
def view_tandem_glycopeptide_database_search_results(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    hypothesis_sample_match_id = id

    state = request.get_json()
    print state
    settings = state['settings']
    context = state['context']

    minimum_score = settings.get('minimum-score', 0.2)

    def filter_context(q):
        return q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            GlycopeptideMatch.ms2_score > minimum_score)

    return render_template(
        "tandem_glycopeptide_search/view_database_search_results.templ",
        hsm=hsm,
        filter_context=filter_context)


@app.route("/view_database_search_results/protein_view/<int:id>", methods=["POST"])
def view_tandem_glycopeptide_protein_results(id):
    parameters = request.get_json()
    print parameters
    print id
    hypothesis_sample_match_id = parameters['context']['hypothesis_sample_match_id']
    protein = g.db.query(Protein).get(id)

    minimum_score = float(parameters['settings'].get("minimum-score", 0.2))

    def filter_context(q):
        return q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            GlycopeptideMatch.ms2_score > minimum_score)

    site_summary = microheterogeneity.GlycoproteinMicroheterogeneitySummary(
        protein, filter_context)

    return render_template(
        "tandem_glycopeptide_search/components/protein_view.templ",
        protein=protein,
        site_summary=site_summary,
        filter_context=filter_context)


@app.route("/view_database_search_results/view_glycopeptide_details/<int:id>")
def view_tandem_glycopeptide_glycopeptide_details(id):
    gpm = g.db.query(GlycopeptideMatch).get(id)
    return render_template(
        "tandem_glycopeptide_search/components/glycopeptide_details.templ", glycopeptide=gpm)


# ----------------------------------------------------------------------
#           View Peak Grouping Database Search Results
# ----------------------------------------------------------------------


# Dispatch through view_database_search_results
def view_composition_database_search_results(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    hypothesis_sample_match_id = id

    def filter_context(q):
        return q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            PeakGroupMatch.ms1_score > 0.2)

    return render_template(
        "peak_group_search/view_database_search_results.templ",
        hsm=hsm,
        filter_context=filter_context)


@app.route("/view_database_search_results/protein_composition_view/<int:id>", methods=["POST"])
def view_composition_glycopeptide_protein_results(id):
    print request.values
    hypothesis_sample_match_id = request.values["hypothesis_sample_match_id"]
    protein = g.db.query(Protein).get(id)

    def filter_context(q):
        return q.filter(
            PeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
            PeakGroupMatch.ms1_score > 0.2)

    return render_template(
        "peak_group_search/components/protein_view.templ",
        protein=protein,
        filter_context=filter_context)


@app.route("/view_database_search_results/glycopeptide_matches_composition_table"
           "/<int:protein_id>/<int:page>", methods=["POST"])
def view_composition_glycopeptide_table_partial(protein_id, page):
    print request.values
    hypothesis_sample_match_id = request.get_json()["context"]["hypothesis_sample_match_id"]
    protein = g.db.query(Protein).get(protein_id)

    def filter_context(q):
        return q.filter(
            PeakGroupMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
            PeakGroupMatch.ms1_score > 0.2)

    paginator = paginate(filter_context(protein.peak_group_matches).order_by(PeakGroupMatch.ms1_score.desc()), page, 50)

    return render_template(
        "peak_group_search/components/glycopeptide_match_table.templ",
        paginator=paginator)


@app.route("/view_database_search_results/view_glycopeptide_composition_details/<int:id>")
def view_peak_grouping_glycopeptide_composition_details(id):
    pgm = g.db.query(PeakGroupMatch).get(id)
    ambiguous_with = g.db.query(PeakGroupMatch).filter(
        PeakGroupMatch.peak_group_id == pgm.peak_group_id,
        PeakGroupMatch.id != pgm.id).all()
    return render_template(
        "peak_group_search/components/glycopeptide_details.templ", pgm=pgm,
        ambiguous_with=ambiguous_with)

# ----------------------------------------
#           CSV Export
# ----------------------------------------


def filterfunc_template(q, model, attr, value, code=1):
    return q.filter(getattr(model, attr) >= value)


@app.route("/view_database_search_results/export_csv/<int:id>", methods=["POST"])
def export_csv_task(id):
    hypothesis_sample_match = g.db.query(HypothesisSampleMatch).get(id)
    tempdir = manager.get_temp_path("glycresoft_export")

    state = request.get_json()
    settings = state["settings"]
    context = state['context']
    minimum_score = settings.get("minimum-score", 0.2)

    if isinstance(hypothesis_sample_match, (MS1GlycopeptideHypothesisSampleMatch,
                                            MS1GlycanHypothesisSampleMatch)):
        filterfunc = functools.partial(filterfunc_template, model=PeakGroupMatch, attr="ms1_score", value=minimum_score)
    elif isinstance(hypothesis_sample_match, MS2GlycopeptideHypothesisSampleMatch):
        filterfunc = functools.partial(filterfunc_template, model=GlycopeptideMatch, attr="ms2_score", value=minimum_score)

    task = ExportCSVTask(manager.path, id, filterfunc, tempdir)

    manager.add_task(task)
    return jsonify(target=hypothesis_sample_match.to_json())


# ----------------------------------------
#           JSON Data API Calls
# ----------------------------------------


@app.route("/api/glycopeptide_matches/<int:id>")
def get_glycopeptide_match_api(id):
    gpm = g.db.query(GlycopeptideMatch).get(id)
    return Response(JSONEncoderType().encode(gpm), mimetype="text/json")


@app.route("/api/tasks")
def api_tasks():
    return jsonify(**{t.id: t.to_json() for t in manager.tasks.values()})


@app.route("/api/hypothesis_sample_matches")
def api_hypothesis_sample_matches():
    hsms = g.db.query(HypothesisSampleMatch).all()
    d = {str(h.id): h.to_json() for h in hsms}
    return jsonify(**d)


@app.route("/api/hypotheses")
def api_hypothesis():
    hypotheses = g.db.query(Hypothesis).all()
    d = {str(h.id): h.to_json() for h in hypotheses}
    return jsonify(**d)


@app.route("/api/samples")
def api_samples():
    samples = manager.samples()
    d = {str(h.name): h.to_json() for h in samples}
    return jsonify(**d)


# ----------------------------------------
#
# ----------------------------------------


@app.route("/ms1_or_ms2_choice")
def branch_ms1_ms2():
    # This would be better done completely client-side, but
    # requires some templating engine/cache on the client
    ms1_choice = request.values.get("ms1_choice")
    ms2_choice = request.values.get("ms2_choice")
    return render_template("components/ms1_or_ms2_choice.templ",
                           ms1_choice=ms1_choice, ms2_choice=ms2_choice)


# ----------------------------------------
#
# ----------------------------------------


@app.route("/hypothesis")
def show_hypotheses():
    return render_template("show_hypotheses.templ", hypotheses=g.db.query(Hypothesis).all())


@app.route("/view_hypothesis/<int:id>")
def view_hypothesis(id):
    return render_template("show_hypotheses.templ", hypotheses=[g.db.query(Hypothesis).get(id)])


# ----------------------------------------
#
# ----------------------------------------


@app.route("/glycan_search_space")
def build_naive_glycan_search():
    return render_template("glycan_search_space.templ")


@app.route("/glycan_search_space", methods=["POST"])
def build_naive_glycan_search_process():
    print request.values
    return jsonify(**dict(request.values))


# ----------------------------------------
#
# ----------------------------------------


@app.route("/glycopeptide_search_space")
def build_naive_glycopeptide_search_space():
    return render_template("glycopeptide_search_space.templ")


@app.route("/glycopeptide_search_space", methods=["POST"])
def build_naive_glycopeptide_search_space_post():
    values = request.values
    constant_modifications = values.getlist("constant_modifications")
    variable_modifications = values.getlist("variable_modifications")
    enzyme = values.get("enzyme")
    hypothesis_name = values.get("hypothesis_name")
    protein_fasta = request.files["protein-fasta-file"]
    site_list = request.files["glycosylation-site-list-file"]
    glycan_file = request.files["glycan-definition-file"]
    glycan_file_type = values.get("glycans-file-format")
    max_missed_cleavages = int(values.get("missed_cleavages"))

    secure_protein_fasta = manager.get_temp_path(secure_filename(protein_fasta.filename))
    secure_site_list_file = manager.get_temp_path(secure_filename(site_list.filename))
    secure_glycan_file = manager.get_temp_path(secure_filename(glycan_file.filename))

    protein_fasta.save(secure_protein_fasta)
    glycan_file.save(secure_glycan_file)
    if site_list.filename != "":
        site_list.save(secure_site_list_file)
    else:
        secure_site_list_file = None

    task = NaiveGlycopeptideHypothesisBuilderTask(
        manager.path, hypothesis_name,
        secure_protein_fasta, secure_site_list_file,
        secure_glycan_file, glycan_file_type, constant_modifications,
        variable_modifications, enzyme, max_missed_cleavages, callback=lambda: 0)
    manager.add_task(task)
    return jsonify(**dict(request.values))


# ----------------------------------------
#
# ----------------------------------------


@app.route("/peak_grouping_match_samples")
def peak_grouping_match_samples():
    return render_template("peak_group_search/peak_grouping_match_samples.templ")


@app.route("/peak_grouping_match_samples", methods=["POST"])
def peak_grouping_match_samples_post():
    '''
    Handle received parameters for /peak_grouping_match_samples and schedule
    all resulting tasks

    Parameters
    ----------
    match_tolerance: float
    grouping_error_tolerance: float
    minimum_scan_count: int
    database_path: str

    hypothesis_choice: str
        Encodes the input type and record id as a comma separated string

    mass_shift_name: list of str
    mass_shift_mass_delta: list of float
    mass_shift_max_count: list of int

    '''
    user_parameters = request.values
    job_parameters = {
        "match_tolerance": float(user_parameters["mass-matching-tolerance"]) * 1e-6,
        "grouping_error_tolerance": float(user_parameters["peak-grouping-tolerance"]) * 1e-6,
        "minimum_scan_count": int(user_parameters.get("minimum-scan-count", 1)),
        "database_path": manager.path,
    }

    input_type, input_id = user_parameters.get("hypothesis_choice").split(",")
    input_id = int(input_id)

    job_parameters["search_type"] = input_type
    job_parameters["hypothesis_id"] = input_id

    # Handle mass shifts --
    mass_shift_names = user_parameters.getlist("mass_shift_name")
    mass_shift_mass_deltas = user_parameters.getlist("mass_shift_mass_delta")
    mass_shift_max_counts = user_parameters.getlist("mass_shift_max_count")
    mass_shifts_tuplets = zip(mass_shift_names, mass_shift_mass_deltas, mass_shift_max_counts)

    mass_shift_map = {}
    for name, mass, max_count in mass_shifts_tuplets:
        if name == "" or mass == "" or max_count == "":
            continue
        mass = float(mass)
        max_count = int(max_count)
        obj = MassShift.get(g.db, name=name, mass=mass)
        g.db.add(obj)
        g.db.commit()
        mass_shift_map[obj.id] = max_count

    job_parameters["mass_shift_map"] = mass_shift_map
    # --

    for sample_name in request.values.getlist('samples'):
        instance_parameters = job_parameters.copy()
        sample_run, sample_manager = manager.find_sample(sample_name)
        instance_parameters["observed_ions_path"] = sample_manager.path
        instance_parameters["sample_run_id"] = sample_run.id
        instance_parameters['callback'] = lambda: 0
        instance_parameters['observed_ions_type'] = 'db'
        task = LCMSSearchTask(**instance_parameters)
        manager.add_task(task)
    return jsonify(**dict(request.values))

# ----------------------------------------
#
# ----------------------------------------


@app.route("/tandem_match_samples")
def tandem_match_samples():
    return render_template("tandem_glycopeptide_search/tandem_match_samples.templ")


@app.route("/tandem_match_samples", methods=["POST"])
def tandem_match_samples_post():
    user_parameters = request.values
    job_parameters = {
        "ms1_tolerance": float(user_parameters["ms1-tolerance"]) * 1e-6,
        "ms2_tolerance": float(user_parameters["ms2-tolerance"]) * 1e-6,
        "database_path": manager.path

    }
    input_type, input_id = user_parameters.get("hypothesis_choice").split(",")

    input_id = int(input_id)
    if input_type == "Hypothesis":
        defer = False
    else:
        defer = True

    db = manager.session()
    if not defer:
        job_parameters["target_hypothesis_id"] = input_id
        target_hypothesis = db.query(Hypothesis).get(input_id)
        decoy_id = target_hypothesis.parameters["decoys"][0]["hypothesis_id"]
        job_parameters['decoy_hypothesis_id'] = decoy_id
        job_parameters['source_hypothesis_sample_match_id'] = None
    else:
        job_parameters['source_hypothesis_sample_match_id'] = input_id
        job_parameters['target_hypothesis_id'] = None
        job_parameters['decoy_hypothesis_id'] = None

    print request.values.__dict__

    for sample_name in request.values.getlist('samples'):
        instance_parameters = job_parameters.copy()
        sample_run, sample_manager = manager.find_sample(sample_name)
        instance_parameters["observed_ions_path"] = sample_manager.path
        instance_parameters["sample_run_id"] = sample_run.id
        instance_parameters['callback'] = lambda: 0
        instance_parameters['observed_ions_type'] = 'db'
        print instance_parameters
        task = TandemMSGlycoproteomicsSearchTask(**instance_parameters)
        manager.add_task(task)
    return jsonify(**dict(request.values))

# ----------------------------------------
# Sample Management
# ----------------------------------------


@app.route("/add_sample", methods=["POST"])
def post_add_sample():
    """Handle an uploaded sample file

    Returns
    -------
    TYPE : Description
    """
    print request.values
    run_name = request.values['sample_name']
    secure_name = secure_filename(run_name)

    secure_name += ".%s" % request.values["file-type"]
    path = manager.get_temp_path(secure_name)
    request.files['observed-ions-file'].save(path)
    dest = manager.get_sample_path(run_name)
    # Construct the task with a callback to add the processed sample
    # to the set of project samples
    callback = functools.partial(manager.add_sample, path=dest)
    task_type = None
    print request.values["file-type"], type(request.values["file-type"])
    if request.values["file-type"] == "decon2ls":
        task_type = Decon2LSIsosParseTask
    elif request.values["file-type"] == "bupid":
        task_type = BUPIDYamlParseTask
    task = task_type(
        manager.path,
        path,
        dest,
        callback)
    manager.add_task(task)
    return redirect("/")


@app.route("/add_sample")
def add_sample():
    return render_template("add_sample_form.templ")


# ----------------------------------------
# Server Shutdown
# ----------------------------------------

def shutdown_server():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()


@app.route('/internal/shutdown', methods=['POST'])
def shutdown():
    shutdown_server()
    manager.stoploop()
    return 'Server shutting down...'

# ----------------------------------------
#
# ----------------------------------------


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
        "PeakGroupMatch": PeakGroupMatch,
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
parser.add_argument("results_database")
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


def main():
    args = parser.parse_args()
    results_database = args.results_database
    global DATABASE, manager, CAN_EXECUTE
    host = None
    if args.external:
        host = "0.0.0.0"
    DATABASE = results_database
    CAN_EXECUTE = not args.no_execute_tasks
    manager = ProjectManager(DATABASE)
    app.debug = DEBUG
    app.secret_key = SECRETKEY
    setup_logging()
    app.run(host=host, use_reloader=False, threaded=True, debug=DEBUG, port=args.port)

if __name__ == "__main__":
    main()
