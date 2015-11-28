from flask import request, g, render_template, Blueprint, jsonify

from glycresoft_sqlalchemy.data_model import Hypothesis
from glycresoft_sqlalchemy.web_app.task.do_ms2_search import TandemMSGlycoproteomicsSearchTask

app = tandem_glycopeptide_search = Blueprint("tandem_glycopeptide_search", __name__)


@app.route("/tandem_match_samples")
def tandem_match_samples():
    return render_template("tandem_glycopeptide_search/tandem_match_samples.templ")


@app.route("/tandem_match_samples", methods=["POST"])
def tandem_match_samples_post():
    user_parameters = request.values
    job_parameters = {
        "ms1_tolerance": float(user_parameters["ms1-tolerance"]) * 1e-6,
        "ms2_tolerance": float(user_parameters["ms2-tolerance"]) * 1e-6,
        "database_path": g.manager.path

    }
    input_type, input_id = user_parameters.get("hypothesis_choice").split(",")

    input_id = int(input_id)
    if input_type == "Hypothesis":
        defer = False
    else:
        defer = True

    db = g.manager.session()
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
        sample_run, sample_manager = g.manager.find_sample(sample_name)
        instance_parameters["observed_ions_path"] = sample_manager.path
        instance_parameters["sample_run_id"] = sample_run.id
        instance_parameters['callback'] = lambda: 0
        instance_parameters['observed_ions_type'] = 'db'
        print instance_parameters
        task = TandemMSGlycoproteomicsSearchTask(**instance_parameters)
        g.manager.add_task(task)
    return jsonify(**dict(request.values))
