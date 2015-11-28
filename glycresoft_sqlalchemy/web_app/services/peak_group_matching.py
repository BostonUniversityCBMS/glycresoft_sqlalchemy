from flask import request, g, render_template, Blueprint, jsonify

from glycresoft_sqlalchemy.data_model import MassShift

from glycresoft_sqlalchemy.web_app.task.do_ms1_peak_group_matching import LCMSSearchTask


app = peak_group_matching = Blueprint("peak_group_matching", __name__)


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
        "database_path": g.manager.path,
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
        sample_run, sample_manager = g.manager.find_sample(sample_name)
        instance_parameters["observed_ions_path"] = sample_manager.path
        instance_parameters["sample_run_id"] = sample_run.id
        instance_parameters['callback'] = lambda: 0
        instance_parameters['observed_ions_type'] = 'db'
        task = LCMSSearchTask(**instance_parameters)
        g.manager.add_task(task)
    return jsonify(**dict(request.values))
