import functools
from werkzeug import secure_filename
from flask import Blueprint, g, request, Response, render_template, redirect

from glycresoft_sqlalchemy.web_app.task.do_decon2ls_parse import Decon2LSIsosParseTask
from glycresoft_sqlalchemy.web_app.task.do_bupid_yaml_parse import BUPIDYamlParseTask


sample_management = Blueprint("sample_management", __name__)


@sample_management.route("/add_sample", methods=["POST"])
def post_add_sample():
    """Handle an uploaded sample file

    Returns
    -------
    TYPE : Description
    """
    run_name = request.values['sample_name']
    if run_name == "":
        run_name = request.files['observed-ions-file'].filename
    secure_name = secure_filename(run_name)

    secure_name += ".%s" % request.values["file-type"]
    path = g.manager.get_temp_path(secure_name)
    request.files['observed-ions-file'].save(path)
    dest = g.manager.get_sample_path(run_name)
    # Construct the task with a callback to add the processed sample
    # to the set of project samples
    callback = functools.partial(g.manager.add_sample, path=dest)
    task_type = None

    if request.values["file-type"] == "decon2ls":
        task_type = Decon2LSIsosParseTask
    elif request.values["file-type"] == "bupid":
        task_type = BUPIDYamlParseTask
    task = task_type(
        g.manager.path,
        path,
        dest,
        callback)
    g.manager.add_task(task)
    return redirect("/")


@sample_management.route("/add_sample")
def add_sample():
    return render_template("add_sample_form.templ")
