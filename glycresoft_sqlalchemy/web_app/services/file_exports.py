import os
import base64
import functools

from flask import Response, g, Blueprint, request, jsonify

from glycresoft_sqlalchemy.data_model import (
    HypothesisSampleMatch, MS1GlycanHypothesisSampleMatch, MS1GlycopeptideHypothesisSampleMatch,
    MS2GlycopeptideHypothesisSampleMatch, GlycopeptideMatch, PeakGroupMatchType)

from glycresoft_sqlalchemy.web_app.task.do_export_csv import ExportCSVTask

file_exports = Blueprint("file_exports", __name__)


@file_exports.route("/internal/file_download/<b64path>")
def download_file(b64path):
    path = base64.b64decode(b64path)
    name = os.path.basename(path)
    if g.manager.temp_dir in path:
        def yielder():
            for line in open(path):
                yield line
        return Response(yielder(), mimetype="application/octet-stream",
                        headers={"Content-Disposition": "attachment; filename=%s;" % name})


def filterfunc_template(q, model, attr, value, code=1):
    return q.filter(getattr(model, attr) >= value)


@file_exports.route("/view_database_search_results/export_csv/<int:id>", methods=["POST"])
def export_csv_task(id):
    hypothesis_sample_match = g.db.query(HypothesisSampleMatch).get(id)
    tempdir = g.manager.get_temp_path("glycresoft_export")

    state = request.get_json()
    settings = state["settings"]
    context = state['context']
    minimum_score = settings.get("minimum_ms2_score", 0.2)

    print hypothesis_sample_match, id

    if isinstance(hypothesis_sample_match, (MS1GlycopeptideHypothesisSampleMatch,
                                            MS1GlycanHypothesisSampleMatch)):
        filterfunc = functools.partial(
            filterfunc_template, model=PeakGroupMatchType, attr="ms1_score", value=minimum_score)
    elif isinstance(hypothesis_sample_match, MS2GlycopeptideHypothesisSampleMatch):
        filterfunc = functools.partial(
            filterfunc_template, model=GlycopeptideMatch, attr="ms2_score", value=minimum_score)

    task = ExportCSVTask(g.manager.path, id, filterfunc, tempdir)

    g.manager.add_task(task)
    return jsonify(target=hypothesis_sample_match.to_json())
