from flask import Response, Blueprint, g, jsonify

from glycresoft_sqlalchemy.data_model import GlycopeptideMatch, Hypothesis, HypothesisSampleMatch, json_type
from glycresoft_sqlalchemy.report import colors


JSONEncoderType = json_type.new_alchemy_encoder()

# ----------------------------------------
#           JSON Data API Calls
# ----------------------------------------

api = Blueprint("api", __name__)


@api.route("/api/glycopeptide_matches/<int:id>")
def get_glycopeptide_match_api(id):
    gpm = g.db.query(GlycopeptideMatch).get(id)
    return Response(JSONEncoderType().encode(gpm), mimetype="text/json")


@api.route("/api/tasks")
def api_tasks():
    return jsonify(**{t.id: t.to_json() for t in g.manager.tasks.values()})


@api.route("/api/hypothesis_sample_matches")
def api_hypothesis_sample_matches():
    hsms = g.db.query(HypothesisSampleMatch).all()
    d = {str(h.id): h.to_json() for h in hsms}
    return jsonify(**d)


@api.route("/api/hypotheses")
def api_hypothesis():
    hypotheses = g.db.query(Hypothesis).all()
    d = {str(h.id): h.to_json() for h in hypotheses}
    return jsonify(**d)


@api.route("/api/samples")
def api_samples():
    samples = g.manager.samples()
    d = {str(h.name): h.to_json() for h in samples}
    return jsonify(**d)


@api.route("/api/colors")
def api_colors():
    return jsonify(**colors.color_dict())
