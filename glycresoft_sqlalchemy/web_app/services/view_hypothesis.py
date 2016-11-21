import logging
from flask import request, g, render_template, Response, Blueprint, jsonify

from glycresoft_sqlalchemy.data_model import (Hypothesis, MS1GlycanHypothesis, MS1GlycopeptideHypothesis,
                                              Protein, MS2GlycopeptideHypothesis)

from glycresoft_sqlalchemy.web_app.utils.state_transfer import request_arguments_and_context
from glycresoft_sqlalchemy.web_app.utils.pagination import paginate

app = view_hypothesis = Blueprint("view_hypothesis", __name__)


@app.route("/view_hypothesis/<int:id>/mass_search", methods=["POST"])
def mass_search(id):
    arguments, state = request_arguments_and_context(request)
    hypothesis = g.db.query(Hypothesis).get(id)
    mass = arguments['mass']
    tolerance = arguments['tolerance']
    output_format = arguments.get("format", 'json')
    hypothesis = g.db.query(Hypothesis).get(id)
    matches = hypothesis.search_by_mass(mass, tolerance)
    if output_format == "json":
        results = [{
            "mass": r.calculated_mass,
            "sequence": r.most_detailed_sequence,
            "type": hypothesis.theoretical_structure_type.__name__,
            "id": r.id
        } for r in matches]
        return jsonify(results)


@app.route("/view_hypothesis/<int:id>", methods=["POST"])
def view_hypothesis_dispatch(id):
    try:
        print id, request.get_json()
        hypothesis = g.db.query(Hypothesis).get(id)
        if isinstance(hypothesis, MS1GlycanHypothesis):
            return view_glycan_composition_hypothesis(id)
        elif isinstance(hypothesis, MS1GlycopeptideHypothesis):
            return view_glycopeptide_composition_hypothesis(id)
        elif isinstance(hypothesis, MS2GlycopeptideHypothesis):
            return view_glycopeptide_hypothesis(id)
    except Exception, e:
        logging.exception("An exception occurred for %r", hypothesis, exc_info=e)
    return Response("<h2>No display method is implemented for %s </h2>" % hypothesis.__class__.__name__)


@app.route("/view_glycan_composition_hypothesis/<int:id>/", methods=["POST"])
def view_glycan_composition_hypothesis(id):
    state = request.get_json()
    settings = state['settings']
    context = state['context']
    hypothesis = g.db.query(Hypothesis).get(id)
    return render_template("view_glycan_hypothesis/container.templ", hypothesis=hypothesis)


@app.route("/view_glycan_composition_hypothesis/<int:id>/<int:page>", methods=["POST"])
def view_glycan_composition_hypothesis_table(id, page=1):
    state = request.get_json()
    settings = state['settings']
    context = state['context']
    hypothesis = g.db.query(Hypothesis).get(id)

    page_size = 50

    def filter_context(q):
        return q.filter_by(
            hypothesis_id=id)
    paginator = paginate(filter_context(hypothesis.glycans), page, page_size)
    return render_template(
        "view_glycan_hypothesis/display_table.templ",
        paginator=paginator, base_index=(page - 1) * page_size)


@app.route("/view_glycopeptide_composition_hypothesis/<int:id>/", methods=["POST"])
def view_glycopeptide_composition_hypothesis(id):
    state = request.get_json()
    settings = state['settings']
    context = state['context']
    hypothesis = g.db.query(Hypothesis).get(id)
    return render_template("view_glycopeptide_composition_hypothesis/container.templ", hypothesis=hypothesis)


@app.route("/view_glycopeptide_composition_hypothesis/protein_view/<int:id>", methods=["POST"])
def view_glycopeptide_composition_hypothesis_protein_view(id):
    protein = g.db.query(Protein).get(id)
    return render_template(
        "view_glycopeptide_composition_hypothesis/components/protein_view.templ",
        protein=protein)


@app.route("/view_glycopeptide_composition_hypothesis/protein_view/<int:id>/<int:page>", methods=["POST"])
def view_glycopeptide_composition_hypothesis_glycopeptide_table(id, page=1):
    state = request.get_json()
    settings = state['settings']
    context = state['context']
    protein = g.db.query(Protein).get(id)

    page_size = 200

    def filter_context(q):
        return q.filter_by(
            protein_id=id)
    paginator = paginate(filter_context(protein.theoretical_glycopeptide_compositions), page, page_size)
    return render_template(
        "view_glycopeptide_composition_hypothesis/components/display_table.templ",
        paginator=paginator, base_index=(page - 1) * page_size)


@app.route("/view_glycopeptide_hypothesis/<int:id>/", methods=["POST"])
def view_glycopeptide_hypothesis(id):
    # state = request.get_json()
    # settings = state['settings']
    # context = state['context']
    hypothesis = g.db.query(Hypothesis).get(id)
    print hypothesis.proteins
    return render_template("view_glycopeptide_hypothesis/container.templ", hypothesis=hypothesis)


@app.route("/view_glycopeptide_hypothesis/protein_view/<int:id>", methods=["POST"])
def view_glycopeptide_hypothesis_protein_view(id):
    protein = g.db.query(Protein).get(id)
    return render_template(
        "view_glycopeptide_hypothesis/components/protein_view.templ",
        protein=protein)


@app.route("/view_glycopeptide_hypothesis/protein_view/<int:id>/<int:page>", methods=["POST"])
def view_glycopeptide_hypothesis_glycopeptide_table(id, page=1):
    state = request.get_json()
    settings = state['settings']
    context = state['context']
    protein = g.db.query(Protein).get(id)

    page_size = 200

    def filter_context(q):
        return q.filter_by(
            protein_id=id)
    paginator = paginate(filter_context(protein.theoretical_glycopeptides), page, page_size)
    return render_template(
        "view_glycopeptide_hypothesis/components/display_table.templ",
        paginator=paginator, base_index=(page - 1) * page_size)
