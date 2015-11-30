from flask import request, g, render_template, Response, Blueprint

from glycresoft_sqlalchemy.data_model import (
    HypothesisSampleMatch, GlycopeptideMatch, PeakGroupMatch, Protein,
    JointPeakGroupMatch)
from glycresoft_sqlalchemy.report import microheterogeneity
from glycresoft_sqlalchemy.web_app.utils.pagination import paginate

view_database_search_results = Blueprint("view_database_search_results", __name__)

app = view_database_search_results

PeakGroupType = PeakGroupMatch


@app.route("/view_database_search_results/<int:id>", methods=["POST"])
def view_database_search_results_dispatch(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    results_type = hsm.results_type
    if PeakGroupType == results_type:
        return view_composition_database_search_results(id)
    elif GlycopeptideMatch == results_type:
        return view_tandem_glycopeptide_database_search_results(id)
    return Response("No Match, " + results_type.__name__)


# ----------------------------------------------------------------------
#           View Tandem Glycopeptide Database Search Results
# ----------------------------------------------------------------------
def view_tandem_glycopeptide_database_search_results(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    hypothesis_sample_match_id = id

    state = request.get_json()
    settings = state['settings']
    context = state['context']

    minimum_score = settings.get('minimum_ms2_score', 0.2)
    monosaccharide_filters = settings.get("monosaccharide_filters", {})

    def filter_context(q):
        return GlycopeptideMatch.glycan_composition_filters(q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            GlycopeptideMatch.ms2_score > minimum_score), monosaccharide_filters)

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

    monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})

    minimum_score = float(parameters['settings'].get("minimum_ms2_score", 0.2))

    def filter_context(q):
        return GlycopeptideMatch.glycan_composition_filters(q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            GlycopeptideMatch.ms2_score > minimum_score), monosaccharide_filters)

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
            PeakGroupType.ms1_score > 0.2)

    return render_template(
        "peak_group_search/view_database_search_results.templ",
        hsm=hsm,
        filter_context=filter_context)


@app.route("/view_database_search_results/protein_composition_view/<int:id>", methods=["POST"])
def view_composition_glycopeptide_protein_results(id):
    print request.values
    hypothesis_sample_match_id = request.values["hypothesis_sample_match_id"]
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
    protein = g.db.query(Protein).get(id)

    parameters = request.get_json()
    print parameters
    print id

    monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})
    minimum_score = float(parameters['settings'].get("minimum_ms1_score", 0.2))

    def filter_context(q):

        return q.filter(
            PeakGroupType.hypothesis_sample_match_id == hypothesis_sample_match_id,
            PeakGroupType.ms1_score > minimum_score)

    return render_template(
        "peak_group_search/components/protein_view.templ",
        protein=protein,
        filter_context=filter_context)


@app.route("/view_database_search_results/glycopeptide_matches_composition_table"
           "/<int:protein_id>/<int:page>", methods=["POST"])
def view_composition_glycopeptide_table_partial(protein_id, page):
    hypothesis_sample_match_id = request.get_json()["context"]["hypothesis_sample_match_id"]
    protein = g.db.query(Protein).get(protein_id)
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

    parameters = request.get_json()
    print parameters
    print id

    monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})
    minimum_score = float(parameters['settings'].get("minimum_ms1_score", 0.2))

    def filter_context(q):
        results_type = hsm.results_type
        theoretical_type, generator = hsm.results().next()
        q = q.join(theoretical_type, results_type.theoretical_match_id == theoretical_type.id)
        q = theoretical_type.glycan_composition_filters(q, monosaccharide_filters)
        return q.filter(
            results_type.hypothesis_sample_match_id == hypothesis_sample_match_id,
            results_type.ms1_score > minimum_score)

    paginator = paginate(filter_context(protein.peak_group_matches).order_by(PeakGroupType.ms1_score.desc()), page, 50)

    return render_template(
        "peak_group_search/components/glycopeptide_match_table.templ",
        paginator=paginator)


@app.route("/view_database_search_results/view_glycopeptide_composition_details/<int:id>")
def view_peak_grouping_glycopeptide_composition_details(id):
    pgm = g.db.query(PeakGroupType).get(id)
    ambiguous_with = g.db.query(PeakGroupType).filter(
        PeakGroupType.peak_group_id == pgm.peak_group_id,
        PeakGroupType.id != pgm.id).all()
    return render_template(
        "peak_group_search/components/glycopeptide_details.templ", pgm=pgm,
        ambiguous_with=ambiguous_with)
