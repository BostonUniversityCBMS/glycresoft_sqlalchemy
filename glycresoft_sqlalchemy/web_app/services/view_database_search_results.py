from flask import request, g, render_template, Response, Blueprint

from glycresoft_sqlalchemy.data_model import (
    HypothesisSampleMatch, GlycopeptideMatch, PeakGroupMatchType, Protein,
    JointPeakGroupMatch, MS1GlycopeptideHypothesisSampleMatch,
    MS1GlycanHypothesisSampleMatch)
from glycresoft_sqlalchemy.report import microheterogeneity
from glycresoft_sqlalchemy.web_app.utils.pagination import paginate

view_database_search_results = Blueprint("view_database_search_results", __name__)

app = view_database_search_results


@app.route("/view_database_search_results/<int:id>", methods=["POST"])
def view_database_search_results_dispatch(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    results_type = hsm.results_type
    if PeakGroupMatchType == results_type:
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

    if parameters['settings'].get("color_palette") == "NGlycanCompositionColorizer":
        palette = microheterogeneity.NGlycanCompositionColorizer
    else:
        palette = microheterogeneity._null_color_chooser

    site_summary = microheterogeneity.GlycoproteinMicroheterogeneitySummary(
        protein, filter_context, color_chooser=palette)

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
            PeakGroupMatchType.ms1_score > 0.2)

    if isinstance(hsm, MS1GlycopeptideHypothesisSampleMatch):
        return render_template(
            "glycopeptide_peak_group_search/view_database_search_results.templ",
            hsm=hsm,
            filter_context=filter_context)
    elif isinstance(hsm, MS1GlycanHypothesisSampleMatch):
        return render_template(
            "view_glycan_peak_group_search_results/view_database_search_results.templ",
            hsm=hsm,
            filter_context=filter_context)


# ----------------------------------------------------------------------
#           Glycopeptide Compostion Peak Grouping Database Search Results
# ----------------------------------------------------------------------


@app.route("/view_database_search_results/protein_composition_view/<int:id>", methods=["POST"])
def view_composition_glycopeptide_protein_results(id):
    print
    print "view_composition_glycopeptide_protein_results"
    print
    print request.values
    hypothesis_sample_match_id = request.values["hypothesis_sample_match_id"]
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
    protein = g.db.query(Protein).get(id)

    # parameters = request.get_json()
    # print parameters
    # print id

    # monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})
    # minimum_score = float(parameters['settings'].get("minimum_ms1_score", 0.2))
    minimum_score = 0.1

    def filter_context(q):

        return q.filter(
            PeakGroupMatchType.hypothesis_sample_match_id == hypothesis_sample_match_id,
            PeakGroupMatchType.ms1_score > minimum_score)

    return render_template(
        "glycopeptide_peak_group_search/components/protein_view.templ",
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
        q = theoretical_type.glycan_composition_filters(q, monosaccharide_filters)
        return q.filter(
            results_type.hypothesis_sample_match_id == hypothesis_sample_match_id,
            results_type.ms1_score > minimum_score)

    paginator = paginate(
        filter_context(
            protein.peak_group_matches).order_by(
            PeakGroupMatchType.ms1_score.desc()), page, 50)

    return render_template(
        "glycopeptide_peak_group_search/components/glycopeptide_match_table.templ",
        paginator=paginator)


@app.route("/view_database_search_results/view_glycopeptide_composition_details/<int:id>")
def view_peak_grouping_glycopeptide_composition_details(id):
    pgm = g.db.query(PeakGroupMatchType).get(id)
    ambiguous_with = g.db.query(PeakGroupMatchType).filter(
        PeakGroupMatchType.fingerprint == pgm.fingerprint,
        PeakGroupMatchType.id != pgm.id).all()
    return render_template(
        "glycopeptide_peak_group_search/components/glycopeptide_details.templ", pgm=pgm,
        ambiguous_with=ambiguous_with)


# ----------------------------------------------------------------------
#           Glycopeptide Compostion Peak Grouping Database Search Results
# ----------------------------------------------------------------------


@app.route("/view_database_search_results/results_view/", methods=["POST"])
def view_glycan_composition_results():

    print request.values
    print request
    parameters = request.get_json()
    print parameters
    hypothesis_sample_match_id = parameters['context']["hypothesis_sample_match_id"]
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

    monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})
    minimum_score = float(parameters['settings'].get("minimum_ms1_score", 0.2))

    def filter_context(q):
        results_type = hsm.results_type
        theoretical_type, generator = hsm.results().next()
        q = q.join(theoretical_type, theoretical_type.id == results_type.theoretical_match_id)
        q = theoretical_type.glycan_composition_filters(q, monosaccharide_filters)
        return q.filter(
            results_type.hypothesis_sample_match_id == hypothesis_sample_match_id,
            results_type.ms1_score > minimum_score)

    GlycoformAbundancePlot = microheterogeneity.GlycoformAbundancePlot
    if parameters['settings'].get("color_palette") == "NGlycanCompositionColorizer":
        palette = microheterogeneity.NGlycanCompositionColorizer
    else:
        palette = microheterogeneity._null_color_chooser

    query = filter_context(hsm.peak_group_matches)
    plot = GlycoformAbundancePlot(query, color_chooser=palette)

    return render_template(
        "view_glycan_peak_group_search_results/components/results_view.templ",
        hypothesis_sample_match=hsm,
        filter_context=filter_context,
        plot=plot)


@app.route("/view_database_search_results/glycan_composition_match_table/"
           "<int:page>", methods=["POST"])
def view_composition_glycan_table_partial(page):
    hypothesis_sample_match_id = request.get_json()["context"]["hypothesis_sample_match_id"]
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

    parameters = request.get_json()
    print parameters
    print page

    monosaccharide_filters = parameters['settings'].get("monosaccharide_filters", {})
    minimum_score = float(parameters['settings'].get("minimum_ms1_score", 0.2))

    def filter_context(q):
        results_type = hsm.results_type
        theoretical_type, generator = hsm.results().next()
        q = q.join(theoretical_type, theoretical_type.id == results_type.theoretical_match_id)
        q = theoretical_type.glycan_composition_filters(q, monosaccharide_filters)
        return q.filter(
            results_type.hypothesis_sample_match_id == hypothesis_sample_match_id,
            results_type.ms1_score > minimum_score)

    paginator = paginate(
        filter_context(
            hsm.peak_group_matches).order_by(
            PeakGroupMatchType.ms1_score.desc()), page, 50)

    return render_template(
        "view_glycan_peak_group_search_results/components/glycan_composition_match_table.templ",
        paginator=paginator)


@app.route("/view_database_search_results/view_glycan_composition_details/<int:id>")
def view_peak_grouping_glycan_composition_details(id):
    pgm = g.db.query(PeakGroupMatchType).get(id)
    ambiguous_with = g.db.query(PeakGroupMatchType).filter(
        PeakGroupMatchType.fingerprint == pgm.fingerprint,
        PeakGroupMatchType.id != pgm.id).all()
    return render_template(
        "view_glycan_peak_group_search_results/components/glycan_composition_details.templ", pgm=pgm,
        ambiguous_with=ambiguous_with)
