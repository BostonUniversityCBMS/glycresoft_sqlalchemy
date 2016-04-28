from flask import request, g, render_template, Response, Blueprint, jsonify

from glycresoft_sqlalchemy.data_model import (
    HypothesisSampleMatch, GlycopeptideMatch, PeakGroupMatchType, Protein,
    MS1GlycopeptideHypothesisSampleMatch, TheoreticalGlycopeptide,
    MS1GlycanHypothesisSampleMatch, TheoreticalGlycopeptideComposition, func)

from glycresoft_sqlalchemy.report import analysis_comparison
from glycresoft_sqlalchemy.report import microheterogeneity
from glycresoft_sqlalchemy.web_app.utils.pagination import paginate
from glycresoft_sqlalchemy.web_app.report import svg_plot
from glycresoft_sqlalchemy.web_app.utils.state_transfer import request_arguments_and_context, MonosaccharideFilterSet
from glycresoft_sqlalchemy.web_app.utils.cache import CachablePartialFunction, ApplicationDataCache

view_database_search_results = Blueprint("view_database_search_results", __name__)

app = view_database_search_results

microheterogeneity_summary_cache = ApplicationDataCache()
protein_table_cache = ApplicationDataCache(1000)


@app.route("/view_database_search_results/<int:id>", methods=["POST"])
def view_database_search_results_dispatch(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    results_type = hsm.results_type
    if PeakGroupMatchType == results_type:
        return view_composition_database_search_results(id)
    elif GlycopeptideMatch == results_type:
        return view_tandem_glycopeptide_database_search_results(id)
    return Response("No Match, " + results_type.__name__)


def get_protein_from_tabulated_table(table, protein_id):
    for row in table:
        if row['protein'].id == protein_id:
            return row


# ----------------------------------------------------------------------
#           View Tandem Glycopeptide Database Search Results
# ----------------------------------------------------------------------

def view_tandem_glycopeptide_protein_results_filter_context(
        q, hypothesis_sample_match_id, minimum_score, monosaccharide_filters):
    return GlycopeptideMatch.glycan_composition_filters(q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            GlycopeptideMatch.ms2_score > minimum_score), monosaccharide_filters)


def protein_table_data_composer_ms2(session, hypothesis_id, hypothesis_sample_match_id, filter_context=lambda q: q):
    table = dict()

    q = session.query(
        Protein, func.count(GlycopeptideMatch.id)).join(
        GlycopeptideMatch).filter(
        GlycopeptideMatch.hypothesis_sample_match_id == hypothesis_sample_match_id,
        Protein.hypothesis_id == hypothesis_id).group_by(
        Protein.id)

    glycopeptide_match_counts = filter_context(q).all()

    theoretical_sequence_counts = session.query(
        Protein, func.count(TheoreticalGlycopeptide.id)).join(
        TheoreticalGlycopeptide).filter(
        Protein.hypothesis_id == hypothesis_id).group_by(
        Protein.id).all()

    for protein, theoretical_count in theoretical_sequence_counts:
        if theoretical_count == 0:
            continue
        table[protein.id] = {
            "protein": protein,
            "theoretical_count": theoretical_count,
            "match_count": 0
        }

    for protein, match_count in glycopeptide_match_counts:
        table[protein.id]['match_count'] = match_count

    return sorted(table.values(), key=lambda x: x['theoretical_count'], reverse=True)


def view_tandem_glycopeptide_database_search_results(id):
    hsm = g.db.query(HypothesisSampleMatch).get(id)
    hypothesis_sample_match_id = id

    state = request.get_json()
    settings = state['settings']
    context = state['context']

    minimum_score = float(settings.get('minimum_ms2_score', 0.2))
    monosaccharide_filters = MonosaccharideFilterSet.fromdict(settings.get("monosaccharide_filters", {}))

    filter_context = CachablePartialFunction(
        view_tandem_glycopeptide_protein_results_filter_context,
        hypothesis_sample_match_id=hypothesis_sample_match_id,
        minimum_score=minimum_score,
        monosaccharide_filters=monosaccharide_filters
        )
    protein_table = protein_table_data_composer_ms2(g.db, hsm.target_hypothesis_id, hsm.id, filter_context)

    protein_table_cache[hsm.id] = protein_table

    template = render_template(
        "tandem_glycopeptide_search/view_database_search_results.templ",
        hsm=hsm,
        filter_context=filter_context,
        protein_table=protein_table)
    return template


@app.route("/view_database_search_results/protein_view/<int:id>", methods=["POST"])
def view_tandem_glycopeptide_protein_results(id):
    parameters = request.get_json()

    hypothesis_sample_match_id = parameters['context']['hypothesis_sample_match_id']
    protein = g.db.query(Protein).get(id)

    monosaccharide_filters = MonosaccharideFilterSet.fromdict(
        parameters['settings'].get("monosaccharide_filters", {}))

    minimum_score = float(parameters['settings'].get("minimum_ms2_score", 0.2))

    filter_context = CachablePartialFunction(
        view_tandem_glycopeptide_protein_results_filter_context,
        hypothesis_sample_match_id=hypothesis_sample_match_id,
        minimum_score=minimum_score,
        monosaccharide_filters=monosaccharide_filters
        )

    if parameters['settings'].get("color_palette") == "NGlycanCompositionColorizer":
        palette = microheterogeneity.NGlycanCompositionColorizer
    else:
        palette = microheterogeneity._null_color_chooser

    site_summary = microheterogeneity_summary_cache.cache(
        protein, microheterogeneity.GlycoproteinMicroheterogeneitySummary,
        protein, filter_context, color_chooser=palette)
    # site_summary = microheterogeneity.GlycoproteinMicroheterogeneitySummary(
    #     protein, filter_context, color_chooser=palette)

    return render_template(
        "tandem_glycopeptide_search/components/protein_view.templ",
        protein=protein,
        site_summary=site_summary,
        filter_context=filter_context)


@app.route("/view_database_search_results/glycopeptide_match_table"
           "/<int:protein_id>/<int:page>", methods=["POST"])
def view_glycopeptide_table_partial(protein_id, page):
    arguments, state = request_arguments_and_context(request)

    hypothesis_sample_match_id = state.hypothesis_sample_match_id
    protein = g.db.query(Protein).get(protein_id)

    monosaccharide_filters = state.monosaccharide_filters
    minimum_score = float(state.minimum_ms2_score)

    filter_context = CachablePartialFunction(
        view_tandem_glycopeptide_protein_results_filter_context,
        hypothesis_sample_match_id=hypothesis_sample_match_id,
        minimum_score=minimum_score,
        monosaccharide_filters=monosaccharide_filters
        )

    paginator = paginate(
        filter_context(
            protein.glycopeptide_matches).order_by(
            GlycopeptideMatch.ms2_score.desc()), page, 50)

    return render_template(
        "tandem_glycopeptide_search/components/glycopeptide_match_table.templ",
        paginator=paginator)


@app.route("/view_database_search_results/protein_view/<int:id>/microheterogeneity_plot_panel", methods=["POST"])
def microheterogeneity_plot_panel(id):
    parameters = request.get_json()

    hypothesis_sample_match_id = parameters['context']['hypothesis_sample_match_id']
    protein = g.db.query(Protein).get(id)

    monosaccharide_filters = MonosaccharideFilterSet.fromdict(
        parameters['settings'].get("monosaccharide_filters", {}))

    minimum_score = float(parameters['settings'].get("minimum_ms2_score", 0.2))

    filter_context = CachablePartialFunction(
        view_tandem_glycopeptide_protein_results_filter_context,
        hypothesis_sample_match_id=hypothesis_sample_match_id,
        minimum_score=minimum_score,
        monosaccharide_filters=monosaccharide_filters
        )

    if parameters['settings'].get("color_palette") == "NGlycanCompositionColorizer":
        palette = microheterogeneity.NGlycanCompositionColorizer
    else:
        palette = microheterogeneity._null_color_chooser

    site_summary = microheterogeneity_summary_cache.cache(
        protein, microheterogeneity.GlycoproteinMicroheterogeneitySummary,
        protein, filter_context, color_chooser=palette)

    # site_summary = microheterogeneity.GlycoproteinMicroheterogeneitySummary(
    #     protein, filter_context, color_chooser=palette)

    return render_template(
        "tandem_glycopeptide_search/components/microheterogeneity_plot_panel.templ",
        site_summary=site_summary, protein=protein)


@app.route("/view_database_search_results/protein_view/<int:id>/protein_overview_panel", methods=["POST"])
def protein_overview_panel(id):
    arguments, state = request_arguments_and_context(request)

    hypothesis_sample_match_id = state.hypothesis_sample_match_id
    protein = g.db.query(Protein).get(id)

    monosaccharide_filters = state.monosaccharide_filters
    minimum_score = float(state.minimum_ms2_score)

    filter_context = CachablePartialFunction(
        view_tandem_glycopeptide_protein_results_filter_context,
        hypothesis_sample_match_id=hypothesis_sample_match_id,
        minimum_score=minimum_score,
        monosaccharide_filters=monosaccharide_filters
        )
    if state.color_palette == "NGlycanCompositionColorizer":
        palette = microheterogeneity.NGlycanCompositionColorizer
    else:
        palette = microheterogeneity._null_color_chooser
    site_summary = microheterogeneity_summary_cache.cache(
        protein, microheterogeneity.GlycoproteinMicroheterogeneitySummary,
        protein, filter_context, color_chooser=palette)
    return render_template(
        "tandem_glycopeptide_search/components/protein_overview_panel.templ",
        site_summary=site_summary, protein=protein, filter_context=filter_context)


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
        return render_glycopeptide_composition_database_search_results(hsm, id)
        # return render_template(
        #     "glycopeptide_peak_group_search/view_database_search_results.templ",
        #     hsm=hsm,
        #     filter_context=filter_context)
    elif isinstance(hsm, MS1GlycanHypothesisSampleMatch):
        return render_template(
            "view_glycan_peak_group_search_results/view_database_search_results.templ",
            hsm=hsm,
            filter_context=filter_context)


# ----------------------------------------------------------------------
#           Glycopeptide Compostion Peak Grouping Database Search Results
# ----------------------------------------------------------------------


def render_glycopeptide_composition_database_search_results(hsm, hypothesis_sample_match_id):
    def filter_context(q):
        return q.filter_by(
            hypothesis_sample_match_id=hypothesis_sample_match_id).filter(
            PeakGroupMatchType.ms1_score > 0.2)

    protein_table = protein_table_data_composer_ms1(g.db, hsm.target_hypothesis_id, hsm.id, filter_context)

    protein_table_cache[hypothesis_sample_match_id] = protein_table

    print protein_table

    return render_template(
        "glycopeptide_peak_group_search/view_database_search_results.templ",
        hsm=hsm,
        filter_context=filter_context,
        protein_table=protein_table)


def protein_table_data_composer_ms1(session, hypothesis_id, hypothesis_sample_match_id, filter_context=lambda q: q):
    table = dict()

    q = session.query(
        Protein, func.count(PeakGroupMatchType.id)).join(
        TheoreticalGlycopeptideComposition).join(
        PeakGroupMatchType,
        PeakGroupMatchType.theoretical_match_id == TheoreticalGlycopeptideComposition.id).filter(
        PeakGroupMatchType.hypothesis_sample_match_id == hypothesis_sample_match_id,
        Protein.hypothesis_id == hypothesis_id).group_by(
        Protein.id)

    glycopeptide_match_counts = filter_context(q).all()

    theoretical_sequence_counts = session.query(
        Protein, func.count(TheoreticalGlycopeptideComposition.id)).join(
        TheoreticalGlycopeptideComposition).filter(
        Protein.hypothesis_id == hypothesis_id).group_by(
        Protein.id).all()

    for protein, theoretical_count in theoretical_sequence_counts:
        if theoretical_count == 0:
            continue
        table[protein.id] = {
            "protein": protein,
            "theoretical_count": theoretical_count,
            "match_count": 0
        }

    for protein, match_count in glycopeptide_match_counts:
        table[protein.id]['match_count'] = match_count

    return sorted(table.values(), key=lambda x: x['theoretical_count'], reverse=True)


@app.route("/view_database_search_results/protein_composition_view/<int:id>", methods=["POST"])
def view_composition_glycopeptide_protein_results(id):
    arguments, context = request_arguments_and_context(request)
    hypothesis_sample_match_id = context.hypothesis_sample_match_id
    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)
    protein = g.db.query(Protein).get(id)
    minimum_score = context.minimum_ms1_score

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
    arguments, context = request_arguments_and_context(request)
    hypothesis_sample_match_id = context.hypothesis_sample_match_id

    hsm = g.db.query(HypothesisSampleMatch).get(hypothesis_sample_match_id)

    monosaccharide_filters = context.monosaccharide_filters
    minimum_score = context.minimum_ms1_score

    try:
        table = protein_table_cache[hypothesis_sample_match_id]
        table_data = get_protein_from_tabulated_table(table, protein_id)
        total_count = table_data['match_count']
    except Exception, e:
        total_count = None

    def filter_context(q):
        results_type = hsm.results_type
        theoretical_type, generator = hsm.results().next()
        q = theoretical_type.glycan_composition_filters(q, monosaccharide_filters)
        return q.filter(
            results_type.hypothesis_sample_match_id == hypothesis_sample_match_id,
            results_type.ms1_score > minimum_score)

    paginator = paginate(
        filter_context(
            g.db.query(PeakGroupMatchType).join(
                TheoreticalGlycopeptideComposition,
                TheoreticalGlycopeptideComposition.id == PeakGroupMatchType.theoretical_match_id).filter(
                TheoreticalGlycopeptideComposition.protein_id == protein_id,
                PeakGroupMatchType.hypothesis_sample_match_id == hypothesis_sample_match_id)).order_by(
            PeakGroupMatchType.ms1_score.desc()), page, 50, total_count)

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
#           Glycan Compostion Peak Grouping Database Search Results
# ----------------------------------------------------------------------


@app.route("/view_database_search_results/results_view/", methods=["POST"])
def view_glycan_composition_results():

    parameters = request.get_json()

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


# ----------------------------------------------------------------------
#           Hypothesis Sample Match Comparison
# ----------------------------------------------------------------------

@app.route("/view_database_search_results/compare")
def select_hypothesis_sample_match_checklist():
    hsms = g.db.query(HypothesisSampleMatch).all()
    return render_template(
        "components/hypothesis_sample_match_comparison_checklist.templ", hsms=hsms)


@app.route("/view_database_search_results/compare/plot", methods=["POST"])
def dispatch_plot():
    params = [int(case.split('-')[-1]) for case in request.values.keys()]
    cases = g.db.query(HypothesisSampleMatch).filter(HypothesisSampleMatch.id.in_(params)).all()
    print cases
    if len({(
            case.target_hypothesis.theoretical_structure_type,
            case.target_hypothesis.ms_level) for case in cases}) > 1:
        return jsonify(message='Invalid: Not all Hypothesis Sample Matches are the same type')
    job = analysis_comparison.HypothesisSampleMatchComparer(g.manager.path, params)
    ax = job.start()
    return svg_plot(ax, width=10, patchless=True)
