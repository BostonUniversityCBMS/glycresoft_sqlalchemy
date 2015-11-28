from flask import request, g, render_template, Blueprint, jsonify
from werkzeug import secure_filename

from glycresoft_sqlalchemy.web_app.common import intify, random_string, remove_empty_rows
from glycresoft_sqlalchemy.search_space_builder.glycan_builder.constrained_combinatorics import CombinatoricCompositionGenerator
from glycresoft_sqlalchemy.web_app.task import do_combinatorial_glycan_hypothesis

app = make_glycan_hypothesis = Blueprint("make_glycan_hypothesis", __name__)


@app.route("/glycan_search_space")
def build_naive_glycan_search():
    return render_template("glycan_search_space.templ")


@app.route("/glycan_search_space", methods=["POST"])
def build_naive_glycan_search_process():
    data = request.values
    custom_reduction_type = data.get("custom_reduction_type")
    custom_derivatization_type = data.get("custom_derivatization_type")

    has_custom_reduction = custom_reduction_type != ""
    has_custom_derivatization = custom_derivatization_type != ""

    reduction_type = data.get("reduction_type")
    derivatization_type = data.get("derivatization_type")

    comb_monosaccharide_name = data.getlist('monosaccharide_name')[:-1]
    comb_compositions = data.getlist('monosaccharide_composition')[:-1]
    comb_lower_bound = map(intify, data.getlist('monosaccharide_lower_bound')[:-1])
    comb_upper_bound = map(intify, data.getlist('monosaccharide_upper_bound')[:-1])

    comb_monosaccharide_name, comb_lower_bound, comb_upper_bound = remove_empty_rows(
        comb_monosaccharide_name, comb_lower_bound, comb_upper_bound)

    include_human_n_glycan = data.get("glycomedb-human-n-glycan")
    include_human_o_glycan = data.get("glycomedb-human-o-glycan")
    include_mammalian_n_glycan = data.get("glycomedb-mammlian-n-glycan")
    include_mammalian_o_glycan = data.get("glycomedb-mammlian-o-glycan")

    constraint_lhs = data.getlist("left_hand_side")[:-1]
    constraint_op = data.getlist("operator")[:-1]
    constraint_rhs = data.getlist("right_hand_side")[:-1]

    constraints = zip(*remove_empty_rows(constraint_lhs, constraint_op, constraint_rhs))

    if any([include_human_n_glycan, include_human_o_glycan,
            include_mammalian_n_glycan, include_mammalian_o_glycan]):
        pass
        # Extract Pre-canned
    else:
        rules_table = CombinatoricCompositionGenerator.build_rules_table(
            comb_monosaccharide_name, comb_lower_bound, comb_upper_bound)
        task = do_combinatorial_glycan_hypothesis.CombinatorialGlycanHypothesisTask(
            g.manager.path, rules_table, constraints, data.get("hypothesis_name", random_string(10)),
            custom_reduction_type if has_custom_reduction else reduction_type,
            custom_derivatization_type if has_custom_derivatization else derivatization_type,
            callback=lambda: 0)
        g.manager.add_task(task)
        # Create Combinatorial
    # else read from text file

    return jsonify(**dict(request.values))
