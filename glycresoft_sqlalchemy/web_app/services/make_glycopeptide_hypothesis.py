from flask import request, g, render_template, Blueprint, jsonify
from werkzeug import secure_filename

from glycresoft_sqlalchemy.web_app.common import intify
from glycresoft_sqlalchemy.web_app.task.do_informed_glycopeptide_hypothesis import (
    IntegratedOmicsGlycopeptideHypothesisBuilderTask)
from glycresoft_sqlalchemy.web_app.task.do_naive_glycopeptide_hypothesis import (
    NaiveGlycopeptideHypothesisBuilderTask)

make_glycopeptide_hypothesis = Blueprint("make_glycopeptide_hypothesis", __name__)

app = make_glycopeptide_hypothesis


@app.route("/glycopeptide_search_space")
def build_glycopeptide_search_space():
    return render_template("glycopeptide_search_space.templ")


@app.route("/glycopeptide_search_space", methods=["POST"])
def build_glycopeptide_search_space_post():
    values = request.values
    print values
    constant_modifications = values.getlist("constant_modifications")
    variable_modifications = values.getlist("variable_modifications")
    enzyme = values.get("enzyme")
    hypothesis_name = values.get("hypothesis_name")

    protein_list = request.files["protein-list-file"]
    protein_list_type = values.get("proteomics-file-type")

    site_list = request.files["glycosylation-site-list-file"]
    glycan_file = request.files.get("glycan-definition-file")
    glycan_database = values.get("glycan-database-source")
    glycan_file_type = values.get("glycans-file-format")

    glycan_options = {}

    max_missed_cleavages = intify(values.get("missed_cleavages"))
    maximum_glycosylation_sites = intify(values.get("max_glycosylation_combinations", 1))

    secure_protein_list = g.manager.get_temp_path(secure_filename(protein_list.filename))
    secure_site_list_file = g.manager.get_temp_path(secure_filename(site_list.filename))

    protein_list.save(secure_protein_list)
    if site_list.filename != "":
        site_list.save(secure_site_list_file)
    else:
        secure_site_list_file = None

    if glycan_database == "" or glycan_database is None:
        glycan_file_type = "txt"
        glycan_options["glycomics_format"] = glycan_file_type

        secure_glycan_file = g.manager.get_temp_path(secure_filename(glycan_file.filename))
        glycan_file.save(secure_glycan_file)

        glycan_options["glycomics_path"] = secure_glycan_file
    else:
        option_type, option_id = glycan_database.split(",", 1)
        option_id = int(option_id)
        glycan_options["glycomics_path"] = g.manager.path

        if option_type == "Hypothesis":
            option_type = "hypothesis"
            glycan_options["source_hypothesis_id"] = option_id
        elif option_type == "HypothesisSampleMatch":
            option_type = "hypothesis-sample-match"
            glycan_options["source_hypothesis_sample_match_id"] = option_id

        glycan_options["glycomics_format"] = option_type

    task = None
    if protein_list_type == "mzIdentML":
        protein_names = values.get("protein_names").split(",")
        task = IntegratedOmicsGlycopeptideHypothesisBuilderTask(
            g.manager.path, hypothesis_name=hypothesis_name, protein_ids=protein_names,
            protein_file=secure_protein_list, site_list_file=secure_site_list_file,
            maximum_glycosylation_sites=maximum_glycosylation_sites,
            callback=lambda: 0, glycan_options=glycan_options)
    else:
        task = NaiveGlycopeptideHypothesisBuilderTask(
            g.manager.path, hypothesis_name=hypothesis_name,
            protein_file=secure_protein_list, site_list_file=secure_site_list_file,
            constant_modifications=constant_modifications, variable_modifications=variable_modifications,
            maximum_glycosylation_sites=maximum_glycosylation_sites,
            enzyme=enzyme, max_missed_cleavages=max_missed_cleavages, callback=lambda: 0,
            glycan_options=glycan_options)
    g.manager.add_task(task)
    return jsonify(**dict(request.values))
