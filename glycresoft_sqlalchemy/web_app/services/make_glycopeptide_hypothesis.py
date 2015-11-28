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
    constant_modifications = values.getlist("constant_modifications")
    variable_modifications = values.getlist("variable_modifications")
    enzyme = values.get("enzyme")
    hypothesis_name = values.get("hypothesis_name")

    protein_list = request.files["protein-list-file"]
    protein_list_type = values.get("proteomics-file-type")

    site_list = request.files["glycosylation-site-list-file"]
    glycan_file = request.files["glycan-definition-file"]
    glycan_file_type = values.get("glycans-file-format")

    max_missed_cleavages = intify(values.get("missed_cleavages"))

    secure_protein_list = g.manager.get_temp_path(secure_filename(protein_list.filename))
    secure_site_list_file = g.manager.get_temp_path(secure_filename(site_list.filename))
    secure_glycan_file = g.manager.get_temp_path(secure_filename(glycan_file.filename))

    protein_list.save(secure_protein_list)
    glycan_file.save(secure_glycan_file)
    if site_list.filename != "":
        site_list.save(secure_site_list_file)
    else:
        secure_site_list_file = None

    task = None
    if protein_list_type == "mzIdentML":
        task = IntegratedOmicsGlycopeptideHypothesisBuilderTask(
            g.manager.path, hypothesis_name=hypothesis_name,
            protein_file=secure_protein_list, site_list_file=secure_site_list_file,
            glycan_file=secure_glycan_file, glycan_file_type=glycan_file_type,
            constant_modifications=constant_modifications,
            callback=lambda: 0)
    else:
        task = NaiveGlycopeptideHypothesisBuilderTask(
            g.manager.path, hypothesis_name,
            secure_protein_list, secure_site_list_file,
            secure_glycan_file, glycan_file_type, constant_modifications,
            variable_modifications, enzyme, max_missed_cleavages, callback=lambda: 0)
    g.manager.add_task(task)
    return jsonify(**dict(request.values))
