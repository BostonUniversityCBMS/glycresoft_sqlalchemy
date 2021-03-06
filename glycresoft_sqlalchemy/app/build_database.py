import argparse
import os
from os import path

from glycresoft_sqlalchemy.proteomics import mzid_sa
from glycresoft_sqlalchemy.search_space_builder import (
    search_space_builder,
    make_decoys,
    exact_search_space_builder,
    naive_glycopeptide_hypothesis,
    integrated_omics)

from glycresoft_sqlalchemy.search_space_builder.glycan_builder import (
    composition_source, constrained_combinatorics, glycomedb_utils
    )

from glycresoft_sqlalchemy.app import let


def build_naive_ms1_glycopeptide(database_path, protein_file, site_list_file, glycan_file,
                                 glycan_file_type, constant_modifications, variable_modifications,
                                 enzyme, max_missed_cleavages=1, **kwargs):
    hypothesis_name = kwargs.pop("hypothesis_name", None)
    if hypothesis_name is None:
        hypothesis_name = "NaiveGlycopeptideHypothesis-%s@%s" % (
            os.path.basename(protein_file), os.path.basename(glycan_file))
    kwargs.setdefault("n_processes", 4)

    job = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesisBuilder(
                 database_path, hypothesis_name, protein_file, site_list_file,
                 glycomics_path=glycan_file, glycomics_format=glycan_file_type,
                 constant_modifications=constant_modifications,
                 variable_modifications=variable_modifications, enzyme=enzyme,
                 max_missed_cleavages=max_missed_cleavages, **kwargs)
    job.start()
    print job.hypothesis_id


def build_informed_ms1_glycopeptide(mzid_path, database_path, site_list_file, glycan_file,
                                    glycan_file_type="txt", protein_name_pattern='.*', hypothesis_id=None,
                                    protein_selector_type='regex', n_processes=4, **kwargs):
    if database_path is None:
        database_path = path.splitext(mzid_path)[0] + '.db'

    if protein_selector_type == "regex":
        protein_ids = list(mzid_sa.protein_names(mzid_path, protein_name_pattern))
    else:
        protein_ids = [name.lstrip().rstrip() for name in protein_name_pattern.split(",")]

    job = integrated_omics.IntegratedOmicsMS1SearchSpaceBuilder(
        database_path, mzid_path=mzid_path, protein_ids=protein_ids,
        glycomics_path=glycan_file, glycomics_format=glycan_file_type,
        n_processes=n_processes, **kwargs)

    job.start()
    print job.hypothesis_id


def build_naive_ms2_glycopeptide(database_path, hypothesis_sample_match_id, **kwargs):
    kwargs.setdefault("n_processes", 4)
    job = search_space_builder.BatchingTheoreticalSearchSpaceBuilder.from_hypothesis(
        database_path, hypothesis_sample_match_id, **kwargs)
    target_hypothesis_id = job.start()
    decoy_job = make_decoys.BatchingDecoySearchSpaceBuilder(
        database_path,
        hypothesis_ids=[target_hypothesis_id],
        prefix_len=kwargs.get("prefix_len", 0),
        suffix_len=kwargs.get("suffix_len", 1),
        decoy_type=kwargs.get("decoy_type", 0),
        n_processes=kwargs.get('n_processes', 4))
    decoy_hypothesis_id, = decoy_job.start()
    print target_hypothesis_id, decoy_hypothesis_id


def build_informed_ms2_glycopeptide(database_path, hypothesis_sample_match_id, **kwargs):
    kwargs.setdefault("n_processes", 4)
    job = exact_search_space_builder.ExactSearchSpaceBuilder.from_hypothesis(
        database_path, hypothesis_sample_match_id, **kwargs)
    target_hypothesis_id = job.start()
    decoy_job = make_decoys.BatchingDecoySearchSpaceBuilder(
        database_path,
        hypothesis_ids=[target_hypothesis_id],
        prefix_len=kwargs.get("prefix_len", 0),
        suffix_len=kwargs.get("suffix_len", 1),
        decoy_type=kwargs.get("decoy_type", 0),
        n_processes=kwargs.get('n_processes', 4))
    decoy_hypothesis_id, = decoy_job.start()
    print target_hypothesis_id, decoy_hypothesis_id


def build_glycan_glycome_db(
        database_path, glycomedb_path=None,
        taxonomy_path=None, taxon_id=None, include_descendent_taxa=False, include_structures=True,
        motif_family=None, reduction=None, derivatization=None, *args, **kwargs):
    print args, kwargs
    job = glycomedb_utils.GlycomeDBHypothesis(
        database_path, hypothesis_id=None, glycomedb_path=glycomedb_path,
        taxonomy_path=taxonomy_path, taxa_ids=taxon_id, include_descendent_taxa=include_descendent_taxa,
        include_structures=include_structures,
        motif_family=motif_family, reduction=reduction, derivatization=derivatization)
    job.start()
    print job.hypothesis_id


def build_glycan_text(database_path, text_file_path, reduction=None, derivatization=None, *args, **kwargs):
    print args, kwargs
    job = composition_source.TextGlycanCompositionHypothesisBuilder(
        database_path=database_path, text_file_path=text_file_path, reduction=reduction,
        derivatization=derivatization, hypothesis_id=None)
    job.start()
    print job.hypothesis_id


def build_glycan_other_hypothesis(database_path, source_hypothesis_id=None, reduction=None,
                                  derivatization=None, *args, **kwargs):
    print args, kwargs
    job = composition_source.OtherGlycanHypothesisGlycanHypothesisBuilder(
        database_path, source_hypothesis_id=source_hypothesis_id, reduction=reduction,
        derivatization=derivatization, *args, **kwargs)
    job.start()
    print job.hypothesis_id


def build_glycan_rules(database_path, rules_file_path, reduction=None, derivatization=None, *args, **kwargs):
    print args, kwargs
    job = constrained_combinatorics.ConstrainedCombinatoricsGlycanHypothesisBuilder(
        database_path, rules_file_path, reduction=reduction, derivatization=derivatization)
    job.start()
    print job.hypothesis_id


def dispatch_glycomics(database_path, glycan_file, glycan_file_type, derivatization, reduction, *args, **kwargs):
    if glycan_file_type == "text":
        build_glycan_text(
            database_path, glycan_file, derivatization=derivatization,
            reduction=reduction, *args, **kwargs)
    elif glycan_file_type == "rules":
        build_glycan_rules(
            database_path, glycan_file, derivatization=derivatization,
            reduction=reduction, *args, **kwargs)
    else:
        raise Exception("Glycan File Type Not Supported")

app = argparse.ArgumentParser("build-database")
app.add_argument("-n", "--n-processes", default=4, type=int, help='Number of processes to run on')
subparsers = app.add_subparsers()


def configure_parser_glycomics_for_glycoproteomics(parser):
    c = parser
    c.add_argument(
        "-g", "--glycan-file",
        help="Path to file of glycans to use. May be either a CSV or a database file formats."
        )
    c.add_argument("-f", "--glycan-file-type", help="The format of --glycan-file")
    c.add_argument("--source-hypothesis-id", required=False, type=int,
                   help="The hypothesis from which to copy all glycan compositions")
    c.add_argument("--source-hypothesis-sample-match-id", required=False, type=int,
                   help="The hypothesis sample match from which to copy all glycan compositions")


def configure_parser_common_glycopeptide_properties(parser):
    c = parser
    c.add_argument("--maximum-glycosylation-sites", type=int, required=False,
                   default=1, help="The maximum number of glycosylation sites to"
                   " populate on a single peptide. This number will exponentially increase"
                   " runtime if greater than 1.")
    c.add_argument(
        "-s", "--site-list-file",
        required=False, help="Path to a file in FASTA format where each sequence is N-glycosylated."
        "Entries are each a space-separated list of positions that start at 1. If omitted, glycosites"
        " will be inferred using canonical sequon patterns.")

    c.add_argument("-n", "--n-processes", default=4, required=False, type=int)


with let(subparsers.add_parser("naive-glycopeptide-ms1")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument(
        "-p", "--protein-file", help="Path to Proteins to digest in Fasta format")

    configure_parser_glycomics_for_glycoproteomics(c)
    configure_parser_common_glycopeptide_properties(c)

    c.add_argument(
        "-c", "--constant-modifications", action='append', help='A text formatted modification specifier'
        ' to be applied at every possible site. Modification specifier of the form \
        "ModificationName (ResiduesAllowed)"'
        " May be specified more than once.")

    c.add_argument(
        "-v", "--variable-modifications", action='append', help='A text formatted modification specifier'
        ' to be applied at every possible site. Modification specifier of the form '
        '"ModificationName (ResiduesAllowed). May be specied more than once.')
    c.add_argument(
        "-e", "--enzyme", required=False, default='trypsin', help="Protease to use for in-silico"
        " digestion of proteins. Defaults to trypsin.")
    c.add_argument(
        "-m", "--missed-cleavages", dest="max_missed_cleavages", type=int, default=2, required=False,
        help="The maximum number of missed cleavages to allow. Defaults to 2")
    c.add_argument("--hypothesis-name", default=None, required=False, help="Name of the hypothesis")
    c.set_defaults(task=build_naive_ms1_glycopeptide)

with let(subparsers.add_parser("informed-glycopeptide-ms1")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument("--hypothesis-name", default=None, required=False, help="Name of the hypothesis")

    configure_parser_glycomics_for_glycoproteomics(c)
    configure_parser_common_glycopeptide_properties(c)

    c.add_argument("-m", "--mzid-path", required=True, help="Path to mzidentML file to read proteomics data from")
    c.add_argument("-b", "--include-baseline-peptides", required=False, action="store_true", default=False,
                   help='Include all peptides spanning glycosylation sites with'
                        ' on constant modifications and no missed cleavages, even if'
                        ' they were not observed in the proteomics data.')
    c.add_argument("-p", '--protein-selection', dest="protein_name_pattern", default='.*', required=False,
                   help='Select the proteins to include ')
    c.add_argument("--protein-selector-type", choices=['regex', 'list'], default='regex', required=False,
                   help="Choose how to select proteins from the mzidentML file, using either a name regex or a list"
                        " of names separated by commas")
    c.set_defaults(task=build_informed_ms1_glycopeptide)


with let(subparsers.add_parser("naive-glycopeptide-ms2")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument("-i", "--hypothesis-sample-match-id", help="The id number of the hypothesis sample match to build from")
    c.set_defaults(task=build_naive_ms2_glycopeptide)

with let(subparsers.add_parser("informed-glycopeptide-ms2")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument("-i", "--hypothesis-sample-match-id", help="The id number of the hypothesis sample match to build from")
    c.set_defaults(task=build_informed_ms2_glycopeptide)

with let(subparsers.add_parser("glycomics")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument(
        "-g", "--glycan-file",
        help="Path to file of glycans to use. May be either a text file containing glycan compositions "
             "(In trivial IUPAC) or glycan structures (in condensed GlycoCT). Alternatively, it may be a "
             "list of constraint rules with which to generate compositions."
        )
    c.add_argument("-f", "--glycan-file-type", choices=('text', 'glycoct', 'rules'),
                   help="The format of --glycan-file")
    c.add_argument("-e", "--derivatization")
    c.add_argument("-r", "--reduction")
    c.set_defaults(task=dispatch_glycomics)

with let(subparsers.add_parser("glycomics-glycomedb")) as c:
    c.add_argument(
        "-d", "--database-path", default=None, required=True,
        help="Path to project database.")
    c.add_argument(
        "-g", "--glycomedb-path", default=None, required=False)
    c.add_argument("-t", "--taxon-id", default=[], type=int, action='append')
    c.add_argument("-s", "--include-structures", required=False, action="store_true")
    c.add_argument("-m", "--motif-family", choices=("N-Linked", "O-Linked"), required=True)
    c.add_argument("-x", "--taxonomy-path")
    c.add_argument("-i", "--include-descendent-taxa", action="store_true")
    c.add_argument("-e", "--derivatization")
    c.add_argument("-r", "--reduction")

    c.set_defaults(task=build_glycan_glycome_db)


def main():
    args = app.parse_args()
    task_fn = args.task
    del args.task
    task_fn(**args.__dict__)

if __name__ == '__main__':
    main()
