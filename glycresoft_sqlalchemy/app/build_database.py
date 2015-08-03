import argparse
import os
from glycresoft_sqlalchemy.search_space_builder import (
    pooling_search_space_builder,
    pooling_make_decoys,
    exact_search_space_builder,
    naive_glycopeptide_hypothesis,
    integrated_omics)
import summarize

commit_checkpoint = os.environ.get("GLYCRESOFT_COMMIT_INTERVAL", 1000)


def build_naive_search_space_simple(glycopeptide_csv, site_list_file, digest_file, db_file_name=None, **kwargs):
    digest = pooling_search_space_builder.parse_digest(digest_file)
    kwargs.update(digest.__dict__)
    kwargs.setdefault("n_processes", 4)
    target_job = pooling_search_space_builder.PoolingTheoreticalSearchSpaceBuilder(
        glycopeptide_csv,
        db_file_name=db_file_name,
        site_list=site_list_file,
        commit_checkpoint=commit_checkpoint,
        **kwargs)
    target = target_job.start()
    print target_job.db_file_name
    print target
    decoy_job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(
        target_job.db_file_name,
        hypothesis_ids=[target],
        prefix_len=kwargs.get("prefix_len", 0),
        suffix_len=kwargs.get("suffix_len", 1),
        decoy_type=kwargs.get("decoy_type", 0),
        n_processes=kwargs.get('n_processes', 4),
        commit_checkpoint=commit_checkpoint)
    decoy, = decoy_job.start()
    print target_job.db_file_name
    summarize.main(target_job.db_file_name)
    return target, decoy


def build_informed_search_space_existing(glycopeptide_csv, hypothesis_id, db_file_name, **kwargs):
    target_job = exact_search_space_builder.ExactSearchSpaceBuilder(
        glycopeptide_csv,
        db_file_name=db_file_name,
        hypothesis_id=hypothesis_id,
        n_processes=kwargs.get('n_processes', 4),
        commit_checkpoint=commit_checkpoint)
    target = target_job.start()
    print target
    decoy_job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(
        db_file_name,
        hypothesis_ids=[hypothesis_id],
        prefix_len=kwargs.get("prefix_len", 0),
        suffix_len=kwargs.get("suffix_len", 1),
        decoy_type=kwargs.get("decoy_type", 0),
        n_processes=kwargs.get('n_processes', 4),
        commit_checkpoint=commit_checkpoint)
    decoy = decoy_job.start()
    return target, decoy


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
                 glycan_file, glycan_file_type, constant_modifications,
                 variable_modifications, enzyme, max_missed_cleavages, **kwargs)
    job.start()


def build_naive_ms2_glycopeptide(database_path, hypothesis_sample_match_id, **kwargs):
    kwargs.setdefault("n_processes", 4)
    job = pooling_search_space_builder.PoolingTheoreticalSearchSpaceBuilder.from_hypothesis(
        database_path, hypothesis_sample_match_id, **kwargs)
    target_hypothesis_id = job.start()
    decoy_job = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(
        database_path,
        hypothesis_ids=[target_hypothesis_id],
        prefix_len=kwargs.get("prefix_len", 0),
        suffix_len=kwargs.get("suffix_len", 1),
        decoy_type=kwargs.get("decoy_type", 0),
        n_processes=kwargs.get('n_processes', 4),
        commit_checkpoint=commit_checkpoint)
    decoy_hypothesis_id, = decoy_job.start()
    print target_hypothesis_id, decoy_hypothesis_id


app = argparse.ArgumentParser("build-database")
app.add_argument("-n", "--n-processes", default=4, type=int, help='Number of processes to run on')
subparsers = app.add_subparsers()

build_naive_simple = subparsers.add_parser("build-naive-simple")

build_naive_simple.add_argument("glycopeptide_csv")
build_naive_simple.add_argument(
    "-g", "--digest-file",
    help="The path to an msdigest XML file for \
    one protein with all the modification and enzyme information")
build_naive_simple.add_argument(
    "-s", "--site-list-file",
    help="The path to a file in FASTA format where the sequence body \
    is a space-separated list of positions that start at 1")
build_naive_simple.add_argument(
    "-d", "--database-path", dest="db_file_name", default=None, required=False,
    help="Path to project database. If omitted, it will be constructed from <glycopeptide_csv>")
build_naive_simple.set_defaults(task=build_naive_search_space_simple)


c = build_naive_ms1_glycopeptide_app = subparsers.add_parser("naive-glycopeptide-ms1")
c.add_argument(
    "-d", "--database-path", default=None, required=True,
    help="Path to project database.")
c.add_argument(
    "-p", "--protein-file", help="Path to Proteins to digest in Fasta format")
c.add_argument(
    "-s", "--site-list-file",
    required=False, help="Path to a file in FASTA format where each sequence is N-glycosylated."
    "Entries are each a space-separated list of positions that start at 1. If omitted, glycosites"
    " will be inferred using canonical sequon patterns.")
c.add_argument(
    "-g", "--glycan-file",
    help="Path to file of glycans to use. May be either a CSV or a database file formats."
    )
c.add_argument("-f", "--glycan-file-type", help="The format of --glycan-file")
c.add_argument(
    "-c", "--constant-modifications", action='append', help='A text formatted modification specifier'
    ' to be applied at every possible site. Modification specifier of the form "ModificationName (ResiduesAllowed)"'
    " May be specified more than once.")

c.add_argument(
    "-v", "--variable-modifications", action='append', help='A text formatted modification specifier'
    ' to be applied at every possible site. Modification specifier of the form "ModificationName (ResiduesAllowed)".'
    " May be specied more than once.")
c.add_argument(
    "-e", "--enzyme", required=False, default='trypsin', help="Protease to use for in-silico"
    " digestion of proteins. Defaults to trypsin.")
c.add_argument(
    "-m", "--missed-cleavages", type=int, default=2, required=False, help="The maximum number of"
    " missed cleavages to allow. Defaults to 2")
c.add_argument("--hypothesis-name", default=None, required=False, help="Name of the hypothesis")
c.set_defaults(task=build_naive_ms1_glycopeptide)


c = build_naive_ms2_glycopeptide_app = subparsers.add_parser("naive-glycopeptide-ms2")
c.add_argument(
    "-d", "--database-path", default=None, required=True,
    help="Path to project database.")
c.add_argument("-i", "--hypothesis-sample-match-id", help="The id number of the hypothesis sample match to build from")
c.set_defaults(task=build_naive_ms2_glycopeptide)


def main():
    args = app.parse_args()
    print args
    task_fn = args.task
    del args.task
    task_fn(**args.__dict__)

if __name__ == '__main__':
    main()
