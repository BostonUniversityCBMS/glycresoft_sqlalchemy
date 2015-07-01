import argparse
import os
from glycresoft_sqlalchemy.search_space_builder import pooling_search_space_builder, pooling_make_decoys, exact_search_space_builder

commit_checkpoint = os.environ.get("GLYCRESOFT_COMMIT_INTERVAL", 2000)
import summarize


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


def main():
    args = app.parse_args()
    print args
    task_fn = args.task
    del args.task
    task_fn(**args.__dict__)

if __name__ == '__main__':
    main()
