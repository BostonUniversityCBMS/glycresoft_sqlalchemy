import argparse
import logging
try:
    logger = logging.getLogger("run_search")
except:
    pass
from glycresoft_sqlalchemy.matching import matching, peak_grouping
from glycresoft_sqlalchemy.matching.glycopeptide.pipeline import GlycopeptideFragmentMatchingPipeline

from glycresoft_sqlalchemy.spectra.bupid_topdown_deconvoluter_sa import BUPIDMSMSYamlParser, process_data_file
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms2.search_space_builder import (
    TheoreticalSearchSpaceBuilder)
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms2.make_decoys import (
    BatchingDecoySearchSpaceBuilder)
from glycresoft_sqlalchemy.spectra.decon2ls_sa import Decon2LSIsosParser
from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, MS2GlycopeptideHypothesisSampleMatch, SampleRun, Hypothesis,
    MassShift, HypothesisSampleMatch, TheoreticalGlycanComposition, TheoreticalGlycopeptideComposition)
from glycresoft_sqlalchemy.utils.database_utils import get_or_create

from glycresoft_sqlalchemy.app import let, fail


ms1_tolerance_default = matching.ms1_tolerance_default
ms2_tolerance_default = matching.ms2_tolerance_default


class ParseMassShiftAction(argparse.Action):
    def __init__(self, option_strings, dest, default=None, **kwargs):
        kwargs['default'] = []
        kwargs['nargs'] = 3
        kwargs["metavar"] = ("NAME", "MASSDELTA", "MAX")
        super(ParseMassShiftAction, self).__init__(option_strings, dest, **kwargs)

    def parse(self, shift):
        return shift[0].replace("\-", '-'), float(shift[1].replace("\-", '-')), int(shift[2])

    def __call__(self, parser, namespace, values, option_string=None):
        getattr(namespace, self.dest).append(self.parse(values))


def index_isos(
        file_path, database_path, grouping_tolerance=8e-5,
        minimum_scan_count=1, max_charge_state=8,
        minimum_abundance_ratio=0.01, minimum_mass=1200., maximum_mass=15000.,
        begin_scan=0, end_scan=float('inf'),
        n_processes=4):
    parser = Decon2LSIsosParser(file_path, database_path)
    database_path = parser.manager.path
    sample_run_id = parser.sample_run.id
    grouper = peak_grouping.Decon2LSPeakGrouper(
        database_path, sample_run_id=sample_run_id, max_charge_state=max_charge_state,
        minimum_abundance_ratio=minimum_abundance_ratio, grouping_error_tolerance=grouping_tolerance,
        minimum_mass=minimum_mass, maximum_mass=maximum_mass, minimum_scan_id=begin_scan,
        maximum_scan_id=end_scan, n_processes=n_processes)
    grouper.start()


def index_bupid(file_path, database_path=None):
    process_data_file(file_path, database_path)


def run_ms2_glycoproteomics_search(
        database_path, observed_ions_path, target_hypothesis_id=None, source_hypothesis_sample_match_id=None,
        decoy_hypothesis_id=None, observed_ions_type='bupid_yaml', sample_run_id=None,
        ms1_tolerance=ms1_tolerance_default,
        ms2_tolerance=ms2_tolerance_default, **kwargs):
    manager = DatabaseManager(database_path)
    manager.initialize()

    n_processes = kwargs.get("n_processes", 4)

    if source_hypothesis_sample_match_id is not None:
        builder = TheoreticalSearchSpaceBuilder.from_hypothesis_sample_match(
            database_path, source_hypothesis_sample_match_id, n_processes=n_processes)
        target_hypothesis_id = builder.start()
        decoy_builder = BatchingDecoySearchSpaceBuilder(
            database_path, hypothesis_ids=[target_hypothesis_id], n_processes=n_processes)
        decoy_hypothesis_id = decoy_builder.start()
        decoy_hypothesis_id = decoy_hypothesis_id[0]
    elif target_hypothesis_id is None:
        fail("A Hypothesis must be provided if an MS1 HypothesisSampleMatch is not specified")
    if observed_ions_type == 'bupid_yaml' and observed_ions_path[-3:] != '.db':
        parser = BUPIDMSMSYamlParser(observed_ions_path, manager.bridge_address())
        observed_ions_path = parser.manager.path
        observed_ions_type = 'db'
    #     sample_name = parser.sample_run_name
    # else:
    #     sample_name = ','.join(x[0] for x in DatabaseManager(
    #     observed_ions_path).session().query(SampleRun.name).all())

    job = GlycopeptideFragmentMatchingPipeline(
        database_path, observed_ions_path, target_hypothesis_id=target_hypothesis_id,
        decoy_hypothesis_id=decoy_hypothesis_id, ms1_tolerance=ms1_tolerance,
        ms2_tolerance=ms2_tolerance, n_processes=n_processes)
    job.start()


def run_ms1_search(
        database_path, observed_ions_path, hypothesis_id=None,
        observed_ions_type='isos', sample_run_id=None,
        grouping_tolerance=8e-5, search_type=None,
        match_tolerance=1e-5, mass_shift=None, n_processes=4,
        minimum_mass=1200, maximum_mass=15000,
        **kwargs):
    search_type = {
        "glycopeptide": "TheoreticalGlycopeptideComposition",
        "glycan": "TheoreticalGlycanComposition"
    }[search_type]

    manager = DatabaseManager(database_path)
    session = manager.session()

    hypothesis = session.query(Hypothesis).get(hypothesis_id)
    guess_search_type = hypothesis.theoretical_structure_type

    if issubclass(guess_search_type, TheoreticalGlycanComposition):
        search_type = "TheoreticalGlycanComposition"
    elif issubclass(guess_search_type, TheoreticalGlycopeptideComposition):
        search_type = "TheoreticalGlycopeptideComposition"
    else:
        raise Exception("Could not infer search type for %s" % guess_search_type)

    mass_shift_map = {}
    for shift_params in mass_shift:
        shift, flag = get_or_create(session, MassShift, name=shift_params[0], mass=shift_params[1])
        session.add(shift)
        mass_shift_map[shift] = shift_params[2]

    session.commit()

    if observed_ions_type == 'isos':
        parser = Decon2LSIsosParser(observed_ions_path, manager.bridge_address())
        observed_ions_path = parser.manager.path

    pipeline = peak_grouping.LCMSPeakClusterSearch(
        database_path, observed_ions_path, hypothesis_id, sample_run_id=1,
        grouping_error_tolerance=grouping_tolerance,
        search_type=search_type, match_tolerance=match_tolerance,
        mass_shift_map=mass_shift_map, minimum_mass=minimum_mass,
        maximum_mass=maximum_mass,
        n_processes=n_processes, **kwargs)
    pipeline.start()


app = argparse.ArgumentParser('database-search')

subparsers = app.add_subparsers()

ms2_glycoproteomics_app = subparsers.add_parser("ms2-glycoproteomics")
with let(ms2_glycoproteomics_app) as c:
    c.add_argument("database_path")

    target_group = c.add_mutually_exclusive_group()
    target_group.add_argument("-a", "--target-hypothesis-id", type=int,
                              help='The identity of the target hypothesis, if it already exists')
    target_group.add_argument(
        "-s", "--source-hypothesis-sample-match-id", type=int,
        help="The identity of the hypothesis sample match from which to build the MS2 Hypothesis")

    c.add_argument("-n", "--n-processes", default=4, required=False, type=int)
    c.add_argument("-i", "--observed-ions-path")
    c.add_argument("-p", "--observed-ions-type", default='bupid_yaml', choices=["bupid_yaml", "db"])
    c.add_argument("-d", "--decoy-hypothesis-id", type=int, default=None, required=False)
    c.add_argument(
        "-t1", "--ms1-tolerance", default=ms1_tolerance_default, required=False, type=float)
    c.add_argument(
        "-t2", "--ms2-tolerance", default=ms2_tolerance_default, required=False, type=float)
    c.set_defaults(task=run_ms2_glycoproteomics_search)


ms1_app = subparsers.add_parser("ms1")
with let(ms1_app) as c:
    c.add_argument("database_path")
    c.add_argument("hypothesis_id")
    c.add_argument("-n", "--n-processes", default=4, required=False, type=int)
    c.add_argument("-i", "--observed-ions-path")
    c.add_argument("-p", "--observed-ions-type", default='isos', choices=["isos", "db"])
    c.add_argument("-t", "--match-tolerance", default=1e-5, required=False, type=float)
    c.add_argument("-g", "--grouping-tolerance", default=2e-5, required=False, type=float)
    c.add_argument("-s", "--search-type", default='glycopeptide', choices=['glycan', 'glycopeptide'])
    c.add_argument("-m", "--mass-shift", action=ParseMassShiftAction)
    c.add_argument('--skip-grouping', action='store_true', required=False)
    c.add_argument('--skip-matching', action='store_true', required=False)
    c.add_argument('--hypothesis-sample-match-id', action='store', default=None, required=False)
    c.add_argument('-l', '--minimum-mass', type=float, required=False, default=None)
    c.add_argument("-u", "--maximum-mass", type=float, required=False, default=None)
    c.set_defaults(task=run_ms1_search)


index_isos_app = subparsers.add_parser("index-isos")
with let(index_isos_app) as c:
    c.add_argument("-n", "--n-processes", default=4, required=False, type=int)
    c.add_argument("file_path", help='path to isos.csv file')
    c.add_argument("-d", "--database-path", required=False, default=None, help='path to output database.'
                   ' Defaults to the same directory as the isos file')
    c.add_argument("-c", "--max-charge-state", type=int, required=False, default=8,
                   help="Only consider peaks with a charge state <= this.")
    c.add_argument("-l", "--minimum-mass", type=float, required=False, default=1200,
                   help="Only consider peaks with a neutral mass > this. Default"
                        " = 1200 Da, reasonable for glycopeptides. For glycans, 500 Da is more appropriate"
                        " depending upon derivatization.")
    c.add_argument("-u", "--maximum-mass", type=float, required=False, default=15000,
                   help="Only consider peaks with a neutral mass < this. Default"
                        " = 15000 Da, reasonable for glycopeptides and glycans.")
    c.add_argument("-g", "--grouping-tolerance", default=2e-5, required=False, type=float)
    c.add_argument("-b", "--begin-scan", type=int, default=0, required=False, help="The scan to start grouping from")
    c.add_argument("-e", "--end-scan", type=int, default=float('inf'), required=False,
                   help="The scan to stop grouping at")
    c.set_defaults(task=index_isos)


index_bupid_app = subparsers.add_parser("index-bupid")
with let(index_bupid_app) as c:
    c.add_argument("file_path", help="path to yaml file")
    c.add_argument("-d", "--database-path", required=False, default=None, help="path to output database."
                   " Defaults to the same directory as the yaml file.")
    c.set_defaults(task=index_bupid)


def main():
    args = app.parse_args()
    logger.debug("Arguments %r", args)
    task = args.task
    del args.task
    task(**args.__dict__)

if __name__ == '__main__':
    main()
