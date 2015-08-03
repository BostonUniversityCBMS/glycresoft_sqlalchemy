import argparse
import sys

from glycresoft_sqlalchemy.data_model import *
from glycresoft_sqlalchemy.report import export_csv
from glycresoft_sqlalchemy.report import plot_glycoforms
from glycresoft_sqlalchemy.report import plot_scores
from . import summarize
from glycresoft_sqlalchemy.web_app.serve import main as webmain

task_map = {
    "export-csv": export_csv.taskmain,
    "plot-glycoforms": plot_glycoforms.taskmain,
    "plot-scores": plot_scores.taskmain,
    "summarize": summarize.taskmain,
    "web": webmain
}


class _NonExitingHelpAction(argparse._HelpAction):

    def __call__(self, parser, namespace, values, option_string=None):
        if namespace.action is None:
            parser.print_help()
            parser.exit()
        else:
            namespace.help = True

help_action = _NonExitingHelpAction(option_strings=["-h", "--help"], dest='help',
                                    help='show this help message and exit')

app = argparse.ArgumentParser("glycresoft-report", add_help=False)
app.add_argument('action', choices=task_map.keys(), help="Subcommand to execute")
app.add_argument("-h", "--help", action=_NonExitingHelpAction)


def taskmain():
    args, argv = app.parse_known_args()
    if hasattr(args, 'help'):
        argv.append("-h")
    sys.argv[1:] = argv
    task_map[args.action]()


if __name__ == '__main__':
    taskmain()
