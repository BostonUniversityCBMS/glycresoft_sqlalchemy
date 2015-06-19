import os
import time

from ..data_model import Experiment, Protein, TheoreticalGlycopeptide, GlycopeptideMatch, DatabaseManager

from jinja2 import Environment, PackageLoader, Undefined, FileSystemLoader
from jinja2 import nodes
from jinja2.ext import Extension


def ms2_score_histogram(session, experiment_id):
    session.query(GlycopeptideMatch).filter(
        GlycopeptideMatch.protein_id == Protein.id,
        Protein.experiment_id == experiment_id)


def group_glycoforms(query):
    q = query.order_by(
            GlycopeptideMatch.q_value.asc()).order_by(
            GlycopeptideMatch.base_peptide_sequence,
            GlycopeptideMatch.glycan_composition_str)
    return q


def q_value_below(query, threshold):
    q = query.filter(GlycopeptideMatch.q_value <= threshold)
    return q


def prepare_environment(env=None):
    try:
        # loader = PackageLoader("glypy", "search/results_template")
        # loader.list_templates()
        raise Exception()
    except:
        loader = FileSystemLoader(os.path.join(os.path.dirname(__file__), 'template'))
    if env is None:
        env = Environment(loader=loader, extensions=[FragmentCacheExtension])
    else:
        env.loader = loader
        env.add_extension(FragmentCacheExtension)
    env.fragment_cache = dict()
    env.filters["group_glycoforms"] = group_glycoforms
    env.filters["q_value_below"] = q_value_below

    return env


def render(database_path, experiment_id):
    env = prepare_environment()
    template = env.get_template("main.tmpl")
    session = DatabaseManager(database_path).session()
    experiment = session.query(Experiment).get(experiment_id)
    return template.render(experiment=experiment)


class FragmentCacheExtension(Extension):
    # a set of names that trigger the extension.
    tags = set(['cache'])

    def __init__(self, environment):
        super(FragmentCacheExtension, self).__init__(environment)

        # add the defaults to the environment
        environment.extend(
            fragment_cache_prefix='',
            fragment_cache=None
        )

    def parse(self, parser):
        # the first token is the token that started the tag.  In our case
        # we only listen to ``'cache'`` so this will be a name token with
        # `cache` as value.  We get the line number so that we can give
        # that line number to the nodes we create by hand.
        lineno = parser.stream.next().lineno

        # now we parse a single expression that is used as cache key.
        args = [parser.parse_expression()]

        # if there is a comma, the user provided a timeout.  If not use
        # None as second parameter.
        if parser.stream.skip_if('comma'):
            args.append(parser.parse_expression())
        else:
            args.append(nodes.Const(None))

        # now we parse the body of the cache block up to `endcache` and
        # drop the needle (which would always be `endcache` in that case)
        body = parser.parse_statements(['name:endcache'], drop_needle=True)

        # now return a `CallBlock` node that calls our _cache_support
        # helper method on this extension.
        return nodes.CallBlock(self.call_method('_cache_support', args),
                               [], [], body).set_lineno(lineno)

    def _cache_support(self, name, timeout, caller):
        """Helper callback."""
        key = self.environment.fragment_cache_prefix + name
        print("{} - In cache for {}".format(time.time(), name))
        print(self.environment.fragment_cache.keys())
        # try to load the block from the cache
        # if there is no fragment in the cache, render it and store
        # it in the cache.
        rv = self.environment.fragment_cache.get(key)
        if rv is not None:
            print "Cache Hit"
            return rv
        print "Cache Miss"
        rv = caller()
        self.environment.fragment_cache[key] = rv
        return rv
