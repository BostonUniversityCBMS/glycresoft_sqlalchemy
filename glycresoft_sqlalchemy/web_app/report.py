import os
import operator
import time
import logging
try:
    logger = logging.getLogger("web_app.report")
except:
    pass
from glycresoft_sqlalchemy.data_model import Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch, DatabaseManager
from glycresoft_sqlalchemy.report.plot_glycoforms import plot_glycoforms_svg
from jinja2 import Environment, PackageLoader, Undefined, FileSystemLoader
from jinja2 import nodes
from jinja2.ext import Extension
try:
    from cStringIO import StringIO
except:
    from io import StringIO

from matplotlib import rcParams as mpl_params
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import urllib


def ms2_score_histogram(session, hypothesis_id):
    session.query(GlycopeptideMatch).filter(
        GlycopeptideMatch.protein_id == Protein.id,
        Protein.hypothesis_id == hypothesis_id)


def q_value_below(query, threshold):
    q = query.filter(GlycopeptideMatch.q_value <= threshold)
    return q


def _string_to_deep_getter(string):
    parts = string.split(".")
    chain = []
    for p in parts:
        if p[-2:] != '()':
            link = operator.attrgetter(p)
        else:
            link = operator.methodcaller(p[:-2])
        chain.append(link)

    def unravel(obj):
        last = obj
        for link in chain:
            last = link(last)
        return last

    return unravel


def fsort(L, f):
    if isinstance(f, basestring):
        f = _string_to_deep_getter(f)
    return sorted(L, key=f, reverse=True)


def png_plot(figure, **kwargs):
    buffer = render_plot(figure, format='png', **kwargs)
    return "<img src='data:image/png;base64,%s'>" % urllib.quote(buffer.getvalue().encode("base64"))


def svg_plot(figure, **kwargs):
    buffer = render_plot(figure, format='svg', **kwargs)
    return buffer.getvalue()


def render_plot(figure, **kwargs):
    if isinstance(figure, Axes):
        figure = figure.get_figure()
    if kwargs.get("bbox_inches") != 'tight':
        figure.patch.set_visible(False)
    buffer = StringIO()
    figure.savefig(buffer, **kwargs)
    return buffer


def mass_histogram(iterable, **kwargs):
    figure = plt.figure()
    plt.hist(iterable, normed=True, **kwargs)
    ax = figure.axes[0]
    ax.set_xlabel("Mass")
    ax.set_ylabel("Frequency")
    return svg_plot(figure)


def plot_glycoforms(protein):
    old_size = mpl_params['figure.figsize']
    mpl_params['figure.figsize'] = 16, 10
    svg = plot_glycoforms_svg(protein)
    mpl_params['figure.figsize'] = old_size
    return svg


def prepare_environment(env=None):
    try:
        raise Exception()
    except:
        loader = FileSystemLoader(os.path.join(os.path.dirname(__file__), 'html'))
    if env is None:
        env = Environment(loader=loader, extensions=[FragmentCacheExtension])
    else:
        env.loader = loader
        env.add_extension(FragmentCacheExtension)
    env.fragment_cache = dict()
    env.filters["q_value_below"] = q_value_below
    env.filters["n_per_row"] = n_per_row
    env.filters['highlight_sequence_site'] = highlight_sequence_site
    env.filters['plot_glycoforms'] = plot_glycoforms_svg
    env.filters['mass_histogram'] = mass_histogram
    env.filters['svg_plot'] = svg_plot
    env.filters['png_plot'] = png_plot
    env.filters['fsort'] = fsort
    env.globals
    return env


def highlight_sequence_site(amino_acid_sequence, site_list, site_type_list):
    if isinstance(site_type_list, basestring):
        site_type_list = [site_type_list for i in site_list]
    sequence = list(amino_acid_sequence)
    for site, site_type in zip(site_list, site_type_list):
        sequence[site] = "<span class='{}'>{}</span>".format(site_type, sequence[site])
    return sequence


def n_per_row(sequence, n=60):
    row_buffer = []
    i = 0
    while i < len(sequence):
        row_buffer.append(
            ''.join(sequence[i:(i + n)])
            )
        i += n
    return '<br>'.join(row_buffer)


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
