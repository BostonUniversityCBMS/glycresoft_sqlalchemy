import os
import operator
import time
import logging
try:
    logger = logging.getLogger("web_app.report")
except:
    pass
from glycresoft_sqlalchemy.data_model import (
    GlycopeptideMatch)
from glycresoft_sqlalchemy.report.plot_glycoforms import plot_glycoforms_svg
from glycresoft_sqlalchemy.report import colors

from glycresoft_sqlalchemy.structure.sequence import Sequence, FrozenGlycanComposition
from glycresoft_sqlalchemy.report.chromatogram import draw_chromatogram
from glycresoft_sqlalchemy.report.sequence_fragment_logo import glycopeptide_match_logo

from jinja2 import Environment, PackageLoader, Undefined, FileSystemLoader, escape
from jinja2 import nodes
from jinja2.ext import Extension
try:
    from cStringIO import StringIO
except:
    from io import StringIO

import matplotlib
from matplotlib import rcParams as mpl_params
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import urllib

mpl_params.update({
    'figure.facecolor': 'white',
    'figure.edgecolor': 'white',
    'font.size': 10,
    # 72 dpi matches SVG
    # this only affects PNG export, as SVG has no dpi setting
    'savefig.dpi': 72,
    # 10pt still needs a little more room on the xlabel:
    'figure.subplot.bottom': .125})


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


def fsort(L, f, *args, **kwargs):
    if isinstance(f, basestring):
        f = _string_to_deep_getter(f)
    return sorted(L, key=lambda x: f(x, *args, **kwargs), reverse=True)


def png_plot(figure, **kwargs):
    buffer = render_plot(figure, format='png', **kwargs)
    return "<img src='data:image/png;base64,%s'>" % urllib.quote(buffer.getvalue().encode("base64"))


def svg_plot(figure, **kwargs):
    buffer = render_plot(figure, format='svg', **kwargs)
    return buffer.getvalue()


def render_plot(figure, **kwargs):
    if isinstance(figure, Axes):
        figure = figure.get_figure()
    if "height" in kwargs:
        figure.set_figheight(kwargs["height"])
    if "width" in kwargs:
        figure.set_figwidth(kwargs['width'])
    if kwargs.get("bbox_inches") != 'tight' or kwargs.get("patchless"):
        figure.patch.set_visible(False)
        figure.axes[0].patch.set_visible(False)
    buffer = StringIO()
    figure.savefig(buffer, **kwargs)
    plt.close(figure)
    return buffer


def plot_glycoforms(protein, filter_context):
    svg = plot_glycoforms_svg(protein, filterfunc=filter_context)
    return svg


def plot_chromatogram(peak_group):
    with matplotlib.rc_context({
            "figure.figsize": (16, 4),
            "axes.edgecolor": 'grey'}):
        ax = draw_chromatogram(peak_group, color='teal', alpha=0.3)
    return ax


def rgbpack(color):
    return "rgba(%d,%d,%d,0.5)" % tuple(i * 255 for i in color)


def glycopeptide_string(sequence, long=False, include_glycan=True):
    sequence = Sequence(sequence)
    parts = []
    template = "(<span class='modification-chip'"\
        " style='background-color:%s;padding-left:1px;padding-right:2px;border-radius:2px;'"\
        " title='%s' data-modification='%s'>%s</span>)"

    n_term_template = template.replace("(", "").replace(")", "") + '-'
    c_term_template = "-" + (template.replace("(", "").replace(")", ""))

    def render(mod, template=template):
        color = colors.get_color(str(mod))
        letter = escape(mod.name if long else mod.name[0])
        name = escape(mod.name)
        parts.append(template % (rgbpack(color), name, name, letter))

    if sequence.n_term != "H":
        render(sequence.n_term, n_term_template)
    for res, mods in sequence:
        parts.append(res.symbol)
        for mod in mods:
            render(mod)
    if sequence.c_term != "OH":
        render(sequence.c_term, c_term_template)
    parts.append((
        ' ' + glycan_composition_string(str(sequence.glycan)) if sequence.glycan is not None else "")
        if include_glycan else "")
    return ''.join(parts)


def formula(composition):
    return ''.join("<b>%s</b><sub>%d</sub>" % (k, v) for k, v in sorted(composition.items()))


def glycan_composition_string(composition):
    composition = FrozenGlycanComposition.parse(composition)
    parts = []
    template = "<span class='monosaccharide-name' style='background-color:%s;padding:2px;border-radius:2px;'>%s %d</span>"
    for k, v in sorted(composition.items(), key=lambda x: x[0].mass()):
        name = str(k)
        color = colors.get_color(str(name))
        parts.append(template % (rgbpack(color), name, v))
    return ' '.join(parts)


def prepare_environment(env=None):
    try:
        loader = PackageLoader("glycresoft_sqlalchemy.web_app", "html")
        loader.list_templates()
    except:
        loader = FileSystemLoader(os.path.join(os.path.dirname(__file__), 'html'))
    if env is None:
        env = Environment(loader=loader, extensions=[FragmentCacheExtension])
    else:
        env.loader = loader
    env.add_extension(FragmentCacheExtension)
    env.fragment_cache = dict()
    env.filters["n_per_row"] = n_per_row
    env.filters['highlight_sequence_site'] = highlight_sequence_site
    env.filters['plot_glycoforms'] = plot_glycoforms
    env.filters['chromatogram'] = plot_chromatogram
    env.filters['svg_plot'] = svg_plot
    env.filters['png_plot'] = png_plot
    env.filters['fsort'] = fsort
    env.filters['glycopeptide_string'] = glycopeptide_string
    env.filters['glycan_composition_string'] = glycan_composition_string
    env.filters["glycopeptide_match_logo"] = glycopeptide_match_logo
    env.filters["formula"] = formula
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
