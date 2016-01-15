from collections import OrderedDict
from itertools import cycle
from glypy.composition.glycan_composition import FrozenGlycanComposition, FrozenMonosaccharideResidue
from matplotlib.colors import cnames, hex2color
from matplotlib import patches as mpatches


from glycresoft_sqlalchemy.utils import simple_repr


def lighten(rgb, factor=0.25):
    '''Given a triplet of rgb values, lighten the color by `factor`%'''
    factor += 1
    return [min(c * factor, 1) for c in rgb]


def darken(rgb, factor=0.25):
    '''Given a triplet of rgb values, darken the color by `factor`%'''
    factor = 1 - factor
    return [(c * factor) for c in rgb]


colors = cycle([hex2color(cnames[name]) for name in (
    "red", "blue", "yellow", "purple", "navy", "grey", "coral", "forestgreen", "limegreen", "maroon", "aqua",
    "lavender", "lightcoral", "mediumorchid")])


class ColorMapper(object):
    colors = [hex2color(cnames[name]) for name in (
              "red", "blue", "yellow", "purple", "navy", "grey", "coral", "forestgreen", "limegreen", "maroon", "aqua",
              "lavender", "lightcoral", "mediumorchid")]

    def __init__(self):
        self.color_name_map = {
            "HexNAc": hex2color(cnames["mediumseagreen"]),
        }
        self.color_generator = cycle(self.colors)

    def get_color(self, name):
        """Given a name, find the color mapped to that name, or
        select the next color from the `colors` generator and assign
        it to the name and return the new color.

        Parameters
        ----------
        name : object
            Any hashable object, usually a string

        Returns
        -------
        tuple: RGB triplet
        """
        try:
            return self.color_name_map[name]
        except KeyError:
            o = self.color_name_map[name] = self.color_generator.next()
            return o

    __getitem__ = get_color

    def __setitem__(self, name, color):
        self.color_name_map[name] = color

    def __repr__(self):
        return repr(self.color_name_map)

    darken = staticmethod(darken)
    lighten = staticmethod(lighten)

    def keys(self):
        return self.color_name_map.keys()

    def items(self):
        return self.color_name_map.items()


_color_mapper = ColorMapper()

color_name_map = _color_mapper.color_name_map
get_color = _color_mapper.get_color


def color_dict():
    return {str(k): v for k, v in _color_mapper.items()}


def _degree_monosaccharide_alteration(x):
    try:
        if not isinstance(x, FrozenMonosaccharideResidue):
            x = FrozenMonosaccharideResidue.from_iupac_lite(x)
        return (len(x.modifications), len(x.substituent_links))
    except:
        print x
        return 0, 0


class GlycanCompositionOrderer(object):
    def __init__(self, priority_residues=None, sorter=None):
        self.priority_residues = priority_residues or []
        self.sorter = sorter or _degree_monosaccharide_alteration

        self.priority_residues = [
            residue if isinstance(
                residue, FrozenMonosaccharideResidue) else FrozenMonosaccharideResidue.from_iupac_lite(
                residue)
            for residue in self.priority_residues
        ]

    def key_order(self, keys):
        keys = list(keys)
        for r in reversed(self.priority_residues):
            try:
                i = keys.index(r)
                keys.pop(i)
                keys = [r] + keys
            except ValueError:
                pass
        return keys

    def __call__(self, a, b):
        if isinstance(a, basestring):
            a = FrozenGlycanComposition.parse(a)
            b = FrozenGlycanComposition.parse(b)
        keys = self.key_order(sorted(set(a) | set(b), key=self.sorter))

        for key in keys:
            if a[key] < b[key]:
                return -1
            elif a[key] > b[key]:
                return 1
            else:
                continue
        return 0

    __repr__ = simple_repr


class CompositionRangeRule(object):
    def __init__(self, name, low=None, high=None, required=True):
        self.name = name
        self.low = low
        self.high = high
        self.required = required

    __repr__ = simple_repr

    def __call__(self, obj):
        try:
            composition = obj.glycan_composition
        except:
            composition = FrozenGlycanComposition.parse(obj)
        if self.name in composition:
            if self.low is None:
                return composition[self.name] <= self.high
            elif self.high is None:
                return self.low <= composition[self.name]
            return self.low <= composition[self.name] <= self.high
        else:
            return not self.required


class CompositionRuleClassifier(object):
    def __init__(self, name, rules):
        self.name = name
        self.rules = rules

    def __iter__(self):
        return iter(self.rules)

    def __call__(self, obj):
        for rule in self:
            if not rule(obj):
                return False
        return True

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.name)

    __repr__ = simple_repr


class allset(frozenset):
    def __contains__(self, k):
        return True


class GlycanCompositionClassifierColorizer(object):
    def __init__(self, rule_color_map=None, default=None):
        self.rule_color_map = rule_color_map or {}
        self.default = default

    def __call__(self, obj):
        for rule, color in self.rule_color_map.items():
            if rule(obj):
                return color
        if self.default:
            return self.default
        raise ValueError("Could not classify %r" % obj)

    def classify(self, obj):
        for rule, color in self.rule_color_map.items():
            if rule(obj):
                return rule.name
        return None

    __repr__ = simple_repr

    def make_legend(self, included=allset(), alpha=0.5):
        return [
            mpatches.Patch(
                label=rule.name, color=color, alpha=0.5) for rule, color in self.rule_color_map.items()
            if rule.name in included
            ]


NGlycanCompositionColorizer = GlycanCompositionClassifierColorizer(OrderedDict([
    (CompositionRuleClassifier("High Mannose", [CompositionRangeRule("HexNAc", 2, 2)]), '#1f77b4'),
    (CompositionRuleClassifier("Hybrid", [CompositionRangeRule("HexNAc", 3, 3)]), '#ff7f0e'),
    (CompositionRuleClassifier("Bi-Antennerary", [CompositionRangeRule("HexNAc", 4, 4)]), '#2ca02c'),
    (CompositionRuleClassifier("Tri-Antennerary", [CompositionRangeRule("HexNAc", 5, 5)]), '#d62728'),
    (CompositionRuleClassifier("Tetra-Antennerary", [CompositionRangeRule("HexNAc", 6, 6)]), '#9467bd'),
    (CompositionRuleClassifier("Penta-Antennerary", [CompositionRangeRule("HexNAc", 7, 7)]), '#8c564b'),
    (CompositionRuleClassifier("Supra-Penta-Antennerary", [CompositionRangeRule("HexNAc", 8)]), 'brown')
]))
NGlycanCompositionOrderer = GlycanCompositionOrderer(["HexNAc", "Hex", "Fucose", "NeuAc"])

_null_color_chooser = GlycanCompositionClassifierColorizer({}, default='blue')


class GlycanLabelTransformer(object):
    def __init__(self, label_series, order_chooser):
        self._input_series = label_series
        self.order_chooser = order_chooser
        self.residues = None
        self._infer_compositions()

    def _infer_compositions(self):
        residues = set()
        for item in self._input_series:
            if isinstance(item, basestring):
                item = FrozenGlycanComposition.parse(item)
            residues.update(item)

        residues = sorted(residues, key=lambda x: (
            len(x.modifications), len(x.substituent_links)))
        self.residues = self.order_chooser.key_order(residues)

    def transform(self):
        for item in self._input_series:
            if isinstance(item, basestring):
                item = FrozenGlycanComposition.parse(item)
            counts = [str(item[r]) for r in self.residues]
            yield '[%s]' % '; '.join(counts)

    __iter__ = transform

    __repr__ = simple_repr

    @property
    def label_key(self):
        return "Key: [%s]" % '; '.join(map(str, self.residues))
