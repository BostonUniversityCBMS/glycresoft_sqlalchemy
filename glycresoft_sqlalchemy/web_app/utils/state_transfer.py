

class MonosaccharideFilterSet(object):
    def __init__(self, constraints=None):
        if constraints is None:
            constraints = tuple()
        self.constraints = constraints

    def __eq__(self, other):
        return self.constraints == other.constraints

    def __hash__(self):
        return hash(self.constraints)

    def items(self):
        for constraint in self.constraints:
            yield constraint.monosaccharide, constraint

    @classmethod
    def fromdict(cls, filters):
        filters = [MonosaccharideFilter(monosaccharide=monosaccharide, **constraint)
                   for monosaccharide, constraint in filters.items()]
        filters = tuple(filters)
        return cls(filters)

    def __repr__(self):
        return "MonosaccharideFilterSet(%r)" % (self.constraints,)


class MonosaccharideFilter(object):
    def __init__(self, monosaccharide, minimum, maximum, include):
        self.monosaccharide = monosaccharide
        self.minimum = minimum
        self.maximum = maximum
        self.include = include

    def __eq__(self, other):
        return self.monosaccharide == other.monosaccharide and\
            self.minimum == other.minimum and\
            self.maximum == other.maximum and\
            self.include == other.include

    def __hash__(self):
        return hash((self.monosaccharide, self.minimum, self.maximum,
                     self.include))

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return "MonosaccharideFilter(monosaccharide=%r, minimum=%d, maximum=%d, include=%r)" % (
            self.monosaccharide, self.minimum, self.maximum, self.include)


def typify(s):
    if isinstance(s, list):
        return typify_list(s)
    elif isinstance(s, dict):
        return typify_dict(s)
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return str(s)


def typify_list(d):
    return [typify(v) for v in d]


def typify_dict(d):
    return {
        k: typify(v) for k, v in d.items()
    }


def request_arguments_and_context(request):
    parameters = request.get_json()
    context = typify_dict(parameters.get("context"))
    settings = typify_dict(parameters.get("settings"))
    arguments = {k: v for k, v in parameters.items() if k not in ("context", "settings")}

    monosaccharide_filters = settings.pop("monosaccharide_filters", {})
    monosaccharide_filters = MonosaccharideFilterSet.fromdict(monosaccharide_filters)

    return arguments, ApplicationState(settings, context, monosaccharide_filters)


class ApplicationState(object):
    def __init__(self, settings, context, monosaccharide_filters):
        self.settings = settings
        self.context = context
        self.monosaccharide_filters = monosaccharide_filters

    def __getattr__(self, key):
        try:
            return self.settings[key]
        except KeyError:
            try:
                return self.context[key]
            except KeyError:
                raise AttributeError(key)

    def __repr__(self):
        return "ApplicationState\n%r \n %r \n%r" % (self.settings, self.context, self.monosaccharide_filters)
