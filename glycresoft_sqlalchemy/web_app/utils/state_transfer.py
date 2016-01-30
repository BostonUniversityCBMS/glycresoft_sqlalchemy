from glycresoft_sqlalchemy.web_app.utils.cache import MonosaccharideFilterSet


def request_arguments_and_context(request):
    parameters = request.get_json()
    context = parameters.get("context")
    settings = parameters.get("settings")
    arguments = {k: v for k, v in parameters.items() if k not in ("context", "settings")}

    monosaccharide_filters = settings.get("monosaccharide_filters", {})
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
        return "ApplicationState\n%r\n%r\n" % (self.settings, self.context, self.monosaccharide_filters)
