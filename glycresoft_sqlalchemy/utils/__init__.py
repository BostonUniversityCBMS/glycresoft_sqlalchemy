try:
    import cPickle as pickle
except:
    import pickle


def simple_repr(self):  # pragma: no cover
    template = "{self.__class__.__name__}({d})"
    d = {"%s=%r" % (k, v) for k, v in sorted(self.__dict__.items(), key=lambda x:x[0]) if not k.startswith("_")}
    return template.format(self=self, d=', '.join(d))
