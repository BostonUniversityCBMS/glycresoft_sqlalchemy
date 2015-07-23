from glycresoft_sqlalchemy.data_model import *
from glycresoft_sqlalchemy.structure import *

import IPython


def copyquery_dict(result_proxy):
    keys = [key.split("_", 1)[1] for key in result_proxy.keys()]
    for row in result_proxy:
        d = {keys[i]: value for i, value in enumerate(row)}
        d.pop('id')
        yield d


IPython.embed()
