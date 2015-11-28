import json
import sqlalchemy.types
import types
from sqlalchemy.ext.declarative import DeclarativeMeta


class JSONType(sqlalchemy.types.PickleType):
    impl = sqlalchemy.types.UnicodeText

    def __init__(self):
        sqlalchemy.types.PickleType.__init__(self, pickler=json)


def tryjson(obj):
    try:
        return obj.to_json()
    except:
        return None


def clean_dict(d):
    result = {}
    for k, v in d.items():
        if not isinstance(k, (str, int)):
            if isinstance(k.__class__, DeclarativeMeta):
                rk = k.id
            else:
                rk = str(k)
        else:
            rk = k
        if isinstance(v, (list, tuple)):
            rv = []
            for x in v:
                if isinstance(x, dict):
                    rv.append(clean_dict(x))
                elif isinstance(x.__class__, DeclarativeMeta):
                    rv.append(tryjson(x))
                else:
                    rv.append(x)
        elif isinstance(v, dict):
            rv = clean_dict(v)
        else:
            rv = v

        result[rk] = rv
    return result


def new_alchemy_encoder():
    class AlchemyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj.__class__, DeclarativeMeta):
                # an SQLAlchemy class
                fields = {}
                for field in [x for x in dir(obj) if not x.startswith('_') and x != 'metadata']:
                    data = obj.__getattribute__(field)
                    try:
                        json.dumps(data)  # this will fail on non-encodable values, like other classes
                        fields[field] = data
                    except TypeError:
                        if isinstance(data, types.MethodType):
                            continue
                        fields[field] = None
                # a json-encodable dict
                return fields

            return json.JSONEncoder.default(self, obj)
    return AlchemyEncoder
