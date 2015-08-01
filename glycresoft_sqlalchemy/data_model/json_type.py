import json
import sqlalchemy.types
from sqlalchemy.ext.declarative import DeclarativeMeta


class JSONType(sqlalchemy.types.PickleType):
    """An example

       >>> class User(declarative_base):
       ...     friends = Column(JSONType())
       ...
       ...     def update_friends(self, access_token):
       ...          graph = facebook.GraphAPI(access_token)
       ...          profile = graph.get_object("me")
       ...          friends = graph.get_connections("me", "friends")
       ...          if "data" in friends:
       ...              self.friends = dict(
       ...                  [ (f["id"], f["name"]) for f in friends["data"] ]
       ...              )

    """
    impl = sqlalchemy.types.UnicodeText

    def __init__(self):
        sqlalchemy.types.PickleType.__init__(self, pickler=json)


def tryjson(obj):
    try:
        return obj.to_json()
    except:
        return None


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
                        fields[field] = None
                # a json-encodable dict
                return fields

            return json.JSONEncoder.default(self, obj)
    return AlchemyEncoder
