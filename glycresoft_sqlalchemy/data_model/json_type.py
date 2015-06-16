import json
import sqlalchemy.types


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
