from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class _Namespace(object):
    def __str__(self):
        return "Namespace(%s)" % self.__dict__

Namespace = _Namespace()
Namespace.initialization_list = []
