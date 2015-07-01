import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .data_model import Base


class ConnectionManager(object):
    echo = False

    def __init__(self, database_uri, database_uri_prefix, connect_args=None):
        self.database_uri = database_uri
        self.database_uri_prefix = database_uri_prefix
        self.connect_args = connect_args or {}

    def connect(self):
        return create_engine(
            "{}{}".format(self.database_uri_prefix, self.database_uri),
            echo=self.echo,
            connect_args=self.connect_args)

    def clear(self):
        pass


class SQLiteConnectionManager(ConnectionManager):
    connect_args = {"timeout": 100}
    database_uri_prefix = "sqlite:///"

    def __init__(self, path, connect_args=None):
        if connect_args is None:
            connect_args = self.connect_args
        super(SQLiteConnectionManager, self).__init__(path, self.database_uri_prefix, connect_args)

    def clear(self):
        try:
            os.remove(self.database_uri)
        except:
            pass


class DatabaseManager(object):
    connection_manager_type = SQLiteConnectionManager

    def __init__(self, path, clear=False):
        self.connection_manager = self.connection_manager_type(path)
        if clear:
            self.connection_manager.clear()
        self.path = path

    def connect(self):
        return self.connection_manager.connect()

    def initialize(self, conn=None):
        if conn is None:
            conn = self.connect()
        Base.metadata.create_all(conn)

    def session(self, connection=None):
        if connection is None:
            connection = self.connect()
        return sessionmaker(bind=connection)()


def initialize(database_path):
    manager = DatabaseManager(database_path, clear=True)
    manager.initialize()
    return manager


def session(database_path):
    manager = DatabaseManager(database_path)
    return manager.session()
