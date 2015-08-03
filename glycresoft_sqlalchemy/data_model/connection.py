import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError

from .data_model import Base
from ..utils import database_utils

import logging

logger = logging.getLogger("database_manager")


class ConnectionManager(object):
    echo = False

    def __init__(self, database_uri, database_uri_prefix, connect_args=None):
        self.database_uri = database_uri
        self.database_uri_prefix = database_uri_prefix
        self.connect_args = connect_args or {}

    def connect(self):
        url = self.construct_url()
        return create_engine(
            url,
            echo=self.echo,
            connect_args=self.connect_args)

    def construct_url(self):
        return "{}{}".format(self.database_uri_prefix, self.database_uri)

    def clear(self):
        pass

    def bridge_address(self):
        '''
        When a ConnectionManager backend is reasonably able to be used for concurrent
        reads and writes (assumed the default), and a new ConnectionManager will be needed
        for another task, provide this interface to transparently share the database connection
        information.

        For backends which do not support concurrent read/write access cleanly (like SQLite),
        override this method to return None and have the default policy on the new task control
        how to create a new storage backend.
        '''
        return self

    def ensure_database(self):
        url = self.construct_url()
        logger.info("Checking %s for database", url)
        if not database_utils.database_exists(url):
            logger.info("No database exists at %s", url)
            database_utils.create_database(url)

    def __repr__(self):
        return '<{} {}{}>'.format(self.__class__.__name__, self.database_uri_prefix, self.database_uri)


class SQLiteConnectionManager(ConnectionManager):
    connect_args = {"timeout": 300}
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

    def bridge_address(self):
        return None

    def connect(self):
        engine = super(SQLiteConnectionManager, self).connect()
        engine.execute("PRAGMA journal_mode=WAL;")
        return engine


class LocalPostgresConnectionManager(ConnectionManager):
    connect_args = {}
    database_uri_prefix = "postgresql+psycopg2://"

    def __init__(self, path, connect_args=None):
        if connect_args is None:
            connect_args = self.connect_args
        super(LocalPostgresConnectionManager, self).__init__(path, self.database_uri_prefix, connect_args)


class DatabaseManager(object):
    connection_manager_type = SQLiteConnectionManager

    def __init__(self, path, clear=False):
        self.connection_manager = self.make_connection_manager(path)
        if clear:
            self.connection_manager.clear()
        self.path = path

    def make_connection_manager(self, path):
        if isinstance(path, ConnectionManager):
            return path
        else:
            return self.connection_manager_type(path)

    def connect(self):
        return self.connection_manager.connect()

    def initialize(self, conn=None):

        if conn is None:
            conn = self.connect()
        try:
            conn.execute("SELECT id FROM Hypothesis LIMIT 1;")
        except Exception:
            logger.exception("Database requires initialization")
            self.connection_manager.ensure_database()
            Base.metadata.create_all(conn)

    def session(self, connection=None):
        if connection is None:
            connection = self.connect()
        return sessionmaker(bind=connection)()

    def bridge_address(self):
        return self.connection_manager.bridge_address()

    def __repr__(self):
        return '<{} {}>'.format(self.__class__.__name__, self.connection_manager)


class LocalPostgresDatabaseManager(DatabaseManager):
    connection_manager_type = LocalPostgresConnectionManager

    def __init__(self, path, clear=False, **kwargs):
        super(LocalPostgresDatabaseManager, self).__init__(path, clear)


def create_remote_manager(driver, user="", password="", host='localhost', port=5432):
    database_uri_prefix = "{driver}://{user}:{password}@{host}:{port}".format(**locals())

    class RemoteConnectionManager(ConnectionManager):
        database_uri_prefix = database_uri_prefix

        def __init__(self, path, connect_args):
            super(RemoteConnectionManager, self).__init__(path, self.database_uri_prefix, connect_args)


def initialize(database_path, clear=False):
    manager = DatabaseManager(database_path, clear=clear)
    manager.initialize()
    return manager


def session(database_path):
    manager = DatabaseManager(database_path)
    return manager.session()
