import os
import warnings
from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker
from sqlalchemy.engine.url import make_url
from sqlalchemy.pool import NullPool, SingletonThreadPool

from .base import Base, Namespace
from ..utils import database_utils

import logging

warnings.filterwarnings("ignore", category=DeprecationWarning)
logger = logging.getLogger("database_manager")
logging.getLogger("sqlalchemy.pool.NullPool").disabled = True


class ConnectionManager(object):
    echo = False

    def __init__(self, database_uri, database_uri_prefix="",
                 connect_args=None):
        self.database_uri = str(database_uri)
        self.database_uri_prefix = database_uri_prefix
        self.connect_args = connect_args or {}

    def clear(self):
        pass

    def connect(self):
        url = self.construct_url()
        engine = create_engine(
            url,
            echo=self.echo,
            connect_args=self.connect_args,
            poolclass=NullPool)
        self._configure_creation(engine)
        return engine

    def construct_url(self):
        return "{}{}".format(self.database_uri_prefix, self.database_uri)

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

    def _configure_creation(self, connection):
        pass

    def __repr__(self):
        return '<{} {}>'.format(
            self.__class__.__name__, self.construct_url())


class SQLiteConnectionManager(ConnectionManager):
    connect_args = {"timeout": 300}
    database_uri_prefix = "sqlite:///"
    _memory_database = "sqlite://"

    def __init__(self, path, connect_args=None):
        if connect_args is None:
            connect_args = self.connect_args
        self._path = path
        super(
            SQLiteConnectionManager,
            self).__init__(
            str(path),
            self.database_uri_prefix,
            connect_args)
        self._is_memory = self._memory_database == self.construct_url()
        self._saved_engine = None

    def connect(self):
        if self._is_memory:
            if self._saved_engine is not None:
                return self._saved_engine

        url = self.construct_url()
        if url == self._memory_database:
            pool = SingletonThreadPool
        else:
            pool = NullPool
        engine = create_engine(
            url,
            echo=self.echo,
            connect_args=self.connect_args,
            poolclass=pool)
        self._configure_creation(engine)
        if self._is_memory:
            self._saved_engine = engine
        return engine

    def bridge_address(self):
        # SQLite doesn't support multiple writing
        # connections, so we should not allow other
        # components to bridge through this connection
        return None

    def clear(self):
        open(self._path, 'wb').write("")
        try:
            os.remove(self._path)
        except:
            pass

    def _configure_creation(self, connection):
        def do_connect(dbapi_connection, connection_record):
            # disable pysqlite's emitting of the BEGIN statement entirely.
            # also stops it from emitting COMMIT before any DDL.
            iso_level = dbapi_connection.isolation_level
            dbapi_connection.isolation_level = None
            try:
                dbapi_connection.execute("PRAGMA page_size = 5120;")
                dbapi_connection.execute("PRAGMA cache_size = 4000;")
                dbapi_connection.execute("PRAGMA foreign_keys = ON;")
            except:
                pass
            dbapi_connection.isolation_level = iso_level

        event.listens_for(connection, "connect")(do_connect)

    def construct_url(self):
        if self.database_uri == self._memory_database or self.database_uri == "":
            return self._memory_database
        else:
            return super(SQLiteConnectionManager, self).construct_url()


class NoConfigSQLiteConnectionManager(SQLiteConnectionManager):
    connect_args = {"timeout": 5}

    def _configure_creation(self, connection):
        pass


class LocalPostgresConnectionManager(ConnectionManager):
    connect_args = {}
    database_uri_prefix = "postgresql+psycopg2://"

    def __init__(self, path, connect_args=None):
        if connect_args is None:
            connect_args = self.connect_args
        super(
            LocalPostgresConnectionManager,
            self).__init__(
            path,
            self.database_uri_prefix,
            connect_args)


class DatabaseManager(object):
    connection_manager_type = SQLiteConnectionManager

    def __init__(self, path, clear=False, connection_manager_type=None):
        if connection_manager_type is not None:
            self.connection_manager_type = connection_manager_type
        self.connection_manager = self.make_connection_manager(path)
        if clear:
            self.clear()
        self.path = path

    def clear(self, reinitialize=False):
        Base.metadata.drop_all(bind=self.connect())
        self.connection_manager.clear()
        if reinitialize:
            self.initialize()

    def make_connection_manager(self, path):
        if isinstance(path, ConnectionManager):
            return path
        try:
            path = str(path)
            if path == "":
                return SQLiteConnectionManager("sqlite://")
            # Ensure that the path is in fact a URI
            make_url(path)
            return ConnectionManager(path)
        except:
            return self.connection_manager_type(path)

    def connect(self):
        return self.connection_manager.connect()

    def initialize(self, conn=None, force=False):

        if conn is None:
            conn = self.connect()
        try:
            conn.execute("SELECT id FROM Hypothesis LIMIT 1;")
            if force:
                raise Exception()
        except:
            logger.info("Database requires initialization")
            self.connection_manager.ensure_database()
            Base.metadata.create_all(conn)
            session = self.session()
            for initializer in Namespace.initialization_list:
                initializer(session)

    def session(self, connection=None):
        if connection is None:
            connection = self.connect()
        return sessionmaker(bind=connection, autoflush=False)()

    __call__ = session

    def query(self, *args, **kwargs):
        return self().query(*args, **kwargs)

    def bridge_address(self):
        return self.connection_manager.bridge_address()

    def __repr__(self):
        return '<{} {}>'.format(
            self.__class__.__name__, self.connection_manager)


class LocalPostgresDatabaseManager(DatabaseManager):
    connection_manager_type = LocalPostgresConnectionManager

    def __init__(self, path, clear=False, **kwargs):
        super(LocalPostgresDatabaseManager, self).__init__(path, clear)


def create_remote_manager(
        driver, user="", password="", host='localhost', port=5432):
    database_uri_prefix = "{driver}://{user}:{password}@{host}:{port}".format(
        **locals())

    class RemoteConnectionManager(ConnectionManager):
        database_uri_prefix = database_uri_prefix

        def __init__(self, path, connect_args):
            super(
                RemoteConnectionManager,
                self).__init__(
                path,
                self.database_uri_prefix,
                connect_args)


def initialize(database_path, clear=False):
    manager = DatabaseManager(database_path, clear=clear)
    manager.initialize()
    return manager


def session(database_path):
    manager = DatabaseManager(database_path)
    return manager.session()
