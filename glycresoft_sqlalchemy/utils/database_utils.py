import os
from copy import copy
from time import time
from sqlalchemy import create_engine, Table
from sqlalchemy.exc import OperationalError, ProgrammingError
from sqlalchemy.sql import ClauseElement
from sqlalchemy.engine.url import make_url
from sqlalchemy.engine.interfaces import Dialect

from sqlalchemy.orm.session import object_session
from sqlalchemy.orm.exc import UnmappedInstanceError


def get_or_create(session, model, defaults=None, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance, False
    else:
        params = dict((k, v) for k, v in kwargs.iteritems() if not isinstance(v, ClauseElement))
        params.update(defaults or {})
        instance = model(**params)
        session.add(instance)
        return instance, True


#: Adopted from sqlalchemy-utils

def get_bind(obj):
    if hasattr(obj, 'bind'):
        conn = obj.bind
    else:
        try:
            conn = object_session(obj).bind
        except UnmappedInstanceError:
            conn = obj

    if not hasattr(conn, 'execute'):
        raise TypeError(
            'This method accepts only Session, Engine, Connection and '
            'declarative model objects.'
        )
    return conn


def quote(mixed, ident):
    """
    Conditionally quote an identifier.
    """
    if isinstance(mixed, Dialect):
        dialect = mixed
    else:
        dialect = get_bind(mixed).dialect
    return dialect.preparer(dialect).quote(ident)


def database_exists(url):
    """Check if a database exists.

    Parameters
    ----------
    url: A SQLAlchemy engine URL.

    Returns
    -------
    bool

    References
    ----------
    http://sqlalchemy-utils.readthedocs.org/en/latest/_modules/sqlalchemy_utils/functions/database.html#database_exists
    """

    url = copy(make_url(url))
    database = url.database
    if url.drivername.startswith('postgresql'):
        url.database = 'template1'
    else:
        url.database = None

    engine = create_engine(url)

    if engine.dialect.name == 'postgresql':
        text = "SELECT 1 FROM pg_database WHERE datname='%s'" % database
        return bool(engine.execute(text).scalar())

    elif engine.dialect.name == 'mysql':
        text = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
                "WHERE SCHEMA_NAME = '%s'" % database)
        return bool(engine.execute(text).scalar())

    elif engine.dialect.name == 'sqlite':
        return database == ':memory:' or os.path.exists(database)

    else:
        text = 'SELECT 1'
        try:
            url.database = database
            engine = create_engine(url)
            engine.execute(text)
            return True

        except (ProgrammingError, OperationalError):
            return False


def create_database(url, encoding='utf8', template=None):
    """Issue the appropriate CREATE DATABASE statement.

    Parameters
    ----------
    url: str
        A SQLAlchemy engine URL.
    encoding: str
        The encoding to create the database as.
    template: str
        The name of the template from which to create the new database.

    References
    ----------
    http://sqlalchemy-utils.readthedocs.org/en/latest/_modules/sqlalchemy_utils/functions/database.html#create_database

    """

    url = copy(make_url(url))

    database = url.database

    if url.drivername.startswith('postgresql'):
        url.database = 'template1'
    elif not url.drivername.startswith('sqlite'):
        url.database = None

    engine = create_engine(url)

    if engine.dialect.name == 'postgresql':
        if engine.driver == 'psycopg2':
            from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
            engine.raw_connection().set_isolation_level(
                ISOLATION_LEVEL_AUTOCOMMIT
            )

        if not template:
            template = 'template0'

        text = "CREATE DATABASE {0} ENCODING '{1}' TEMPLATE {2}".format(
            quote(engine, database),
            encoding,
            quote(engine, template)
        )
        engine.execute(text)

    elif engine.dialect.name == 'mysql':
        text = "CREATE DATABASE {0} CHARACTER SET = '{1}'".format(
            quote(engine, database),
            encoding
        )
        engine.execute(text)

    elif engine.dialect.name == 'sqlite' and database != ':memory:':
        open(database, 'w').close()

    else:
        text = 'CREATE DATABASE {0}'.format(quote(engine, database))
        engine.execute(text)


def drop_database(url):
    """Issue the appropriate DROP DATABASE statement.

    Parameters
    ----------
    url: str

    References
    ----------
    http://sqlalchemy-utils.readthedocs.org/en/latest/_modules/sqlalchemy_utils/functions/database.html#drop_database
    """

    url = copy(make_url(url))

    database = url.database

    if url.drivername.startswith('postgresql'):
        url.database = 'template1'
    elif not url.drivername.startswith('sqlite'):
        url.database = None

    engine = create_engine(url)

    if engine.dialect.name == 'sqlite' and url.database != ':memory:':
        os.remove(url.database)

    elif engine.dialect.name == 'postgresql' and engine.driver == 'psycopg2':
        from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
        engine.raw_connection().set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)

        # Disconnect all users from the database we are dropping.
        version = list(
            map(
                int,
                engine.execute('SHOW server_version').first()[0].split('.')
            )
        )
        pid_column = (
            'pid' if (version[0] >= 9 and version[1] >= 2) else 'procpid'
        )
        text = '''
        SELECT pg_terminate_backend(pg_stat_activity.%(pid_column)s)
        FROM pg_stat_activity
        WHERE pg_stat_activity.datname = '%(database)s'
          AND %(pid_column)s <> pg_backend_pid();
        ''' % {'pid_column': pid_column, 'database': database}
        engine.execute(text)

        # Drop the database.
        text = 'DROP DATABASE {0}'.format(quote(engine, database))
        engine.execute(text)

    else:
        text = 'DROP DATABASE {0}'.format(quote(engine, database))
        engine.execute(text)


class toggle_indices(object):

    def __init__(self, session, table):
        if isinstance(table, type):
            table = table.__table__
        self.table = table
        self.session = session

    def drop(self):
        conn = self.session.connection()
        for index in self.table.indexes:
            index.drop(conn)
        try:
            self.session.commit()
        except:
            pass

    def create(self):
        conn = self.session.connection()
        for index in self.table.indexes:
            index.create(conn)
        try:
            self.session.commit()
        except:
            pass


def temp_table(table, metdata=None):
    if metdata is None:
        if hasattr(table, "metadata"):
            metadata = table.metadata
        else:
            from glycresoft_sqlalchemy.data_model import Base
            metadata = Base.metadata
    if hasattr(table, "__table__"):
        table = table.__table__
    cols = [c.copy() for c in table.columns]
    constraints = [c.copy() for c in table.constraints]
    return Table("TempTable_%s" % hex(int(time())), metadata, *(cols + constraints))
