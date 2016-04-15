import os
from functools import partial

from .vendor import sqlitedict, appdir

from glycresoft_sqlalchemy.report import colors
from glycresoft_sqlalchemy.structure.data import unimod

dirs = appdir.AppDirs("GlycReSoft", "Zaia Lab", "1.0", roaming=True)

pjoin = os.path.join

data_directory = dirs.user_data_dir
cache_directory = dirs.user_cache_dir

if not os.path.exists(data_directory):
    os.makedirs(data_directory)

try:
    invalidation_errors = [OSError, WindowsError]
except:
    invalidation_errors = [OSError]


class ResourcePath(str):
    valid = True

    def invalidate(self):
        self.valid = False

    def validate(self):
        if not self.valid:
            if self.exists:
                self.remove()

    def remove(self):
        try:
            os.remove(self)
        except invalidation_errors:
            pass

    @property
    def exists(self):
        return os.path.exists(self)


class Resource(object):
    def __init__(self, name, path, **kwargs):
        self.name = name
        self.path = ResourcePath(path)
        self.held = kwargs.get('held', False)
        self.owners = kwargs.get('owners', set())
        self.ready = kwargs.get("ready", False)

    def __str__(self):
        return self.path

    def __repr__(self):
        return "Resource(name=%r, path=%r)"

    def acquired(self, owner):
        if owner not in self.owners:
            self.owners.add(owner)

    def release(self, owner):
        if owner not in self.owners:
            raise ValueError("%r is not a valid owner" % owner)
        self.owners.remove(owner)
        if len(self.owners) == 0:
            self.held = False


display_store = ResourcePath(pjoin(data_directory, "display_store.db"))

unimod_store = ResourcePath(pjoin(data_directory, "unimod.db"))

glycomedb_store = ResourcePath(pjoin(data_directory, "glycome-db.db"))
glycomedb_download_cache = ResourcePath(pjoin(data_directory, "glycome-db-download-cache"))

taxonomylite_store = ResourcePath(pjoin(data_directory, "taxonomylite.db"))


def make_absolute_sqlite_sqlalchemy_uri(path):
    return "sqlite:///%s" % path


def configure_color_store():
    '''Use a disk-based data-store to persist color assignments
    '''
    color_map = colors._color_mapper.color_name_map
    cmap = sqlitedict.open(display_store, "colors", autocommit=True)
    cmap.update(color_map)
    colors._color_mapper.color_name_map = cmap

configure_color_store()
unimod.load = partial(unimod.load, make_absolute_sqlite_sqlalchemy_uri(unimod_store))
