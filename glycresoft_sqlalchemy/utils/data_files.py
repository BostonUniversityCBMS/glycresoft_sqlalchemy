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

display_store = pjoin(data_directory, "display_store.db")
unimod_store = pjoin(data_directory, "unimod.db")
glycomedb_store = pjoin(data_directory, "glycome-db.db")


def make_absolute_sqlite_sqlalchemy_uri(path):
    return "sqlite:///%s" % path


def configure_color_store():
    '''Use a disk-basd data-store to persist color assignments
    '''
    color_map = colors._color_mapper.color_name_map
    cmap = sqlitedict.open(display_store, "colors", autocommit=True)
    cmap.update(color_map)
    colors._color_mapper.color_name_map = cmap

configure_color_store()
unimod.load = partial(unimod.load, make_absolute_sqlite_sqlalchemy_uri(unimod_store))
