import os

from .vendor import sqlitedict, appdir

from glycresoft_sqlalchemy.report import colors

dirs = appdir.AppDirs("GlycReSoft", "Zaia Lab", "1.0", roaming=True)

pjoin = os.path.join

data_directory = dirs.user_data_dir
cache_directory = dirs.user_cache_dir

if not os.path.exists(data_directory):
    os.makedirs(data_directory)

display_store = pjoin(data_directory, "display_store.db")


def configure_color_store():
    '''Use a disk-basd data-store to persist color assignments
    '''
    color_map = colors._color_mapper.color_name_map
    cmap = sqlitedict.open(display_store, "colors", autocommit=True)
    cmap.update(color_map)
    colors._color_mapper.color_name_map = cmap

configure_color_store()
