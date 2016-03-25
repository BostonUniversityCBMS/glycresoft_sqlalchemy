from . import uid
import os
import tempfile
import shutil


class TempFileManager(object):
    def __init__(self, base_directory=None):
        if base_directory is None:
            base_directory = tempfile.mkdtemp("glycresoft"[::-1])
        self.base_directory = base_directory
        self.cache = {}

    def get(self, key=None):
        if key is None:
            _key = ""
        else:
            _key = key
        name = "%s_%x" % (_key, uid())
        path = os.path.join(self.base_directory, name)
        self.cache[key] = path
        return path

    def clear(self):
        shutil.rmtree(self.base_directory)

    def __repr__(self):
        return "TempFileManager(%s)" % self.base_directory
