import sys
from setuptools import setup, find_packages, Extension

# With gratitude to the SqlAlchemy setup.py authors

from distutils.command.build_ext import build_ext
from distutils.errors import (CCompilerError, DistutilsExecError,
                              DistutilsPlatformError)

ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)
if sys.platform == 'win32':
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    ext_errors += (IOError,)

c_ext = "pyx"
try:
    from Cython.Build import cythonize
except:
    c_ext = "c"

extensions = [
    Extension("glycresoft_sqlalchemy.structure.composition.ccomposition",
              ["glycresoft_sqlalchemy/structure/composition/ccomposition." + c_ext])
]

if c_ext == "pyx":
    extensions = cythonize(extensions, annotate=True)

cmdclass = {}


class BuildFailed(Exception):

    def __init__(self):
        self.cause = sys.exc_info()[1]  # work around py 2/3 different syntax

    def __str__(self):
        return str(self.cause)


class ve_build_ext(build_ext):
    # This class allows C extension building to fail.

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            raise BuildFailed()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except ext_errors:
            raise BuildFailed()
        except ValueError:
            # this can happen on Windows 64 bit, see Python issue 7511
            if "'path'" in str(sys.exc_info()[1]):  # works with both py 2/3
                raise BuildFailed()
            raise

cmdclass['build_ext'] = ve_build_ext


required = []
with open('requirements.txt') as f:
    required = f.read().splitlines()


def status_msgs(*msgs):
    print('*' * 75)
    for msg in msgs:
        print(msg)
    print('*' * 75)


def run_setup(include_cext=True):
    setup(
      name='glycresoft_sqlalchemy',
      version='0.0.5',
      packages=find_packages(),
      install_requires=required,
      entry_points={
                'console_scripts': [
                    "glycresoft-build-database = glycresoft_sqlalchemy.app.build_database:main",
                    "glycresoft-database-search = glycresoft_sqlalchemy.app.run_search:main",
                    "glycresoft-report = glycresoft_sqlalchemy.app.reporting:taskmain"
                ],
            },
      cmdclass=cmdclass,
      zip_safe=False,
      include_package_data=True,
      package_data={
        'glycresoft_sqlalchemy': ["*.csv", "*.xml", "*.json", "data/*.csv"],
        'glycresoft_sqlalchemy.structure': ["structure/data/*.csv", "structure/data/*.json"]
        },
      ext_modules=extensions if include_cext else None,

    )

try:
    run_setup(True)
except Exception as exc:
    status_msgs(
        str(exc),
        "WARNING: The C extension could not be compiled, " +
        "speedups are not enabled.",
        "Failure information, if any, is above.",
        "Retrying the build without the C extension now."
    )

    run_setup(False)

    status_msgs(
        "WARNING: The C extension could not be compiled, " +
        "speedups are not enabled.",
        "Plain-Python build succeeded."
    )
