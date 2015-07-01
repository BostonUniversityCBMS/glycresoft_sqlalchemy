from setuptools import setup, find_packages

setup(
  name='glycresoft_sqlalchemy',
  version='0.0.5',
  packages=find_packages(),
  entry_points={
            'console_scripts': [
                "glycresoft-build-database = glycresoft_sqlalchemy.app.build_database:main",
                "glycresoft-database-search = glycresoft_sqlalchemy.app.run_search:main"
            ],
        },
  include_package_data=True,
)
