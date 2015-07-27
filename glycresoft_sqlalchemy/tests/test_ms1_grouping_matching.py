import os

ec = os.system(
    r"glycresoft-database-search ms1 -n 5 -i datafiles\20140918_01_isos.db datafiles\naive_glycopeptide.db 1"
    r" -p db -g 2e-5 --skip-grouping")
assert ec == 0
