import os

ec = os.system(
    r"glycresoft-database-search ms1 -n 4 -i datafiles\20140918_01_isos.db datafiles\hybrid_ms1.db 1"
    r" -p db -g 2e-5")
assert ec == 0
