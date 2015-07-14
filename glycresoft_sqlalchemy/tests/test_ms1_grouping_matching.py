import os

ec = os.system(
    r"glycresoft-database-search ms1 -n 4 -i datafiles\20141128_05_isos_AGP.db datafiles\integrated_omics_simple.db 1"
    r" -p db --skip-grouping --skip-matching --hypothesis-sample-match-id 2")
assert ec == 0
