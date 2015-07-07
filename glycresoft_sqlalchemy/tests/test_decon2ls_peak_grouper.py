import logging
try:
    logging.basicConfig(level=logging.DEBUG, filemode='w',
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from glycresoft_sqlalchemy.matching import peak_grouping
from glycresoft_sqlalchemy.spectra import decon2ls_sa
import os


def test_main():
    try:
        os.remove("datafiles/20141128_05_isos_AGP.db")
    except:
        pass
    decon2ls_sa.Decon2LSIsosParser("datafiles/20141128_05_isos_AGP.csv")
    job = peak_grouping.Decon2LSPeakGrouper("datafiles/20141128_05_isos_AGP.db", n_processes=6)
    job.start()


if __name__ == '__main__':
    test_main()
