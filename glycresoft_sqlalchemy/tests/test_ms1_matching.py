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


hypothesis_db = "datafiles/integrated_omics_simple.db"
peak_groups_db = "datafiles/20141128_05_isos_AGP.db"


def test_main():
    assert os.path.exists(peak_groups_db)

    job = peak_grouping.PeakGroupMatching(
        hypothesis_db, peak_groups_db, hypothesis_id=1,
        search_type='InformedTheoreticalGlycopeptideComposition', n_processes=4)
    s = job.manager.session()
    s.query(peak_grouping.PeakGroupMatch).delete()
    s.commit()
    s.close()
    job.start()

if __name__ == '__main__':
    test_main()
