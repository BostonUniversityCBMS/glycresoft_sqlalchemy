import csv
import os

from glycresoft_sqlalchemy.data_model.observed_ions import SampleRun, MSScan, Decon2LSPeak, DatabaseManager

isos_to_db_map = {
    "abundance": "intensity",
    "average_mw": "average_mass",
    "charge": "charge",
    "fwhm": "full_width_half_max",
    "mono_abundance": "monoisotopic_intensity",
    "mono_plus2_abundance": "monoisotopic_plus_2_intensity",
    "monoisotopic_mw": "monoisotopic_mass",
    "mostabundant_mw": "most_abundant_mass",
    "scan_num": "scan_id",
    "signal_noise": "signal_to_noise"
}

TDecon2LSPeak = Decon2LSPeak.__table__


def parse_decon2ls(isos_path, database_path=None):
    if database_path is None:
        database_path = os.path.splitext(isos_path)[0] + '.db'
    dbm = DatabaseManager(database_path)
    dbm.initialize()
    session = dbm.session()
    sample_run = SampleRun(name=os.path.basename(isos_path), parameters={"deconvoluted_by": 'decon2ls'})
    session.add(sample_run)
    session.commit()
    interval = 100000
    last = 0
    conn = session.connection()
    last_scan_id = -1
    for i, row in enumerate(csv.DictReader(open(isos_path))):
        remap = {v: row[k] for k, v in isos_to_db_map.items()}
        if remap['scan_id'] != last_scan_id:
            session.add(MSScan(id=remap['scan_id'], sample_run_id=sample_run.id))
            last_scan_id = remap['scan_id']
            if last + interval == i:
                session.commit()
                print "Commit!"
                last = i
                conn = session.connection()
        conn.execute(TDecon2LSPeak.insert().values(**remap))
    session.commit()
    return dbm
