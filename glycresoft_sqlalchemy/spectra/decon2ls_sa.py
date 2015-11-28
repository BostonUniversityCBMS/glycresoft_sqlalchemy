import csv
import os
import uuid

from glycresoft_sqlalchemy.data_model import PipelineModule
from glycresoft_sqlalchemy.data_model.observed_ions import (
    Decon2LSLCMSSampleRun, MSScan, Decon2LSPeak, DatabaseManager)

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


class Decon2LSIsosParser(PipelineModule):
    def __init__(self, file_path, database_path=None, interval=100000):
        self.file_path = file_path
        if database_path is None:
            database_path = os.path.splitext(file_path)[0] + '.db'
        self.interval = interval
        self.manager = self.manager_type(database_path)
        self.start()

    def run(self):
        self.manager.initialize()
        session = self.manager.session()
        sample_run = Decon2LSLCMSSampleRun(
            name=os.path.basename(self.file_path), parameters={"deconvoluted_by": 'decon2ls'},
            uuid=uuid.uuid4().hex)
        session.add(sample_run)
        session.commit()
        sample_run_id = sample_run.id
        self.sample_run = sample_run

        last = 0
        conn = session.connection()
        last_scan_id = -1
        last_key = -1
        for i, row in enumerate(csv.DictReader(open(self.file_path))):
            remap = {v: row[k] for k, v in isos_to_db_map.items()}
            remap['scan_peak_index'] = i
            if remap['scan_id'] != last_scan_id:
                scan = MSScan(time=remap['scan_id'], sample_run_id=sample_run_id)
                session.add(scan)
                session.flush()
                last_key = scan.id
                last_scan_id = remap['scan_id']
                if last + self.interval == i:
                    session.commit()
                    print "Commit!"
                    last = i
                    conn = session.connection()
            remap['scan_id'] = last_key
            conn.execute(TDecon2LSPeak.insert().values(**remap))
        session.commit()


# def parse_decon2ls(isos_path, database_path=None):
#     if database_path is None:
#         database_path = os.path.splitext(isos_path)[0] + '.db'
#     dbm = DatabaseManager(database_path)
#     dbm.initialize()
#     session = dbm.session()
#     sample_run = Decon2LSLCMSSampleRun(name=os.path.basename(isos_path), parameters={"deconvoluted_by": 'decon2ls'})
#     session.add(sample_run)
#     session.commit()
#     interval = 100000
#     last = 0
#     conn = session.connection()
#     last_scan_id = -1
#     for i, row in enumerate(csv.DictReader(open(isos_path))):
#         remap = {v: row[k] for k, v in isos_to_db_map.items()}
#         if remap['scan_id'] != last_scan_id:
#             session.add(MSScan(id=remap['scan_id'], sample_run_id=sample_run.id))
#             last_scan_id = remap['scan_id']
#             if last + interval == i:
#                 session.commit()
#                 print "Commit!"
#                 last = i
#                 conn = session.connection()
#         conn.execute(TDecon2LSPeak.insert().values(**remap))
#     session.commit()
#     return dbm
