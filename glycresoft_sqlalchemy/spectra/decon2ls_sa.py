import csv
import os
import uuid

from glycresoft_sqlalchemy.data_model import PipelineModule
from glycresoft_sqlalchemy.data_model.observed_ions import (
    Decon2LSLCMSSampleRun, MSScan, Decon2LSPeak)

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
    """
    Parse an "_isos.csv" file produced by `Decon2LS` into :class:`MSScan` and
    :class:`Decon2LSPeak` records in the database given by :attr:`manager`. These
    records are owned by a newly created :class:`Decon2LSLCMSSampleRun` record associated
    with the file.

    Attributes
    ----------
    file_path : str
        Path to "_isos.csv" file to be parsed
    interval : int
        # of records to produce between commits.
    manager : :class:`DatabaseManager`
        Broker database session
    sample_run : :class:`Decon2LSLCMSSampleRun`
        The `SampleRun` instance associated with this experimental data.

    On creation, if the `database_path` argument is `None`, a new database will be constructed
    extrapolated from `file_path` with the '.db' extension, and the :attr:`manager` will initialize
    it.

    .. warning::
        This `PipelineModule` will start immediately upon instantiation. This behavior is deprecated
        and will change at the earliest opportunity.
    """
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
                    last = i
                    conn = session.connection()
            remap['scan_id'] = last_key
            conn.execute(TDecon2LSPeak.insert().values(**remap))
        session.commit()
