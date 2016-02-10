import os

from pyteomics import mgf

from glycresoft_sqlalchemy.data_model import PipelineModule, BUPIDDeconvolutedLCMSMSSampleRun
from .spectrum_model import mass_charge_ratio


def scan_to_ions(scan, **kwargs):
    header = {}
    tandem_queue = []
    mz_array = []
    intensity_array = []
    charge_array = []

    for peak in scan:
        if peak.charge is not None:
            mz_array.append(mass_charge_ratio(peak.neutral_mass, peak.charge))
            charge_array.append(peak.charge)
        else:
            mz_array.append(peak.neutral_mass)
        intensity_array.append(peak.intensity)
        if len(peak.tandem_data) > 0:
            tandem_queue.append(peak)
    header.update(kwargs)
    ions = {
        "m/z array": mz_array,
        "intensity array": intensity_array,
        "charge array": charge_array,
        "params": header
    }
    if all(e is None for e in ions['charge array']):
        ions.pop('charge array')
    return ions, tandem_queue


def precursor_to_tandem_ions(tandem_scan):
    header = {}
    if tandem_scan.precursor_charge_state is not None:
        header['precursor_mass'] = tandem_scan.precursor_neutral_mass
        header['precursor_mz'] = mass_charge_ratio(
            tandem_scan.precursor_neutral_mass, tandem_scan.precursor_charge_state)
        header['precursor_charge'] = tandem_scan.precursor_charge_state
    else:
        header['precursor_mz'] = tandem_scan.precursor_neutral_mass
    header['precursor_id'] = tandem_scan.id
    header['scan_time'] = tandem_scan.time
    mz_array = []
    intensity_array = []
    charge_array = []

    for peak in tandem_scan.tandem_data:
        if peak.charge is not None:
            mz_array.append(mass_charge_ratio(peak.neutral_mass, peak.charge))
            charge_array.append(peak.charge)
        else:
            mz_array.append(peak.neutral_mass)
        intensity_array.append(peak.intensity)
    ions = {
        "m/z array": mz_array,
        "intensity array": intensity_array,
        "charge array": charge_array,
        "params": header
    }
    if all(e is None for e in ions['charge array']):
        ions.pop('charge array')
    return ions


def process_scan(scan, **kwargs):
    return precursor_to_tandem_ions(scan)


def process_database(sample_run, **kwargs):
    for tandem_scan in sample_run.tandem_scans.yield_per(100):
        yield process_scan(tandem_scan, scan_id=tandem_scan.id)


def to_mgf(sample_run, output_path=None, header=None, **kwargs):
    if header is None:
        header = {}
    return mgf.write(process_database(sample_run, **kwargs), output_path, header)


class BUPIDToMGFConverter(PipelineModule):
    def __init__(self, database_path, sample_run_id, output_path=None):
        self.manager = self.manager_type(database_path)
        self.sample_run_id = sample_run_id
        self.output_path = output_path

    def run(self):
        session = self.manager()
        sample_run = session.query(BUPIDDeconvolutedLCMSMSSampleRun).get(self.sample_run_id)

        if sample_run is None:
            raise ValueError("sample_run_id %d did not refer to a BUPID Sample Run" % self.sample_run_id)

        header = {
            "deconvoluted_with": "BUPID Top-Down"
        }

        if self.output_path is None:
            self.output_path = os.path.join(os.path.dirname(
                self.database_path), os.path.splitext(sample_run.name)[0] + '.mgf')

        to_mgf(sample_run, self.output_path, header)
