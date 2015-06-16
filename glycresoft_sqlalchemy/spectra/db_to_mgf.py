from glypy.search.spectra.spectrum_model import MSMSSqlDB, mass_charge_ratio
from pyteomics import mgf


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
        intensity_array.append(peak.get('intensity', 1))
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


def precursor_to_tandem_ions(precursor):
    header = {}
    if precursor.charge is not None:
        header['precursor_mass'] = precursor.neutral_mass
        header['precursor_mz'] = mass_charge_ratio(precursor.neutral_mass, precursor.charge)
        header['precursor_charge'] = precursor.charge
    else:
        header['precursor_mz'] = precursor.neutral_mass
    header['precursor_id'] = precursor.id

    mz_array = []
    intensity_array = []
    charge_array = []

    for peak in precursor.tandem_data:
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
    scan_entry, tandem_queue = scan_to_ions(scan, **kwargs)
    yield scan_entry
    for precursor in tandem_queue:
        yield precursor_to_tandem_ions(precursor)


def process_database(msmsdb, **kwargs):
    for scan_id in msmsdb.scan_ids():
        for spectra in process_scan(msmsdb[scan_id], scan_id=scan_id):
            yield spectra


def to_mgf(msmsdb, output_path=None, header=None, **kwargs):
    if header is None:
        header = {}
    return mgf.write(process_database(msmsdb, **kwargs), output_path, header)


if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    database = MSMSSqlDB(args[0])
    try:
        output_path = args[1]
        if output_path == '-':
            output_path = None
    except:
        output_path = None
    to_mgf(database, output_path)
