import os
import yaml
import itertools
import logging

from ..data_model import DatabaseManager
from ..data_model.observed_ions import SampleRun, MSScan, TandemScan, Peak
from . import neutral_mass
from .constants import constants as ms_constants

try:
    logger = logging.getLogger("bupid_topdown_deconvoluter")
except:
    pass

try:
    raise Exception()
    from yaml import CLoader as Loader
except:
    from yaml import Loader

ScalarEvent = yaml.events.ScalarEvent
MappingStartEvent = yaml.events.MappingStartEvent
MappingEndEvent = yaml.events.MappingEndEvent
SequenceStartEvent = yaml.events.SequenceStartEvent
SequenceEndEvent = yaml.events.SequenceEndEvent

ScalarNode = yaml.nodes.ScalarNode

BEGIN = 0
SCANS = 1
PEAKS = 2


class NullDict(dict):
    def __setitem__(self, k, v):
        if isinstance(v, int):
            dict.__setitem__(self, k, v)


class BUPIDMSMSYamlParser(object):
    manager_type = DatabaseManager

    def __init__(self, file_path, database_path=None):
        if database_path is None:
            database_path = os.path.splitext(file_path)[0] + '.db'
        self.file_path = file_path
        self.manager = self.manager_type(database_path)
        self.file_handle = open(file_path)
        self.loader = Loader(self.file_handle)
        self.loader.constructed_objects = NullDict()

        self.state = BEGIN
        self.manager.initialize()
        if not self._is_in_database():
            self.sample_run_name = os.path.basename(file_path)
            return
        session = self.manager.session()
        self.sample_run = SampleRun(name=os.path.basename(file_path), parameters={
            "file_path": file_path,
            "deconvoluted_with": "BUPID Top Down Deconvoluter"
            })
        self.sample_run_name = self.sample_run.name
        session.add(self.sample_run)
        session.commit()
        self.parse()

    def _is_in_database(self):
        """Check if this sample name is already present in the target database file.

        Returns
        -------
        bool : Whether the sample name is present in the database already
        """
        session = self.manager.session()
        result = session.query(
            SampleRun.id).filter(
            SampleRun.name == os.path.basename(
                self.file_path)).count() == 0
        session.close()
        return result

    def seek_scans(self):
        """Advance the parse stream to the beginning of the scans section
        """
        seeking = True
        loader = self.loader
        while seeking:
            if not hasattr(loader.peek_event(), 'value') or loader.peek_event().value != "scan":
                loader.get_event()
            else:
                seeking = False
        loader.get_event()
        self.state = SCANS
        return loader

    def process_scans(self):
        """Step the parser through the scans section, eliminating as many
        object creation steps as possible to minimize memory consumption
        """
        logger.info("Process Scans (State: %r)", self.state)
        if self.state == BEGIN:
            loader = self.seek_scans()
        elif self.state == SCANS:
            loader = self.loader
        else:
            raise Exception("Incorrect state: %s" % self.state)
        has_more_scans = True
        session = self.manager.session()
        while has_more_scans:
            next_event = loader.peek_event()
            if isinstance(next_event, SequenceStartEvent):
                loader.get_event()
            elif isinstance(next_event, MappingEndEvent):
                loader.get_event()
            elif isinstance(next_event, MappingStartEvent):
                node = loader.compose_mapping_node(next_event.anchor)
                mapping = loader.construct_mapping(node)
                # The only purpose of this inner loop is to register the
                # scan references with the YAML parser.

                # peak_id = mapping['id']
                # scan = MSScan()
                # session.add(scan)
                # session.commit()
                # z = mapping['z']
                # precursor_peak = Peak(
                #     id=peak_id, charge=z,
                #     neutral_mass=neutral_mass(mapping['mz'], z), scan_id=scan.id)
                # session.add(precursor_peak)
                # yield scan, precursor_peak

            elif isinstance(next_event, ScalarEvent):
                next_event = loader.get_event()
            elif isinstance(next_event, SequenceEndEvent):
                has_more_scans = False
            else:
                raise TypeError("Unexpected Event: %r" % loader.get_event())
        session.commit()
        session.close()

    def seek_peaks(self):
        seeking = True
        loader = self.loader
        while seeking:
            if not hasattr(loader.peek_event(), 'value') or loader.peek_event().value != "peaks":
                loader.get_event()
            else:
                seeking = False
        loader.get_event()
        self.state = PEAKS
        return loader

    def process_peaks(self):
        logger.info("Process Peaks (State: %r)", self.state)
        if self.state == BEGIN:
            self.process_scans()
        if self.state == SCANS:
            loader = self.seek_peaks()
        elif self.state == PEAKS:
            loader = self.loader
        else:
            raise Exception("Incorrect state: %s" % self.state)

        has_more_peaks = True
        while has_more_peaks:
            next_event = loader.peek_event()
            if isinstance(next_event, SequenceStartEvent):
                loader.get_event()
            elif isinstance(next_event, MappingEndEvent):
                loader.get_event()
            elif isinstance(next_event, MappingStartEvent):
                node = loader.compose_mapping_node(next_event.anchor)
                loader.anchors.pop(next_event.anchor)
                self.resolve_peak(node)
            elif isinstance(next_event, SequenceEndEvent):
                has_more_peaks = False

    def resolve_peak(self, node):
        loader = self.loader
        mass = []
        charge = []
        intensity = []
        scans = []
        for key, value in node.value:
            if key.value == 'id':
                scan_id = loader.construct_scalar(value)
            elif key.value == 'num':
                continue
            elif key.value == 'scans':
                if isinstance(value, ScalarNode):
                    scans = [loader.construct_scalar(value)]
                else:
                    scans = map(loader.construct_mapping, value.value)
            elif key.value in {'mass', 'z', 'intensity'}:
                if isinstance(value, ScalarNode):
                    payload = [loader.construct_scalar(value)]
                else:
                    payload = loader.construct_sequence(value)
                if key.value == 'mass':
                    mass = payload
                elif key.value == 'z':
                    charge = payload
                elif key.value == 'intensity':
                    intensity = payload
                else:
                    raise TypeError("Unknown value series %r" % key)

        precursor_neutral_mass = neutral_mass(scans[0]['mz'], scans[0]['z'])
        scan = TandemScan(
            id=scan_id, precursor_neutral_mass=precursor_neutral_mass, sample_run_id=self.sample_run.id)

        tandem_peaks = [None] * len(mass)
        for i in range(len(mass)):
            tandem_peaks[i] = Peak(neutral_mass=mass[i], charge=charge[i], intensity=intensity[i], scan_id=scan_id)

        session = self.manager.session()
        session.add(scan)
        session.bulk_save_objects(tandem_peaks)
        session.commit()
        session.close()

    def parse(self):
        self.process_scans()
        self.process_peaks()
        self.loader.dispose()
        self.file_handle.close()

    def to_db(self):
        return self.manager
