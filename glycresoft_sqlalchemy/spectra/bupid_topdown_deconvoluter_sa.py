import os
import yaml
import itertools
import re
import logging
import uuid

from ..data_model import DatabaseManager, PipelineModule
from ..data_model.observed_ions import BUPIDDeconvolutedLCMSMSSampleRun, TandemScan, Peak
from . import neutral_mass
from ..utils import sqlitedict
from .constants import constants as ms_constants

try:
    logger = logging.getLogger("bupid_topdown_deconvoluter")
except:
    print("Could not create logger for bupid_topdown_deconvoluter")

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
AliasEvent = yaml.events.AliasEvent

ScalarNode = yaml.nodes.ScalarNode

BEGIN = 0
SCANS = 1
PEAKS = 2
IN_PEAKS = 3
IN_PEAK_SCANS = 4
DONE = 5

states = {v: k for k, v in dict(
    BEGIN=0,
    SCANS=1,
    PEAKS=2,
    IN_PEAKS=3,
    IN_PEAK_SCANS=4,
    DONE=5
).items()}


int_pattern = re.compile(r":int$")
float_pattern = re.compile(r":float$")


class StreamingYAMLPusher(object):
    def __init__(self, file_handle):
        self.file_handle = file_handle
        self.loader = Loader(file_handle)
        self.state = BEGIN
        self.scan_counter = 0

    def clean(self, event):
        del event.start_mark
        del event.end_mark
        return event

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
        if self.state == BEGIN:
            loader = self.seek_scans()
        elif self.state == SCANS:
            loader = self.loader
        else:
            raise Exception("Incorrect state: %s" % self.state)
        has_more_scans = True
        while has_more_scans:
            next_event = loader.peek_event()
            if isinstance(next_event, SequenceStartEvent):
                loader.get_event()
            elif isinstance(next_event, MappingEndEvent):
                loader.get_event()
                yield self.clean(next_event)
            elif isinstance(next_event, MappingStartEvent):
                loader.get_event()
                yield self.clean(next_event)
                self.scan_counter += 1

            elif isinstance(next_event, ScalarEvent):
                loader.get_event()
                yield self.clean(next_event)

            elif isinstance(next_event, SequenceEndEvent):
                has_more_scans = False
            else:
                raise TypeError("Unexpected Event: %r" % loader.get_event())

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
            next_event = loader.get_event()
            # logger.debug("Next event = %r, State = %r", next_event, states[self.state])
            # print 2, next_event, states[self.state]
            if next_event is None:
                has_more_peaks = False
                continue
            if isinstance(next_event, SequenceStartEvent):
                yield self.clean(next_event)
                if self.state is PEAKS:
                    self.state = IN_PEAKS
            elif isinstance(next_event, MappingEndEvent):
                yield self.clean(next_event)

                if self.state is IN_PEAK_SCANS:
                    self.state = IN_PEAKS
                if self.state is PEAKS:
                    has_more_peaks = False
                    logger.info("No more peaks")
            elif isinstance(next_event, MappingStartEvent):
                yield self.clean(next_event)
                if self.state is PEAKS:
                    self.state = IN_PEAKS
                elif self.state is IN_PEAKS:
                    self.state = IN_PEAK_SCANS
            elif isinstance(next_event, SequenceEndEvent):
                yield self.clean(next_event)

                if self.state is IN_PEAKS:
                    self.state = PEAKS
                # elif self.state is PEAKS:
                #     has_more_peaks = False
                #     logger.info("No more peaks")
            elif isinstance(next_event, ScalarEvent):
                if self.state is PEAKS:
                    has_more_peaks = False
                    continue
                yield self.clean(next_event)

            elif isinstance(next_event, AliasEvent):
                yield self.clean(next_event)
            else:
                raise Exception("Unexpected Event %r" % next_event)

    def parse(self):
        try:
            yield (SCANS)
            for event in self.process_scans():
                yield event
            yield (PEAKS)
            for event in self.process_peaks():
                yield event
            yield (DONE)
        except Exception, e:
            print e, type(e)
            raise
        self.loader.dispose()
        self.file_handle.close()


class StreamingYAMLRenderer(object):
    def __init__(self, database_path, source, sample_run_id=None):
        self.database_path = database_path
        self.manager = DatabaseManager(database_path)
        self.source = source
        self.anchors = sqlitedict.open()
        self.cache = dict()
        self.sample_run_id = sample_run_id

    def run(self):
        source = self.source
        state = BEGIN
        counter = 0
        scan_time = None
        scans = []
        tandem_peaks = []
        mass = []
        charge = []
        intensity = []

        current_map = None
        current_list = None
        current_key = None
        current_anchor = None
        session = self.manager.session()
        for event in source:
            # logger.debug("Next event = %r", event)
            if state is BEGIN:
                if event == SCANS:
                    state = SCANS
                    logger.info("State -> SCANS")
            if state is SCANS:
                if event == PEAKS:
                    state = PEAKS
                    logger.info("State -> PEAKS")

                elif isinstance(event, MappingStartEvent):
                    current_anchor = event.anchor
                    current_map = {}
                elif isinstance(event, MappingEndEvent):
                    self.anchors[current_anchor] = current_map
                    current_map = None
                    current_anchor = None
                elif isinstance(event, ScalarEvent):
                    if event.implicit[1]:
                        current_key = event.value
                    else:
                        if current_list is not None:
                            current_list.append(event.value)
                        elif current_key is not None:
                            # This section contains only ints or floats
                            if int_pattern.match(event.tag):
                                value = int(event.value)
                            else:
                                value = float(event.value)
                            current_map[current_key] = value
            elif state is PEAKS or state is IN_PEAKS:
                if event == DONE:
                    logger.info("Done!")
                    break
                elif isinstance(event, MappingStartEvent):
                    current_map = {}
                    current_anchor = event.anchor
                    state = IN_PEAKS
                elif isinstance(event, MappingEndEvent):
                    if len(scans) == 0:
                        continue
                    scan_time = scans[0]['id']
                    precursor_charge_state = scans[0]['z']
                    precursor_neutral_mass = neutral_mass(scans[0]['mz'], precursor_charge_state)
                    scan = TandemScan(
                        time=scan_time, precursor_charge_state=precursor_charge_state,
                        precursor_neutral_mass=precursor_neutral_mass, sample_run_id=self.sample_run_id)
                    session.add(scan)
                    session.flush()
                    scan_id = scan.id
                    tandem_peaks = []
                    for i in range(len(mass)):
                        if charge[i] > precursor_charge_state:
                            continue
                        tandem_peaks.append(dict(
                            neutral_mass=mass[i], charge=charge[i],
                            intensity=intensity[i], scan_id=scan_id,
                            scan_peak_index=i))
                    session.bulk_insert_mappings(Peak, tandem_peaks)
                    mass = []
                    charge = []
                    intensity = []
                    scans = []
                    counter += 1
                    if counter % 1000 == 0:
                        logger.info("%d scans done.", counter)
                    state = PEAKS
                elif isinstance(event, ScalarEvent):
                    if event.implicit[1]:
                        current_key = event.value
                        if current_key == 'z':
                            current_list = charge.append
                        elif current_key == 'mass':
                            current_list = mass.append
                        elif current_key == 'intensity':
                            current_list = intensity.append
                        elif current_key == 'scans':
                            current_list = scans.append
                    else:
                        if int_pattern.match(event.tag):
                            value = int(event.value)
                        else:
                            value = float(event.value)
                        if current_key == 'id':
                            pass
                        elif current_list is not None:
                            current_list(value)

                elif isinstance(event, AliasEvent):
                    if current_list is not None:
                        current_list(self.anchors[event.anchor])
                    else:
                        raise Exception("Expected an active sequence when recovering anchor")

        session.commit()
        session.close()


class BUPIDMSMSYamlParser(PipelineModule):
    def __init__(self, file_path, database_path=None):
        if database_path is None:
            database_path = os.path.splitext(file_path)[0] + '.db'
        self.file_path = file_path
        self.producer = None
        self.consumer = None
        self.queue = None
        self.manager = DatabaseManager(database_path)
        self.manager.initialize()
        if not self._is_in_database():
            self.sample_run_name = os.path.basename(file_path)
            return
        session = self.manager.session()
        self.sample_run = BUPIDDeconvolutedLCMSMSSampleRun(name=os.path.basename(file_path), parameters={
            "file_path": file_path,
            "deconvoluted_with": "BUPID Top Down Deconvoluter"
            }, uuid=uuid.uuid4().hex)
        self.sample_run_name = self.sample_run.name
        session.add(self.sample_run)
        session.commit()
        self.sample_run_id = self.sample_run.id
        self.start()

    def _is_in_database(self):
        """Check if this sample name is already present in the target database file.

        Returns
        -------
        bool : Whether the sample name is present in the database already
        """
        session = self.manager.session()
        result = session.query(
            BUPIDDeconvolutedLCMSMSSampleRun.id).filter(
            BUPIDDeconvolutedLCMSMSSampleRun.name == os.path.basename(
                self.file_path)).count() == 0
        session.close()
        return result

    def run(self):
        self.producer = StreamingYAMLPusher(open(self.file_path, 'rb'))
        event_stream = self.producer.parse()
        self.consumer = StreamingYAMLRenderer(self.database_path, event_stream, self.sample_run_id)
        self.consumer.run()


def process_data_file(file_path, database_path=None):
    job = BUPIDMSMSYamlParser(file_path, database_path)
