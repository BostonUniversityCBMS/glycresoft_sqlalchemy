from .data_model import DatabaseManager, Experiment

from contextlib import contextmanager
from multiprocessing import Process, Manager
from Queue import Empty as QueueEmpty

SETUP = 0
READY = 1
WAIT = 2
WORKING = 3
DIE = 4


class CommitPooler(object):
    def __init__(self, database_path):
        self.database_path = database_path
        self.manager = Manager()
        self.process = None
        self.control_queue = self.manager.Queue()
        self.work_queue = None

    def queue(self):
        self.work_queue = self.manager.Queue()
        return self.work_queue

    def start(self):
        queues = [self.queue()]
        self.process = CommitPoolerProcess(self.database_path, self.control_queue, queues)
        self.process.start()
        return queues[0]

    def terminate(self):
        self.process.terminate()

    def __del__(self):
        try:
            self.process.terminate()
        except:
            pass

@contextmanager
def sessioncontrol(session):
    print "in"
    yield session
    print "out"
    session.commit()


class CommitPoolerProcess(Process):
    def __init__(self, database_path, control_queue, queues=None):
        Process.__init__(self)
        if queues is None:
            queues = []
        self.control_queue = control_queue
        self.queues = queues
        self.state = READY
        self.database_path = database_path
        self.session = None

    def run(self):
        counter = 0
        self.session = DatabaseManager(self.database_path).session()
        while self.state == READY:
            for queue in self.queues:
                try:
                    if queue.empty():
                        continue
                    item = queue.get(True, 1)
                    self.session.add(item)
                    counter += 1

                    if counter % 100000 == 0:
                        self.session.commit()
                except QueueEmpty:
                    print "No Data"
                    self.session.commit()
                except EOFError:
                    print "No Data"
                    self.session.commit()
                except IOError:
                    self.clean_up()
                    break

        print "DONE"

    def clean_up(self):
        print "Cleaning Up"
        self.session.commit()
        self.state = DIE


def tester(queue):
    queue.put(Experiment())
    print("Sent Data")
