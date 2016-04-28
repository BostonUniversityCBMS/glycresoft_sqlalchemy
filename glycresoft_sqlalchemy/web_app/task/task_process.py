import logging
import traceback
from os import path
from uuid import uuid4
from multiprocessing import Process, Pipe
from threading import Event, Thread, RLock
from Queue import Queue, Empty as QueueEmptyException


logger = logging.getLogger("task_process")
logger.setLevel("ERROR")

NEW = intern('new')
RUNNING = intern('running')
ERROR = intern('error')
FINISHED = intern('finished')


def noop():
    pass


def printop(*args, **kwargs):
    print(args, kwargs)


def configure_log(log_file_path, callable, args):
    import logging
    logger = logging.getLogger()
    handler = logging.FileHandler(log_file_path)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
        "%H:%M:%S")
    handler.setFormatter(formatter)
    logging.captureWarnings(True)
    warner = logging.getLogger('py.warnings')
    warner.setLevel("CRITICAL")

    logger.handlers = []
    logger.addHandler(handler)
    logger.setLevel("DEBUG")
    logger.propagate = False
    return callable(*args)


class CallInterval(object):
    """Call a function every `interval` seconds from
    a separate thread.

    Attributes
    ----------
    stopped: threading.Event
        A semaphore lock that controls when to run `call_target`
    call_target: callable
        The thing to call every `interval` seconds
    args: iterable
        Arguments for `call_target`
    interval: number
        Time between calls to `call_target`
    """

    def __init__(self, interval, call_target, *args):
        self.stopped = Event()
        self.interval = interval
        self.call_target = call_target
        self.args = args
        self.thread = Thread(target=self.mainloop)
        self.thread.daemon = True

    def mainloop(self):
        while not self.stopped.wait(self.interval):
            try:
                self.call_target(*self.args)
            except Exception, e:
                logger.exception("An error occurred in %r", self, exc_info=e)

    def start(self):
        self.thread.start()

    def stop(self):
        self.stopped.set()


class Task(object):
    """
    Represents a separate process that is performing an long running operation against
    the database with a distinct endpoint. This process is executing a series of function
    launched from within Python as opposed to simply running a 3rd-party executable, and is
    able to set up two-way communication between the main process and the "task" process.


    Attributes
    ----------
    args : iterable
        The list of arguments to the task's main function
    callback : callable
        A function to be called when the task's main function completes
    id : str
        A unique identifier for the task. Usually a UUID
    log_file_path : str
        The path for the task write all logging messages to
    message_buffer : list
        An accumulator for outbound messages from this task
    name : str
        A more human-readable description of this "task" process which may
        be displayed to the user. Defaults to :attr:`id`.
    process : multiprocessing.Process
        The actual Process object doing the work
    state : str
        One of several constants describing whether task is new, has started,
        has finished with an error, or has completed successfully
    task_fn : callable
        The function to call in the "task" process
    """
    def __init__(self, task_fn, args, callback=printop, **kwargs):
        self.id = str(uuid4())
        self.task_fn = task_fn
        self.pipe, child_conn = Pipe(True)
        self.state = NEW
        self.process = None
        self.args = list(args)
        self.args.append(child_conn)
        self.callback = callback
        self.log_file_path = kwargs.get("log_file_path", "%s.log" % self.id)
        self.name = kwargs.get('name', self.id)
        self.message_buffer = []

    def start(self):
        self.process = Process(target=configure_log, args=(self.log_file_path, self.task_fn, self.args))
        self.state = RUNNING
        self.process.start()

    def get_message(self):
        if len(self.message_buffer) > 0:
            message = self.message_buffer.pop(0)
            message.source = self
            return message
        if self.pipe.poll():
            message = self.pipe.recv()
            message.source = self
            return message
        return None

    def add_message(self, message):
        self.message_buffer.append(message)

    def update(self):
        if self.process is not None:
            result = self.process.is_alive()
            if not result and self.state == RUNNING:
                exitcode = self.process.exitcode
                if exitcode == 0:
                    self.state = FINISHED
                    self.on_complete()
                elif exitcode is None:
                    self.state = ERROR
                else:
                    self.state = ERROR
            return result

        return False

    def __getstate__(self):
        return {
            "id": self.id,
            "state": self.state,
            "task_fn": self.task_fn,
            "args": self.args[:-1],
            "callback": self.callback
        }

    def __setstate__(self, state):
        self.id = state['id']
        self.task_fn = state['task_fn']
        self.state = state['state']
        self.args = state['args']
        self.callback = state['callback']
        self.pipe, child_conn = Pipe(True)
        self.args.append(child_conn)
        self.process = None

    def to_json(self):
        return dict(id=self.id, name=self.name, status=self.state)

    def messages(self):
        message = self.get_message()
        while message is not None:
            yield message
            message = self.get_message()

    def __repr__(self):
        return "<Task {} {}>".format(self.id, self.state)

    def on_complete(self):
        self.callback()


class NullPipe(object):
    """
    A class to stub out multiprocessing.Pipe
    """
    def send(self, *args, **kwargs):
        logger.info(*args, **kwargs)

    def recv(self):
        return ""

    def poll(self, timeout=None):
        return False


class Message(object):
    """
    Represent a message sent between a Task and the main process.

    Attributes
    ----------
    message : object
        Anything, but preveriably something JSON serializeable
    source : object
        Description
    type : str
        A constant similar to logging levels. Options in use include
        ("info", "error", "update")
    """
    def __init__(self, message, type="info", source=None):
        self.message = message
        self.source = source
        self.type = type

    def __str__(self):
        return "%s:%s - %r" % (self.source, self.type, self.message)

    @classmethod
    def traceback(cls):
        return cls(traceback.format_exc(), "error")


class TaskManager(object):
    """Track and schedule `Task` objects and associated processes.

    Attributes
    ----------
    task_dir: str
        file system directory path for writing task-specific information
    tasks: dict
        A task id -> `Task` object mapping for all tasks, running or otherwise
    task_queue: Queue.Queue
        A `Queue.Queue` for holding `Task` objects currently waiting to be ran
    currently_running: dict
        A task id -> `Task` object mapping for all tasks currently running
    n_running: int
        The number of tasks currently running
    max_running: int
        The maximum number of tasks allowed to run at once
    timer: CallInterval
        A `CallInterval` object who schedules :meth:`TaskManager.tick`
    messages: Queue
        A `Queue` for holding task messages to be read by clients
    """
    interval = 5

    def __init__(self, task_dir=None, max_running=1):
        if task_dir is None:
            task_dir = "./"
        self.task_dir = task_dir
        self.tasks = {}
        # self.read_tasks()
        self.n_running = 0
        self.max_running = max_running
        self.task_queue = Queue()
        self.currently_running = {}
        self.completed_tasks = set()
        self.timer = CallInterval(self.interval, self.tick)
        self.timer.start()
        self.messages = Queue()
        self.running_lock = RLock()
        self.halting = False

    def add_task(self, task):
        """Add a `Task` object to the set of all tasks being managed
        by this instance.

        Once a task is added, it will be checked during each `tick`

        Parameters
        ----------
        task : Task
            The task to be scheduled
        """
        self.tasks[task.id] = task
        self.messages.put(Message({"id": task.id, "name": task.name}, "task-queued"))

    def get_task_log_path(self, task):
        return path.join(getattr(self, "task_dir", ""), task.id + '.log')

    def startloop(self):
        self.timer.start()

    def stoploop(self):
        self.timer.stop()

    def terminate(self):
        self.stoploop()
        # Clean up here

    def tick(self):
        """Check each managed task for status updates, schedule new tasks
        when space becomes available, remove finished tasks and handle errors.

        This method is called every :attr:`TaskManager.interval` seconds.
        """
        try:
            self.check_state()
            self.launch_new_tasks()
        except Exception, e:
            logger.exception("an error occurred in `tick`", exc_info=e)

    def check_state(self):
        """Iterate over all tasks in :attr:`TaskManager.tasks` and check their status.
        """
        logger.debug(
            "Checking task manager state:\n %d tasks running\nRunning: %r\n%r",
            self.n_running, self.currently_running, self.tasks)
        for task_id, task in list(self.tasks.items()):
            running = task.update()
            logger.debug("Checking %r", task)
            for message in task.messages():
                self.messages.put(message)

            if task.state == NEW:
                self.task_queue.put(task)
            elif task.state == FINISHED:
                if self.running_lock.acquire(0):
                    self.currently_running.pop(task.id)
                    self.tasks.pop(task.id)
                    self.n_running -= 1
                    task.callback()
                    self.messages.put(Message({"id": task.id, "name": task.name}, "task-complete"))
                    self.running_lock.release()
                    self.completed_tasks.add(task.id)
            elif task.state == ERROR:
                if task.id in self.currently_running:
                    if self.running_lock.acquire(0):
                        self.currently_running.pop(task.id)
                        self.messages.put(Message({"id": task.id, "name": task.name}, "task-error"))
                        self.tasks.pop(task.id)
                        self.n_running -= 1
                        self.running_lock.release()
            elif task.state == RUNNING:
                if running:
                    continue
                else:
                    print task.id, running

        # for task_id, task in list(self.currently_running.items()):
        #     running = task.update()
        #     logger.debug("Checking %r", task)
        #     for message in task.messages():
        #         self.messages.put(message)

        #     if task.state == FINISHED:
        #         self.currently_running.pop(task.id)
        #         self.completed_tasks.add(task.id)
        #         try:
        #             self.tasks.pop(task.id)
        #         except KeyError:
        #             pass
        #         self.n_running -= 1
        #         task.callback()
        #         self.messages.put(Message({"id": task.id, "name": task.name}, "task-complete"))
        #     elif not running:
        #         print task.id, "not running", task.state

    def launch_new_tasks(self):
        while((self.n_running < self.max_running) and (self.task_queue.qsize() > 0)):
            try:
                task = self.task_queue.get(False)
                if task.state != NEW or task.id in self.completed_tasks:
                    continue
                self.run_task(task)
            except QueueEmptyException:
                break

    def run_task(self, task):
        """Start running a task if space is available

        Parameters
        ----------
        task : Task
        """
        if self.running_lock.acquire(0):
            self.currently_running[task.id] = task
            self.running_lock.release()
            task.log_file_path = self.get_task_log_path(task)
            self.n_running += 1
            task.start()
            self.add_message(Message({"id": task.id, "name": task.name}, 'task-start'))

    def add_message(self, message):
        self.messages.put(message)

    def task_list(self):
        tasks = []
        for t_id, task in self.tasks.items():
            tasks.append({
                "name": task.name,
                "id": task.id,
                "status": task.state
            })
        return tasks
