from .task_process import Task, Message
import time


def echo(*args, **kwargs):
    time.sleep(15)
    args[-1].send(Message("Echo.... %s" % [args[:-1]], 'update'))


class DummyTask(Task):
    def __init__(self, *args, **kwargs):
        Task.__init__(self, echo, args, **kwargs)
