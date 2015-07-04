from .task_process import Task, Message


def echo(*args, **kwargs):
    args[-1].send(Message("Echo.... %s" % [args[:-1]], 'update'))


class DummyTask(Task):
    def __init__(self, *args, **kwargs):
        Task.__init__(self, echo, args, **kwargs)
