import time
import logging
try:
    logger = logging.getLogger("worker_utils")
except:
    pass

import multiprocessing

try:
    range = xrange
except:
    pass


def async_worker_pool(worker_pool, work_stream, task_fn, result_callback=None,
                      logger=logger, update_window=30, initial_load=1500,
                      maxload=10000):
    if result_callback is None:
        result_callback = ResultCounter()
    logger.info("async_worker_pool starting")
    work_queue = []
    work_loader = (work_queue.append(worker_pool.apply_async(task_fn, [item])) for item in work_stream)
    work_left = True
    for i in range(initial_load):
        try:
            work_loader.next()
        except StopIteration:
            work_left = False

    last_length = len(work_queue)
    begin_time = timer = time.time()
    elapsed = False

    while len(work_queue) > 0 or work_left:
        next_round = []
        for i in range(len(work_queue)):
            task = work_queue[i]
            task.wait(0.1)
            if task.ready():
                try:
                    result_callback(task.get())
                except Exception, e:
                    logger.exception("An error occurred", exc_info=e)
            else:
                next_round.append(task)
        work_queue = next_round
        elapsed = time.time() - timer > update_window
        if last_length != len(work_queue) or elapsed:
            last_length = len(work_queue)
            logger.info("%d tasks remaining", last_length)
            timer = time.time()
        if len(work_queue) < maxload:
            for i in range(initial_load):
                try:
                    work_loader.next()
                except StopIteration:
                    work_left = False
    end_time = time.time()
    logger.info("Time elapsed: %0.2fs", (end_time - begin_time))


class ResultCounter(object):
    def __init__(self, n=0, logger=None):
        self.n = n
        self.logger = logger

    def __call__(self, i):
        self.n += i
        if self.logger is not None:
            self.logger.info("%d completed", self.n)
