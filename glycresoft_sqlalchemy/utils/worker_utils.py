import os
import uuid
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

try:
    from profilehooks import profile
except:
    # Stub in profiling decorator when not installed.
    def profile(f, *args, **kwargs):
        return f
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

    job_completion_counter = 0

    while len(work_queue) > 0 or work_left:
        next_round = []
        for i in range(len(work_queue)):
            task = work_queue[i]
            task.wait(0.1)
            if task.ready():
                try:
                    if task.ready():
                        result_callback(task.get())
                        job_completion_counter += 1
                except Exception, e:
                    logger.exception("An error occurred", exc_info=e)
            else:
                next_round.append(task)
        work_queue = next_round
        elapsed = (time.time() - timer) > update_window
        if last_length != len(work_queue) or elapsed:
            last_length = len(work_queue)
            logger.info("%d tasks remaining (%d completed)", last_length, job_completion_counter)
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


def profile_task(f, directory=".", **kwargs):
    try:
        name = f.func_name
    except:
        try:
            name = f.func.func_name
        except:
            name = ""

    def produce_new_name():
        return "%s-%s" % (name, str(uuid.uuid4()))

    def profiled_fn(*args, **kwargs):
        fname = os.path.join(directory, produce_new_name())
        wrapped = profile(f, filename=fname)(*args, **kwargs)
        return wrapped

    return profiled_fn
