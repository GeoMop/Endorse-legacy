from functools import wraps
import logging
import time

def report(fn):
    @wraps(fn)
    def do_report(*args, **kwargs):
        init_time = time.perf_counter()
        logging.info(f"{fn.__module__}.{fn.__name__} ... ")
        result = fn(*args, **kwargs)
        duration = time.perf_counter() - init_time
        logging.info(f"        ... {duration}")
        return result
    return do_report
