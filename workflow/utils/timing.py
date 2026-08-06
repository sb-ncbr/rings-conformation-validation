from datetime import timedelta
import time
import logging
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class ExecutionTimer:
    def __init__(self):
        self.times = {}

    def add(self, name, seconds):
        self.times[name] = seconds

    def report(self):
        logger.info("Execution time summary:")
        for name, seconds in self.times.items():
            formatted = str(timedelta(seconds=seconds))
            logger.info(
                "  %-30s %s",
                name,
                formatted
            )


@contextmanager
def timed_step(timer, name):
    start = time.perf_counter()

    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        formatted = str(timedelta(seconds=elapsed))
        timer.add(name, elapsed)
        logger.info(
            "Finished %s in %s",
            name,
            formatted
        )