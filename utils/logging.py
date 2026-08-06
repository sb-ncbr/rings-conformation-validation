import logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir: Path):
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    class IgnoreSHGOFilter(logging.Filter):
        def filter(self, record):
            return not (
                record.pathname.endswith("scipy/optimize/_shgo.py")
                and record.levelno == logging.INFO
            )
    
    shgo_filter = IgnoreSHGOFilter()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )


    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    workflow_handler = logging.FileHandler(
        log_dir / f"workflow_{timestamp}.log"
    )
    workflow_handler.setLevel(logging.INFO)
    workflow_handler.setFormatter(formatter)
    workflow_handler.addFilter(shgo_filter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    console_handler.addFilter(shgo_filter)

    root_logger.addHandler(workflow_handler)
    root_logger.addHandler(console_handler)


    validation_logger = logging.getLogger("validation")
    validation_logger.setLevel(logging.WARNING)
    validation_logger.propagate = False

    validation_handler = logging.FileHandler(
        log_dir / f"validation_{timestamp}.log"
    )
    validation_handler.setLevel(logging.WARNING)
    validation_handler.setFormatter(formatter)

    validation_logger.addHandler(validation_handler)
    validation_logger.addHandler(console_handler)
