import json
import logging
import subprocess
import os
import shutil
import sys
import re
from pathlib import Path
from multiprocessing import cpu_count
from typing import List
from HelperModule.helper_functions import unzip_file


CPU_COUNT = cpu_count()
PQ_CONFIG = "config.json"
PQ_CMD = Path("pqnext") / "WebChemistry.Queries.Service"

logger = logging.getLogger(__name__)


def create_config_for_pq(path_to_main_output: Path, path_to_pdb_local: str) -> None:
    logger.info("Creating configuration file for Pattern Query...")
    config = {
        "InputFolders": [path_to_pdb_local],
        "Queries": [{"Id": "RingsInHetResidues",
                    "QueryString": "Rings().Inside(HetResidues())"}],
        "StatisticsOnly": False,
        "MaxParallelism": CPU_COUNT,
        "OutputMMCIF": True
    }

    try:
        with open(path_to_main_output / PQ_CONFIG, "w") as outfile:
            json.dump(config, outfile)
        logger.info("Configuration file successfully created.")
    except OSError as e:
        logger.error(f"Error writing to {PQ_CONFIG}: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        sys.exit(1)


def start_program(results_folder: Path, pq_cmd):
    command = [
        pq_cmd,
        results_folder,
        str(results_folder / PQ_CONFIG)
    ]

    logger.info(f"Running Pattern Query on CPU count: {CPU_COUNT}...")

    pq_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, universal_newlines=True)

    for line in pq_process.stdout:
        error_pattern = r"^\[.*?\] Error:"
        if re.match(error_pattern, line):
            logger.error(f"Error while running Pattern Query {line}")
            sys.exit(1)
        print(line, end='')


def get_results(src: Path, dst: Path):
    # Unzipping the results from Pattern Query
    logger.info('Unzipping the results from Pattern Query...')

    try:
        unzip_file(src, dst)
    except FileNotFoundError as e:
        logger.error(f"File with PQ results was not found: {e}.")
        sys.exit(1)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)


def extract_rings(pdb_dir: Path, main_dir: Path):

    create_config_for_pq(main_dir, str(pdb_dir))
    start_program(main_dir, pq_cmd=PQ_CMD)
    get_results(main_dir / 'result' / 'result.zip', main_dir / 'result')


