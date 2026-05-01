import json
import logging
import subprocess
import os
import sys
import re
from argparse import ArgumentParser
from pathlib import Path
from multiprocessing import cpu_count
from HelperModule.helper_functions import (unzip_file, is_mono_installed,
                                           is_valid_directory, file_exists)
from HelperModule.constants import DEFAULT_DICT_NAME, MAIN_DIR

CPU_COUNT = cpu_count()
PQ_CONFIG = "config.json"
PDB_DIR = "pdb_copy_local"
PQ_CMD = Path("PatternQuery_1.1.25.8.19") / "WebChemistry.Queries.Service.exe"


def create_config_for_pq(path_to_main_output: Path, path_to_pdb_local: str) -> None:
    logging.info("Creating configuration file for Pattern Query...")
    config = {
        "InputFolders": [path_to_pdb_local],
        "Queries": [{"Id": "RingsInHetResidues",
                    "QueryString": "Rings().Inside(HetResidues())"}],
        "StatisticsOnly": False,
        "MaxParallelism": CPU_COUNT
    }

    try:
        with open(path_to_main_output / PQ_CONFIG, "w") as outfile:
            json.dump(config, outfile)
        logging.info("Configuration file successfully created.")
    except OSError as e:
        logging.error(f"Error writing to {PQ_CONFIG}: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        sys.exit(1)


def start_program(results_folder: Path, pq_cmd):
    commands = {
        'posix': ['mono', pq_cmd],
        'nt': [pq_cmd]
    }

    command = commands.get(os.name, [])
    command.extend([results_folder, str(results_folder / PQ_CONFIG)])

    logging.info(f"Running Pattern Query on CPU count: {CPU_COUNT}...")

    pq_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, universal_newlines=True)

    for line in pq_process.stdout:
        error_pattern = r"^\[.*?\] Error:"
        if re.match(error_pattern, line):
            logging.error(f"Error while running Pattern Query {line}")
            sys.exit(1)
        print(line, end='')


def prerequisites_are_met(input_dir: str, output_dir: str) -> bool:

    input_path = Path(input_dir).resolve()
    if not is_valid_directory(input_dir):
        return False

    #TODO 
    # if not is_valid_directory(input_path / CCP4_DIR):
    #     return False

    if not is_valid_directory(input_path / PDB_DIR):
        return False

    if not file_exists(input_path / DEFAULT_DICT_NAME):
        return False

    #TODO generate this file form data.csv
    # if not file_exists(input_path / PDB_INFO_FILE):
    #     return False

    if not file_exists(PQ_CMD):
        return False

    if os.name == 'posix' and not is_mono_installed():
        return False

    output_path = Path(output_dir).resolve()
    main_workflow_output_dir = output_path / MAIN_DIR
    try:
        os.makedirs(main_workflow_output_dir)
    except FileExistsError:
        logging.error(f"The directory {main_workflow_output_dir} already exists.")
        return False
    except PermissionError:
        logging.error(f"The directory {main_workflow_output_dir} is not writable.")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        return False

    return True


def get_results(src: Path, dst: Path):
    # Unzipping the results from Pattern Query
    logging.info('Unzipping the results from Pattern Query...')

    try:
        unzip_file(src, dst)
    except FileNotFoundError as e:
        logging.error(f"File with PQ results was not found: {e}.")
        sys.exit(1)
    except Exception as e:
        logging.error(str(e))
        sys.exit(1)

#TODO use later with onedata
# def unzip_all(path_to_archives: Path) -> None:
#     lst = path_to_archives.glob('*.zip')
#     for zip_ in lst:
#         try:
#             unzip_file(zip_, path_to_archives)
#         except Exception as e:
#             logging.error(str(e))
#             sys.exit(1)


# def preprocess_data(data_path: Path) -> None:
#     unzip_all(data_path / PDB)
#     unzip_all(data_path / CCP4_DIR)


def main(input_path: str, output_path: str):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    logging.info('Starting PrepareDataset...')
    if not prerequisites_are_met(input_path, output_path):
        sys.exit(1)

    #TODO
    # preprocess_data(Path(input_path).resolve())

    path_to_local_pdb = Path(input_path).resolve() / PDB_DIR
    main_workflow_output_dir = Path(output_path).resolve() / MAIN_DIR

    create_config_for_pq(main_workflow_output_dir, str(path_to_local_pdb))

    start_program(main_workflow_output_dir, pq_cmd=PQ_CMD)

    get_results(main_workflow_output_dir / 'result' / 'result.zip', main_workflow_output_dir / 'result')

    logging.info('PrepareDataset has completed successfully')


if __name__ == "__main__":
    parser = ArgumentParser(description="Get the dataset of rings using the PatternQuery")
    required = parser.add_argument_group('required named arguments')

    required.add_argument('-i', '--input', type=str, required=True,
                          help='Path to the directory with input data (local pdb, ccp4 files, etc.)')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='Path to the output directory')

    args = parser.parse_args()
    main(args.input, args.output)