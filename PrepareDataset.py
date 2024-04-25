import re
from argparse import ArgumentParser
from HelperModule.Ring import Ring
from HelperModule.helper_functions import (unzip_file, is_mono_installed,
                                           read_component_dictionary, is_valid_directory, file_exists)
from HelperModule.constants import *
import logging
from gemmi import cif
from multiprocessing import cpu_count
from typing import Dict, List
import json
import subprocess
import os
import sys

CPU_COUNT = cpu_count()


def is_target_ring_in_name(ring: Ring, all_names: list[str]) -> bool:
    if ring is Ring.BENZENE:

        # we need to check both 'benz' and 'phen'
        substrings = ring.name_substring.split(';')

        for substring in substrings:
            if any(substring in e for e in all_names):
                return True

    return any(ring.name_substring in e for e in all_names)


def extract_ligand_names(doc: cif.Document) -> Dict[Ring, List[str]]:
    logging.info("Extracting ligand names...")
    extracted_names = {ring: [] for ring in Ring}
    for i, ligand_block in enumerate(doc):
        compound_name = ligand_block.find_value('_chem_comp.name')  # chemical name
        try:
            other_names = ligand_block.find(['_pdbx_chem_comp_identifier.identifier'])
            synonyms = ligand_block.find_value('_chem_comp.pdbx_synonyms')

            all_names = [compound_name.lower()] + [list(e)[0].lower() for e in list(other_names)] + \
                        [synonyms.lower()]

            ligand_name = ligand_block.find_value('_chem_comp.id')

            if len(ligand_name) > 3:
                logging.warning(f"The ligand name is longer than allowed three characters. Skipping {ligand_name}...")
                continue

            if str(ligand_name) in ('PHE', 'TYR', 'TRP'):
                continue

            for ring in Ring:
                if is_target_ring_in_name(ring, all_names):
                    extracted_names[ring].append(ligand_name)

        except Exception:
            logging.warning(f'Error while extracting ligand names from block with name {compound_name}. Skipping...')
            continue

    if not all(extracted_names.values()):
        logging.warning('No rings found. Exiting...')
        sys.exit(1)

    return extracted_names


def create_config_for_pq(path_to_main_output: Path, path_to_pdb_local: str, ligands_dict: Dict[Ring, List[str]]) -> None:
    logging.info("Creating configuration file for Pattern Query...")
    config = {
        "InputFolders": [path_to_pdb_local],
        "Queries": [],
        "StatisticsOnly": False,
        "MaxParallelism": CPU_COUNT
    }

    for ring, ligands in ligands_dict.items():
        config["Queries"].append(create_query(ring.name.lower(), ring.pattern_query +
                                              f".Inside(Residues({ligands}))"))

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


def create_query(query_id: str, query_string: str) -> Dict[str, str]:
    return {
        "Id": query_id,
        "QueryString": query_string
    }


def start_program(results_folder: Path, pq_cmd):
    commands = {
        'posix': ['mono', pq_cmd],
        'nt': [pq_cmd]
    }

    command = commands.get(os.name, [])
    command.extend([results_folder, str(results_folder / PQ_CONFIG)])

    logging.info(f"Running Pattern Query on CPU count: {CPU_COUNT}...")

    pq_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = [stream.strip() for stream in pq_process.communicate()]

    error_pattern = r"^\[.*?\] Error:"
    for line in stdout.splitlines():
        if re.match(error_pattern, line):
            logging.error(f"Error while running Pattern Query {line}")
            sys.exit(1)

    if stdout:
        print(stdout)


def prerequisites_are_met(input_dir: str, output_dir: str) -> bool:

    input_path = Path(input_dir).resolve()
    if not is_valid_directory(input_dir):
        return False

    if not is_valid_directory(input_path / CCP4_DIR):
        return False

    if not is_valid_directory(input_path / PDB):
        return False

    if not file_exists(input_path / DEFAULT_DICT_NAME):
        return False

    if not file_exists(input_path / PDB_INFO_FILE):
        return False

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


def main(input_path: str, output_path: str):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    logging.info('Starting PrepareDataset...')
    if not prerequisites_are_met(input_path, output_path):
        sys.exit(1)

    document = read_component_dictionary(Path(input_path).resolve() / DEFAULT_DICT_NAME)

    path_to_local_pdb = Path(input_path).resolve() / PDB
    main_workflow_output_dir = Path(output_path).resolve() / MAIN_DIR

    ligands_dict = extract_ligand_names(document)
    create_config_for_pq(main_workflow_output_dir, str(path_to_local_pdb), ligands_dict)

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
