from argparse import ArgumentParser
from HelperModule.Ring import Ring
from HelperModule.helper_functions import unzip_file, is_mono_installed
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
                logging.warning(f"The ligand name is longer than allowed three charachters. Skipping {ligand_name}...")
                continue

            if str(ligand_name) in ('PHE', 'TYR', 'TRP'):
                continue

            for ring in Ring:
                if is_target_ring_in_name(ring, all_names):
                    extracted_names[ring].append(ligand_name)

        except Exception:
            logging.debug(f'Error while extracting ligand names from block with name {compound_name}. Skipping...')
            continue

    if not all(extracted_names.values()):
        logging.info('No rings found. Exiting...')
        sys.exit(1)

    return extracted_names


def create_config_for_pq(path_to_pdb_local: str, ligands_dict: Dict[Ring, List[str]]) -> None:
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
        with open(PQ_CONFIG, "w") as outfile:
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


def start_program(results_folder, pq_cmd):

    commands = {
        'posix': ['mono', pq_cmd],
        'nt': [pq_cmd]
    }

    command = commands.get(os.name, [])
    command.extend([results_folder, PQ_CONFIG])

    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stderr:
            logging.error(stderr.decode())

    except subprocess.SubprocessError as e:
        logging.error(f"Error executing command: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


def main(path_to_pdb_local: str, output_path: str):

    if not os.path.exists(path_to_pdb_local):
        sys.exit(f"The directory {os.path.abspath(path_to_pdb_local)} does not exist.")

    if not os.path.exists(PQ_CMD):
        sys.exit(f"The file {os.path.abspath(PQ_CMD)} was not found")

    if os.name == 'posix' and not is_mono_installed():
        sys.exit("Mono is not installed")

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    main_workflow_output_dir = os.path.join(output_path, MAIN_DIR)
    try:
        os.makedirs(main_workflow_output_dir)
    except FileExistsError:
        sys.exit(f"The directory {os.path.abspath(main_workflow_output_dir)} already exists.")
    except PermissionError:
        sys.exit(
            f"You do not have the necessary permissions to create the directory"
            f"{os.path.abspath(main_workflow_output_dir)}.")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

    logging.info('Reading components dictionary...')
    try:
        document = cif.read(DEFAULT_DICT_NAME)
    except FileNotFoundError:
        sys.exit(f"The file {os.path.abspath(DEFAULT_DICT_NAME)} does not exist.")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

    logging.info("Extracting ligand names...")
    ligands_dict = extract_ligand_names(document)

    # Create text files with ligand names
    #
    # for ring in ligands_dict.keys():
    #     with open(f"{ring.name}_ligands.txt", 'w') as output:
    #         for ligand in ligands_dict[ring]:
    #             output.write(ligand + '\n')

    logging.info("Creating configuration file for Pattern Query...")
    create_config_for_pq(path_to_pdb_local, ligands_dict)

    logging.info("Running Pattern Query...")
    start_program(MAIN_DIR, pq_cmd=PQ_CMD)

    logging.info('Unzipping the results from Pattern Query...')
    path_to_results = os.path.join(MAIN_DIR, 'result')
    unzip_file(os.path.join(path_to_results, 'result.zip'), path_to_results)


if __name__ == "__main__":
    parser = ArgumentParser(description="Get the dataset of rings using the PatternQuery")
    required = parser.add_argument_group('required named arguments')

    required.add_argument('-d', '--pdb_local', type=str, required=True,
                          help='Path to the local PDB')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='Path to the output directory')

    args = parser.parse_args()
    main(args.pdb_local, args.output)
