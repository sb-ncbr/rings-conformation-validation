from argparse import ArgumentParser
from HelperModule.Ring import Ring
from HelperModule.helper_functions import unzip_file
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

    if not all(extracted_names.values()):
        logging.info('No rings found. Exiting...')
        sys.exit()

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
        
    with open(PQ_CONFIG, "w") as outfile:
        json.dump(config, outfile)


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

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stderr:
        logging.error(stderr.decode())


def main(path_to_pdb_local: str, pq_cmd: str):
    
    try:
        os.makedirs(MAIN_DIR, exist_ok=True)
    except PermissionError:
        sys.exit(f"Error: You do not have the necessary permissions to create the directory {os.path.abspath(MAIN_DIR)}.")
    
    logging.basicConfig(level=logging.DEBUG, 
					    format='%(asctime)s - %(levelname)s - %(message)s') 
    

    logging.info('Reading components dictionary...')
    document = cif.read(DEFAULT_DICT_NAME)
    logging.info("Extracting ligand names...")
    ligands_dict = extract_ligand_names(document)
    logging.info("Creating configuration file for Pattern Query...")
    create_config_for_pq(path_to_pdb_local, ligands_dict)
    logging.info("Running Pattern Query...")
    start_program(MAIN_DIR, pq_cmd=PQ_CMD_DIR)
    
    logging.info('Unzipping the results from Pattern Query...')
    
    path_to_results = os.path.join(MAIN_DIR, 'result')
    unzip_file(os.path.join(path_to_results, 'result.zip'), path_to_results)


if __name__ == "__main__":
    parser = ArgumentParser(description="Get the dataset of rings using the PatternQuery")
    required = parser.add_argument_group('required named arguments')
    
    required.add_argument('-d', '--pdb_local', type=str,
                          help='Path to the local PDB')
    
    args = parser.parse_args()
    main(args.pdb_local)