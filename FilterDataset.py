import logging
import os
import sys
import shutil
from argparse import ArgumentParser
from gemmi import cif
from HelperModule.Ring import Ring
from HelperModule.getter_functions import get_data_from_pdb, get_bonds_from_cif
from HelperModule.helper_functions import are_bonds_correct
from HelperModule.constants import *


def run_filter(input_dir: str, ring: Ring, output_dir: str, document: cif.Document) -> None:
    total = 0
    target_count = 0

    for root, _, files in os.walk(input_dir):
        for file in files:
            if not file.endswith('.pdb'):
                continue
            total += 1
            filepath = os.path.join(root, file)
            ligand, atom_names = get_data_from_pdb(Path(filepath), ring)
            ligand_block = document.find_block(ligand)

            if ligand_block is None:
                logging.warning(f"Ligand_block is None for {ligand}")
                continue

            atom_bonds = get_bonds_from_cif(ligand_block)

            if are_bonds_correct(atom_names, atom_bonds, ring):
                output_pdb_dir = os.path.join(output_dir, ligand, 'patterns')
                os.makedirs(output_pdb_dir, exist_ok=True)

                # ligand name can be up to 3 chars
                if os.path.basename(filepath)[len(ligand)] != '_':
                    # assigns a new unique name for the purposes of following analysis
                    new_name_path = os.path.join(output_pdb_dir,
                                                 ligand + '_' + os.path.basename(filepath))
                else:
                    new_name_path = os.path.join(output_pdb_dir, os.path.basename(filepath))
                logging.info(f"Copying {os.path.basename(filepath)}")
                shutil.copy(filepath, new_name_path)
                target_count += 1

    logging.info(f"[{ring.name.capitalize()}]: {target_count} patterns were found.")


def main(ring: str, output_path: str, input_path: str):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s'
                        )

    ring = ring.upper()
    if ring not in Ring.__members__.keys():
        logging.error(f"Ring {ring} is not a valid Ring. Currently supported: {[e.name for e in Ring]} Exiting...")
        sys.exit(1)

    logging.info(f"[{ring.capitalize()}]: Starting FilterDataset...")

    main_workflow_output_dir = os.path.join(output_path, MAIN_DIR)
    if not os.path.exists(main_workflow_output_dir):
        logging.error(f"The directory {os.path.abspath(main_workflow_output_dir)} does not exist. Exiting...")
        sys.exit(1)

    path_to_comp_dict = os.path.join(input_path, DEFAULT_DICT_NAME)
    logging.info('Reading components dictionary...')
    if not os.path.exists(path_to_comp_dict):
        logging.error(f"The file {path_to_comp_dict} does not exist. Exiting...")
        sys.exit(1)

    current_ring_path = os.path.join(main_workflow_output_dir, ring.lower())

    # that is the output dir from the previous script (previous step)
    dir_with_patterns = os.path.join(main_workflow_output_dir, "result", ring.lower())

    if not os.path.exists(dir_with_patterns) or not os.listdir(dir_with_patterns):
        logging.error(f'The directory "{dir_with_patterns}" does not exist or is empty')
        sys.exit(1)

    dir_for_filtered_patterns = os.path.join(current_ring_path, 'filtered_ligands')

    document = cif.read(path_to_comp_dict)
    run_filter(dir_with_patterns, Ring[ring], dir_for_filtered_patterns, document)
    logging.info(f'[{ring.capitalize()}]: FilterDataset has completed successfully')


if __name__ == '__main__':
    parser = ArgumentParser(description="Filter out the wrong patterns, which do not possess "
                                        "the required atomic bonds")
    required = parser.add_argument_group('required named arguments')

    required.add_argument('-r', '--ring', required=True, type=str,
                          help=f'Choose the target ring type. Currently supported:'
                               f' {[e.name for e in Ring]}')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='Path to the output directory. Should be the same as in the previous step.')
    required.add_argument('-i', '--input', type=str, required=True,
                          help='Path to the directory with input data (local pdb, ccp4 files, etc.)')

    args = parser.parse_args()

    main(args.ring, args.output, args.input)
