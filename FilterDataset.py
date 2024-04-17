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


def run_filter(input_dir: str, ring: Ring, output_dir: str) -> None:
    total = 0
    target_count = 0
    document = cif.read(DEFAULT_DICT_NAME)

    for root, _, files in os.walk(input_dir):
        for file in files:
            if not file.endswith('.pdb'):
                continue
            total += 1
            filepath = os.path.join(root, file)
            ligand, atom_names = get_data_from_pdb(filepath, ring)
            ligand_block = document.find_block(ligand)

            if ligand_block is None:
                logging.warning(f"ligand_block is None for {ligand}")
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

    logging.info(f"{target_count} {ring.name.lower()} patterns were found.")


def main(input_ring: str, output_path: str):

    main_workflow_output_dir = os.path.join(output_path, MAIN_DIR)
    if not os.path.exists(main_workflow_output_dir):
        sys.exit(f"The directory {os.path.abspath(main_workflow_output_dir)} does not exist.")

    if not os.path.exists(DEFAULT_DICT_NAME):
        sys.exit(f"No component dictionary <{DEFAULT_DICT_NAME}> was found in the current directory. Exiting...")

    ring = input_ring.upper()
    if ring not in Ring.__members__.keys():
        logging.error(f"Ring {ring} is not a valid Ring.")
        sys.exit(f'Check ring name! Currently supported: {[e.name for e in Ring]}')

    current_ring_path = os.path.join(main_workflow_output_dir, input_ring)

    input_dir = os.path.join(main_workflow_output_dir, "result", input_ring)

    if not os.path.exists(input_dir) or not os.listdir(input_dir):
        sys.exit(f'The input directory "{input_dir}" does not exist or is empty')

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info(f"Running FilterDataset for the ring {ring.lower()}")

    output_dir = os.path.join(current_ring_path, 'filtered_ligands')
    run_filter(input_dir, Ring[ring], output_dir)


if __name__ == '__main__':
    parser = ArgumentParser(description="Filter out the wrong patterns, which do not possess "
                                        "the required atomic bonds")
    required = parser.add_argument_group('required named arguments')

    required.add_argument('-r', '--ring', required=True, type=str,
                          help=f'Choose the target ring type. Currently supported:'
                               f' {[e.name for e in Ring]}')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='Path to the output directory. Should be the same as in the previous step.')

    args = parser.parse_args()

    main(args.ring, args.output)
