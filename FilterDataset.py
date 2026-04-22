import logging
import os
import sys
import shutil
import pandas as pd
import pickle
from pathlib import Path
from argparse import ArgumentParser
from gemmi import cif
from HelperModule.Ring import Ring
from HelperModule.getter_functions import (
    get_bonds_from_cif,
    get_atoms_from_pdb,
)
from HelperModule.helper_functions import are_bonds_correct
from HelperModule.constants import *


def process_correct_rings(output_dir, ligand, filepath):
    output_pdb_dir = output_dir / ligand / "patterns"

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    filename = filepath.name
    if not filename.startswith(f"{ligand}_"):
        filename = f"{ligand}_{filename}"

    new_name_path = output_pdb_dir / filename
    shutil.copy(filepath, new_name_path)


def run_filter(
    input_path: Path, ring: Ring, output_dir: Path, document: cif.Document
) -> None:
    patterns_df = pd.read_csv(input_path / "patterns.csv")
    pkl_path = Path(f"{ring.name.lower()}_ligand_ring_map.pkl")
    processed_data_dict = {}

    if pkl_path.is_file():
        with pkl_path.open("rb") as f:
            processed_data_dict = pickle.load(f)

    target_count = 0

    for row in patterns_df.itertuples(index=False):
        ligand = row.Residues.split()[0]

        filepath = input_path / "patterns" / (row.Id + ".pdb")
        atom_names = get_atoms_from_pdb(filepath, ring)

        key = (ligand, frozenset(atom_names))

        if key in processed_data_dict:
            # skipping because already filtered out this ring as wrong
            if not processed_data_dict[key]:
                continue

            process_correct_rings(output_dir, ligand, filepath)
            target_count += 1
            continue

        ligand_block = document.find_block(ligand)
        if ligand_block is None:
            logging.warning(f"Ligand_block is None for {ligand}")
            continue

        atom_bonds = get_bonds_from_cif(ligand_block)
        is_correct = are_bonds_correct(atom_names, atom_bonds, ring)
        processed_data_dict[key] = is_correct

        if is_correct:
            process_correct_rings(output_dir, ligand, filepath)
            target_count += 1

    with pkl_path.open("wb") as f:
        pickle.dump(processed_data_dict, f)

    logging.info(f"[{ring.name.capitalize()}]: {target_count} patterns were found.")


def main(ring: str, output_path: str, input_path: str):
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    ring = ring.upper()
    if ring not in Ring.__members__.keys():
        logging.error(
            f"Ring {ring} is not a valid Ring. Currently supported: {[e.name for e in Ring]} Exiting..."
        )
        sys.exit(1)

    logging.info(f"[{ring.capitalize()}]: Starting FilterDataset...")

    main_workflow_output_dir = Path(output_path) / MAIN_DIR
    if not os.path.exists(main_workflow_output_dir):
        logging.error(
            f"The directory {main_workflow_output_dir.resolve()} does not exist. Exiting..."
        )
        sys.exit(1)

    path_to_comp_dict = Path(input_path) / DEFAULT_DICT_NAME

    if not path_to_comp_dict.exists():
        logging.error(f"The file {path_to_comp_dict} does not exist. Exiting...")
        sys.exit(1)

    # that is the output dir from the previous script (previous step)
    dir_with_patterns = main_workflow_output_dir / "result" / ring.lower()

    if not dir_with_patterns.exists() or not any(dir_with_patterns.iterdir()):
        logging.error(f'The directory "{dir_with_patterns}" does not exist or is empty')
        sys.exit(1)

    current_ring_path = main_workflow_output_dir / ring.lower()
    dir_for_filtered_patterns = current_ring_path / "filtered_ligands"

    logging.info("Reading components dictionary...")
    document = cif.read(str(path_to_comp_dict))
    run_filter(dir_with_patterns, Ring[ring], dir_for_filtered_patterns, document)
    logging.info(f"[{ring.capitalize()}]: FilterDataset has completed successfully")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Filter out the wrong patterns, which do not possess "
        "the required atomic bonds"
    )
    required = parser.add_argument_group("required named arguments")

    required.add_argument(
        "-r",
        "--ring",
        required=True,
        type=str,
        help=f"Choose the target ring type. Currently supported:"
        f" {[e.name for e in Ring]}",
    )
    required.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output directory. Should be the same as in the previous step.",
    )
    required.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to the directory with input data (local pdb, ccp4 files, etc.)",
    )

    args = parser.parse_args()

    main(args.ring, args.output, args.input)
