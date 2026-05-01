from datetime import timedelta
import logging
import os
import sys
import shutil
import time
import pandas as pd
from pathlib import Path
from argparse import ArgumentParser
from gemmi import cif
from HelperModule.Ring import Ring
from HelperModule.getter_functions import (
    get_data_from_cif,
    get_atoms_from_pdb,
)
from HelperModule.helper_functions import are_bonds_correct
from HelperModule.constants import *


def process_correct_rings(output_dir, ligand, filepath):
    output_pdb_dir = output_dir / "filtered_ligands" / ligand / "patterns"

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    filename = filepath.name
    if not filename.startswith(f"{ligand}_"):
        filename = f"{ligand}_{filename}"

    new_name_path = output_pdb_dir / filename
    shutil.copy(filepath, new_name_path)


def filter_by_bond_type(
    patterns_df, dir_with_patterns: Path, ring: Ring, output_dir: Path, document: cif.Document
) -> None:
    processed_data_dict = {}
    target_count = 0
    target_ring_rows = []

    for row in patterns_df.itertuples(index=False):
        ligand = row.Residues.split()[0]

        filepath = dir_with_patterns / "patterns" / (row.Id + ".pdb")

        # only for debug/manual adding of rings, etc
        if not filepath.exists():
            logging.warning(f"{str(filepath)} does not exist.")
            continue

        atom_names = get_atoms_from_pdb(filepath, ring)

        key = (ligand, frozenset(atom_names))

        if key in processed_data_dict:
            # skipping because already filtered out this ring as wrong
            if not processed_data_dict[key]:
                continue

            process_correct_rings(output_dir, ligand, filepath)
            target_ring_rows.append(row._asdict())
            target_count += 1
            continue

        ligand_block = document.find_block(ligand)
        if ligand_block is None:
            logging.warning(f"Ligand_block is None for {ligand}")
            continue

        bond_df, atom_df = get_data_from_cif(ligand_block)
        is_correct = are_bonds_correct(
            atom_names, bond_df, atom_df, ring, filepath, ligand
        )
        processed_data_dict[key] = is_correct

        if is_correct:
            process_correct_rings(output_dir, ligand, filepath)
            target_ring_rows.append(row._asdict())
            target_count += 1


    res_df = pd.DataFrame(target_ring_rows)
    res_df.to_csv(output_dir / f"filtered_patterns_{ring.name.lower()}.csv")
    return res_df


def main(output_path: str, input_path: str):
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    logging.info(f"Starting FilterDataset...")

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
    dir_with_patterns = main_workflow_output_dir / "result" / "RingsInHetResidues"

    if not dir_with_patterns.exists() or not any(dir_with_patterns.iterdir()):
        logging.error(f'The directory "{dir_with_patterns}" does not exist or is empty')
        sys.exit(1)

    logging.info("Creating rings statistics...")
    df = pd.read_csv(dir_with_patterns / "patterns.csv")

    stats_df = (
        df.groupby("Atoms", as_index=False)
        .agg(
            TotalRingCount=("Atoms", "size"),
            UniqueLigands=("Signature", "nunique"),
            UniquePdbs=("ParentId", "nunique"),
        )
        .sort_values("TotalRingCount", ascending=False)
        .reset_index(drop=True)
    )

    true_totals = {
        "TotalRingCount": len(df),
        "UniqueLigands": df["Signature"].nunique(),
        "UniquePdbs": df["ParentId"].nunique(),
    }

    totals_row = pd.DataFrame([true_totals])
    totals_row["Atoms"] = "TOTAL"

    stats_df_with_total = pd.concat([stats_df, totals_row], ignore_index=True)
    csv_stats = main_workflow_output_dir / "all_rings_stats_nofilter.csv"
    stats_df_with_total.to_csv(csv_stats)
    logging.info(f"Done. Exported to {csv_stats}")

    target_atoms = [ring.atoms for ring in Ring]

    grouped = dict(tuple(df.groupby("Atoms")))

    dfs_by_ring = {atom: grouped[atom].copy() for atom in target_atoms if atom in grouped}

    logging.info("Reading components dictionary...")
    document = cif.read(str(path_to_comp_dict))

    for ring in Ring:
        logging.info(f"Processing {ring.name.lower()}...")
        output_path = main_workflow_output_dir / ring.name.lower()
        output_path.mkdir(parents=True, exist_ok=True)

        df = filter_by_bond_type(
            dfs_by_ring[ring.atoms],
            dir_with_patterns,
            ring,
            output_path,
            document
        )

        logging.info(f"{ring.name} | Rings={len(df)} | Unique ligands={df["Signature"].nunique()} | Unique PDBs={df["ParentId"].nunique()}")

    logging.info("FilterDataset has completed successfully.")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Filter out the wrong patterns, which do not possess "
        "the required atomic bonds"
    )
    required = parser.add_argument_group("required named arguments")

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
    start = time.perf_counter()
    
    main(args.output, args.input)

    elapsed = time.perf_counter() - start
    formatted = str(timedelta(seconds=elapsed))
    logging.info(f"Total time: {formatted}")
