import logging
import os
import sys
from pathlib import Path
from zipfile import ZipFile
from typing import Set
from gemmi import cif
from HelperModule.Ring import Ring


def are_bonds_correct(
    atom_names: Set[str], bond_df, atom_df, ring: Ring, filepath, ligand
):
    atom_map = atom_df.set_index("atom_id")["type_symbol"]

    bond_df["type_atom_1"] = bond_df["atom_id_1"].map(atom_map)
    bond_df["type_atom_2"] = bond_df["atom_id_2"].map(atom_map)

    cols = ["atom_id_1", "atom_id_2"]
    bond_df[cols] = bond_df[cols].apply(lambda col: col.str.strip('"'))

    mask = bond_df["atom_id_1"].isin(atom_names) & bond_df["atom_id_2"].isin(atom_names)
    current_ring_df = bond_df[mask]

    metal_atoms = ["FE", "MN", "CO", "RU", "TI", "ZR", "NI", "CR", "RH"]

    if ring is Ring.BENZENE:
        return (current_ring_df["aromatic"] == "Y").all()
    if ring in (Ring.CYCLOHEXANE, Ring.OXANE, Ring.OXOLANE):
        return (current_ring_df["value_order"] == "SING").all()
    if ring is Ring.CYCLOPENTANE:
        if not (current_ring_df["value_order"] == "SING").all():
            return False
        metal_bonds = bond_df[
            (
                bond_df["atom_id_1"].isin(atom_names)
                & bond_df["type_atom_2"].isin(metal_atoms)
            )
            | (
                bond_df["atom_id_2"].isin(atom_names)
                & bond_df["type_atom_1"].isin(metal_atoms)
            )
        ]

        # TODO: simplify after debugging/analysis
        if metal_bonds.empty:
            return True
        if len(metal_bonds) > 1:
            logging.debug(
                f"More that one metal atom is connected to cyclopentane {filepath}"
            )
        return len(metal_bonds) != ring.atom_number

    return False


def unzip_file(src: Path, dst: Path) -> None:
    try:
        if not src.exists():
            raise FileNotFoundError(f"Source file for unzipping not found: {str(src)}")

        with ZipFile(src, "r") as zip_obj:
            zip_obj.extractall(dst)
        os.remove(src)
    except Exception as e:
        logging.error(f"An error occurred during extraction: {e}")


def is_mono_installed():
    # Check if 'mono' executable exists in any of the directories in the PATH environment variable
    for path in os.environ.get("PATH", "").split(os.pathsep):
        mono_executable = os.path.join(path, "mono")
        if os.path.exists(mono_executable):
            return True
    logging.error(f"The Mono package is not installed.")
    return False


def read_component_dictionary(path_to_comp_dict: Path) -> cif.Document:
    logging.info("Reading components dictionary...")
    try:
        document = cif.read(str(path_to_comp_dict))
        return document
    except FileNotFoundError:
        logging.error(
            f"File {str(path_to_comp_dict)} not found. Please check the file path."
        )
        sys.exit("Exiting...")
    except PermissionError:
        logging.error(
            f"Permission denied. Make sure you have the necessary permissions "
            f"to access the file {str(path_to_comp_dict)}."
        )
        sys.exit("Exiting...")
    except Exception as e:
        logging.error(
            f"An error occurred while trying to read {str(path_to_comp_dict)}: {e}"
        )
        sys.exit("Exiting...")


def is_valid_directory(directory: str | Path) -> bool:
    directory_path = Path(directory).resolve()
    if not directory_path.exists():
        logging.error(f"The directory {str(directory_path)} does not exist.")
        return False
    return True


def file_exists(input_file: str | Path) -> bool:
    input_path = Path(input_file).resolve()
    if not input_path.exists():
        logging.error(f"The file {str(input_path)} was not found.")
        return False
    return True
