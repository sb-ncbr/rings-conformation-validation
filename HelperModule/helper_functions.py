import logging
import os
import sys
from pathlib import Path
from zipfile import ZipFile
from typing import Set
from gemmi import cif
from HelperModule.Ring import Ring


def is_oxane(current_ring_df, all_single, filepath):
    if not all_single:
        return False
    
    aromatic_count = (current_ring_df["aromatic"].str.upper() == "Y").sum()
    if aromatic_count > 0:
        logging.warning(f"[INCONSISTEND LABELLING oxane]: Skipping a ring with all bonds labelled as single AND at {aromatic_count} bond(s) labelled as aromatic: {filepath}")
        return False
    return True


def is_oxolane(current_ring_df, all_single, filepath):
    if not all_single:
        return False
    
    aromatic_count = (current_ring_df["aromatic"].str.upper() == "Y").sum()
    if aromatic_count > 0:
        logging.warning(f"[INCONSISTEND LABELLING oxolane]: Skipping a ring with all bonds labelled as single AND {aromatic_count} bond(s) labelled as aromatic: {filepath}")
        return False
    return True


def is_cyclopentane(current_ring_df, bond_df, atom_names, metal_atoms, all_single, filepath):
    if not all_single:
        return False
    
    aromatic_count = (current_ring_df["aromatic"].str.upper() == "Y").sum()
    if aromatic_count > 0:
        logging.warning(f"[INCONSISTEND LABELLING cyclopentane]: Skipping a ring with all bonds labelled as single AND {aromatic_count} bond(s) labelled as aromatic: {filepath}")
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
    return len(metal_bonds) != Ring.CYCLOPENTANE.atom_number


def is_cyclohexane(current_ring_df, all_single, filepath):
    if not all_single:
        return False
    
    aromatic_count = (current_ring_df["aromatic"].str.upper() == "Y").sum()
    # condition <aromatic_count> == 6 is checked in is_benzene function
    if 0 < aromatic_count < 6:
        logging.warning(f"[INCONSISTEND LABELLING cyclohexane]: Skipping a ring with all bonds labelled as single AND {aromatic_count} bond(s) labelled as aromatic: {filepath}")
        return False
    # double check
    if aromatic_count == 0:
        return True


def is_benzene(current_ring_df, all_single, filepath):
    double_bonds_count = (current_ring_df["value_order"].str.upper() == "DOUB").sum()

    if (current_ring_df["aromatic"].str.upper() == "Y").all():
        if all_single:
            logging.warning(f"[INCONSISTEND LABELLING benzene]: Skipping a ring with all bonds labelled as aromatic AND single: {filepath}")
            return False
        
        if double_bonds_count == 3:
            logging.info(f"Benzene | {filepath} has THREE bonds labelled as double and all aromatic")
        elif double_bonds_count == 2:
            logging.info(f"Benzene | {filepath} has TWO bonds labelled as double and all aromatic")
        else:
            logging.warning(f"[FOR INVESTIGATION benzene]: Skipping a ring with {double_bonds_count} bonds labelled as double and all aromatic: {filepath}")
            return False
    
        return True
    
    if (current_ring_df["aromatic"].str.upper() == "N").all():
        
        if double_bonds_count == 3:
            logging.info(f"Possible Benzene included | {filepath} has THREE bonds labelled as double BUT all bonds as NOT aromatic")
            return True
        
    if (~current_ring_df["aromatic"].str.upper().isin(["N", "Y"])).all():
        if double_bonds_count == 3:
            logging.info(f"Possible Benzene included | {filepath} has THREE bonds labelled as double BUT all bonds are WITHOUT or INVALID aromatic flag")
            return True

    return False
    

def classify_ring(
    atoms_shape: str, atom_names: Set[str], bond_df, atom_df, filepath, ligand
) -> Ring | None:
    atom_map = atom_df.set_index("atom_id")["type_symbol"]

    bond_df["type_atom_1"] = bond_df["atom_id_1"].map(atom_map)
    bond_df["type_atom_2"] = bond_df["atom_id_2"].map(atom_map)

    cols = ["atom_id_1", "atom_id_2"]
    bond_df[cols] = bond_df[cols].apply(lambda col: col.str.strip('"'))

    mask = bond_df["atom_id_1"].isin(atom_names) & bond_df["atom_id_2"].isin(atom_names)
    current_ring_df = bond_df[mask]
    all_single = (current_ring_df["value_order"].str.upper() == "SING").all()
    match atoms_shape:
        case "C*6":
            if is_benzene(current_ring_df, all_single, filepath):
                return Ring.BENZENE
            if is_cyclohexane(current_ring_df, all_single, filepath):
                return Ring.CYCLOHEXANE
            return None
        case "C*5":
            metal_atoms = ["FE", "MN", "CO", "RU", "TI", "ZR", "NI", "CR", "RH", "IR", "RE", "OS"]
            if is_cyclopentane(current_ring_df, bond_df, atom_names, metal_atoms, all_single, filepath):
                return Ring.CYCLOPENTANE
            return None
        case "C*5-O*1":
            if is_oxane(current_ring_df, all_single, filepath):
                return Ring.OXANE
            return None
        case "C*4-O*1":
            if is_oxolane(current_ring_df, all_single, filepath):
                return Ring.OXOLANE
            return None
        case _:
            return None


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
    logging.error("The Mono package is not installed.")
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
