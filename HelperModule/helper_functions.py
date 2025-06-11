import datetime
import logging
import os
import shutil
import sys
from pathlib import Path
from zipfile import ZipFile

from gemmi import cif

from HelperModule.Ring import Ring


def are_bonds_correct(atom_names, bonds, ring: Ring):
    names_set = set(atom_names)
    metal_atoms = {"FE": 0, "MN": 0, "CO": 0, "RU": 0, "TI": 0, "ZR": 0, "NI": 0, "CR": 0, "RH": 0}
    count = 0
    max_count = ring.atom_number
    double_count = 0
    for bond in bonds:
        atom_1, atom_2 = bond[0].strip('"'), bond[1].strip('"')
                
        if atom_1 in metal_atoms.keys() and atom_2 in names_set:
            metal_atoms[atom_1] += 1
            
        elif atom_2 in metal_atoms.keys() and atom_1 in names_set:
            metal_atoms[atom_2] += 1
            
        elif (atom_1 not in names_set) or (atom_2 not in names_set):
            continue
        else:
            count += 1
        
            
        if ring is Ring.CYCLOPENTANE and any(v == 5 for v in metal_atoms.values()):
            return False
            
            
        if ring is Ring.BENZENE:
            if bond[2] == 'DOUB':
                double_count += 1
            if double_count == 3:
                return True

        # for cyclohexanes/cyclopentanes
        else:
            if bond[2] != 'SING':
                return False
            if count == max_count:
                return True
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
    for path in os.environ.get('PATH', '').split(os.pathsep):
        mono_executable = os.path.join(path, 'mono')
        if os.path.exists(mono_executable):
            return True
    logging.error(f"The Mono package is not installed.")
    return False


def read_component_dictionary(path_to_comp_dict: Path) -> cif.Document:
    logging.info('Reading components dictionary...')
    try:
        document = cif.read(str(path_to_comp_dict))
        return document
    except FileNotFoundError:
        logging.error(f"File {str(path_to_comp_dict)} not found. Please check the file path.")
        sys.exit('Exiting...')
    except PermissionError:
        logging.error(f"Permission denied. Make sure you have the necessary permissions "
                      f"to access the file {str(path_to_comp_dict)}.")
        sys.exit('Exiting...')
    except Exception as e:
        logging.error(f"An error occurred while trying to read {str(path_to_comp_dict)}: {e}")
        sys.exit('Exiting...')


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


