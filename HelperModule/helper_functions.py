import datetime
import logging
import os
import shutil
from zipfile import ZipFile
from HelperModule.Ring import Ring


def are_bonds_correct(atom_names, bonds, ring: Ring):
    names_set = set(atom_names)
    metal_atoms = {"FE": 0, "MN": 0, "CO": 0, "RU": 0, "TI": 0, "ZR": 0, "NI": 0, "CR": 0}
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


def unzip_file(src: str, dst: str) -> None:
    try:
        with ZipFile(src, "r") as zip_obj:
            zip_obj.extractall(dst)
        os.remove(src)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


