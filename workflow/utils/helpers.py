import logging
import os
from pathlib import Path
import re
from zipfile import ZipFile
import gemmi
import pandas as pd
from workflow.models.Ring import Ring

logger = logging.getLogger(__name__)


def extract_extended_pdb_code(string: str):
    return string.split('_', maxsplit=1)[1].rsplit('_', maxsplit=1)[0]


def get_old_pdb_id(extended_id: str) -> str | None:
    """
    Extract old 4-character PDB ID from extended PDB ID.

    Example:
        00002xyz -> 2xyz
        abcd1234 -> None
    """
    extended_id = extended_id.lower()

    if not extended_id.split('_')[1].startswith("0000"):
        logger.warning(f'Skipping structure {extended_id}')
        return None

    old_id = extended_id[-4:]

    return old_id


def unzip_file(src: Path, dst: Path) -> None:
    try:
        if not src.exists():
            raise FileNotFoundError(f"Source file for unzipping not found: {str(src)}")

        with ZipFile(src, "r") as zip_obj:
            zip_obj.extractall(dst)
        os.remove(src)
    except Exception as e:
        logging.error(f"An error occurred during extraction: {e}")


def count_local_files(directory):
    if not os.path.exists(directory):
        return 0
    return len([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))])


def save_to_file(filename, file_list):
    with open(filename, "w") as f:
        for item in file_list:
            f.write(f"{item}\n")


def get_pdb_names_filtered(main_dir: Path):
    csv_files = []
    for ring in Ring:
        csv_files.append(os.path.join(str(main_dir), ring.name.lower(), f"filtered_patterns_{ring.name.lower()}.csv"))

    unique_ids = set()

    for f in csv_files:
        for chunk in pd.read_csv(f, usecols=["ParentId"], chunksize=100_000):
            # unique_ids.update(chunk["ParentId"].dropna())
            unique_ids.update(chunk["ParentId"].dropna().astype(str).str[-4:])

    logger.info(f"Unique pdb ids across all filtered rings: {len(unique_ids)}")
    return unique_ids


def get_data_from_cif(ligand_block: gemmi.cif.Block):

    bond_table = ligand_block.find(
        [
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order",
            "_chem_comp_bond.pdbx_aromatic_flag",
        ]
    )

    atom_table = ligand_block.find(
        ["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"]
    )

    bond_df = pd.DataFrame(
        bond_table, columns=["atom_id_1", "atom_id_2", "value_order", "aromatic"]
    )

    atom_df = pd.DataFrame(atom_table, columns=["atom_id", "type_symbol"])

    return bond_df, atom_df


def get_atom_names(cif_file):
    with open(cif_file, encoding="utf-8-sig") as f:
        data = f.read()
    
    ring_structure = gemmi.read_structure_string(data, format=gemmi.CoorFormat.Mmcif,)
    model = ring_structure[0]
    chain = model[0]
    res = chain[0]
    atom_names = set()
    for atom in res:
        atom_names.add(atom.name)
    return atom_names


def get_atoms_count_from_shape(s: str) -> int:
    return sum(map(int, re.findall(r"\*(\d+)", s)))
