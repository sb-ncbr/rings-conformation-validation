import pandas as pd
from pathlib import Path
from gemmi import cif
from typing import Dict
from HelperModule.Ring import Ring


def get_atoms_from_pdb(path_to_pdb: Path, ring: Ring) -> list[str]:
    atom_names = []
    with open(path_to_pdb, "r") as file:
        lines = file.read().splitlines()

        for line in lines[1 : ring.atom_number + 1]:
            atom_names.append(line[12:16].strip())
    return atom_names


def get_bonds_from_cif(ligand_block: cif.Block) -> list[list[str]]:

    table = ligand_block.find(
        [
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order",
            "_chem_comp_bond.pdbx_aromatic_flag",
        ]
    )

    return [list(x) for x in list(table)]


def get_data_from_cif(ligand_block: cif.Block) -> Dict[str, list[list[str]]]:

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
        bond_table, columns=["atom_id_1", "atom_id_2", "order", "aromatic"]
    )

    atom_df = pd.DataFrame(atom_table, columns=["atom_id", "type_symbol"])

    return bond_df, atom_df
