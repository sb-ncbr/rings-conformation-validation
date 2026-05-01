import pandas as pd
from pathlib import Path
from gemmi import cif
from typing import Set
from HelperModule.Ring import Ring


def get_data_from_cif(ligand_block: cif.Block):

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
