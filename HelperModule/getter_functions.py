from pathlib import Path
from gemmi import cif
from HelperModule.Ring import Ring


def get_data_from_pdb(path_to_pdb: Path, ring: Ring) -> tuple[str, list[str]]:
    atom_names = []
    with open(path_to_pdb, "r") as file:
        lines = file.read().splitlines()
        ligand = lines[1][17:20].strip()

        for line in lines[1 : ring.atom_number + 1]:
            atom_names.append(line[12:16].strip())
    return ligand, atom_names


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
