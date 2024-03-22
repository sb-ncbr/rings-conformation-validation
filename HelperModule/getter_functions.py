from gemmi import cif
from HelperModule.Ring import Ring


def get_data_from_pdb(path_to_pdb: str, ring: Ring) -> tuple[str, list[str]]:
    atom_names = []
    with open(path_to_pdb, 'r') as file:
        lines = file.read().splitlines()
        ligand = lines[1][17:20].strip()

        for line in lines[1:ring.atom_number + 1]:
            atom_names.append(line[12:16].strip())
    return ligand, atom_names


def get_ligands_from_file(filename: str, ligands_set: set) -> None:
    with open(filename, 'r') as ligands:
        curr_ligand = ligands.readline().strip()
        while curr_ligand:
            # exclude PHE and TYR at this step for separate analysis later
            if curr_ligand not in ('PHE', 'TYR'):
                ligands_set.add(curr_ligand)
            curr_ligand = ligands.readline().strip()


def get_bonds_from_cif(ligand_block: cif.Block) -> list[list[str]]:

    table = ligand_block.find(['_chem_comp_bond.atom_id_1', '_chem_comp_bond.atom_id_2',
                               '_chem_comp_bond.value_order', '_chem_comp_bond.pdbx_aromatic_flag'])
    
    return [list(x) for x in list(table)]

