import logging
from typing import Set
from workflow.models.Ring import Ring


logger = logging.getLogger(__name__)
validation_logger = logging.getLogger("validation")


def is_oxane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
    if not all_single:
        return False
    
    if aromatic_sum > 0:
        validation_logger.info(f"[Skipped oxane in ligand {ligand} atoms: {sorted_atoms}]: all bonds labelled as single AND at {aromatic_sum} bond(s) labelled as aromatic: {filepath}")
        return False
    return True


def is_oxolane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
    if not all_single:
        return False
    
    if aromatic_sum > 0:
        validation_logger.info(f"[Skipped oxolane in ligand {ligand} atoms: {sorted_atoms}] all bonds labelled as single AND {aromatic_sum} bond(s) labelled as aromatic: {filepath}")
        return False
    return True


def is_cyclopentane(bond_df, atom_names, metal_atoms, all_single, aromatic_sum, filepath, ligand):
    if not all_single:
        return False
    
    if aromatic_sum > 0:
        validation_logger.info(f"[Skipped cyclopentane in ligand {ligand} atoms: {sorted(atom_names)}] all bonds labelled as single AND {aromatic_sum} bond(s) labelled as aromatic: {filepath}")
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


def is_cyclohexane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
    if not all_single:
        return False
    
    # condition <aromatic_count> == 6 is checked in is_benzene function
    if 0 < aromatic_sum < Ring.CYCLOHEXANE.atom_number:
        validation_logger.info(f"[Skipped cyclohexane in ligand {ligand} atoms: {sorted_atoms}] all bonds labelled as single AND {aromatic_sum} bond(s) labelled as aromatic: {filepath}")
        return False
    # double check
    if aromatic_sum == 0:
        return True


def is_benzene(current_ring_df, all_single, aromatic_sum, filepath, ligand, sorted_atoms):
    double_bonds_count = (current_ring_df["value_order"].str.upper() == "DOUB").sum()

    if aromatic_sum == Ring.BENZENE.atom_number:
        if all_single:
            validation_logger.info(f"[Skipped benzene in ligand {ligand} atoms: {sorted_atoms}] all bonds labelled as aromatic AND single: {filepath}")
            return False
        if double_bonds_count not in (2, 3):
            validation_logger.info(f"[Skipped benzene in ligand {ligand} atoms: {sorted_atoms}]: Skipping a ring with {double_bonds_count} bonds labelled as double and all aromatic: {filepath}")
            return False
    
        return True
    
    if aromatic_sum == 0:
        
        if double_bonds_count == 3:
            validation_logger.info(f"[Possible Benzene included in ligand {ligand} atoms: {sorted_atoms}] has THREE bonds labelled as double BUT all bonds as NOT aromatic: {filepath}")
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
    aromatic_sum = (current_ring_df["aromatic"].str.upper() == "Y").sum()
    sorted_atoms = sorted(atom_names)
    match atoms_shape:
        case "C*6":
            if is_benzene(current_ring_df, all_single, aromatic_sum, filepath, ligand, sorted_atoms):
                return Ring.BENZENE
            if is_cyclohexane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
                return Ring.CYCLOHEXANE
            return None
        case "C*5":
            metal_atoms = ["FE", "MN", "CO", "RU", "TI", "ZR", "NI", "CR", "RH", "IR", "RE", "OS"]
            if is_cyclopentane(bond_df, atom_names, metal_atoms, all_single, aromatic_sum, filepath, ligand):
                return Ring.CYCLOPENTANE
            return None
        case "C*5-O*1":
            if is_oxane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
                return Ring.OXANE
            return None
        case "C*4-O*1":
            if is_oxolane(all_single, aromatic_sum, filepath, ligand, sorted_atoms):
                return Ring.OXOLANE
            return None
        case _:
            return None