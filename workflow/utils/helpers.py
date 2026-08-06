import logging
import os
from pathlib import Path
import re
import gemmi
import pandas as pd
from workflow.models.Ring import Ring

logger = logging.getLogger(__name__)


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

#TEMP
def add_missing_fields(path):
    lines = path.read_text().splitlines()

    # fast check if already fixed
    if any(line.strip() == "_atom_site.auth_seq_id" for line in lines):
        return  # already contains the fields

    output = []
    in_atom_loop = False

    for line in lines:
        if line.startswith("_atom_site."):
            if line == "_atom_site.pdbx_PDB_model_num":
                output.append("_atom_site.auth_seq_id")
                output.append("_atom_site.auth_asym_id")
            output.append(line)
            continue

        if line.startswith(("ATOM", "HETATM")):
            fields = line.split()

            # label_asym_id and label_seq_id
            label_asym = fields[6]
            label_seq = fields[8]

            # insert before model number
            fields.insert(-1, label_seq)
            fields.insert(-1, label_asym)

            output.append(" ".join(fields))
        else:
            output.append(line)

    path.write_text("\n".join(output) + "\n")


def get_atom_names(cif_file):
    # print(cif_file)
    
    # doc = gemmi.cif.read(str(cif_file))
    # ring_structure = gemmi.make_structure_from_block(doc.sole_block())
    # print(len(ring_structure))
    # # print(doc.sole_block().name)

    ring_structure = gemmi.read_structure(str(cif_file))
    model = ring_structure[0]
    chain = model[0]
    res = chain[0]
    atom_names = set()
    for atom in res:
        atom_names.add(atom.name)
    return atom_names


def get_atoms_count_from_shape(s: str) -> int:
    return sum(map(int, re.findall(r"\*(\d+)", s)))

# # FOR DEBUG: from all files in local pdb
# def initialize_pending_queue(pdb_dir: Path, pending_dir: Path):
#     pending_dir.mkdir(parents=True, exist_ok=True)

#     created = 0

#     for cif_file in pdb_dir.glob("*.cif.gz"):
#         link = pending_dir / cif_file.name

#         if not link.exists():
#             os.symlink(
#                 cif_file.resolve(),
#                 link
#             )
#             created += 1

#     return created

# # from file generated during rsync
# def create_pending_links(pdb_dir: Path, pending_file: Path, pending_dir: Path):

#     pending_dir.mkdir(parents=True, exist_ok=True)

#     with open(pending_file) as f:
#         for filename in f:
#             filename = filename.strip()

#             src = pdb_dir / filename
#             dst = pending_dir / filename

#             if src.exists() and not dst.exists():
#                 dst.symlink_to(src.resolve())