from collections import defaultdict
import logging
import shutil
from typing import Dict
import gemmi
import pandas as pd
from gemmi import cif
from pathlib import Path
from workflow.models.Ring import Ring

from workflow.preprocessing.rings_classification import classify_ring
from workflow.utils.helpers import add_missing_fields, get_atom_names, get_data_from_cif


logger = logging.getLogger(__name__)


def process_correct_rings(output_dir, ligand, filepath):
    output_pdb_dir = output_dir / "filtered_ligands" / ligand / "patterns"

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    filename = filepath.name
    if not filename.startswith(f"{ligand}_"):
        filename = f"{ligand}_{filename}"

    new_name_path = output_pdb_dir / filename
    shutil.copy(filepath, new_name_path)


def filter_by_bond_type(
    patterns_df, dir_with_patterns: Path, atoms_shape: str, main_workflow_output_dir: Path, document: gemmi.cif.Document
) -> Dict[str, str] | None:
    if patterns_df is None:
        return
    
    processed_data_dict = {}
    target_ring_rows = defaultdict(list)
    # atoms_number = get_atoms_count_from_shape(atoms_shape)

    for row in patterns_df.itertuples(index=False):
        ligand = row.Residues.split()[0]

        # skipping unknown ligand
        if ligand == "UNL":
            continue

        # 0a_pdb_000010ad_structures__pdb_000010ad_0.cif
        cif_filepath = dir_with_patterns / "patterns" / (row.Id + ".cif")

        # remove BOM for subsequent gemmi processing
        data = cif_filepath.read_bytes()
        if data.startswith(b"\xef\xbb\xbf"):
            cif_filepath.write_bytes(data[3:])

        add_missing_fields(cif_filepath)
        atom_names = get_atom_names(cif_filepath)

        key = (ligand, frozenset(atom_names))

        if key in processed_data_dict:
            # skipping because already filtered out this ring as wrong
            if processed_data_dict[key] is None:
                continue
            
            ring_name = processed_data_dict[key].name.lower()
            output_dir = main_workflow_output_dir / ring_name
            process_correct_rings(output_dir, ligand, cif_filepath)
            target_ring_rows[ring_name].append(row._asdict())
            continue

        ligand_block = document.find_block(ligand)
        if ligand_block is None:
            logger.warning(f"Ligand_block is None for {ligand}")
            continue

        bond_df, atom_df = get_data_from_cif(ligand_block)
        determined_ring_type = classify_ring(
            atoms_shape, atom_names, bond_df, atom_df, cif_filepath, ligand
        )
        processed_data_dict[key] = determined_ring_type

        if determined_ring_type is None:
            continue
        
        ring_name = determined_ring_type.name.lower()
        output_dir = main_workflow_output_dir / ring_name
        process_correct_rings(output_dir, ligand, cif_filepath)
        target_ring_rows[ring_name].append(row._asdict())

    return target_ring_rows


def generate_stats(df, main_dir: Path):
    stats_df = (
        df.groupby("Atoms", as_index=False)
        .agg(
            TotalRingCount=("Atoms", "size"),
            UniqueLigands=("Signature", "nunique"),
            UniquePdbs=("ParentId", "nunique"),
        )
        .sort_values("TotalRingCount", ascending=False)
        .reset_index(drop=True)
    )

    true_totals = {
        "TotalRingCount": len(df),
        "UniqueLigands": df["Signature"].nunique(),
        "UniquePdbs": df["ParentId"].nunique(),
    }

    totals_row = pd.DataFrame([true_totals])
    totals_row["Atoms"] = "TOTAL"

    stats_df_with_total = pd.concat([stats_df, totals_row], ignore_index=True)
    csv_stats = main_dir / "all_rings_stats_nofilter.csv"
    stats_df_with_total.to_csv(csv_stats)
    logger.info(f"Rings statistics exported to {csv_stats}")


def filter_rings(main_workflow_output_dir: Path, path_to_comp_dict: Path, dir_with_patterns: Path):
    if not (dir_with_patterns / "patterns").exists():
        logger.error(f"Folder with patterns (PQ output) does not exist")
        raise RuntimeError
    
    document = cif.read(str(path_to_comp_dict))
    df = pd.read_csv(dir_with_patterns / "patterns.csv")
    generate_stats(df, main_workflow_output_dir)

    target_atoms = [ring.atoms for ring in Ring]
    grouped = dict(tuple(df.groupby("Atoms")))
    dfs_by_ring_structure = {atom: grouped[atom].copy() for atom in target_atoms if atom in grouped}

    result_dict = {}
    atoms_shape_groups = set()
    for ring in Ring:
        atoms_shape_groups.add(ring.atoms)
        output_path = main_workflow_output_dir / ring.name.lower()
        output_path.mkdir(parents=True, exist_ok=True)

    for atoms_shape in atoms_shape_groups:
        logger.info(f"Processing ring(s) with {atoms_shape} shape...")
        
        target_ring_rows = filter_by_bond_type(
            dfs_by_ring_structure.get(atoms_shape),
            dir_with_patterns,
            atoms_shape,
            main_workflow_output_dir,
            document
        )
        if target_ring_rows is not None:
            result_dict |= target_ring_rows

    for ring_name, rows in result_dict.items():
        res_df = pd.DataFrame(rows)
        res_df.to_csv(main_workflow_output_dir / ring_name / f"filtered_patterns_{ring_name}.csv")
        logger.info(
            "%-15s | Rings=%-10s | Unique ligands=%-12s | Unique PDBs=%-10s",
            ring_name,
            f"{len(res_df):,}".replace(",", " "),
            f"{res_df['Signature'].nunique():,}".replace(",", " "),
            f"{res_df['ParentId'].nunique():,}".replace(",", " ")
        )
        logger.info(f"{ring_name} | Rings={len(res_df)} | Unique ligands={res_df["Signature"].nunique()} | Unique PDBs={res_df["ParentId"].nunique()}")
        

