from collections import defaultdict
import json
import logging
import shutil
from typing import Dict, Set
import gemmi
import pandas as pd
from gemmi import cif
from pathlib import Path
from workflow.models.Ring import Ring
from workflow.preprocessing.rings_classification import classify_ring
from workflow.utils.helpers import add_missing_fields, get_atom_names, get_data_from_cif


logger = logging.getLogger(__name__)

RINGS_CLASSES = "rings_classified.json"
CACHE = Path("cache") / "filter_rings"
PROCESSED_CACHE = CACHE / "processed.json"


def process_correct_rings(output_dir, ligand, filepath, ring_id_with_ext_pdb):
    output_pdb_dir = output_dir / "filtered_ligands" / ligand / "patterns"

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    filename = ring_id_with_ext_pdb + ".cif"
    if not filename.startswith(f"{ligand}_"):
        filename = f"{ligand}_{filename}"

    new_name_path = output_pdb_dir / filename
    shutil.copy(filepath, new_name_path)


def filter_by_bond_type(
    patterns_df, dir_with_patterns: Path, atoms_shape: str, main_workflow_output_dir: Path, document: gemmi.cif.Document, processed_data_dict, processed: Set[str], to_process: Set[str]
) -> Dict[str, str] | None:
    if patterns_df is None:
        return
    
    target_ring_rows = defaultdict(list)

    for pq_ring_id, parent_id, residue, signature in patterns_df[["Id", "ParentId", "Residues", "Signature"]].itertuples(index=False, name=None):
        ligand = residue.split()[0]

        # skipping unknown ligand
        if ligand == "UNL":
            continue

        # path to file: 0a_pdb_000010ad_structures__pdb_000010ad_0.cif
        # row.Id: pdb_000010ad_0
        prefix = parent_id[9:11]
        path_stem = "_".join([prefix, parent_id, "structures", "", pq_ring_id])
        cif_filepath = dir_with_patterns / "patterns" / (path_stem + ".cif")

        if str(cif_filepath) not in to_process:
            continue

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
                processed.add(str(cif_filepath))
                continue
            
            ring_name = processed_data_dict[key].name.lower()
            output_dir = main_workflow_output_dir / ring_name
            process_correct_rings(output_dir, ligand, cif_filepath, pq_ring_id)
            target_ring_rows[ring_name].append({"Id": pq_ring_id, "ParentId": parent_id, "Residues": residue, "Signature": signature, "AtomNames": "-".join(sorted(atom_names))})
            processed.add(str(cif_filepath))
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
            processed.add(str(cif_filepath))
            continue
        
        ring_name = determined_ring_type.name.lower()
        output_dir = main_workflow_output_dir / ring_name
        
        process_correct_rings(output_dir, ligand, cif_filepath, pq_ring_id)
        processed.add(str(cif_filepath))
        target_ring_rows[ring_name].append({"Id": pq_ring_id, "ParentId": parent_id, "Residues": residue,  "Signature": signature, "AtomNames": "-".join(sorted(atom_names))})

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


def filter_rings(main_workflow_output_dir: Path, path_to_comp_dict: Path, dir_with_patterns: Path, state_dir: Path):
    pq_output = dir_with_patterns / "patterns"
    if not pq_output.exists():
        logger.error(f"Folder with patterns (PQ output) does not exist")
        raise RuntimeError


    all_files = {str(p) for p in pq_output.rglob("*.cif")}

    processed = set()
    if PROCESSED_CACHE.exists():
        with open(PROCESSED_CACHE, 'r') as f:
            processed = set(json.load(f))

    to_process = all_files - processed

    ring_dirs = [main_workflow_output_dir / ring.name.lower() / "filtered_ligands" for ring in Ring]
    filtered_files = {
        p.stem
        for ring_dir in ring_dirs
        if ring_dir.exists()
        for p in ring_dir.rglob("*.cif")
    }
    logger.info(f"Available input files: {len(all_files)}")
    logger.info(f"Filtered files: {len(filtered_files)}")

    if len(to_process) == 0:
        logger.info("No unprocessed files left.")
        return

    logger.info(f"Need to be processed: {len(to_process)}")
 
    document = cif.read(str(path_to_comp_dict))
    df = pd.read_csv(dir_with_patterns / "patterns.csv")
    df["Id"] = (df["Id"]
        .str.split("__")
        .str[1]
    )
    df["ParentId"] = (df["ParentId"]
            .str.split("__")
            .str[1]
        )

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

    rings_classes_json = state_dir / RINGS_CLASSES

    processed_data_dict = {}
    if rings_classes_json.exists():
        with open(rings_classes_json) as f:
            cache = json.load(f)

            processed_data_dict = {
                (item["ligand"], frozenset(item["atoms"])): (
                    Ring[item["ring_type"]] if item["ring_type"] else None
                )
                for item in cache
            }


    for atoms_shape in atoms_shape_groups:
        logger.info(f"Processing ring(s) with {atoms_shape} shape...")
        
        target_ring_rows = filter_by_bond_type(
            dfs_by_ring_structure.get(atoms_shape),
            dir_with_patterns,
            atoms_shape,
            main_workflow_output_dir,
            document,
            processed_data_dict,
            processed,
            to_process
        )

        if target_ring_rows is not None:
            result_dict |= target_ring_rows

    rings_classes = [{
        "ligand": ligand,
        "atoms": sorted(atom_names),
        "ring_type": ring_type.name if ring_type else None,
        } for (ligand, atom_names), ring_type in processed_data_dict.items()]

    with open(rings_classes_json, "w") as f:
        json.dump(rings_classes, f, indent=2)

    stats = {}
    for ring_name, rows in result_dict.items():
        res_df = pd.DataFrame(rows)
        res_df.to_csv(main_workflow_output_dir / ring_name / f"filtered_patterns_{ring_name}.csv")
        
        stats[ring_name] = {
                "rings": len(res_df),
                "unique_ligands": int(res_df["Signature"].nunique()),
                "unique_pdbs": int(res_df["ParentId"].nunique()),
            }

    with open(main_workflow_output_dir / "filtered_ring_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    PROCESSED_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(PROCESSED_CACHE, 'w') as f:
        json.dump(sorted(processed), f, indent=2)
        

