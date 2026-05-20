from argparse import ArgumentParser
from collections import defaultdict
from datetime import timedelta
from pathlib import Path
import pickle
import time
import numpy as np
import pandas as pd
import json
import logging
from multiprocessing import Pool, cpu_count
from HelperModule.constants import EL_DENSITY_OUTPUT_DIR, MAIN_DIR, PDB_DIR
from HelperModule.Ring import Ring
import gemmi

CPU_COUNT = cpu_count()

logging.basicConfig(
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def get_atom_names_as_string(cif_file):
    ring_structure = gemmi.read_structure(str(cif_file))
    model = ring_structure[0]
    chain = model[0]
    res = chain[0]
    atom_names = []
    for atom in res:
        atom_names.append(atom.name)
    atom_names_str = "-".join(atom_names)
    return atom_names_str


# for oxane, oxolane
def get_conf_info(conf_name: str):
    match conf_name:
        # Chair
        case name if "C" in name:
            return "Chair"

        # Boat
        case name if "B" in name:
            return "Boat"

        # Envelope
        case name if "E" in name:
            return "Envelope"

        # Half-chair
        case name if "H" in name:
            return "Half chair"

        # Skew
        case name if "S" in name:
            return "Skew"
        
        # Twist
        case name if "T" in name:
            return "Twist"
        
        # Flat
        case name if name == 'Flat':
            return "Flat"

        case _:
            return "Unknown"
        

def is_favorable_conf(ring: Ring, conf_name):
    return conf_name in (ring.favourable_confs or [])


def generate_ring_data(df: pd.DataFrame, ring: Ring):
    """
    Build JSON-like dict for a single ring.
    """
    group = df[df["ring_type"] == ring.name.capitalize()]
    if group.empty:
        return None

    ring_data = {
        "name": ring.name.capitalize(),
        "ccd": int(group["ccd_id"].nunique()),   # unique CCD count
        "pdb": int(group["pdb_id"].nunique()),    # unique PDB count
        "total": len(group),                        # total rows
        "confs": []
    }

    conf_counts = group["conformation"].value_counts()

    if ring in [Ring.OXANE, Ring.OXOLANE]:
        conf_classes= defaultdict(lambda: {"subtypes": [], "total": 0})

        for conf_name, count in conf_counts.items():
            conf_class = get_conf_info(conf_name)
            conf_classes[conf_class]["subtypes"].append({
                                                    "name": conf_name,
                                                    "count": count
                                                })
            conf_classes[conf_class]["total"] += count

        for conf_class, data in conf_classes.items():
            conf_data = {
                    "name": conf_class,
                    "count": int(data["total"]),
                    "favorable": is_favorable_conf(ring, conf_class),
                    "subtypes": data["subtypes"] if conf_class != "Flat" else []
                }
            ring_data["confs"].append(conf_data)
                
  
    else:
        for conf_name, count in conf_counts.items():
            conf_data = {
                        "name": conf_name,
                        "count": int(count),
                        "favorable": is_favorable_conf(ring, conf_name),
                        "subtypes": []
                    }

            ring_data["confs"].append(conf_data)

    return ring_data


def format_to_percentage(value):
    if value == "N/A":
        return "N/A"
    try:
        numerator, denominator = map(int, value.split(';'))

        percentage = (numerator / denominator) * 100
        return f"{percentage:.0f}%"
    except (ValueError, ZeroDivisionError):
        return None


def format_to_float(value):
    if value == "N/A":
        return None
    try:
        numerator, denominator = map(int, value.split(';'))
        return numerator / denominator
    except (ValueError, ZeroDivisionError):
        return None



def format_ring_df(df):

    df.rename(columns={"Coverage": "Ring Coverage"}, inplace=True)
    df['Ring Coverage'] = df['Ring Coverage'].replace(np.nan, 'N/A')
    df['Ring Coverage Float'] = df['Ring Coverage'].apply(format_to_float)
    df['Ring Coverage (%)'] = df['Ring Coverage'].apply(format_to_percentage)

    return df


def build_filepath_to_cif(row, ring_name: str, path_to_data):
    return path_to_data / ring_name / "filtered_ligands" / row['Ligand ID'] / "patterns" / str(row["Ring_ID"] + ".cif")


def get_atoms_from_row(row, ring_name: str, path_to_data):
    path = build_filepath_to_cif(row, ring_name, path_to_data)
    return get_atom_names_as_string(path)


def create_json(df, stats_json_path):
    logging.info(f"Creating JSON for statistics")
    rings_json = []
    for ring in Ring:
        ring_data = generate_ring_data(df, ring)
        if ring_data is None:
            continue
        rings_json.append(ring_data)


    counts = df['experimental_method'].dropna().value_counts()

    exp_methods = [
        {
            "method": method,
            "count": int(count)
        }
        for method, count in counts.items()
    ]

    stats = {
        "summary": {
            "rings": int(len(df)),
            "ligands": int(df['ccd_id'].nunique()),
            "pdbEntries": int(df['pdb_id'].nunique())
        },
        "experimentalMethods": exp_methods,
        "rings": rings_json
    }

    with open(stats_json_path, "w") as f:
        json.dump(stats, f, indent=2)

    logging.info(f"Generated stats: {len(rings_json)} rings, {len(exp_methods)} methods")
    


def update_ring_columns(df, ring: Ring, path_to_data):
    df['Ring Type'] = ring.name.capitalize()
    df.rename(columns={"Ligand_name": "Ligand ID"}, inplace=True)
    if ring not in [Ring.OXANE, Ring.OXOLANE]:
        df['Conformation'] = (df['Conformation'].
                              str.replace('_', ' ').str.capitalize())

    df['Conformation'] = df['Conformation'].replace({
        'Tw boat': 'Twist boat',
        'P': 'Flat'  # for oxolane/oxane change P as Planar to Flat
    })

    df["Residue ID"] = (
        df["Residues"]
        .str.split()
        .str[1]
    )

    df["Chain ID"] = (
        df["Residues"]
        .str.split()
        .str[2]
    )

    df["PDB_ins_code"] = (
        df["Residues"]
        .str.split()
        .str[3]
        .fillna("")
    )

    df["atom_names"] = df.apply(
        lambda row: get_atoms_from_row(row, ring.name.lower(), path_to_data),
        axis=1
    )

    df["id"] = (
    df["PDB ID"].astype(str) + "_" +
    df["Chain ID"].astype(str) + "_" +
    df["Ligand ID"].astype(str) + "_" +
    df["Residue ID"].astype(str) +
    df["PDB_ins_code"].apply(lambda x: f"_{x}" if x else "") + "_" +
    df["atom_names"].astype(str)
)

    df.drop(columns=["Ring_ID"], inplace=True)
    df.drop(columns=["pq_id"], inplace=True)
    df.drop(columns=["atom_names"], inplace=True)
    df.drop(columns=["Residues"], inplace=True)


def process_ring(ring: Ring, path_to_data):
    ring_lower = ring.name.lower()
    

    hr_results_filepath = path_to_data / ring_lower / "hr_analysis_output" / "result_conf_chart.csv"
    if not hr_results_filepath.exists():
        return None
    logging.info(f"Processing {ring_lower}")
    conf_df = pd.read_csv(hr_results_filepath, delimiter=';')

    patterns_df = pd.read_csv(
        path_to_data / ring_lower / f"filtered_patterns_{ring_lower}.csv",
        usecols=["Id", "ParentId", "Residues"],
        dtype=str
    )
    patterns_df.rename(columns={"ParentId": "PDB ID"}, inplace=True)
    conf_df["pq_id"] = conf_df["Ring_ID"].str.split("_", n=1).str[1]

    merged_with_pq_df = conf_df.merge(
        patterns_df,
        left_on="pq_id",
        right_on="Id",
        how="left"  # keep all rows from conf_df, add matches from patterns_df
    ).drop(columns=["Id"])

    el_density_path = path_to_data / ring_lower / EL_DENSITY_OUTPUT_DIR / f"{ring_lower}_params__analysis_output.csv"
    el_density_df = pd.read_csv(el_density_path, delimiter=',')
    merged_with_coverage_df = pd.merge(merged_with_pq_df, el_density_df[['Id', 'Coverage']],
                               left_on="Ring_ID",
                               right_on="Id",
                               how='left').drop(columns=["Id"])

    update_ring_columns(merged_with_coverage_df, ring, path_to_data)

    return merged_with_coverage_df


def format_single_method(method: str) -> str:
    SPECIAL_MAP = {
        "X-RAY": "X-ray",
        "NMR": "NMR",
    }
    words = method.upper().split()

    formatted = []
    for i, w in enumerate(words):
        if w in SPECIAL_MAP:
            token = SPECIAL_MAP[w]
        else:
            token = w.lower()

        if i == 0 and w not in SPECIAL_MAP:
            token = token.capitalize()

        formatted.append(token)

    return " ".join(formatted)


def format_methods(df):
    df["experimental_method"] = df["experimental_method"].apply(
        lambda x: ", ".join(
            format_single_method(part.strip()) for part in x.split(",")
        ) if pd.notna(x) else x
    )
    counts = df['experimental_method'].value_counts()
    print(counts)


def check_for_duplicates(df, output_dir):
    dup_mask = df.index.duplicated(keep=False)

    n_duplicates = dup_mask.sum()

    if n_duplicates > 0:
        logging.error(f"Duplicates found: {n_duplicates} rows affected")

        # select all duplicated rows
        duplicated_rows = df[dup_mask]
        duplicated_rows.to_csv(output_dir / "duplicate_rows.csv")

        df.index[dup_mask].to_series().to_csv(
            output_dir / "duplicate_indexes.csv", index=False
        )

        df_clean = df[~dup_mask]   # removes ALL duplicates
        return df_clean

    else:
        logging.info("No duplicates among IDs")
        return df


def combine_with_valtrends_data(df, path_to_valtrends_data):
    logging.info(f"Combining with valtrends data")
    valtrends_df = pd.read_csv(
        path_to_valtrends_data,
        sep=";",
        usecols=["PDB ID", "averageLigandRSCC", "averageLigandRSR"]
    )

    valtrends_df["PDB ID"] = (
    valtrends_df["PDB ID"]
    .astype(str)
    .str.replace('"', '')
    .str.strip()
    .str.lower()
)

    combined_df = df.merge(
        valtrends_df,
        on="PDB ID",
        how="left"  # keep all rows from all_rings_df, add matches from valtrends_df
    )

    return combined_df


def reformat_final_data(df):
    logging.info(f"Reformatting final data")
    df.set_index('id', inplace=True)
    df = df.sort_index()

    df['averageLigandRSCC'] = df['averageLigandRSCC'].round(3)
    df['averageLigandRSR'] = df['averageLigandRSR'].round(2)

    df['averageLigandRSCC'] = df['averageLigandRSCC'].fillna('N/A').astype(str)
    df['averageLigandRSR'] = df['averageLigandRSR'].fillna('N/A').astype(str)

    df["Resolution"] = (
        pd.to_numeric(df["Resolution"].replace(".", np.nan), errors="coerce")
        .round(2)
        .astype("string")
        .fillna("N/A")
    )

    new_column_order = ['PDB ID', 'Chain ID', 'Ligand ID', 'Residue ID', 'PDB_ins_code', 'Ring Type', 'Resolution', 'Ring Coverage',
                        'Ring Coverage Float', 'Ring Coverage (%)', 'Conformation', 'conf_main_type', 'Experimental Method',
                        'averageLigandRSR', 'averageLigandRSCC']
    df = df[new_column_order]
    new_column_names = ['pdb_id', 'chain_id', 'ccd_id', 'residue_id', 'PDB_ins_code', 'ring_type', 'resolution', 'ring_coverage_counts',
                        'ring_coverage_float', 'ring_coverage', 'conformation', 'conf_main_type', 'experimental_method', 'rsr',
                        'rscc']
    df.columns = new_column_names

    format_methods(df)
    return df


def extract_metadata(cif_filepath: Path):
    try:
        doc = gemmi.cif.read(str(cif_filepath))
        block = doc.sole_block()

        method = block.find_value("_exptl.method")
        if method:
            method = method.strip("'")
        else:
            # loop = block.find_loop("_exptl.method")
            methods = [m.strip("'") for m in block.find_loop("_exptl.method")]
            method = ", ".join(methods)
        
        resolution = block.find_value("_refine.ls_d_res_high")
        return method, resolution, cif_filepath.stem.removesuffix('.cif')
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    

def add_metadata_from_cif(df, input_dir):
    logging.info("Extracting info about experimental methods and resolution...")
    cache = {}
    cache_file = Path("./cache") / "methods_and_resolution.pkl"
    if cache_file.exists():
        logging.info("Cache file loaded.")
        with cache_file.open('rb') as f:
            cache = pickle.load(f)
        df[["Experimental Method", "Resolution"]] = df["PDB ID"].map(cache).apply(pd.Series) 
        return df
    
    total = len(df["PDB ID"].unique())
    paths_list = [Path(input_dir) / PDB_DIR / f"{pdb_id}.cif.gz" for pdb_id in sorted(df["PDB ID"].unique())]

    with Pool(int(CPU_COUNT)) as p:
        for i, output in enumerate(p.imap(extract_metadata, paths_list),1):
            method, res, pdb_id = output
            logging.info(f"[{i}/{total}] | pdb_id: {pdb_id} method: {method} | res: {res}")
            cache[pdb_id] = (method, res)

    logging.info("Saving to cache file...")
    with cache_file.open('wb') as f:
        pickle.dump(cache, f)

    logging.info("Enriching the dataframe with exp. method and resolution data...")
    df[["Experimental Method", "Resolution"]] = df["PDB ID"].map(cache).apply(pd.Series)   
    return df  
            

def main(output_dir, input_dir):
    final_output_path = Path(output_dir) / "web"
    final_output_path.mkdir(parents=True, exist_ok=False)
    stats_json_path = final_output_path / "stats.json"
    valtrends_data_path = Path(input_dir) / "data.csv"

    all_rings: list[pd.DataFrame] = []

    for ring in Ring:
        ring_df = process_ring(ring, Path(output_dir) / MAIN_DIR)
        if ring_df is None:
            continue

        if ring in [Ring.OXANE, Ring.OXOLANE]:
            ring_df["conf_main_type"] = ring_df["Conformation"].apply(get_conf_info)
        else:
            ring_df["conf_main_type"] = ring_df["Conformation"]
        all_rings.append(format_ring_df(ring_df))

    all_rings_df = pd.concat(all_rings, ignore_index=False)
    combined_df = combine_with_valtrends_data(all_rings_df, valtrends_data_path)
    complete_df = add_metadata_from_cif(combined_df, input_dir)
    reformatted_df = reformat_final_data(complete_df)
    final_df = check_for_duplicates(reformatted_df, final_output_path)
    create_json(final_df, stats_json_path)
    logging.info(f"Final data are being written to {final_output_path}")

    excel = final_output_path / "rings.xlsx"
    final_df.to_excel(excel, index=True)

    data = pd.read_excel(excel, engine="openpyxl", index_col=0)
    data.to_parquet(final_output_path / "rings.parquet")


if __name__ == "__main__":
    parser = ArgumentParser()
    required = parser.add_argument_group("required named arguments")

    required.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the workflow output directory (should contain validation_data folder).",
    )
    required.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to the directory with input data (local pdb, ccp4 files, etc.)",
    )

    args = parser.parse_args()
    start = time.perf_counter()
    
    main(args.output, args.input)

    elapsed = time.perf_counter() - start
    formatted = str(timedelta(seconds=elapsed))
    logging.info(f"Total time: {formatted}")

