from argparse import ArgumentParser
from collections import defaultdict
from datetime import timedelta
from pathlib import Path
import time
import numpy as np
import pandas as pd
import json
import logging
from HelperModule.constants import MAIN_DIR
import swifter
import gemmi



# here should be listed all analyzed rings
RINGS_FAVORABLE_CONFS_MAP = {
    "cyclohexane": ["Chair"],
    "cyclopentane": ["Envelope", "Half chair"],
    "benzene": ["Flat"],
    "oxane": ["Chair"],
    "oxolane": "ALL_EXCEPT_UNFAVORABLE"
}

RINGS_UNFAVORABLE_CONFS_MAP = {
    "oxolane": ["Flat"]
}
ATOMS_COUNT = {
    "cyclohexane": 6,
    "cyclopentane": 5,
    "benzene": 6,
    "oxane": 6,
    "oxolane": 5
}


logging.basicConfig(
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def get_atom_names(cif_file):
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
        

def is_favorable_conf(ring_name, conf_name):
    favorable = RINGS_FAVORABLE_CONFS_MAP.get(ring_name)

    if favorable == "ALL_EXCEPT_UNFAVORABLE":
        return conf_name not in RINGS_UNFAVORABLE_CONFS_MAP.get(ring_name, [])

    return conf_name in (favorable or [])


def generate_ring_data(df: pd.DataFrame, ring_name):
    """
    Build JSON-like dict for a single ring.
    """
    group = df[df["ring_type"] == ring_name.capitalize()]

    ring_data = {
        "name": ring_name.capitalize(),
        "ccd": int(group["ccd_id"].nunique()),   # unique CCD count
        "pdb": int(group["pdb_id"].nunique()),    # unique PDB count
        "total": len(group),                        # total rows
        "confs": []
    }

    conf_counts = group["conformation"].value_counts()

    if ring_name.lower() in ["oxane", "oxolane"]:
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
                    "favorable": is_favorable_conf(ring_name, conf_class),
                    "subtypes": data["subtypes"] if conf_class != "Flat" else []
                }
            ring_data["confs"].append(conf_data)
                
  
    else:
        for conf_name, count in conf_counts.items():
            conf_data = {
                        "name": conf_name,
                        "count": int(count),
                        "favorable": is_favorable_conf(ring_name, conf_name),
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


def format_coverage(value, ring_type):
    if value == "N/A":
        return "N/A"

    # get integer part of value
    score = int(float(value))

    # get atom count (default None if unknown ring type)
    atoms = ATOMS_COUNT.get(ring_type)

    if atoms is None:
        raise ValueError(f"Unknown ring type: {ring_type}")

    return f"{score};{atoms}"


def format_ring_df(df, ring):

    df.rename(columns={"Coverage": "Ring Coverage"}, inplace=True)
    df['Ring Coverage'] = df['Ring Coverage'].replace(np.nan, 'N/A')
    # Apply row-wise
    df["Ring Coverage"] = df.swifter.apply(
        lambda row: format_coverage(row["Ring Coverage"], ring),
        axis=1
    )
    df['Ring Coverage Float'] = df['Ring Coverage'].swifter.apply(format_to_float)
    df['Ring Coverage (%)'] = df['Ring Coverage'].swifter.apply(format_to_percentage)

    df['Resolution (A)'] = df['Resolution (A)'].round(2).astype(str)
    df['Resolution (A)'] = df['Resolution (A)'].replace('nan', 'N/A')

    return df


def build_filepath_to_cif(row, ring, path_to_data):
    return path_to_data / ring / "filtered_ligands" / row['Ligand ID'] / "patterns" / str(row["Ring_ID"] + ".cif")


def get_atoms_from_row(row, ring, path_to_data):
    path = build_filepath_to_cif(row, ring, path_to_data)
    return get_atom_names(path)


def create_json(df, stats_json_path):
    logging.info(f"Creating JSON for statistics")
    rings_json = []
    for ring in RINGS_FAVORABLE_CONFS_MAP.keys():
        ring_data = generate_ring_data(df, ring)
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
            "ligands": int(df['residue_id'].nunique()),
            "pdbEntries": int(df['pdb_id'].nunique())
        },
        "experimentalMethods": exp_methods,
        "rings": rings_json
    }

    with open(stats_json_path, "w") as f:
        json.dump(stats, f, indent=2)

    logging.info(f"Generated stats: {len(rings_json)} rings, {len(exp_methods)} methods")
    


def update_ring_columns(df, ring, path_to_data):
    df['Ring Type'] = ring.capitalize()
    df.rename(columns={"Ligand_name": "Ligand ID"}, inplace=True)

    if ring not in ["oxane", "oxolane"]:
        df['Conformation'] = (df['Conformation'].
                              str.replace('_', ' ').str.capitalize())

    df['Conformation'] = df['Conformation'].replace({
        'Tw boat': 'Twist boat',
        'P': 'Flat'  # for oxolane change P as Planar to Flat
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
        lambda row: get_atoms_from_row(row, ring, path_to_data),
        axis=1
    )

    df["id"] = (
    df["Entry ID"].astype(str) + "_" +
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


def process_ring(ring, path_to_data):
    logging.info(f"Processing {ring}")

    ring_data_path = path_to_data / ring / "final_results" / "result_summary.xlsx"
    ring_df = pd.read_excel(ring_data_path)
    patterns_df = pd.read_csv(
        path_to_data / ring / f"filtered_patterns_{ring}.csv",
        usecols=["Id", "Residues"],
        dtype=str
    )
    ring_df["pq_id"] = ring_df["Ring_ID"].str.split("_", n=1).str[1]

    pq_rings_combined_df = ring_df.merge(
        patterns_df,
        left_on="pq_id",
        right_on="Id",
        how="left"  # keep all rows from ring_df, add matches from patterns_df
    ).drop(columns=["Id"])

    update_ring_columns(pq_rings_combined_df, ring, path_to_data)

    return pq_rings_combined_df


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
    .str.upper()
)

    combined_df = df.merge(
        valtrends_df,
        left_on="Entry ID",
        right_on="PDB ID",
        how="left"  # keep all rows from all_rings_df, add matches from valtrends_df
    ).drop(columns=["PDB ID"])

    return combined_df


def reformat_final_data(df):
    logging.info(f"Reformatting final data")
    df.set_index('id', inplace=True)
    df = df.sort_index()

    df['averageLigandRSCC'] = df['averageLigandRSCC'].round(3)
    df['averageLigandRSR'] = df['averageLigandRSR'].round(2)

    df['averageLigandRSCC'] = df['averageLigandRSCC'].fillna('N/A').astype(str)
    df['averageLigandRSR'] = df['averageLigandRSR'].fillna('N/A').astype(str)

    new_column_order = ['Entry ID', 'Chain ID', 'Ligand ID', 'Residue ID', 'PDB_ins_code', 'Ring Type', 'Resolution (A)', 'Ring Coverage',
                        'Ring Coverage Float', 'Ring Coverage (%)', 'Conformation', 'Experimental Method',
                        'averageLigandRSR', 'averageLigandRSCC']
    df = df[new_column_order]
    new_column_names = ['pdb_id', 'chain_id', 'ccd_id', 'residue_id', 'PDB_ins_code', 'ring_type', 'resolution', 'ring_coverage_counts',
                        'ring_coverage_float', 'ring_coverage', 'conformation', 'experimental_method', 'rsr',
                        'rscc']
    df.columns = new_column_names

    format_methods(df)
    return df


def main(output_dir, input_dir):
    final_output_path = Path(output_dir) / "web"
    final_output_path.mkdir(parents=True, exist_ok=False)
    stats_json_path = final_output_path / "stats.json"
    valtrends_data_path = Path(input_dir) / "data.csv"

    all_rings: list[pd.DataFrame] = []

    for ring in RINGS_FAVORABLE_CONFS_MAP.keys():
        ring_df = process_ring(ring, Path(output_dir) / MAIN_DIR)
        all_rings.append(format_ring_df(ring_df, ring))

    all_rings_df = pd.concat(all_rings, ignore_index=False)
    combined_df = combine_with_valtrends_data(all_rings_df, valtrends_data_path)
    reformatted_df = reformat_final_data(combined_df)
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

