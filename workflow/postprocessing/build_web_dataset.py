from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import json
import logging
from multiprocessing import Pool, cpu_count
from workflow.models.Ring import Ring
import gemmi

CPU_COUNT = cpu_count()
EL_DENSITY_OUTPUT_DIR = "el-density-output"

logger = logging.getLogger(__name__)


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


def create_json(df, stats_json_path):
    logger.info(f"Creating JSON for statistics")
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

    logger.info(f"Generated stats: {len(rings_json)} rings, {len(exp_methods)} methods")
    


def update_ring_columns(df, ring: Ring):
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

    df["id"] = (
    df["extended_PDB_ID"].astype(str) + "_" +
    df["Chain ID"].astype(str) + "_" +
    df["Ligand ID"].astype(str) + "_" +
    df["Residue ID"].astype(str) +
    df["PDB_ins_code"].apply(lambda x: f"_{x}" if x else "") + "_" +
    df["AtomNames"].astype(str)
)

    df.drop(columns=["Ring_ID"], inplace=True)
    df.drop(columns=["pq_id"], inplace=True)
    df.drop(columns=["AtomNames"], inplace=True)
    df.drop(columns=["Residues"], inplace=True)


def get_old_pdb_id(extended_id):
    pdb_id = extended_id.split("_")[1]
    return pdb_id[4:] if pdb_id[:4] == "0000" else None


def process_ring(ring: Ring, path_to_data):
    ring_lower = ring.name.lower()
    

    hr_results_filepath = path_to_data / ring_lower / "hr_analysis_output" / "result_conf_chart.csv"
    if not hr_results_filepath.exists():
        return None
    logger.info(f"Processing {ring_lower}")
    conf_df = pd.read_csv(hr_results_filepath, delimiter=';')

    patterns_df = pd.read_csv(
        path_to_data / ring_lower / f"filtered_patterns_{ring_lower}.csv",
        usecols=["Id", "ParentId", "Residues", "AtomNames"],
        dtype=str
    )
    patterns_df.rename(columns={"ParentId": "extended_PDB_ID"}, inplace=True)
    patterns_df["PDB ID"] = patterns_df["extended_PDB_ID"].apply(get_old_pdb_id)
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

    update_ring_columns(merged_with_coverage_df, ring)

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
    logger.info(counts)


def check_for_duplicates(df, output_dir):
    dup_mask = df.index.duplicated(keep=False)

    n_duplicates = dup_mask.sum()

    if n_duplicates > 0:
        logger.error(f"Duplicates found: {n_duplicates} rows affected")

        # select all duplicated rows
        duplicated_rows = df[dup_mask]
        duplicated_rows.to_csv(output_dir / "duplicate_rows.csv")

        df.index[dup_mask].to_series().to_csv(
            output_dir / "duplicate_indexes.csv", index=False
        )

        df_clean = df[~dup_mask]   # removes ALL duplicates
        return df_clean

    else:
        logger.info("No duplicates among IDs")
        return df


def reformat_final_data(df):
    logger.info(f"Reformatting final data")
    df.set_index('id', inplace=True)
    df = df.sort_index()

    df["Resolution"] = (
        pd.to_numeric(df["Resolution"].replace(".", np.nan), errors="coerce")
        .round(2)
        .astype("string")
        .fillna("N/A")
    )

    new_column_order = ['PDB ID', 'extended_PDB_ID', 'Chain ID', 'Ligand ID', 'Residue ID', 'PDB_ins_code', 'Ring Type', 'Resolution', 'Ring Coverage',
                        'Ring Coverage Float', 'Ring Coverage (%)', 'Conformation', 'conf_main_type', 'Experimental Method'
                        ]
    df = df[new_column_order]
    new_column_names = ['pdb_id', 'extended_PDB_ID', 'chain_id', 'ccd_id', 'residue_id', 'PDB_ins_code', 'ring_type', 'resolution', 'ring_coverage_counts',
                        'ring_coverage_float', 'ring_coverage', 'conformation', 'conf_main_type', 'experimental_method'
                        ]
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
        logger.error(e, stack_info=True, exc_info=True)
    

def add_metadata_from_cif(df, pdb_dir, methods_info: Path):

    logger.info("Extracting info about experimental methods and resolution...")
    methods_info.parent.mkdir(exist_ok=True)

    metadata = {}

    # load existing metadata
    if methods_info.exists():
        old = pd.read_csv(methods_info, sep="\t")
        metadata = dict(
            zip(
                old["extended_PDB_ID"],
                zip(old["Experimental Method"], old["Resolution"])
            )
        )

    # extract missing entries
    missing = set(df["extended_PDB_ID"]) - set(metadata)

    if missing:
        missing_metadata = {}
        paths = [
            pdb_dir / pdb_id[9:11] / pdb_id / "structures" / f"{pdb_id}.cif.gz"
            for pdb_id in missing
        ]

        with Pool(int(CPU_COUNT)) as p:
            for method, res, pdb_id in p.imap(extract_metadata, paths):
                logger.info(
                    "%s | method: %s | resolution: %s",
                    pdb_id, method, res
                )
                missing_metadata[pdb_id] = (method, res)

        metadata.update(missing_metadata)

        # overwrite metadata file with updated content
        metadata_df = pd.DataFrame(
            [
                (pdb_id, method, res)
                for pdb_id, (method, res) in metadata.items()
            ],
            columns=[
                "extended_PDB_ID",
                "Experimental Method",
                "Resolution",
            ],
        )

        metadata_df.to_csv(
            methods_info,
            sep="\t",
            index=False,
        )

    else:
        logger.info("All metadata already available.")

    df[["Experimental Method", "Resolution"]] = (
        df["extended_PDB_ID"].map(metadata).apply(pd.Series)
    )

    return df
            

def create_data_for_web(output_dir: Path, main_dir: Path, pdb_dir: Path, methods_info: Path):
    final_output_path = output_dir / "web"
    final_output_path.mkdir(parents=True, exist_ok=True)
    stats_json_path = final_output_path / "stats.json"

    all_rings: list[pd.DataFrame] = []

    for ring in Ring:
        ring_df = process_ring(ring, main_dir)
        if ring_df is None:
            continue

        if ring in [Ring.OXANE, Ring.OXOLANE]:
            ring_df["conf_main_type"] = ring_df["Conformation"].apply(get_conf_info)
        else:
            ring_df["conf_main_type"] = ring_df["Conformation"]
        all_rings.append(format_ring_df(ring_df))

    all_rings_df = pd.concat(all_rings, ignore_index=False)
    complete_df = add_metadata_from_cif(all_rings_df, pdb_dir, methods_info)
    reformatted_df = reformat_final_data(complete_df)
    final_df = check_for_duplicates(reformatted_df, final_output_path)
    create_json(final_df, stats_json_path)
    logger.info(f"Final data are being written to {final_output_path}")

    excel = final_output_path / "rings.xlsx"
    final_df.to_excel(excel, index=True)

    data = pd.read_excel(excel, engine="openpyxl", index_col=0)
    data.to_parquet(final_output_path / "rings.parquet")


