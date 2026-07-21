from datetime import timedelta
from pathlib import Path
import logging
import os
import sys
from typing import List
import requests
import subprocess
import argparse
import time

logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )


def create_rsync_filters(
    pdb_filenames: List[str],
    output_dir: str,
    chunk_size: int = 10000
):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    files = []
    current_chunk = []

    for pdb_id in pdb_filenames:
        current_chunk.append(f"+ /{pdb_id}\n")

        if len(current_chunk) >= chunk_size:
            current_chunk.append("- /*\n")

            filter_file = output_dir / f"rsync_filter_{len(files):03d}.txt"
            filter_file.write_text("".join(current_chunk))

            files.append(filter_file)
            current_chunk = []

    if current_chunk:
        current_chunk.append("- /*\n")

        filter_file = output_dir / f"rsync_filter_{len(files):03d}.txt"
        filter_file.write_text("".join(current_chunk))

        files.append(filter_file)

    return files


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--id_file", 
        type=str, 
        help="Path to text file with PDB ID (one per line). If exists, only those structures will be downloaded."
    )
    parser.add_argument(
        "-o", "--output_dir", 
        type=str, 
        required=True,
        help="Path for downloaded mmcif structures"
    )
    return parser.parse_args()


def run_full_sync(target_dir: str):
    cmd = [
        "rsync", "-rLtz", "--delete", 
        '--out-format=%o %f', 
        "--port=33444", 
        "rsync.rcsb.org::ftp_data/structures/all/mmCIF/", 
        target_dir
    ]
    return execute_and_parse_rsync(cmd)


def count_local_files(directory):
    if not os.path.exists(directory):
        return 0
    return len([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))])


def run_selective_sync(target_dir: str, pdb_filenames: List[str]):

    filter_files = create_rsync_filters(
        pdb_filenames,
        "rsync_filters",
        chunk_size=10000
    )

    for filter_file in filter_files:
        cmd = [
            "rsync",
            "-rLtz",
            "--out-format=%o %f",
            "--port=33444",
            f"--filter=. {filter_file}",
            "rsync.rcsb.org::ftp_data/structures/all/mmCIF/",
            target_dir
        ]
    
    updated, _ = execute_and_parse_rsync(cmd)
            
    return updated


def execute_and_parse_rsync(cmd):
    updated_files = []
    deleted_files = []

    with subprocess.Popen(
    cmd, 
    stdout=subprocess.PIPE, 
    stderr=subprocess.STDOUT, 
    text=True
    ) as process:
        
        for line in process.stdout:
            clean_line = line.strip()
        
            logging.info(f"Rsync: {clean_line}")
            
            parts = clean_line.split(maxsplit=1)
            if len(parts) == 2:
                operation, filename = parts
                if operation == "recv" and filename != '.':
                    updated_files.append(filename)
                elif operation == "del":
                    deleted_files.append(filename)

    return_code = process.wait()
    if return_code != 0:
        logging.error(f"Rsync exited with error (Code: {return_code})")
        sys.exit(return_code)
        
    return updated_files, deleted_files


def save_to_file(filename, file_list):
    with open(filename, "w") as f:
        for item in file_list:
            f.write(f"{item}\n")


def get_obsolete():
    url = "https://files.wwpdb.org/pub/pdb/data/status/obsolete.dat"

    text = requests.get(url).text

    obsolete_map = {}

    for line in text.splitlines():
        if line.startswith("OBSLTE"):
            old_id = line[20:24].strip().lower()
            replacements = line[31:].split() # can be more than 1 

            if replacements:
                obsolete_map[old_id] = [x.lower() for x in replacements]
    return obsolete_map


def main():
    start = time.perf_counter()
    args = parse_arguments()
    output_dir = args.output_dir

    target_dir = Path(output_dir) / 'input_data' / 'pdb_copy_local'
    if not target_dir.exists():
        target_dir.mkdir(parents=True, exist_ok=True)
    
    initial_count = count_local_files(target_dir)
    
    if args.id_file:
        
        id_file_path = args.id_file
        if not os.path.exists(id_file_path):
            logging.error(f"{id_file_path} does not exist")
            sys.exit(1)
        
        pdb_files = []
        with open(id_file_path, "r") as f:
            for line in f:
                pdb_files.append(line.strip().lower() + '.cif.gz')
                
        logging.info(f"Number of PDB IDs in file: {len(pdb_files)}")
                
        updated = run_selective_sync(target_dir, pdb_files)
        diff = set(pdb_files).difference(os.listdir(target_dir))
        obsolete_map = get_obsolete()

        new_ids = []
        for old_id in diff:
            if old_id.lower() in obsolete_map:
                new_ids.extend(obsolete_map[old_id.lower()])
                
        updated_obs = run_selective_sync(target_dir, new_ids)
        logging.info(f"Updated or new structures from obsolete:   {len(updated_obs)}")

    else:
        updated, deleted = run_full_sync(target_dir)

    final_count = count_local_files(target_dir)
        
    save_to_file("updated_files.txt", updated)
    save_to_file("deleted_files.txt", deleted)
    
    
    logging.info(f"Initial number of PDB structures stored locally: {initial_count}")
    logging.info(f"Updated or new structures:   {len(updated)}")
    
    logging.info(f"Deleted structures:     {len(deleted)}")
    logging.info(f"Current number of PDB structures stored locally: {final_count}")
    logging.info(f"Difference among IDs in file and downloaded files:     {set(pdb_files).difference(os.listdir(target_dir))} ")

    elapsed = time.perf_counter() - start
    formatted = str(timedelta(seconds=elapsed))
    logging.info(f"Total time: {formatted}")


if __name__ == "__main__":
    main()
