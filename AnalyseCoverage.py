import csv
from datetime import timedelta
from itertools import islice
import logging
import argparse
import numpy as np
import pickle
import shutil
import time
import pandas as pd
from typing import List, Set
from collections import defaultdict
import gemmi
from multiprocessing import Pool, cpu_count
from pathlib import Path
from HelperModule.Ring import Ring

CPU_COUNT = cpu_count()
NAME = "Ring Coverage"


def run_exe_wrapper(params):
    return run_exe(*params)


def run_exe(path_to_ccp4file: Path, rings_paths: List[Path], more_or_equal: bool, closest_voxel: bool):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    try:
        # data["result"] of run_calculation is in this format: { benzene: [{ABC_2xyz_0: 3;5}, {...}], oxane: {...} }
        output = run_calculation(path_to_ccp4file, rings_paths, more_or_equal, closest_voxel)
        results_list = []
        for ring_type, ring_records in output["data"].items():
            for inner_dict in ring_records:
                for ring_id, coverage in inner_dict.items():
                    results_list.append((ring_id, ring_type, ring_id.split('_')[0], coverage))

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return {"data": results_list, "metadata": output["metadata"]} # results: List of tuples (ring_id, ring_type, ligand, coverage)


def map_pdb_to_rings_filepaths(rootdir: Path, ccp4_dir: Path, rings: Set[str]):
    try:
        res = defaultdict(set)
        all_ccp4_files = ccp4_dir.glob('**/*')
        pdb_ids_for_which_ccp4_is_available = {x.stem.removesuffix('.ccp4') for x in all_ccp4_files}

        if len(pdb_ids_for_which_ccp4_is_available) == 0:
            return None
        base = Path(rootdir) / "validation_data"
        for ring_type in rings:
            for f in (base / ring_type / "filtered_ligands").rglob("*"):
                if f.is_file():
                    pdb_id = f.stem.split('_')[1]
                    if pdb_id in pdb_ids_for_which_ccp4_is_available:
                        res[pdb_id].add(f)
        logging.info(f"[{NAME}]: There are {len(res)} pdb structures and {sum(len(v) for v in res.values())} rings with corresponding CCP4 file "
                     f"available.")

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return res


# compare intensity, corresponding to the given position, to the threshold for isosurface (MORE vs MORE OR EQUAL)
def determine_atom_coverage(pos, map, sigma_lvl, more_or_equal, closest_voxel):
    try:
        if more_or_equal:
            return get_intensity(pos, map, closest_voxel) >= sigma_lvl
        return get_intensity(pos, map, closest_voxel) > sigma_lvl
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


# get intensity corresponding to the given position (trilinear interpolation vs itensity of the closest voxel)
def get_intensity(pos, map, closest_voxel):
    try:
        if closest_voxel:
            return map.grid.get_nearest_point(pos).value
        return map.grid.interpolate_value(pos)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


# only of pdb files with 5chars ligand name, coords are shifted 2 columns to the right
def get_coverage_from_nonstd_pdb(dens_map, file_path, atom_count, sigma_lvl, more_or_equal, closest_voxel):
    covered_atoms_count = 0
    with open(file_path, "r") as f:
        next(f)
        lines = list(islice(f, atom_count))
        
        for line in lines:
            x = float(line[32:40])  # 31+2 to 38+2
            y = float(line[40:48])  # 39+2 to 46+2
            z = float(line[48:56])  # 47+2 to 54+2
 
            pos = gemmi.Position(x, y, z)
            if determine_atom_coverage(pos, dens_map, sigma_lvl, more_or_equal, closest_voxel):
                covered_atoms_count += 1
                            
    return covered_atoms_count, atom_count


def get_coverage_using_gemmi(dens_map, ring_path: Path, input_density_ccp4: Path, sigma_lvl, more_or_equal, closest_voxel):
    # we process only one ring in pdb format, so there is only one model, one chain and one residue
    ring_pdbfile = gemmi.read_pdb(str(ring_path.resolve()))
    model = ring_pdbfile[0]
    chain = model[0]
    res = chain[0]
    covered_atoms_count = 0
    total_atom_count = 0
    for atom in res:
        total_atom_count += 1
        if determine_atom_coverage(atom.pos, dens_map, sigma_lvl, more_or_equal, closest_voxel):
            covered_atoms_count += 1
    return covered_atoms_count, total_atom_count


def run_calculation(input_density_ccp4: Path, rings_paths: List[Path], more_or_equal, closest_voxel):
    try:
        start = time.perf_counter()
        output = defaultdict(list) # { benzene: [{ABC_2xyz_0: 3;5}, {...}], oxane: {...} }
        dens_map = gemmi.read_ccp4_map(str(input_density_ccp4))
        dens_map.setup(float('nan'))

        arr = np.array(dens_map.grid, copy=False)
        arr = arr[~np.isnan(arr)]
        std = arr.std()
        sigma_lvl = 1.5 * std

        for ring_path in rings_paths:
            
            ring_type = ring_path.parents[3].name
            atom_count = Ring[ring_type.upper()].atom_number

            ligand = ring_path.parents[1].name
            if len(ligand) == 5:
                covered_atoms_count, total_atom_count = get_coverage_from_nonstd_pdb(dens_map, ring_path,
                                                                                     atom_count, sigma_lvl,
                                                                                     more_or_equal, closest_voxel)
            else:
                covered_atoms_count, total_atom_count = get_coverage_using_gemmi(dens_map, ring_path,
                                                                                 input_density_ccp4,
                                                                                 sigma_lvl,
                                                                                 more_or_equal,
                                                                                 closest_voxel)

            coverage = f'{covered_atoms_count};{total_atom_count}'
            ring_id = ring_path.name.split(".")[0]
            curr_result_record = {ring_id: coverage}
            output[ring_type].append(curr_result_record)
        total_time = time.perf_counter() - start

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return {"data": output,
            "metadata": {
                "ccp4_name": input_density_ccp4.stem.removesuffix('.ccp4'),
                "n_rings": len(rings_paths),
                "total_time": total_time,
            }}


def split_into_subsets(parent_csv_path, root_dir, filename_stem):
    df = pd.read_csv(parent_csv_path)
    for ring_type, subdf in df.groupby("Ring"):
        result_dir = Path(root_dir).resolve() / "validation_data" / ring_type / "el-density-output"
        result_dir.mkdir(parents=True, exist_ok=True)
        res_path = result_dir / f"{ring_type}{filename_stem}.csv"
        subdf.to_csv(res_path, index=False)
        logging.info(f"[{NAME}]: Exporting to {res_path}")


def main(root_dir: str, inputdir: str, more_or_equal: bool, closest_voxel: bool):
    logging.info(f"[{NAME}]: Starting...")
    try:
        params = ''
        if closest_voxel:
            params += "c"
        if more_or_equal:
            params += "m"

        rings: Set[str] = {ring.name.lower() for ring in Ring}
        output_path = Path(root_dir).resolve() / "validation_data" / "el-density-output"
        output_path.mkdir(parents=True, exist_ok=True)
        ccp4_dir = Path(inputdir).resolve()
        saves_path = Path("cache") / "el_density_saves"
        saves_path.mkdir(parents=True, exist_ok=True)

        filename_stem = f"_params_{params}_analysis_output"
        pkl_path = saves_path / f"{filename_stem}.pkl"
        csv_path = output_path / f"{filename_stem}.csv"

        processed_data_dict = {} # "2xyz": { 'benzene': { "ABC_2xyz_0": "3;6", "ABC_2xyz_1": "6;6" }}, oxane: {...}}
        if pkl_path.is_file():
            with pkl_path.open("rb") as f:
                processed_data_dict = pickle.load(f)

        output_path.mkdir(parents=True, exist_ok=True)
        
        logging.info(f"[{NAME}]: Scanning the directory with ccp4 files...")
        pdb_to_ring_paths_map = map_pdb_to_rings_filepaths(Path(root_dir), ccp4_dir, rings)

        if pdb_to_ring_paths_map is None:
            logging.info(f"[{NAME}]: No files for analysis were found.")
            return
        
        ccp4_filestems_to_process = []
        precomputed_rows = []
        for pdb_id in pdb_to_ring_paths_map: # 2xyz -> {Path(.../benzene/filtered_ligands/ABC/patterns/ABC_2xyz_0.pdb), Path(...), ...}
            if pdb_id in processed_data_dict:
                for ring_type, ring_ids in processed_data_dict[pdb_id].items():
                    for ring_id, coverage in ring_ids.items():
                        ligand = ring_id.split('_')[0]
                        precomputed_rows.append((ring_id, ring_type, ligand, coverage))
            else:    
                ccp4_filestems_to_process.append(pdb_id)

        
        modified_filepaths = [(ccp4_dir / f"{filestem}.ccp4.gz",
                               pdb_to_ring_paths_map[filestem],
                               more_or_equal,
                               closest_voxel) for filestem in ccp4_filestems_to_process]

        with open(csv_path, mode='w', newline='', buffering=1) as f:
            w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            header = ('Id', 'Ring', 'Ligand', 'Coverage')
            w.writerow(header)

            if len(precomputed_rows) != 0:
                logging.info(f"[{NAME}]: Writing precomputed data for {len(precomputed_rows)} rings.")
                w.writerows(precomputed_rows)
                logging.info(f"[{NAME}]: Done.")

            if len(ccp4_filestems_to_process) == 0:
                logging.info(f"[{NAME}]: No files to process. Computation will not start.")
                split_into_subsets(csv_path, root_dir, filename_stem)
                return

            logging.info(f"[{NAME}]: Running electron density coverage analysis on CPU count: {CPU_COUNT}")
            with Pool(int(CPU_COUNT)) as p:
                logging.info(f"[{NAME}]: Starting analysis for {len(ccp4_filestems_to_process)} ccp4 files...")
                total = len(modified_filepaths)
                # List of tuples (ring_id, ring_type, ligand, coverage)
                for i, output in enumerate(p.imap_unordered(run_exe_wrapper, modified_filepaths),1):
                    logging.info(f"[{NAME}]: {i}/{total} | {output['metadata']['ccp4_name']} | rings: {output['metadata']['n_rings']} | {output['metadata']['total_time']:.2f}s")
                    for ring_id, ring_type, ligand, coverage in output["data"]:
                        w.writerow((ring_id, ring_type, ligand, coverage))
                        pdb_id = ring_id.split('_')[1]
                        processed_data_dict.setdefault(pdb_id, {}).setdefault(ring_type, {})[ring_id] = coverage

                logging.info(f"[{NAME}]: Finished analysis for {len(ccp4_filestems_to_process)} ccp4 files.")

        split_into_subsets(csv_path, root_dir, filename_stem)

        with pkl_path.open("wb") as f:
            pickle.dump(processed_data_dict, f)

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


if __name__ == '__main__':
    start = time.perf_counter()
    parser = argparse.ArgumentParser(description='ED coverage analysis. Output is two numbers: first is the number of '
                                                  'covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('rootdir', type=str,
                        help='Root directory of the result data (<ROOTDIR>/validation_data/etc)')
    parser.add_argument('input_dir',
                        type=str, help='Directory with ccp4 files')
    
    parser.add_argument('-m', '--more_or_equal',
                        action='store_true', help='Atom is considered to be covered by the electron density when the '
                                                  'corresponding intensity is MORE OR EQUAL to the threshold for the '
                                                  'isosurface')
    parser.add_argument('-c', '--closest_voxel',
                        action='store_true', help='Instead of trilinear interpolation, the intensity of the closest '
                                                  'voxel is used')
    
    args = parser.parse_args()
        
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    main(args.rootdir, args.input_dir, args.more_or_equal, args.closest_voxel)

    elapsed = time.perf_counter() - start
    formatted = str(timedelta(seconds=elapsed))
    logging.info(f"[{NAME}]: Total time: {formatted}")