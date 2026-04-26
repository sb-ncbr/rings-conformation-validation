import csv
import logging
import argparse
import pickle
import shutil
import time
import pandas as pd
from typing import List, Set
from collections import defaultdict
import gemmi
import statistics as st
import math
from multiprocessing import Pool, cpu_count
from pathlib import Path
from HelperModule.Ring import Ring

CPU_COUNT = cpu_count()
NAME = "Ring Coverage"


def _create_output_folder(output_folder: Path):
    try:
        if output_folder.exists():
            shutil.rmtree(str(output_folder.resolve()))

        output_folder.mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output_folder


def run_exe_wrapper(params):
    return run_exe(*params)


def run_exe(path_to_ccp4file: Path, rings_paths: List[str], more_or_equal: bool, closest_voxel: bool):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    try:
        # output of run_calculation is in this format: { benzene: [{ABC_2xyz_0: 3;5}, {...}], oxane: {...} }
        output = run_calculation(path_to_ccp4file, rings_paths, more_or_equal, closest_voxel)
        results_list = []
        for ring_type in output:
            for inner_dict in output[ring_type]:
                for ring_id, coverage in inner_dict.items():
                    results_list.append((ring_id, ring_type, ring_id.split('_')[0], coverage))

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return results_list # List of tuples (ring_id, ring_type, ligand, coverage)


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


def run_calculation(input_density_ccp4: Path, rings_paths: List[str], more_or_equal, closest_voxel):
    try:
        output = defaultdict(list) # { benzene: [{ABC_2xyz_0: 3;5}, {...}], oxane: {...} }
        map = gemmi.read_ccp4_map(str(input_density_ccp4))
        map.setup(float('nan'))

        # calculate the sigma values
        grid_values = []
        for point in map.grid:
            if not math.isnan(point.value):
                grid_values.append(point.value)

        std = st.pstdev(grid_values)
        sigma_lvl = 1.5 * std

        for ring_path in rings_paths:
            ring_pdbfile = gemmi.read_pdb(str(Path(ring_path).resolve()))
            total_atom_count = 0
            covered_atoms_count = 0
            for model in ring_pdbfile:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            total_atom_count = total_atom_count + 1
                            if determine_atom_coverage(atom.pos, map, sigma_lvl, more_or_equal, closest_voxel):
                                covered_atoms_count = covered_atoms_count + 1

            coverage = f'{covered_atoms_count};{total_atom_count}'
            ring_type = ring_path.parents[3].name
            ring_id = ring_path.name.split(".")[0]
            curr_result_record = {ring_id: coverage}
            output[ring_type].append(curr_result_record)

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output


def split_into_subsets(parent_csv_path, root_dir, filename_stem):
    df = pd.read_csv(parent_csv_path)
    for ring_type, subdf in df.groupby("Ring"):
        result_dir = Path(root_dir).resolve() / "validation_data" / ring_type / "el-density-output"
        result_dir.mkdir(parents=True, exist_ok=True)
        res_path = result_dir / f"{ring_type}{filename_stem}.csv"
        subdf.to_csv(res_path, index=False)
        logging.info(f"[{NAME}]: CSV with results has been created at {res_path}")


def main(root_dir: str, input_dir: str, more_or_equal: bool, closest_voxel: bool):

    try:
        params = ''
        if closest_voxel:
            params += "c"
        if more_or_equal:
            params += "m"

        rings: Set[str] = {ring.name.lower() for ring in Ring}
        ccp4_dir = Path(input_dir).resolve() / "ccp4"
        output_path = Path(root_dir).resolve() / "validation_data" / "el-density-output"
        output_path.mkdir(parents=True, exist_ok=True)

        saves_path = Path("cache") / "el_density_saves"
        saves_path.mkdir(parents=True, exist_ok=True)

        filename_stem = f"_params_{params}_analysis_output"
        pkl_path = saves_path / f"{filename_stem}.pkl"
        csv_path = output_path / f"{filename_stem}.csv"

        processed_data_dict = {} # "2xyz": { 'benzene': { "ABC_2xyz_0": "3;6", "ABC_2xyz_1": "6;6" }}, oxane: {...}}
        if pkl_path.is_file():
            with pkl_path.open("rb") as f:
                processed_data_dict = pickle.load(f)

        _create_output_folder(output_path)

        pdb_to_ring_paths_map = map_pdb_to_rings_filepaths(Path(root_dir), ccp4_dir, rings)

        if pdb_to_ring_paths_map is None:
            logging.info(f"[{NAME}]: No files for analysis were found.")
            return
        
        ccp4_filestems_to_process = []
        precomputed_rows = []
        for pdb_id in pdb_to_ring_paths_map: # 2xyz -> {.../benzene/filtered_ligands/ABC/patterns/ABC_2xyz_0.pdb, ..., ...}
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
                for i, result_for_ccp4 in enumerate(p.imap_unordered(run_exe_wrapper, modified_filepaths),1):
                    logging.info(f"[{NAME}]: {i}/{total}")
                    for ring_id, ring_type, ligand, coverage in result_for_ccp4:
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
                        type=str, help='Directory with input files, containing folder ccp4')
    
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
    logging.info(f"[{NAME}]: Total time: {time.perf_counter() - start:.2f}s")