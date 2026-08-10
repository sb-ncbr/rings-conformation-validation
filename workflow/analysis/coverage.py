import csv
import logging
import numpy as np
import time
import pandas as pd
from typing import List, Set
from collections import defaultdict
import gemmi
from multiprocessing import Pool, cpu_count
from pathlib import Path
from workflow.models.Ring import Ring
from workflow.utils.helpers import extract_extended_pdb_code, get_old_pdb_id

EL_DENSITY_OUTPUT_DIR = "el-density-output"
CPU_COUNT = cpu_count()

logger = logging.getLogger(__name__)


def run_exe_wrapper(params):
    return run_exe(*params)


def run_exe(path_to_ccp4file: Path, rings_paths: List[Path], more_or_equal: bool, closest_voxel: bool):
    try:
        # data["result"] of run_calculation is in this format: { benzene: [{ABC_pdb_00002xyz_0: 3;5}, {...}], oxane: {...} }
        output = run_calculation(path_to_ccp4file, rings_paths, more_or_equal, closest_voxel)
        if output is None:
            return {"data": None, "metadata": path_to_ccp4file}
        results_list = []
        for ring_type, ring_records in output["data"].items():
            for inner_dict in ring_records:
                for ring_id, coverage in inner_dict.items():
                    results_list.append((ring_id, ring_type, ring_id.split('_')[0], coverage))

    except Exception as e:
        logger.error(e, stack_info=True, exc_info=True)
    return {"data": results_list, "metadata": output["metadata"]} # results: List of tuples (ring_id, ring_type, ligand, coverage)


def map_pdb_to_rings_filepaths(main_dir: Path, ccp4_dir: Path, rings: Set[str]):
    try:
        res = defaultdict(set)
        all_ccp4_files = ccp4_dir.glob('**/*')
        pdb_ids_for_which_ccp4_is_available = {x.stem.removesuffix('.ccp4') for x in all_ccp4_files}

        if len(pdb_ids_for_which_ccp4_is_available) == 0:
            return None
        for ring_type in rings:
            for f in (main_dir / ring_type / "filtered_ligands").rglob("*"):
                if f.is_file():
                    # ABC_pdb_00002xyz_0.cif -> pdb_00002xyz
        
                    extended_pdb = extract_extended_pdb_code(f.stem)
                    pdb_id = get_old_pdb_id(extended_pdb)
                    
                    if pdb_id in pdb_ids_for_which_ccp4_is_available:
                        res[extended_pdb].add(f)
        logger.info(f"There are {len(res)} structures and {sum(len(v) for v in res.values())} rings with corresponding CCP4 file "
                     f"available.")

    except Exception as e:
        logger.error(e, stack_info=True, exc_info=True)
    return res


# compare intensity, corresponding to the given position, to the threshold for isosurface (MORE vs MORE OR EQUAL)
def determine_atom_coverage(pos, map, sigma_lvl, more_or_equal, closest_voxel):
    try:
        if more_or_equal:
            return get_intensity(pos, map, closest_voxel) >= sigma_lvl
        return get_intensity(pos, map, closest_voxel) > sigma_lvl
    except Exception as e:
        logger.error(e, stack_info=True, exc_info=True)


# get intensity corresponding to the given position (trilinear interpolation vs itensity of the closest voxel)
def get_intensity(pos, map, closest_voxel):
    try:
        if closest_voxel:
            return map.grid.get_nearest_point(pos).value
        return map.grid.interpolate_value(pos)
    except Exception as e:
        logger.error(e, stack_info=True, exc_info=True)


def get_coverage(dens_map, ring_path: Path, sigma_lvl, more_or_equal, closest_voxel):
    # we process only one ring in pdb format, so there is only one model, one chain and one residue
    ring_structure = gemmi.read_structure(str(ring_path.resolve()))
    model = ring_structure[0]
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
        output = defaultdict(list) # { benzene: [{ABC_pdb_00002xyz_0: 3;5}, {...}], oxane: {...} }
        dens_map = gemmi.read_ccp4_map(str(input_density_ccp4))
        dens_map.setup(float('nan'))

        arr = np.array(dens_map.grid, copy=False)
        arr = arr[~np.isnan(arr)]
        std = arr.std()
        sigma_lvl = 1.5 * std

        for ring_path in rings_paths:
            
            ring_type = ring_path.parents[3].name
            covered_atoms_count, total_atom_count = get_coverage(dens_map, ring_path,
                                                                 sigma_lvl,
                                                                 more_or_equal,
                                                                 closest_voxel)

            coverage = f'{covered_atoms_count};{total_atom_count}'
            ring_id = ring_path.name.split(".")[0]
            curr_result_record = {ring_id: coverage}
            output[ring_type].append(curr_result_record)
        total_time = time.perf_counter() - start

    except Exception as e:
        logger.exception(f"Failed processing {input_density_ccp4}")
        with open("failed_density_maps.txt", "a") as f:
            f.write(f"{input_density_ccp4}\n")
        input_density_ccp4.unlink()
        return None

    return {"data": output,
            "metadata": {
                "ccp4_name": input_density_ccp4.stem.removesuffix('.ccp4'),
                "n_rings": len(rings_paths),
                "total_time": total_time,
            }}


def split_into_subsets(parent_csv_path, main_dir, filename_stem):
    df = pd.read_csv(parent_csv_path)
    for ring_type, subdf in df.groupby("Ring"):
        result_dir = main_dir/ ring_type / EL_DENSITY_OUTPUT_DIR
        result_dir.mkdir(parents=True, exist_ok=True)
        res_path = result_dir / f"{ring_type}{filename_stem}.csv"
        subdf.to_csv(res_path, index=False)
        logger.info(f"Exporting to {res_path}")


def analyse_coverage(main_dir: str, ccp4_dir: str, more_or_equal: bool, closest_voxel: bool):
    logger.info(f"Starting...")
    try:
        params = ''
        if closest_voxel:
            params += "c"
        if more_or_equal:
            params += "m"

        rings: Set[str] = {ring.name.lower() for ring in Ring}
        output_path = main_dir / EL_DENSITY_OUTPUT_DIR
        output_path.mkdir(parents=True, exist_ok=True)

        filename_stem = f"_params_{params}_analysis_output"
        csv_path = output_path / f"{filename_stem}.csv"


        output_path.mkdir(parents=True, exist_ok=True)
        processed_data_dict = defaultdict(lambda: defaultdict(dict)) # "pdb_00002xyz": { 'benzene': { "ABC_pdb_00002xyz_0": "3;6", "ABC_pdb_00002xyz_1": "6;6" }}, oxane: {...}}

        if csv_path.exists():
            df = pd.read_csv(csv_path, header=0)

            for ring_id, ring_type, ligand, coverage in df.itertuples(index=False, name=None):

                # ABC_pdb_00002xyz_0 -> pdb_00002xyz
                extended_pdb_code = ring_id[len(ligand) + 1 :].rsplit("_", 1)[0]
                processed_data_dict[extended_pdb_code][ring_type][ring_id] = coverage
        
        logger.info(f"Scanning the directory with ccp4 files...")
        pdb_to_ring_paths_map = map_pdb_to_rings_filepaths(main_dir, ccp4_dir, rings)

        if pdb_to_ring_paths_map is None:
            logger.info(f"No files for analysis were found.")
            return
        ccp4_filestems_to_process = []
        precomputed_rows = []
        for ext_pdb_id in pdb_to_ring_paths_map: # pdb_00002xyz -> {Path(.../benzene/filtered_ligands/ABC/patterns/ABC_pdb_00002xyz_0.pdb), Path(...), ...}
            if ext_pdb_id in processed_data_dict:
                for ring_type, ring_ids in processed_data_dict[ext_pdb_id].items():
                    for ring_id, coverage in ring_ids.items():
                        ligand = ring_id.split('_')[0]
                        precomputed_rows.append((ring_id, ring_type, ligand, coverage))
            else:
                ccp4_filestems_to_process.append(get_old_pdb_id(ext_pdb_id))

        
        modified_filepaths = [(ccp4_dir / f"{filestem}.ccp4.gz",
                               pdb_to_ring_paths_map[f"pdb_0000{filestem}"],
                               more_or_equal,
                               closest_voxel) for filestem in ccp4_filestems_to_process]

        with open(csv_path, mode='w', newline='', buffering=1) as f:
            w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            header = ('Id', 'Ring', 'Ligand', 'Coverage')
            w.writerow(header)

            if len(precomputed_rows) != 0:
                logger.info(f"Writing precomputed data for {len(precomputed_rows)} rings.")
                w.writerows(precomputed_rows)
                logger.info(f"Done.")

            if len(ccp4_filestems_to_process) == 0:
                logger.info(f"No files to process. Computation will not start.")
                split_into_subsets(csv_path, main_dir, filename_stem)
                return

            logger.info(f"Running electron density coverage analysis on CPU count: {CPU_COUNT}")
            with Pool(int(CPU_COUNT)) as p:
                logger.info(f"Starting analysis for {len(ccp4_filestems_to_process)} ccp4 files...")
                total = len(modified_filepaths)
                # List of tuples (ring_id, ring_type, ligand, coverage)
                for i, output in enumerate(p.imap_unordered(run_exe_wrapper, modified_filepaths),1):
                    if output["data"] is None:
                        logger.info(f"{i}/{total} {output["metadata"]}")
                        continue
                    logger.info(f"{i}/{total} | {output['metadata']['ccp4_name']} | rings: {output['metadata']['n_rings']} | {output['metadata']['total_time']:.2f}s")
                    for ring_id, ring_type, ligand, coverage in output["data"]:
                        w.writerow((ring_id, ring_type, ligand, coverage))

                        # ring id is of type A1CS3_pdb_00007ibg_0
                        extended_pdb_code = extract_extended_pdb_code(ring_id)
      
                        processed_data_dict.setdefault(extended_pdb_code, {}).setdefault(ring_type, {})[ring_id] = coverage

                logger.info(f"Finished analysis for {len(ccp4_filestems_to_process)} ccp4 files.")

        split_into_subsets(csv_path, main_dir, filename_stem)

    except Exception as e:
        logger.error(e, stack_info=True, exc_info=True)


    
#     parser.add_argument('-m', '--more_or_equal',
#                         action='store_true', help='Atom is considered to be covered by the electron density when the '
#                                                   'corresponding intensity is MORE OR EQUAL to the threshold for the '
#                                                   'isosurface')
#     parser.add_argument('-c', '--closest_voxel',
#                         action='store_true', help='Instead of trilinear interpolation, the intensity of the closest '
#                                                   'voxel is used')
    