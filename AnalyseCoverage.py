import csv
import logging
import argparse
import pickle
import shutil
import sys
from typing import List
import gemmi
import statistics as st
import math
from multiprocessing import Pool, cpu_count
from pathlib import Path
from HelperModule.Ring import Ring

CPU_COUNT = cpu_count()


def _create_output_folder(output_folder: Path):
    try:
        if output_folder.exists():
            shutil.rmtree(str(output_folder.resolve()))

        output_folder.mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output_folder


def run_exe(ring_path: Path, ccp4_dir_path: Path, more_or_equal: bool, closest_voxel: bool):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    try:
        pq_pdb_name = ring_path.name.split(".")[0]
        pdb_id = pq_pdb_name.split('_')[1]
        ligand_id = ring_path.parent.parent.name
        input_density_ccp4 = str((ccp4_dir_path / (pdb_id + '.ccp4.gz')).resolve())
        input_cycle_pdb = str(ring_path.resolve())

        logging.info(f"Analysing file: {input_cycle_pdb}...")
        output = run_calculation(input_density_ccp4, input_cycle_pdb, more_or_equal, closest_voxel)
        result = (pq_pdb_name, ligand_id, output)

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return result


def get_filepaths(rootdir: Path, ccp4_dir: Path, ring_type: str):
    try:
        l = []

        all_ccp4_files = ccp4_dir.glob('**/*')
        pdb_ids_for_which_ccp4_is_available = [x.stem.removesuffix('.ccp4') for x in all_ccp4_files]
        for f in Path(rootdir / 'validation_data' / ring_type / 'filtered_ligands').rglob("*"):

            if f.is_file():
                # get path
                stem = f.stem
                pdb_id = stem.split('_')[1]
                if pdb_id in pdb_ids_for_which_ccp4_is_available:
                    l.append(f)
        logging.info(f"[{ring_type.capitalize()}]: There are {len(l)} PDB structures with corresponding CCP4 file "
                     f"available.")
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return l


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


def run_calculation(input_density_ccp4, input_cycle_pdb, more_or_equal, closest_voxel):
    try:
        output = None
        str = gemmi.read_pdb(input_cycle_pdb)
        map = gemmi.read_ccp4_map(input_density_ccp4)
        map.setup(float('nan'))

        # calculate the sigma values
        grid_values = []
        for point in map.grid:
            if not math.isnan(point.value):
                grid_values.append(point.value)

        std = st.pstdev(grid_values)
        sigma_lvl = 1.5 * std

        total_atom_count = 0
        covered_atoms_count = 0
        for model in str:
            for chain in model:
                for res in chain:
                    for atom in res:
                        total_atom_count = total_atom_count + 1
                        if determine_atom_coverage(atom.pos, map, sigma_lvl, more_or_equal, closest_voxel):
                            covered_atoms_count = covered_atoms_count + 1

        output = f'{covered_atoms_count};{total_atom_count}'

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output


def main(output_dir: str, input_dir: str, rings: List[str], more_or_equal: bool, closest_voxel: bool):

    try:
        params = ''
        if closest_voxel:
            params += "c"
        if more_or_equal:
            params += "m"

        ccp4_dir = Path(input_dir) / "ccp4"

        for ring_type in rings:
            path_to_output = Path(output_dir).resolve() / "validation_data" / ring_type / "el-density-output"

            saves_path = Path(input_dir) / "el_density_saves"
            saves_path.mkdir(parents=True, exist_ok=True)

            filename_stem = f"{ring_type}_params_{params}_analysis_output"
            pkl_path = saves_path / f"{filename_stem}.pkl"
            csv_path = path_to_output / f"{filename_stem}.csv"

            # e.g. {"CVM_4iut_0": "4;6"}
            processed_data_dict = {}
            if pkl_path.is_file():
                with pkl_path.open("rb") as f:
                    processed_data_dict = pickle.load(f)

            _create_output_folder(path_to_output)

            filepaths = get_filepaths(Path(output_dir), ccp4_dir, ring_type)
            if len(filepaths) == 0:
                logging.info(f"No files for analysis found for ring {ring_type}")
                continue
            
            files_to_process = []
            precomputed_rows = []
            for f in filepaths:
                key = f.stem
                if key in processed_data_dict:
                    precomputed_rows.append((key, *processed_data_dict[key]))
                    continue
                    
                files_to_process.append(f)

            
            modified_filepaths = [(f, ccp4_dir, more_or_equal, closest_voxel) for f in files_to_process]

            with open(csv_path, mode='w', newline='', buffering=1) as f:
                w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

                if len(precomputed_rows) != 0:
                    logging.info(f"[{ring_type.capitalize()}]: Writing precomputed data for {len(precomputed_rows)} rings.")
                    w.writerows(precomputed_rows)
                logging.info(f"[{ring_type.capitalize()}]: Done.")

                if len(files_to_process) == 0:
                    logging.info(f"[{ring_type.capitalize()}]: No files to process. Computation will not start.")
                    continue

                with Pool(int(CPU_COUNT)) as p:
                    logging.info(f"[{ring_type.capitalize()}]: Starting analysis for {len(files_to_process)} files...")
                    for ring_id, ligand, result in p.starmap(run_exe, modified_filepaths):
                        w.writerow((ring_id, ligand, result))
                        processed_data_dict[ring_id] = (ligand, result)

                    logging.info(f"[{ring_type.capitalize()}]: Finished analysis for {len(files_to_process)} files.")

            with pkl_path.open("wb") as f:
                pickle.dump(processed_data_dict, f)

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ED coverage analysis. Output is two numbers: first is the number of '
                                                  'covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('rootdir', type=str,
                        help='Root directory of the result data (<ROOTDIR>/validation_data/etc)')
    parser.add_argument('input_dir',
                        type=str, help='Directory with input files, containing folder ccp4')

    parser.add_argument("-r", "--rings", nargs="+", type=str,
                        choices=[r.name for r in Ring], help="Choose ring type(s)")
    
    parser.add_argument('-m', '--more_or_equal',
                        action='store_true', help='Atom is considered to be covered by the electron density when the '
                                                  'corresponding intensity is MORE OR EQUAL to the threshold for the '
                                                  'isosurface')
    parser.add_argument('-c', '--closest_voxel',
                        action='store_true', help='Instead of trilinear interpolation, the intensity of the closest '
                                                  'voxel is used')
    
    args = parser.parse_args()
    if args.rings is None:
        selected_rings = [r.name for r in Ring]
    else:
        selected_rings = args.rings

    for ring in selected_rings:
        if ring.upper() not in Ring.__members__.keys():
            logging.error(
                f"Ring {ring} is not a valid Ring. Currently supported: {[e.name for e in Ring]} Exiting..."
            )
            sys.exit(1)
        
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    logging.info(f"Running electron density coverage analysis on CPU count: {CPU_COUNT}")
    main(args.rootdir, args.input_dir, selected_rings, args.more_or_equal, args.closest_voxel)