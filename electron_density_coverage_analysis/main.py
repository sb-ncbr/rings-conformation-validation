import csv
import logging
from multiprocessing import Pool, cpu_count
import argparse
from pathlib import Path
import pickle
import shutil
import gemmi
import statistics as st
import math

CPU_COUNT = cpu_count()


def _create_output_folder(output_folder: Path):
    try:
        if output_folder.exists():
            shutil.rmtree(str(output_folder.resolve()))

        output_folder.mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output_folder


def process_args(args: argparse.Namespace):
    try:
        a = argparse.Namespace()
        a.s = True
        a.d = False
        a.closest_voxel = False
        a.more_or_equal = False
        if args.closest_voxel:
            a.closest_voxel = True
        if args.more_or_equal:
            a.more_or_equal = True
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return a


def run_exe(ring_path: Path, ccp4_dir_path: Path, arguments: argparse.Namespace):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    try:
        pq_pdb_name = ring_path.name.split(".")[0]
        pdb_id = pq_pdb_name.split('_')[1]
        ligand_id = ring_path.parent.parent.name
        ccp4_filepath = (ccp4_dir_path / (pdb_id + '.ccp4.gz')).resolve()
        arguments.input_cycle_pdb = str(ring_path.resolve())
        arguments.input_density_ccp4 = str(ccp4_filepath)

        logging.info(f"Analysing file: {arguments.input_cycle_pdb}...")
        output = run_as_function(arguments)
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


def run_analysis(args: argparse.Namespace):

    try:
        arguments = process_args(args)
        params = ''
        if arguments.closest_voxel:
            params = params + "c"
        if arguments.more_or_equal:
            params = params + "m"

        ring_types = ['cyclohexane', 'cyclopentane', 'benzene', 'oxane', 'oxolane']
        ccp4_dir = Path(args.input_dir) / "ccp4"

        for ring_type in ring_types:
            path_to_output = Path(args.rootdir).resolve() / "validation_data" / ring_type / "el-density-output"

            saves_path = Path(args.input_dir) / "el_density_saves"
            saves_path.mkdir(parents=True, exist_ok=True)

            filename_stem = f"{ring_type}_params_{params}_analysis_output"
            pkl_path = saves_path / f"{filename_stem}.pkl"
            csv_path = path_to_output / f"{filename_stem}.csv"

            # e.g. {"CVM_4iut_0": "4;6"} for simple mode
            processed_data_dict = {}
            if pkl_path.is_file():
                with pkl_path.open("rb") as f:
                    processed_data_dict = pickle.load(f)

            _create_output_folder(path_to_output)

            filepaths = get_filepaths(Path(args.rootdir), ccp4_dir, ring_type)
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
            
            modified_filepaths = [(f, ccp4_dir, arguments) for f in files_to_process]

            with open(csv_path, mode='w', newline='', buffering=1) as f:
                w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

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


def run_as_function(args: argparse.Namespace):
    return run_calculation(args)


# compare intensity, corresponding to the given position, to the threshold for isosurface (MORE vs MORE OR EQUAL)
def determine_atom_coverage(pos, map, sigma_lvl, args: argparse.Namespace):
    try:
        if args.more_or_equal:
            return get_intensity(pos, map, args) >= sigma_lvl
        return get_intensity(pos, map, args) > sigma_lvl
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


# get intensity corresponding to the given position (trilinear interpolation vs itensity of the closest voxel)
def get_intensity(pos, map, args: argparse.Namespace):
    try:
        if args.closest_voxel:
            return map.grid.get_nearest_point(pos).value
        return map.grid.interpolate_value(pos)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


def run_calculation(args: argparse.Namespace):
    try:
        output = None
        str = gemmi.read_pdb(args.input_cycle_pdb)
        map = gemmi.read_ccp4_map(args.input_density_ccp4)
        map.setup(float('nan'))

        # calculate the sigma values
        grid_values = []
        for point in map.grid:
            if not math.isnan(point.value):
                grid_values.append(point.value)

        std = st.pstdev(grid_values)
        sigma_lvl = 1.5 * std

        if args.s:
            total_atom_count = 0
            covered_atoms_count = 0
            for model in str:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            total_atom_count = total_atom_count + 1
                            if determine_atom_coverage(atom.pos, map, sigma_lvl, args):
                                covered_atoms_count = covered_atoms_count + 1

            output = f'{covered_atoms_count};{total_atom_count}'

        if args.d:
            output = []
            for model in str:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            if determine_atom_coverage(atom.pos, map, sigma_lvl, args):
                                output.append(f'{atom.serial};y;')
                            else:
                                output.append(f'{atom.serial};n;')
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output


def main():
    parser = argparse.ArgumentParser(description='ED coverage analysis')
    parser.add_argument('rootdir', type=str,
                        help='Root directory of the result data (<ROOTDIR>/validation_data/etc)')
    parser.add_argument('input_dir',
                        type=str, help='Directory with input files, containing folder ccp4')

    parser.add_argument('-s',
                        action='store_true', help='Simple mode - output is two numbers: first is the number of '
                                                  'covered atoms, the second is the total number of atoms in a cycle')
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
    logging.info(f"Running electron density coverage analysis on CPU count: {CPU_COUNT}")
    run_analysis(args)


if __name__ == '__main__':
    main()