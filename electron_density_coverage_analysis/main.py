import csv
import logging
from multiprocessing import Pool, cpu_count
import argparse
from pathlib import Path
import shutil
from electron_density_coverage_analysis import run_as_function

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


def run_exe(ligand_filepath: Path, ccp4_dir_path: Path, arguments: argparse.Namespace):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        )
    try:
        pq_pdb_name = ligand_filepath.name.split(".")[0]
        pdb_id = pq_pdb_name.split('_')[1]
        residue_id = ligand_filepath.parent.parent.name
        ccp4_filepath = (ccp4_dir_path / (pdb_id + '.ccp4.gz')).resolve()
        arguments.input_cycle_pdb = str(ligand_filepath.resolve())
        arguments.input_density_ccp4 = str(ccp4_filepath)

        logging.info(f"Analysing file: {arguments.input_cycle_pdb}...")
        output = run_as_function(arguments)
        result = (pq_pdb_name, residue_id, output)

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

        ring_types = ['cyclohexane', 'cyclopentane', 'benzene']

        for ring_type in ring_types:
            path_to_output = Path(args.rootdir).resolve() / "validation_data" / ring_type / "el-density-output"

            _create_output_folder(path_to_output)

            filepaths = get_filepaths(Path(args.rootdir), Path(args.ccp4_dir), ring_type)
            if len(filepaths) == 0:
                logging.info(f"No files for analysis found for ring {ring_type}")
                continue

            cvs_filename = ring_type + '_params_' + params + '_analysis_output.csv'
            with open(path_to_output / cvs_filename, mode='w', newline='') as f:
                w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

                with Pool(int(CPU_COUNT)) as p:
                    modified_filepaths = [(f, Path(args.ccp4_dir), arguments) for f in filepaths]
                    logging.info(f"[{ring_type.capitalize()}]: Starting analysis for {len(filepaths)} files...")
                    rows = p.starmap(run_exe, modified_filepaths)
                    w.writerows(rows)
                    logging.info(f"[{ring_type.capitalize()}]: Finished analysis for {len(filepaths)} files.")

    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


def main():
    parser = argparse.ArgumentParser(description='ED coverage analysis')
    parser.add_argument('rootdir', type=str,
                        help='Root directory of the result data (<ROOTDIR>/validation_data/etc)')
    parser.add_argument('ccp4_dir',
                        type=str, help='Directory with ccp4 files')

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
