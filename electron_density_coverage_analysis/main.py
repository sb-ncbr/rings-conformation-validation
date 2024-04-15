import csv
import logging
from multiprocessing import Pool, cpu_count
import os
import urllib.request as r
import argparse
from pathlib import Path
import shutil
from electron_density_coverage_analysis import run_as_function

####
# CHANGE PATH
CCP4_DIR = Path('D:\\ccp4')
####
EXE = Path('./electron_density_coverage_analysis.py')

OUTPUT_DIR = Path('./output')
# NO_CCP4_AVAILABLE_FILE = Path(OUTPUT_DIR / 'no_ccp4_pdb_ids.txt')
CPU_COUNT = cpu_count() / 2


#TODO: dangerous?
def _create_output_folder():
    try:
        if OUTPUT_DIR.exists():
            shutil.rmtree(str(OUTPUT_DIR.resolve()))

        OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    
    return OUTPUT_DIR


def _download_single_file(data: tuple):
    link, path, pdb_id = data
    try:
        r.urlretrieve(link, path)
    except:
        return pdb_id


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


def run_exe(ligand_filepath: Path, arguments: argparse.Namespace):
    try:
        PQ_pdb_name = ligand_filepath.name.split(".")[0]
        pdb_id = PQ_pdb_name.split('_')[1]
        residue_id = ligand_filepath.parent.parent.name
        ccp4_filepath = str(Path(CCP4_DIR / f'{pdb_id}.ccp4').resolve())
        
        arguments.input_cycle_pdb = str(ligand_filepath.resolve())
        arguments.input_density_ccp4 = ccp4_filepath
        logging.info(f"Running exe for {pdb_id}, ligand {residue_id}")
        output = run_as_function(arguments)
        result = (PQ_pdb_name, residue_id, output)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return result


def get_filepathes(rootdir: Path, ring_type: str):
    try:
        l = []
        logging.info(f"Reading path for ccp4 files...")
        all_ccp4_files = CCP4_DIR.glob('**/*')
        logging.info(f"Getting pdb files for which ccp4 is available...")
        pdb_ids_for_which_ccp4_is_available = [x.stem for x in all_ccp4_files]
        for f in Path(rootdir / ring_type / 'filtered_ligands').rglob("*"):
            if f.is_file():
                # logging.info(f"Getting file {f.stem}")
                # get path
                stem = f.stem
                pdb_id = stem.split('_')[1]
                if pdb_id in pdb_ids_for_which_ccp4_is_available:
                    l.append(f)
                    # logging.info(f"Length of list of files is {len(l)}")
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

        # ring_types = ['cyclohexane', 'cyclopentane', 'benzene']
        for ring_type in ['cyclopentane']:
            logging.info(f"Analysing ring {ring_type}...")
            cvs_filename = ring_type + '_params_' + params + '_analysis_output.csv'
            with open(str(Path(OUTPUT_DIR / cvs_filename).resolve()), mode='w', newline='') as f:
                rows = []
                w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                filepathes = get_filepathes(Path(args.rootdir), ring_type)
                # print(filepathes)
                with Pool(int(CPU_COUNT)) as p:
                    logging.info("Creating modified paths...")
                    modified_filepathes = [(f, arguments) for f in filepathes]
                    logging.info(f"Found {len(modified_filepathes)} modified filepathes")
                    # logging.info("Writing rows to the CSV file in real-time...")
                    # for filepath, arguments in modified_filepathes:
                    #     row = run_exe(filepath, arguments)
                    #     w.writerow(row)
                    logging.info(f"getting rows...")
                    rows = p.starmap(run_exe, modified_filepathes)
                    logging.info(f"writing rows...")
                    w.writerows(rows)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)


def main():
    parser = argparse.ArgumentParser(description='ED coverage analysis')
    parser.add_argument('rootdir', type=str, help='Root directory (e.g. "validation_data")')
    
    parser.add_argument('-s', action='store_true', help='Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('-m', '--more_or_equal', action='store_true', help='Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface')
    parser.add_argument('-c', '--closest_voxel', action='store_true', help='Instead of trilinear interpolation, the intensity of the closest voxel is used')
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG, 
					    format='%(asctime)s - %(levelname)s - %(message)s',
                        filename='density-cyclopent-FINAL.log') 
    logging.info(f"Cpu count = {CPU_COUNT}")
    _create_output_folder()
    run_analysis(args)

if __name__ == '__main__':
    main()
