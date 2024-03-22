import csv
from multiprocessing import Pool, cpu_count
import os
import subprocess
import urllib.request as r
import argparse
from pathlib import Path
import shutil

CCP4_DIR = Path('./ccp4')
EXE = Path('./electron_density_coverage_analysis.py')
params = '-s'
OUTPUT_DIR = Path('./output')
NO_CCP4_AVAILABLE_FILE = Path(OUTPUT_DIR / 'no_ccp4_pdb_ids.txt')
CPU_COUNT = cpu_count() / 2

def _create_ccp4_folder():
    # create ccp4 folder
    if CCP4_DIR.exists():
        shutil.rmtree(str(CCP4_DIR.resolve()))

    CCP4_DIR.mkdir(parents=True, exist_ok=False)
    return CCP4_DIR

def _create_output_folder():
    if OUTPUT_DIR.exists():
        shutil.rmtree(str(OUTPUT_DIR.resolve()))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    return OUTPUT_DIR


def _download_single_file(data: tuple):
    link, path, pdb_id = data
    try:
        r.urlretrieve(link, path)
    except:
        return pdb_id

def download_ccp4_multithread(args: argparse.Namespace):
    _create_ccp4_folder()
    _create_output_folder()
    
    input_dir = Path(args.rootdir)
    pdbe_base_link = 'https://www.ebi.ac.uk/pdbe/coordinates/files/'
    no_ccp4_available_file = str(NO_CCP4_AVAILABLE_FILE.resolve())

    l = []

    for dirname, dirs, files in os.walk(str(input_dir.resolve())):
        for f in files:
            if f.endswith(".pdb"):
                pdb_id = f.split("_")[1]
                link = pdbe_base_link + pdb_id + ".ccp4"
                path = str(Path(CCP4_DIR / f'{pdb_id}.ccp4'))

                l.append((link, path, pdb_id))

    # removing duplicates
    l = set(l)
    l = list(l)

    pool = Pool(int(CPU_COUNT))
    results = pool.map(_download_single_file, l)

    results_without_none = [x for x in results if x is not None]

    pool.close()
    pool.join()

    with open(no_ccp4_available_file, 'w') as w:
        for i in results_without_none:
            w.write(i + '\n')
            # print(i)

def process_args(args: argparse.Namespace):
    s = "-s"
    if args.closest_voxel:
        s = s + "c"
    if args.more_or_equal:
        s = s + "m"

    global params
    params = s

# ['tests/example_input/validation_data/cyclopentane\\filtered_ligands\\06L\\patterns\\06L_8gxp_0.pdb']
def run_exe(ligand_filepath: Path):
    PQ_pdb_name = ligand_filepath.name.split(".")[0]
    pdb_id = PQ_pdb_name.split('_')[1]
    residue_id = ligand_filepath.parent.parent.name
    ccp4_filepath = str(Path(CCP4_DIR / f'{pdb_id}.ccp4').resolve())
    
    args = [
        'python',
        str(EXE.resolve()),
        # PARAMS should be global
        params,
        str(ligand_filepath.resolve()),
        ccp4_filepath]
    output = subprocess.run(args,
        capture_output=True,
        text=True)

    result = (PQ_pdb_name, residue_id, output.stdout)

    return result

def get_filepathes(rootdir: Path, ligand_type: str):
    l = []
    no_ccp4_available_file = str(NO_CCP4_AVAILABLE_FILE.resolve())
    with open(no_ccp4_available_file) as f:
        no_ccp4 = f.read().splitlines()

    for f in Path(rootdir / ligand_type).rglob("*"):
        if f.is_file():
            # get path
            stem = f.stem
            pdb_id = stem.split('_')[1]
            if pdb_id not in no_ccp4:
                l.append(f)
    return l

def run_analysis(args: argparse.Namespace):
    process_args(args)
    # args.rootdir - contains benzene, cyclohexane, cyclopentane
    list_of_dirs = os.listdir(args.rootdir)
    for ligand_type in list_of_dirs:
        cvs_filename = ligand_type + '_params' + params + '_analysis_output.csv'
        with open(str(Path(OUTPUT_DIR / cvs_filename).resolve()), mode='w', newline='') as f:
            rows = []
            w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            filepathes = get_filepathes(Path(args.rootdir), ligand_type)
            with Pool(int(CPU_COUNT)) as p:
                # runs programm
                rows = p.map(run_exe, filepathes)
                w.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description='ED coverage analysis')
    parser.add_argument('rootdir', type=str, help='Root directory (e.g. "validation_data")')
    
    parser.add_argument('-s', action='store_true', help='Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('-m', '--more_or_equal', action='store_true', help='Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface')
    parser.add_argument('-c', '--closest_voxel', action='store_true', help='Instead of trilinear interpolation, the intensity of the closest voxel is used')
    
    args = parser.parse_args()
    download_ccp4_multithread(args)
    
    run_analysis(args)

if __name__ == '__main__':
    main()
