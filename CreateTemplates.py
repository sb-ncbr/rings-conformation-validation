import math
import os
import sys
import shutil
import tempfile
import logging
from csv import reader, writer
from pathlib import Path

from ConfAnalyser import ConfAnalyser
from subprocess import Popen, PIPE
from argparse import ArgumentParser
from gemmi import cif
from HelperModule.Ring import Ring
from HelperModule.getter_functions import get_data_from_pdb, get_bonds_from_cif
from HelperModule.constants import *


class SiteBinderException(Exception):
    """Exception raised when SiteBinder has failed."""

    def __init__(self, message="Error in SiteBinder."):
        self.message = message
        super().__init__(self.message)


class NoInputFilesException(Exception):
    """Exception raised when no input files are found for creating templates."""

    def __init__(self, message="No input files for creating templates were found."):
        self.message = message
        super().__init__(self.message)


class NoConformationsException(Exception):
    """Exception raised when no conformations are found for creating templates."""

    def __init__(self, message="No conformations for creating templates were found."):
        self.message = message
        super().__init__(self.message)


def extract_rmsd(tmp: Path):
    # get only names and rmsd values from the output of SB, change the delimiters

    with open(tmp / 'rmsd.csv', 'r') as f, open(tmp / 'extracted_rmsd.csv', 'w') as out:
        for line in f.readlines():
            name, rmsd, *_ = line.split(',')
            new_line = ';'.join([name, rmsd])
            out.write("".join([new_line, '\n']))


def sort_rmsd(tmp: Path):
    with open(tmp / 'extracted_rmsd.csv', 'r', newline='\n') as f, \
            open(tmp / 'rmsd_sorted.csv', 'w', newline='\n') as out:
        csv1 = reader(f, delimiter=';', doublequote=False)
        next(csv1, None)
        csv_writer = writer(out, delimiter=';', escapechar='\\', doublequote=False)
        csv_writer.writerows(sorted(csv1, key=lambda x: x[1], reverse=True))


def remove_outliers(lines_count, tmp: Path):
    n = math.ceil(0.1 * lines_count)
    logging.info(f"Searching outliers among {lines_count} pdb files.")
    file_path = tmp / 'rmsd_sorted.csv'

    with open(file_path, 'r+') as file:
        for _ in range(n):
            file.readline()
        remaining_content = file.read()
        file.seek(0)
        file.write(remaining_content)
        file.truncate()


def get_paths(ring: Ring, output_file: Path, path_to_tmp_dir: Path, main_workflow_output_dir: Path):
    """Reads 'rmsd_sorted.csv' input file and creates (updates) output file,
    containing paths to pdb files from the read file.
    Example of line in the input file:
    CVM_2xqu_52;0.004
    Example of resulting line in the output file:
    /this/is/path/to/the/CVM_2xqu_52.pdb

    Args:
        ring (Ring): ring type
        output_file (Path): name of the output file
        path_to_tmp_dir (Path): path to temp directory
        main_workflow_output_dir (Path): path to output directory
    """

    input_file = path_to_tmp_dir / 'rmsd_sorted.csv'
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out:

        for line in in_file.readlines():
            name, *_ = line.split(';')
            ligand_name = extract_letters_before_underscore(name)
            if ligand_name is None:
                logging.error(f"Wrong filename format of {name}. Skipping...")
                continue
            p = main_workflow_output_dir / ring.name.lower() / FILTERED_DATA / ligand_name / "patterns" / (name + '.pdb')
            out.write(str(p) + '\n')


def run_site_binder(sb_input: Path, site_binder: Path, tmp: Path):
    commands = [site_binder,
                sb_input,
                tmp / 'rmsd.csv',
                tmp / 'pairing.csv',
                tmp / 'info.txt',
                tmp / 'avg.pdb',
                tmp / 'superimposed_folder']

    args_sb = []
    if os.name == 'posix':
        args_sb.extend(['mono'])
    args_sb.extend(commands)

    sb_process = Popen(args_sb, stdout=PIPE)
    out, err = sb_process.communicate()
    if out:
        logging.error(out)
        raise SiteBinderException
    if err:
        logging.error(err)
        raise SiteBinderException


def validate_atoms_order(ligand_block: cif.Block, input_atoms: tuple[str, ...], atoms_count: int) -> \
        tuple[str, ...] | None:
    atom_bonds_all = get_bonds_from_cif(ligand_block)
    atom_bonds_filtered = set()

    for bond in atom_bonds_all:
        a_1, a_2, *_ = bond
        atom_1, atom_2 = a_1.strip('"'), a_2.strip('"')
        if (atom_1 not in input_atoms) or (atom_2 not in input_atoms):
            continue
        atom_bonds_filtered.add((atom_1, atom_2))

    if len(atom_bonds_filtered) != atoms_count:
        logging.error(
            f"Received: {input_atoms}\nFound only: {atom_bonds_filtered}\n{ligand_block.find_value('_chem_comp.id')}")
        return None

    for i in range(len(input_atoms) - 1):
        if ((input_atoms[i], input_atoms[i + 1]) not in atom_bonds_filtered) \
                and ((input_atoms[i + 1], input_atoms[i]) not in atom_bonds_filtered):
            correctly_ordered_atoms = correct_atoms_order(atom_bonds_filtered)
            return correctly_ordered_atoms

    return input_atoms


def correct_atoms_order(bonds: set[tuple[str, str]]) -> tuple[str, ...]:
    result = []
    bonds_copy = bonds.copy()
    first, last = bonds_copy.pop()
    result.append(first)

    curr = first
    while curr != last:
        for bond in bonds_copy:
            if bond[0] == curr:
                result.append(bond[1])
                curr = bond[1]
            elif bond[1] == curr:
                result.append(bond[0])
                curr = bond[0]
            else:
                continue

            bonds_copy.remove(bond)
            break

    return tuple(result)


def create_conf_analyser_input(ring: Ring, paths_to_pdbs_filename: Path,
                               atom_names_filename: Path, ligands_dir: Path, document: cif.Document) -> None:
    unique_lines: set[tuple[str, tuple[str, ...]]] = set()
    with open(paths_to_pdbs_filename, 'w', encoding='utf-8') as paths:

        for root, _, files in os.walk(ligands_dir):
            for file in files:

                if not file.endswith('.pdb'):
                    continue

                path_to_pdb = Path(root) / file
                paths.write(str(path_to_pdb) + '\n')

                ligand, atom_names = get_data_from_pdb(path_to_pdb, ring)
                unique_lines.add((ligand, tuple(atom_names)))

        with open(atom_names_filename, 'w', encoding='utf-8') as atoms:
            for ligand_name, atom_names in unique_lines:
                ligand_block = document.find_block(ligand_name)
                if ligand_block is None:
                    logging.warning(
                        f"The block {ligand_name} was not found in {DEFAULT_DICT_NAME}. Skipping...\n")
                    continue

                ordered_atoms = validate_atoms_order(ligand_block, atom_names, ring.atom_number)
                if ordered_atoms is None:
                    continue

                atoms.write("".join([ligand_name, ' ']))
                for atom_name in ordered_atoms:
                    atoms.write("".join([atom_name.strip(), ' ']))
                atoms.write('\n')

    if len(unique_lines) == 0:
        raise NoInputFilesException()


def map_ring_to_molecule_type(ring: Ring) -> ConfAnalyser.MoleculeType:
    match ring:
        case Ring.CYCLOHEXANE:
            return ConfAnalyser.MoleculeType.Cyclohexane
        case Ring.CYCLOPENTANE:
            return ConfAnalyser.MoleculeType.Cyclopentane
        case Ring.BENZENE:
            return ConfAnalyser.MoleculeType.Benzene
        case _:
            return ConfAnalyser.MoleculeType.Undefined


def extract_letters_before_underscore(string: str):
    parts = string.split('_')

    if len(parts) == 3 and len(parts[1]) == 4:
        return parts[0]

    return None


def run_conf_analyser(ring: Ring, main_workflow_output_dir: Path, document: cif.Document) -> None:
    # Prepare input files
    logging.info("Creating input files for ConfAnalyser...")
    with tempfile.TemporaryDirectory() as tmp_str:
        # print(f'TMP: {type(tmp_str)} {tmp_str}')
        tmp = Path(tmp_str)
        paths_to_pdbs_path = tmp / 'paths_to_pdbs.txt'
        atom_names_path = tmp / 'atom_names.txt'

        ligands_dir = main_workflow_output_dir / ring.name.lower() / FILTERED_DATA
        try:
            create_conf_analyser_input(ring, paths_to_pdbs_path, atom_names_path, ligands_dir, document)
        except NoInputFilesException as e:
            logging.error(str(e))
            sys.exit('Exiting...')

        result_dict = {}
        mol_type = map_ring_to_molecule_type(ring)

        try:
            logging.info(f"Starting ConfAnalyser...")
            result_dict = ConfAnalyser.ConfAnalyser(paths_file=str(paths_to_pdbs_path), names_file=str(atom_names_path),
                                                    molecule_type=mol_type).result()

        except Exception as e:
            logging.error(str(e))
            sys.exit('Exiting...')

        conformations_found = False
        for conf in result_dict.keys():

            if conf in ('unanalysed', 'undefined') or not result_dict[conf]:
                continue
            conformations_found = True

        if not conformations_found:
            raise NoConformationsException

        # Process the results
        for conf in result_dict.keys():

            if conf in ('unanalysed', 'undefined') or not result_dict[conf]:
                continue

            logging.info(f'Processing conformation {conf}...')
            path_to_conf_template_dir: Path = main_workflow_output_dir / ring.name.lower() / TEMPLATES_DIR / conf
            os.makedirs(path_to_conf_template_dir, exist_ok=True)

            paths_to_pdbs_for_curr_conf: Path = tmp / (conf + '.txt')
            with open(paths_to_pdbs_for_curr_conf, 'w') as file:
                for pdb in result_dict[conf]:
                    ligand_name = extract_letters_before_underscore(str(Path(pdb).name))
                    if ligand_name is None:
                        logging.error(f"Wrong filename format of {pdb}. Skipping...")
                        continue
                    p = main_workflow_output_dir / ring.name.lower() / FILTERED_DATA / ligand_name / "patterns" / pdb
                    file.write(str(p) + '\n')

            lines_count = len(result_dict[conf])

            while lines_count > 1:

                try:
                    run_site_binder(paths_to_pdbs_for_curr_conf, SB_AVG, tmp)
                except Exception as e:
                    logging.error(str(e))
                    sys.exit("Exiting...")
                extract_rmsd(tmp)
                sort_rmsd(tmp)
                remove_outliers(lines_count, tmp)

                # file <paths_to_pdbs_for_curr_conf> is being updated by get_paths()
                get_paths(ring, paths_to_pdbs_for_curr_conf, tmp, main_workflow_output_dir)
                with open(paths_to_pdbs_for_curr_conf, 'r') as f:
                    lines_count = len(f.readlines())

            # only one pdb file (the found template) is left in <paths_to_pdbs_for_curr_conf>
            with open(paths_to_pdbs_for_curr_conf, 'r') as f:
                path_to_template_pdb = f.readline().strip()
                new_path: Path = path_to_conf_template_dir / (conf.lower() + '.pdb')
                shutil.copy(Path(path_to_template_pdb), new_path)

                logging.info(f"Found template {path_to_template_pdb} --> {conf.lower()}.pdb")


def read_component_dictionary(path_to_comp_dict: Path) -> cif.Document:
    logging.info('Reading components dictionary...')
    try:
        document = cif.read(str(path_to_comp_dict))
        return document
    except FileNotFoundError:
        logging.error(f"File {str(path_to_comp_dict)} not found. Please check the file path.")
        sys.exit('Exiting...')
    except PermissionError:
        logging.error(f"Permission denied. Make sure you have the necessary permissions "
                      f"to access the file {str(path_to_comp_dict)}.")
        sys.exit('Exiting...')
    except Exception as e:
        logging.error(f"An error occurred while trying to read {str(path_to_comp_dict)}: {e}")
        sys.exit('Exiting...')


def is_valid_directory(directory: str) -> bool:
    directory_path = Path(directory).resolve()
    if not directory_path.exists():
        logging.error(f"The directory {str(directory_path)} does not exist.")
        return False
    return True


def validate_args(ring: str, output_dir: str, input_dir: str) -> bool:

    # Validate the ring argument
    if ring.upper() not in Ring.__members__.keys():
        logging.error(f"Ring {ring} is not a valid Ring. Currently supported: {[e.name for e in Ring]}")
        return False

    # Validate the input directory
    if not is_valid_directory(input_dir):
        return False

    # Validate the output directory
    if not is_valid_directory(output_dir):
        return False

    # Validate the main workflow output directory
    main_workflow_output_dir = Path(output_dir) / MAIN_DIR
    if not is_valid_directory(main_workflow_output_dir):
        return False

    return True


def prepare_templates_dir(templates_dir: Path):
    if templates_dir.exists():
        shutil.rmtree(str(templates_dir.resolve()))
        logging.info(f"Template directory {str(templates_dir)} was not empty, cleaning it up...")

    os.makedirs(templates_dir)


def main(ring: str, output_dir: str, input_dir: str):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s'
                        )
    # Validate input arguments
    if not validate_args(ring, output_dir, input_dir):
        sys.exit(1)

    logging.info(f"[{ring.capitalize()}]: Starting CreateTemplates...")
    main_workflow_output_dir = (Path(output_dir) / MAIN_DIR).resolve()
    path_to_comp_dict: Path = (Path(input_dir) / DEFAULT_DICT_NAME).resolve()
    prepare_templates_dir(main_workflow_output_dir / ring.lower() / TEMPLATES_DIR)
    document = read_component_dictionary(path_to_comp_dict)

    try:
        run_conf_analyser(Ring[ring.upper()], main_workflow_output_dir, document)

    except NoConformationsException as e:
        logging.error(str(e))
        sys.exit('Exiting...')
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
        sys.exit('Exiting...')

    logging.info(f"[{ring.capitalize()}]: CreateTemplates has completed successfully")


if __name__ == "__main__":
    parser = ArgumentParser(description=f"Creates template molecules for all the target ring's "
                                        f"conformations")
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--ring', required=True, type=str,
                          help=f'Choose the target ring type. Currently supported:{os.linesep}'
                               f'{[e.name for e in Ring]}')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='Path to the output directory. Should be the same as in the previous step.')
    required.add_argument('-i', '--input', type=str, required=True,
                          help='Path to the directory with input data (local pdb, ccp4 files, etc.)')
    args = parser.parse_args()

    main(args.ring, args.output, args.input)
