import math
import os
import sys
import shutil
import tempfile
import logging
from csv import reader, writer
from ConfAnalyser import ConfAnalyser
from subprocess import Popen, PIPE
from argparse import ArgumentParser
from gemmi import cif
from HelperModule.Ring import Ring
from HelperModule.getter_functions import get_data_from_pdb, get_bonds_from_cif
from HelperModule.constants import *


def extract_rmsd(tmp: str):
    # get only names and rmsd values from the output of SB, change the delimiters

    with open(os.path.join(tmp, 'rmsd.csv'), 'r') as f, open(os.path.join(tmp, 'extracted_rmsd.csv'), 'w') as out:
        for line in f.readlines():
            name, rmsd, *_ = line.split(',')
            new_line = ';'.join([name, rmsd])
            out.write("".join([new_line, os.linesep]))


def sort_rmsd(tmp: str):
    with open(os.path.join(tmp, 'extracted_rmsd.csv'), 'r', newline=os.linesep) as f, \
            open(os.path.join(tmp, 'rmsd_sorted.csv'), 'w', newline=os.linesep) as out:
        csv1 = reader(f, delimiter=';', doublequote=False)
        next(csv1, None)
        csv_writer = writer(out, delimiter=';', escapechar='\\', doublequote=False)
        csv_writer.writerows(sorted(csv1, key=lambda x: x[1], reverse=True))


def remove_outliers(lines_count, tmp: str):
    n = math.ceil(0.1 * lines_count)
    logging.info(f"Searching outliers among {lines_count} pdb files.")
    file_path = os.path.join(tmp, 'rmsd_sorted.csv')

    with open(file_path, 'r+') as file:
        for _ in range(n):
            file.readline()
        remaining_content = file.read()
        file.seek(0)
        file.write(remaining_content)
        file.truncate()


def get_paths(ring: Ring, output_file: str, path_to_tmp_dir: str):
    """Reads 'rmsd_sorted.csv' input file and creates (updates) output file, containing paths to pdb files from the read file.
    Example of line in the input file:
    CVM_2xqu_52;0.004
    Example of resulting line in the output file:
    /this/is/path/to/the/CVM_2xqu_52.pdb

    Args:
        ring (Ring): ring type
        output_file (str): name of the output file
        path_to_tmp_dir (str): path to temp directory
    """
    input_file = os.path.join(path_to_tmp_dir, 'rmsd_sorted.csv')
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out:

        for line in in_file.readlines():
            name, *_ = line.split(';')
            ligand_name = extract_letters_before_underscore(name)
            if ligand_name is None:
                logging.error(f"Wrong filename format of {name}. Skipping...")
                continue
            out.write("".join([os.path.join(MAIN_DIR, ring.name.lower(), FILTERED_DATA,
                                            ligand_name, 'patterns', name + '.pdb'), os.linesep]))


def run_site_binder(sb_input: str, site_binder: str, tmp: str):
    args_sb = [site_binder,
               sb_input,
               os.path.join(tmp, 'rmsd.csv'),
               os.path.join(tmp, 'pairing.csv'),
               os.path.join(tmp, 'info.txt'),
               os.path.join(tmp, 'avg.pdb'),
               os.path.join(tmp, 'superimposed_folder')]

    sb_process = Popen(args_sb, stdout=PIPE)
    out, err = sb_process.communicate()
    if out:
        logging.error(out)
        raise Exception
    if err:
        logging.error(err)
        raise Exception


def validate_atoms_order(ligand_block: cif.Block, input_atoms: tuple[str], atoms_count: int) -> \
        tuple[str] | None:
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
            f"Received: {input_atoms}{os.linesep}Found only: {atom_bonds_filtered}{os.linesep}{ligand_block.find_value('_chem_comp.id')}")
        return None

    for i in range(len(input_atoms) - 1):
        if ((input_atoms[i], input_atoms[i + 1]) not in atom_bonds_filtered) \
                and ((input_atoms[i + 1], input_atoms[i]) not in atom_bonds_filtered):
            correctly_ordered_atoms = correct_atoms_order(atom_bonds_filtered)
            return correctly_ordered_atoms

    return input_atoms


def correct_atoms_order(bonds: set[tuple[str, str]]) -> tuple[str]:
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


def create_conf_analyser_input(ring: Ring, paths_to_pdbs_filename: str, atom_names_filename: str) -> None | int:
    unique_lines: set[tuple[str, tuple[str]]] = set()
    with open(paths_to_pdbs_filename, 'w') as paths:
        ligands_dir = os.path.join(MAIN_DIR, ring.name.lower(), FILTERED_DATA)
        for root, _, files in os.walk(ligands_dir):
            for file in files:

                if not file.endswith('.pdb'):
                    continue

                path_to_pdb = os.path.join(root, file)
                paths.write("".join([path_to_pdb, os.linesep]))

                ligand, atom_names = get_data_from_pdb(path_to_pdb, ring)
                unique_lines.add((ligand, tuple(atom_names)))

        document = cif.read(DEFAULT_DICT_NAME)

        with open(atom_names_filename, 'w') as atoms:
            for ligand_name, atom_names in unique_lines:
                ligand_block = document.find_block(ligand_name)
                if ligand_block is None:
                    logging.error(
                        f"The block {ligand_name} was not found in {DEFAULT_DICT_NAME}. Skipping...{os.linesep}")

                ordered_atoms = validate_atoms_order(ligand_block, atom_names, ring.atom_number)
                if ordered_atoms is None:
                    continue

                atoms.write("".join([ligand_name, ' ']))
                for atom_name in ordered_atoms:
                    atoms.write("".join([atom_name, ' ']))
                atoms.write(os.linesep)

    if len(unique_lines) == 0:
        logging.error("No input files for creating templates were found.")
        return -1


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


def run_conf_analyser(ring: Ring) -> None:
    # Prepare input files
    with tempfile.TemporaryDirectory() as tmp:
        paths_to_pdbs_path = os.path.join(tmp, 'paths_to_pdbs.txt')
        atom_names_path = os.path.join(tmp, 'atom_names.txt')
        return_code = create_conf_analyser_input(ring, paths_to_pdbs_path, atom_names_path)

        # no input data
        if return_code == -1:
            return

        result_dict = {}
        mol_type = map_ring_to_molecule_type(ring)

        try:
            result_dict = ConfAnalyser.ConfAnalyser(paths_file=paths_to_pdbs_path, names_file=atom_names_path,
                                                    molecule_type=mol_type).result()

        except Exception as e:
            logging.exception(f"An exception in ConfAnalyser: {e}")
            raise Exception

        # Process the results
        for conf in result_dict.keys():

            if conf in ('unanalysed', 'undefined') or not result_dict[conf]:
                continue

            logging.info(f'Processing conformation {conf}...')
            path_to_conf_template_dir = os.path.join(MAIN_DIR, ring.name.lower(), TEMPLATES_DIR, conf)
            os.makedirs(path_to_conf_template_dir, exist_ok=True)

            paths_to_pdbs_for_curr_conf = os.path.join(tmp, f"{conf}.txt")
            with open(paths_to_pdbs_for_curr_conf, 'w') as file:
                for pdb in result_dict[conf]:
                    ligand_name = extract_letters_before_underscore(pdb)
                    if ligand_name is None:
                        logging.error(f"Wrong filename format of {pdb}. Skipping...")
                        continue

                    file.write("".join([os.path.join(MAIN_DIR, ring.name.lower(), FILTERED_DATA,
                                                     ligand_name, 'patterns', pdb), os.linesep]))

            lines_count = len(result_dict[conf])

            while lines_count > 1:

                try:
                    run_site_binder(paths_to_pdbs_for_curr_conf, SB_AVG, tmp)
                except Exception:
                    sys.exit("Error in SiteBinder occurred. Exiting...")
                extract_rmsd(tmp)
                sort_rmsd(tmp)
                remove_outliers(lines_count, tmp)

                # file <paths_to_pdbs_for_curr_conf> is beeing updated by get_paths()
                get_paths(ring, paths_to_pdbs_for_curr_conf, tmp)
                with open(paths_to_pdbs_for_curr_conf, 'r') as f:
                    lines_count = len(f.readlines())

            # only one pdb file (the found template) is left in <paths_to_pdbs_for_curr_conf>
            with open(paths_to_pdbs_for_curr_conf, 'r') as f:
                path_to_template_pdb = f.readline().strip()
                new_path = os.path.join(path_to_conf_template_dir, f'{conf.lower()}.pdb')
                shutil.copy(path_to_template_pdb, new_path)

                logging.info(f"Found template {path_to_template_pdb} --> {conf.lower()}.pdb")


def main(ring: str):
    if ring.upper() not in Ring.__members__.keys():
        sys.exit(f"Supported rings are: {[e.name for e in Ring]}")

    ring: Ring = Ring[ring.upper()]

    os.makedirs(os.path.join(MAIN_DIR, ring.name.lower(), TEMPLATES_DIR), exist_ok=True)

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    run_conf_analyser(ring)


if __name__ == "__main__":
    parser = ArgumentParser(description=f"Creates template molecules for all the target ring's "
                                        f"conformations")
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--ring', required=True, type=str,
                          help=f'Choose the target ring type. Currently supported:{os.linesep}'
                               f'{[e.name for e in Ring]}')
    args = parser.parse_args()

    main(args.ring)
