import argparse
import logging
import gemmi
import statistics as st
import math

def run_as_function(args: argparse.Namespace):
    return run_analysis(args)

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

def run_analysis(args: argparse.Namespace):
    try:
        output = None
        str = gemmi.read_pdb(args.input_cycle_pdb)
        map = gemmi.read_ccp4_map(args.input_density_ccp4)
        # map.setup()
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
            # print(covered_atoms_count, total_atom_count, sep=";", end='')

        if args.d:
            output = []
            for model in str:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            if determine_atom_coverage(atom.pos, map, sigma_lvl, args):
                                output.append(f'{atom.serial};y;')
                                # print(atom.serial, "y", sep=";", end=';')
                            else:
                                output.append(f'{atom.serial};n;')
                                # print(atom.serial, "n", sep=";", end=';')
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return output

def main():
    parser = argparse.ArgumentParser(description="Determines electron density coverage of a cycle. By default, uses trilinear interpolation to infer the electron density value, and considers an atom to be covered by the electron density when the corresponding intensity is MORE than the threshold for the isosurface (1.5 sigma)")
    parser.add_argument('-s', action='store_true', help='Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('-d', action='store_true', help='Detailed mode - output a sequence of paried values, where the first member of each pair is a serial atom number, and the second member is "y" - if atom is covered, and "n" - if it is not.')
    parser.add_argument('-m', '--more_or_equal', action='store_true', help='Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface')
    parser.add_argument('-c', '--closest_voxel', action='store_true', help='Instead of trilinear interpolation, the intensity of the closest voxel is used')

    parser.add_argument('input_cycle_pdb', type=str, help='Input PDB file with coordinates of a cycle, produced by PatternQuery')
    parser.add_argument('input_density_ccp4', type=str, help='Input electron density file for the corresponding protein structure from the PDB in CCP4 format')

    args = parser.parse_args()
    # print(args)
    output = run_analysis(args)
    return output

if __name__ == '__main__':
    output = main()
    if isinstance(output, str):
        print(output)
    elif isinstance(output, list):
        for i in output:
            print(i)