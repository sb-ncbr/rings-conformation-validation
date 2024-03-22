import argparse
import gemmi
import statistics as st
import math

# parsing args
parser = argparse.ArgumentParser(description="Determines electron density coverage of a cycle. By default, uses trilinear interpolation to infer the electron density value, and considers an atom to be covered by the electron density when the corresponding intensity is MORE than the threshold for the isosurface (1.5 sigma)")
parser.add_argument('-s', action='store_true', help='Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle')
parser.add_argument('-d', action='store_true', help='Detailed mode - output a sequence of paried values, where the first member of each pair is a serial atom number, and the second member is "y" - if atom is covered, and "n" - if it is not.')
parser.add_argument('-m', '--more_or_equal', action='store_true', help='Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface')
parser.add_argument('-c', '--closest_voxel', action='store_true', help='Instead of trilinear interpolation, the intensity of the closest voxel is used')

parser.add_argument('input_cycle_pdb', type=str, help='Input PDB file with coordinates of a cycle, produced by PatternQuery')
parser.add_argument('input_density_ccp4', type=str, help='Input electron density file for the corresponding protein structure from the PDB in CCP4 format')

args = parser.parse_args()
str = gemmi.read_pdb(args.input_cycle_pdb)
map = gemmi.read_ccp4_map(args.input_density_ccp4)
map.setup()

# calculate the sigma values
grid_values = []
for point in map.grid:
    if not math.isnan(point.value):
        grid_values.append(point.value)

std = st.pstdev(grid_values)
sigma_lvl = 1.5 * std

# compare intensity, corresponding to the given position, to the threshold for isosurface (MORE vs MORE OR EQUAL)
def determine_atom_coverage(pos):
    if args.more_or_equal:
        return get_intensity(pos) >= sigma_lvl
    return get_intensity(pos) > sigma_lvl

# get intensity corresponding to the given position (trilinear interpolation vs itensity of the closest voxel)
def get_intensity(pos):
    if args.closest_voxel:
        return map.grid.get_nearest_point(pos).value
    return map.grid.interpolate_value(pos)


if args.s:
    total_atom_count = 0
    covered_atoms_count = 0
    for model in str:
        for chain in model:
            for res in chain:
                for atom in res:
                    total_atom_count = total_atom_count + 1
                    if determine_atom_coverage(atom.pos):
                        covered_atoms_count = covered_atoms_count + 1

    print(covered_atoms_count, total_atom_count, sep=";", end='')

if args.d:
    for model in str:
        for chain in model:
            for res in chain:
                for atom in res:
                    if determine_atom_coverage(atom.pos):
                        print(atom.serial, "y", sep=";", end=';')
                    else:
                        print(atom.serial, "n", sep=";", end=';')