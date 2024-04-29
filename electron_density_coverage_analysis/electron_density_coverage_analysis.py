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
