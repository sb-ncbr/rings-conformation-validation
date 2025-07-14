from Bio import PDB
import numpy as np
from glob import glob
import math
import json
from pathlib import Path
from numba import jit
from scipy.optimize import shgo
import os
import logging
from argparse import ArgumentParser


def cross(a, b):
    return np.array([a[1]*b[2] - a[2]*b[1],
                     a[2]*b[0] - a[0]*b[2],
                     a[0]*b[1] - a[1]*b[0]], dtype='float64')


def norm(a):
    return math.sqrt(np.sum(a*a))


@jit(nopython=True, cache=True, fastmath=True)
def distance_from_plane(point, plane_point, normal_vector):
    return np.dot(point - plane_point, normal_vector) / np.linalg.norm(normal_vector)


@jit(nopython=True, cache=True, fastmath=True)
def distances_from_plane(coordinates, plane_point, normal_vector):
    return np.sum(np.array([(np.dot(coordinate - plane_point, normal_vector) / np.linalg.norm(normal_vector)) ** 2 for coordinate in coordinates]))


def obj_fun(normal_vector, plane_point, coordinates):
    try:
        return distances_from_plane(coordinates, plane_point, normal_vector)
    except ZeroDivisionError:
        return 1000


def calculate_hr_angles_from_pdb(file):
    atoms = [atom for atom in PDB.PDBParser(QUIET=True).get_structure("structure", file)[0].get_atoms() if atom.element != "H"]

    if len(atoms) < 5:
        raise ValueError("Too few atoms to form ring")

   
    kdtree = PDB.NeighborSearch(atoms)
    sorted_atoms = [atoms[0],
                    kdtree.search(center=atoms[0].coord, radius=1.8, level="A")[1]]
    for x in range(1, len(atoms) - 1):
        nearest_atoms = kdtree.search(center=sorted_atoms[x].coord, radius=1.8, level="A")
        non_sorted = [atom for atom in nearest_atoms if atom not in sorted_atoms]
        if len(non_sorted) != 1:
            raise ValueError("Ambiguous ring connectivity")
        sorted_atoms.append(non_sorted[0])

    coordinates = np.array([atom.coord for atom in sorted_atoms], dtype=np.float64)
    center = np.mean(coordinates, axis=0)
    
    opt = shgo(obj_fun, bounds=[[-1, 1]] * 3, args=(center, coordinates))
    distances = [distance_from_plane(c, center, opt.x) for c in coordinates]

    
    max_central_value = 0
    central_index = None
    for i, d in enumerate(distances):
        if d > 0 and d > distances[(i + 1) % len(distances)] and d > distances[(i - 1) % len(distances)] or \
           d < 0 and d < distances[(i + 1) % len(distances)] and d < distances[(i - 1) % len(distances)]:
            central_value = abs(distances[(i - 1) % len(distances)] - d) + abs(d - distances[(i + 1) % len(distances)])
            if central_value > max_central_value:
                max_central_value = central_value
                central_index = i

    if central_index is None:
        raise ValueError("Could not determine central atom")

    for _ in range(central_index):
        distances.append(distances.pop(0))
        sorted_atoms.append(sorted_atoms.pop(0))

    if distances[0] < 0:
        distances = [-x for x in distances]

    if distances[1] > distances[-1]:
        distances = [distances[0]] + distances[1:][::-1]
        sorted_atoms = [sorted_atoms[0]] + sorted_atoms[1:][::-1]

    #HR angle calc
    coords = np.array([atom.coord for atom in sorted_atoms[1:] + sorted_atoms[:3]], dtype='float64')
    coords -= np.mean(coords[:len(sorted_atoms)], axis=0)

    r = np.array([coords[i+1] - coords[i] for i in range(len(sorted_atoms))])
    a = np.array([coords[2*(i+1)] - coords[2*i] for i in range(3)])
    p = np.array([cross(r[i-1], r[i]) if i > 0 else [0,0,0] for i in range(len(sorted_atoms))])
    q = np.array([cross(a[i], p[2*i+1]) for i in range(len(sorted_atoms)-3)])
    n = cross(a[1], a[0])

    thetas = []
    for v in q:
        angle = 90 - math.degrees(math.acos(np.dot(v, n) / (norm(v) * norm(n))))
        thetas.append(angle)

    if len(thetas) == 2:
        tr = -1 if thetas[0] > thetas[1] else 1
        return {"theta1": thetas[0]*tr, "theta2": thetas[1]*tr, "theta3": None}
    elif len(thetas) == 3:
        tr = -1 if thetas[2] < 0 else 1
        return {"theta1": thetas[0]*tr, "theta2": thetas[1]*tr, "theta3": thetas[2]*tr}
    else:
        raise ValueError("Invalid number of HR angles")

def process_all_ligands(input_dir, output_json):
    result = {}
    excluded = []
    for file in glob(f"{input_dir}/*/*/*.pdb"):
        try:
            ligand_id = Path(file).stem
            hr = calculate_hr_angles_from_pdb(file)
            result[ligand_id] = hr
        except Exception as e:
            excluded.append({"file": file, "reason": str(e)})
    
    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(result, f, indent=4)

    print(f"HR angle extraction complete. {len(result)} succeeded, {len(excluded)} failed.")
    if excluded:
        print("Excluded files:")
        for e in excluded:
            print(f"{e['file']} - {e['reason']}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python CalculateHR.py <filtered_ligands_path> <output.json>")
        exit(1)

    process_all_ligands(sys.argv[1], sys.argv[2])