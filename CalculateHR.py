from Bio.PDB import NeighborSearch, MMCIFParser
import numpy as np
from glob import glob
import math
import json
from pathlib import Path
from numba import jit
from scipy.optimize import shgo
import sys
from HelperModule.Ring import Ring
from HelperModule.constants import MAIN_DIR

HOMOCYCLES = {Ring.CYCLOHEXANE, Ring.CYCLOPENTANE, Ring.BENZENE}
HETEROCYCLES = {Ring.OXANE, Ring.OXOLANE}


def cross(a, b):
    return np.array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ], dtype=float)

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


#shared HR calc function 
#DOI: https://doi.org/10.1186/s13321-026-01154-0
def calculate_HR(coords, N, apply_tr):
    r = np.array([coords[i+1] - coords[i] for i in range(N)])
    a = np.array([coords[2*(i+1)] - coords[2*i] for i in range(3)])
    p = np.array([cross(r[i-1], r[i]) if i > 0 else [0,0,0] for i in range(N)])
    q = np.array([cross(a[i], p[2*i+1]) for i in range(N-3)])
    n = cross(a[1], a[0])

    thetas = []
    for v in q:
        angle = 90 - math.degrees(math.acos(np.dot(v, n) / (norm(v) * norm(n))))
        thetas.append(angle)
    if apply_tr:
        if len(thetas) == 2:           
            tr = -1 if thetas[0] > thetas[1] else 1
            thetas = [t * tr for t in thetas]
        elif len(thetas) == 3:         
            tr = -1 if thetas[2] < 0 else 1
            thetas = [t * tr for t in thetas]
            
    if len(thetas) == 2:
        return {"theta1": thetas[0], "theta2": thetas[1], "theta3": None}
    elif len(thetas) == 3:
        return {"theta1": thetas[0], "theta2": thetas[1], "theta3": thetas[2]}
    else:
        raise ValueError("Invalid number of HR angles")


#for rings with no heteroatom 
#DOI: https://doi.org/10.1186/s13321-026-01154-0
def calculate_HR_homocycles(atoms):
    kdtree = NeighborSearch(atoms)
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
    return calculate_HR(coords, len(sorted_atoms), apply_tr=True)


#Ordering atoms starting with heteroatom
def order_ring_by_heteroatom(atoms, hetero_targets):
    if isinstance(hetero_targets, str):
        hetero_targets = [hetero_targets]
    

    hets = [
        a for a in atoms 
        if a.element in hetero_targets or a.get_name().strip()[0] in hetero_targets
    ]

    if len(hets) != 1:
        raise ValueError(f"Expected exactly one heteroatom in {hetero_targets}, found {len(hets)}")
    
    hetero_atom = hets[0]

    adjacency = {a: [] for a in atoms}
    for i, a in enumerate(atoms):
        for j, b in enumerate(atoms):
            if i < j:
                if np.linalg.norm(a.coord - b.coord) <= 1.8:
                    adjacency[a].append(b)
                    adjacency[b].append(a)

    if len(adjacency[hetero_atom]) != 2:
        raise ValueError("Heteroatom does not have two ring neighbors")

    c2 = adjacency[hetero_atom][0]
    ring = [hetero_atom, c2]

    prev = hetero_atom
    curr = c2
    while True:
        nxts = [n for n in adjacency[curr] if n is not prev]
        if len(nxts) != 1:
            raise ValueError("Ambiguous connectivity detected")
        nxt = nxts[0]
        if nxt is hetero_atom:
            break
        ring.append(nxt)
        prev, curr = curr, nxt

    return ring


#for cycles containing heteroatom(O) 
def calculate_HR_heterocycles(atoms):
    targets = ["O"]
    sorted_atoms = order_ring_by_heteroatom(atoms, hetero_targets=targets)
    N = len(sorted_atoms)
    if N not in [5, 6]:
        raise ValueError(f"Ring size {N} not supported")


    #shift the sorted atoms to the right position 
    sorted_atoms = sorted_atoms[-2:] + sorted_atoms[:-2]
    #ensures cyclic wrapping
    working_atoms = sorted_atoms + sorted_atoms[:3]

    xs = np.array([atom.coord for atom in working_atoms], dtype='float64')
    xs -= xs.mean(axis=0)

    return calculate_HR(xs, N, apply_tr=False)


#main logic
def calculate_hr_angles_from_cif(file, ring):
    #generalized input handling for ring w and w/o heteroatoms
    atoms = [atom for atom in MMCIFParser(QUIET=True).get_structure("structure", file)[0].get_atoms() if atom.element != "H"]
    if len(atoms) < 5:
        raise ValueError("Too few atoms to form ring")

    if ring in HOMOCYCLES:
        return calculate_HR_homocycles(atoms)
    return calculate_HR_heterocycles(atoms)


def process_all_ligands(ring: Ring, input_dir: str):
    result = {}
    excluded = []

    basedir = Path(input_dir) / MAIN_DIR/ ring.name.lower() / "filtered_ligands"
    for file in glob(f"{basedir}/*/*/*.cif"):
        try:
            ring_id = Path(file).stem
            hr = calculate_hr_angles_from_cif(file, ring)
            result[ring_id] = hr
        except Exception as e:
            excluded.append({"file": file, "reason": str(e)})

    outputpath = Path(input_dir) / MAIN_DIR / ring.name.lower() / "hr_analysis_output"
    Path(outputpath).mkdir(parents=True, exist_ok=True)
    with open(outputpath / "output_HR.json", "w") as f:
        json.dump(result, f, indent=4)

    print(f"HR angle extraction complete. {len(result)} succeeded, {len(excluded)} failed.")
    if excluded:
        print("Excluded files:")
        for e in excluded:
            print(f"{e['file']} - {e['reason']}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python CalculateHR.py <workflow_output_path> containing validation_data folder")
        exit(1)
    input_dir = sys.argv[1]

    for ring in Ring:
        process_all_ligands(ring, input_dir)