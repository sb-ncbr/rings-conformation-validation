from Bio import PDB
import numpy as np
from glob import glob
from scipy.optimize import shgo
from numba import jit
from sys import argv
import math


def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],'g')

def norm(a):
    return math.sqrt(np.sum(a*a))

def tofloat(a):
    b=[]
    for i in a:
        b.append(float(i))
    return b

@jit(nopython=True, cache=True, fastmath=True)
def distance_from_plane(point, plane_point, normal_vektor):
    distance = np.dot(point - plane_point, normal_vektor) / np.linalg.norm(normal_vektor)
    return distance

@jit(nopython=True, cache=True, fastmath=True)
def distances_from_plane(coordinates, plane_point, normal_vector):
    return np.sum(np.array([(np.dot(coordinate - plane_point, normal_vector) / np.linalg.norm(normal_vector)) ** 2 for coordinate in coordinates]))

def obj_fun(normal_vector, plane_point, coordinates):
    try:
        return distances_from_plane(coordinates, plane_point, normal_vector)# + distances_from_plane([(coordinates[x] + coordinates[(x+1)%len(coordinates)])/2 for x in range(len(coordinates))], plane_point, normal_vector)
    except ZeroDivisionError:
        return 1000

class Cycle:
    def __init__(self,
                 file):

        atoms = [atom for atom in PDB.PDBParser(QUIET=True).get_structure("structure", file)[0].get_atoms() if atom.element != "H"]

        # Atoms of cycles can be stored in file without right ordering.
        kdtree = PDB.NeighborSearch(atoms)
        sorted_atoms = [atoms[0],
                        kdtree.search(center=atoms[0].coord,
                                      radius=1.8,
                                      level="A")[1]]
        for x in range(1, len(atoms) - 1):
            nearest_atoms = kdtree.search(center=sorted_atoms[x].coord,
                                          radius=1.8,
                                          level="A")
            non_sorted_nearest_atoms = [atom for atom in nearest_atoms if atom not in sorted_atoms]
            if len(non_sorted_nearest_atoms) != 1:
                exit("FATAL ERROR 1")
            sorted_atoms.append(non_sorted_nearest_atoms[0])

        # find best plane
        coordinates = np.array([atom.coord for atom in sorted_atoms], dtype=np.float64)
        cycle_coordinates_mean = np.array([np.mean([c[0] for c in coordinates]),
                                           np.mean([c[1] for c in coordinates]),
                                           np.mean([c[2] for c in coordinates])],
                                          dtype=np.float64)
        optimisation = shgo(obj_fun,
                            bounds=[[-1, 1] for _ in range(3)],
                            args=(cycle_coordinates_mean, coordinates))
        distances = [distance_from_plane(c, cycle_coordinates_mean, optimisation.x) for c in coordinates]

        max_central_atom_value = 0
        central_atom_index = None
        for x, distance in enumerate(distances):
            if distance > 0 and distance > distances[(x + 1) % len(distances)] and distance > distances[
                (x - 1) % len(distances)] or distance < 0 and distance < distances[
                (x + 1) % len(distances)] and distance < distances[(x - 1) % len(distances)]:
                central_atom_value = abs(distances[(x - 1) % len(distances)] - distance) + abs(
                    distance - distances[(x + 1) % len(distances)])
            else:
                continue
            if central_atom_value > max_central_atom_value:
                max_central_atom_value = central_atom_value
                central_atom_index = x
        for i in range(central_atom_index):
            distances.append(distances.pop(0))
            sorted_atoms.append(sorted_atoms.pop(0))

        if distances[0] < 0:
            distances = [x * -1 for x in distances]

        if distances[1] > distances[-1]:
            distances = [distances[0]] + distances[1:][::-1]
            sorted_atoms = [sorted_atoms[0]] + sorted_atoms[1:][::-1]

        self.distances = distances
        self.file = file
        self.atoms = sorted_atoms

        # calculate Hill-Reilly angles of puckering
        coordinates_flappucker = np.array([atom.coord for atom in self.atoms[1:] + self.atoms[:3]])
        r=np.zeros((len(self.atoms),3),dtype='float64')
        a=np.zeros((3,3),dtype='float64')
        p=np.zeros((len(self.atoms),3),dtype='float64')
        q=np.zeros((len(self.atoms)-3,3),dtype='float64')
        center = np.add.reduce(coordinates_flappucker[0:len(self.atoms)])/len(self.atoms)
        coordinates_flappucker = coordinates_flappucker - center
        for i in range(0,len(self.atoms)):
            r[i]=coordinates_flappucker[i+1] - coordinates_flappucker[i]
        for i in range(0,3):
            a[i] = coordinates_flappucker[2*(i+1)] - coordinates_flappucker[2*i]
        for i in range(1,len(self.atoms)):
            p[i]=cross(r[i-1],r[i])
        for i in range(0,len(self.atoms)-3):
            q[i] = cross(a[i],p[2*i+1])
        n=cross(a[1],a[0])
        theta=[]
        for i in q:
            theta.append(90 - math.degrees(math.acos(np.dot(i,n) / (norm(i)*norm(n)))))
        if len(theta) == 2:
            tr = 1
            if theta[0] > theta[1]:
                tr = -1
            self.theta1 = theta[0] * tr
            self.theta2 = theta[1] * tr
        if len(theta) == 3:
            tr = 1
            if theta[2] < 0:
                tr = -1
            self.theta1 = theta[0] * tr
            self.theta2 = theta[1] * tr
            self.theta3 = theta[2] * tr


def superimpose(superimposer, ref_atoms, atoms):
    superimposer.set_atoms(ref_atoms, atoms)
    rmsd1 = superimposer.rms
    superimposer.apply(atoms)

    for atom in atoms:
        atom.set_coord(atom.get_coord() * -1)
    superimposer.set_atoms(ref_atoms, atoms)
    rmsd2 = superimposer.rms
    if rmsd2 < rmsd1:
        superimposer.apply(atoms)
        rmsd = rmsd2
    else:
        for atom in atoms:
            atom.set_coord(atom.get_coord() * -1)
        rmsd = rmsd1
    return rmsd



def load_cycles(type_of_cycle):

    cycles = []
    for cycle_file in glob(f"data/{type_of_cycle}/filtered_ligands/*/*/*.pdb"):
        cycles.append(Cycle(cycle_file))
    return cycles


def load_templates(files):
    templates = {}
    for template_file in files:
        templates[template_file.split("/")[-1][:-4]] = Cycle(template_file)
    return templates


type_of_cycle = argv[1]
print(f"Selection of conformation for {type_of_cycle} cycles. ")
filtered_ligands_path = argv[2]
output_dir = argv[3]

sup = PDB.Superimposer()

cycles = []
for cycle_file in glob(f"{filtered_ligands_path}/*/*/*.pdb"):
    cycles.append(Cycle(cycle_file))
QM_templates = load_templates(glob(f"QM_optimised_templates/{type_of_cycle}/*.pdb"))
for cycle in cycles:
    item1 = cycle.file.split("/")[-1].split("_")[0]
    item2 = cycle.file.split("/")[-1].split(".")[0]
    print(f"Selection of conformation for {item1}, {item2}")
    cycle.rmsds = {}
    best_achieved_hr_distance = 1000 # hill-reilly
    for conformation in QM_templates.keys():
        QM_template = QM_templates[conformation]
        if type_of_cycle in ["cyclohexane", "benzene"]:
            hr_distance = math.dist((QM_template.theta1, QM_template.theta2, QM_template.theta3), (cycle.theta1, cycle.theta2, cycle.theta3))
        elif type_of_cycle == "cyclopentane":
            hr_distance = math.dist((QM_template.theta1, QM_template.theta2), (cycle.theta1, cycle.theta2))
        if hr_distance < best_achieved_hr_distance:
            best_achieved_hr_distance = hr_distance
            cycle.conformation = conformation
        cycle.rmsds[conformation] = superimpose(sup, QM_template.atoms, cycle.atoms)
with open(f"{output_dir}/results_individual_molecules.txt", "w") as output_file:
    for cycle in cycles:
        item1 = cycle.file.split("/")[-1].split("_")[0]
        item2 = cycle.file.split("/")[-1].split(".")[0]
        output_file.write(f"[{item1}] {item2}: {cycle.conformation.upper()}\n")
with open(f"{output_dir}/result_rmsd_chart.csv", "w") as output_file_rmsd:
    output_file_rmsd.write("Ligand_name;Ring_ID;" + ";".join([conformation.upper() for conformation in sorted(QM_templates.keys())]) + "\n")
    for cycle in cycles:
        item1 = cycle.file.split("/")[-1].split("_")[0]
        item2 = cycle.file.split("/")[-1].split(".")[0]
        output_file_rmsd.write(f"{item1};{item2};{';'.join([str(round(float(cycle.rmsds[conformation]), 3)) for conformation in sorted(QM_templates.keys())])}\n")
print(f"Selection of conformation for {type_of_cycle} cycles has completed successfully.")