from Bio import PDB
import numpy as np
from glob import glob
from scipy.optimize import shgo
from numba import jit
from sys import argv
import math
from pathlib import Path
import json

import json
import math
from sys import argv

#compares HR angles of ligands to standard HR for cyclohexane, cyclopentane, benzene and oxane. Chooses the best conformation for each ligand
def rmsd_hr_angles(hr1, hr2, cycle_type):
    try:
        if cycle_type in ["cyclohexane", "benzene", "oxane"]:
            #get theta1 and theta2, fail if missing
            t1_1 = hr1["theta1"]
            t2_1 = hr1["theta2"]
            t1_2 = hr2["theta1"]
            t2_2 = hr2["theta2"]
            
            #theta3 might be None or missing, so use 0 if None or missing
            t3_1 = hr1.get("theta3") or 0
            t3_2 = hr2.get("theta3") or 0

            a1 = [float(t1_1), float(t2_1), float(t3_1)]
            a2 = [float(t1_2), float(t2_2), float(t3_2)]

        elif cycle_type == "cyclopentane":
            a1 = [float(hr1["theta1"]), float(hr1["theta2"])]
            a2 = [float(hr2["theta1"]), float(hr2["theta2"])]
        else:
            raise ValueError(f"Unsupported cycle type: {cycle_type}")

    except KeyError as e:
        raise ValueError(f"Missing angle key: {e}")
    except TypeError as e:
        raise ValueError(f"Missing or null angle value: {e}")
    except ValueError as e:
        raise ValueError(f"Invalid angle value: {e}")

    sqdiff = [(x - y) ** 2 for x, y in zip(a1, a2)]
    return math.sqrt(sum(sqdiff) / len(sqdiff))



if __name__ == "__main__":
    type_of_cycle = argv[1]
    ligands_hr_json_path = argv[2]
    output_csv_path = argv[3]

    #Load standard conformations HR angles
    with open(f"standardHR/{type_of_cycle}.json") as f:
        standard_HRs = json.load(f)

    #Load all ligand HR data 
    with open(ligands_hr_json_path) as f:
        ligands_hr = json.load(f)

    with open(output_csv_path, "w") as out_f:
        header = "compound;ring;" + ";".join([f"RMSD_{c}" for c in sorted(standard_HRs.keys())]) + ";BEST_CONFORMATION;theta1;theta2;theta3\n"
        out_f.write(header)

        for ligand_id, hr_angles in ligands_hr.items():
            best_rmsd = float("inf")
            best_conf = None
            rmsds = {}

            for conf_name, std_hr in standard_HRs.items():
                rmsd = rmsd_hr_angles(std_hr, hr_angles, type_of_cycle)
                rmsds[conf_name] = rmsd
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_conf = conf_name

            item1 = ligand_id.split("_")[0]
            item2 = ligand_id

            rmsd_values_str = ";".join([f"{rmsds[c]:.3f}" for c in sorted(standard_HRs.keys())])
            if type_of_cycle == "cyclopentane":
                line = f"{item1};{item2};{rmsd_values_str};{best_conf.upper()};{hr_angles['theta1']};{hr_angles['theta2']};\n"
            else:
                theta3 = hr_angles.get('theta3')
                theta3_str = "" if theta3 is None else theta3
                line = f"{item1};{item2};{rmsd_values_str};{best_conf.upper()};{hr_angles['theta1']};{hr_angles['theta2']};{theta3_str}\n"
            out_f.write(line)
