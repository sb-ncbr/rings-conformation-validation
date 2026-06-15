from sys import argv
import math
from pathlib import Path
import json
from HelperModule.Ring import Ring
from HelperModule.constants import MAIN_DIR

SIX_MEMBERED_RINGS = {Ring.CYCLOHEXANE, Ring.BENZENE, Ring.OXANE}
FIVE_MEMBERED_RINGS = (Ring.CYCLOPENTANE, Ring.OXOLANE)

def dist_hr_angles(hr1, hr2, ring):
    try:
        if ring in SIX_MEMBERED_RINGS:
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
            return math.dist(a1, a2)

        # the rest: FIVE_MEMBERED_RINGS
        a1 = [float(hr1["theta1"]), float(hr1["theta2"])]
        a2 = [float(hr2["theta1"]), float(hr2["theta2"])]
        return math.dist(a1, a2)

    except KeyError as e:
        raise ValueError(f"Missing angle key: {e}")
    except TypeError as e:
        raise ValueError(f"Missing or null angle value: {e}")
    except ValueError as e:
        raise ValueError(f"Invalid angle value: {e}")


def main(input_dir):

    for ring in Ring:
        #Load standard conformations HR angles
        with open(f"standardHR/{ring.name.lower()}.json") as f:
            standard_HRs = json.load(f)

        #Load all ligand HR data
        hr_analysis_output = Path(input_dir) / MAIN_DIR / ring.name.lower() / "hr_analysis_output"
        with open(hr_analysis_output / "output_HR.json") as f:
            rings_hr = json.load(f)

        with open(hr_analysis_output / "result_conf_chart.csv", "w") as out_f:
            header = "Ligand_name;Ring_ID;Conformation;theta1;theta2;theta3;dist_to_ideal;dist_to_second_best;fit_ratio\n"
            out_f.write(header)

            for ring_id, hr_angles in rings_hr.items():
                best_dist = float("inf")
                best_conf = None
                hr_dists = {}

                for conf_name, std_hr in standard_HRs.items():
                    hr_dist = dist_hr_angles(std_hr, hr_angles, ring)
                    hr_dists[conf_name] = hr_dist

                    if hr_dist < best_dist:
                        best_dist = hr_dist
                        best_conf = conf_name

                sorted_dists = sorted(hr_dists.values())
                second_dist = sorted_dists[1]
                fit_ratio = best_dist / second_dist if second_dist != 0 else ""

                ligand_id = ring_id.split("_")[0]

                if ring in FIVE_MEMBERED_RINGS:
                    line = f"{ligand_id};{ring_id};{best_conf.upper()};{hr_angles['theta1']};{hr_angles['theta2']};{best_dist};{second_dist};{fit_ratio}\n"
                else:
                    theta3 = hr_angles.get('theta3')
                    theta3_str = "" if theta3 is None else theta3
                    line = f"{ligand_id};{ring_id};{best_conf.upper()};{hr_angles['theta1']};{hr_angles['theta2']};{theta3_str};{best_dist};{second_dist};{fit_ratio}\n"
                out_f.write(line)


if __name__ == "__main__":
    workflow_output = argv[1]
    main(workflow_output)

