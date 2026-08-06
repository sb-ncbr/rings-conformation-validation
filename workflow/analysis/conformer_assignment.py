from workflow.analysis.calculate_HR import calc_HR_for_all_rings
from workflow.analysis.compare_HR import compare_HR


def assign_conformations(main_dir):
    calc_HR_for_all_rings(main_dir)
    compare_HR(main_dir)