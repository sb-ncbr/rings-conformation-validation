import pandas as pd

# RING = "cyclohexane"
RING = "cyclopentane"

def statistic_RMSD(ring):
    file_path = 'result_rmsd_chart.csv'
    output_file = f'{ring}_minimum.csv'

    #Read and clean data
    RMSD_data = pd.read_csv(file_path, delimiter=';')
    RMSD_data = RMSD_data.dropna(axis=1)

    #Add header
    new_headers = []
    if ring == "cyclopentane":
        new_headers = ['Ligand ID', 'Ring ID', 'HALF_CHAIR', 'FLAT', 'ENVELOPE'] #cyclopentane
    elif ring == "cyclohexane":
        new_headers = ['Ligand ID', 'Ring ID', 'HALF_CHAIR', 'FLAT', 'CHAIR', 'TW_BOAT_L','TW_BOAT_R', 'BOAT'] #cyclohexane
    else:
        print("Error: Wrong ring")
        return
    
    RMSD_data.columns = new_headers

    #find minimum
    values_columns = RMSD_data.columns[2:]
    RMSD_data['MinValue'] = RMSD_data[values_columns].min(axis=1)
    RMSD_data['Conformation'] = RMSD_data[values_columns].idxmin(axis=1)

    #Output
    RMSD_data.to_csv(output_file, sep=';', index=False)
    return RMSD_data
    
    
print(statistic_RMSD(ring=RING))


