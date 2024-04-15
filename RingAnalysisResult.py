import pandas as pd
import argparse
import os

def statistic_RMSD(ring_type):
    file_path = f'{ring_type}_rmsd_chart.csv'
    output_file = 'CH_minimum.csv' #TODO: correct this name

    RMSD_data = pd.read_csv(file_path, delimiter=';')
    RMSD_data = RMSD_data.dropna(axis=1)

    # Add header
    if ring_type == 'cyclopentane':
        new_headers = ['Ligand ID', 'Ring ID', 'HALF_CHAIR', 'FLAT', 'ENVELOPE']
    elif ring_type == 'cyclohexane' or 'benzene':
        new_headers = ['Ligand ID', 'Ring ID', 'HALF_CHAIR', 'FLAT', 'CHAIR', 'TW_BOAT_L','TW_BOAT_R', 'BOAT'] 
    else:
        raise ValueError("Invalid ring type. Please specify cyclopentane, cyclohexane or benzene.")
    RMSD_data.columns = new_headers

    # Find minimum
    values_columns = RMSD_data.columns[2:]
    RMSD_data['MinValue'] = RMSD_data[values_columns].min(axis=1)
    RMSD_data['Conformation'] = RMSD_data[values_columns].idxmin(axis=1)

    # Output
    RMSD_data.to_csv(output_file, sep=';', index=False)
    return RMSD_data

def addResolution(RMSD_data):
    file_informatin_path = 'CH_PDB_INFORMATION.csv'
    output_file = 'CH_RMSD_resolution.csv'
    trash_file = 'trash_resolution_CH.csv'

    file2 = pd.read_csv(file_informatin_path,  sep=',', quoting=3, quotechar='"', error_bad_lines=False, warn_bad_lines=True)

    # Extract Entry ID from Ring ID and clean column
    RMSD_data['Entry ID_x'] = RMSD_data['Ring ID'].str.extract(r'_(\w+)_\d+')
    RMSD_data['Entry ID'] = RMSD_data['Entry ID_x'].str.upper()
    def strip_quoting(value):
        if isinstance(value, str):
            return value.strip('"')
        return value
    file2 = file2.applymap(strip_quoting)
    file2['Resolution (Å)'] = file2['Resolution (Å)'].str.extract('(\d+\.\d+)').astype(float)

    # Merge the two DataFrames on the common column "Entry ID"
    merged_resolution = pd.merge(RMSD_data, file2[['Entry ID','Experimental Method', 'Deposition Date', 'Resolution (Å)']],
                     how='left', left_on='Entry ID', right_on='Entry ID')

    # Check for entries in file1 that do not have a match in file2
    unmatched_entries = RMSD_data[~RMSD_data['Entry ID'].isin(file2['Entry ID'])]

    # Save unmatched entries to the trash file
    unmatched_entries.to_csv(trash_file, sep=';', index=False)

    # Delete the 'Entry ID_x' column
    merged_resolution.drop(['Entry ID_x'], axis=1, inplace=True)

    # Output
    merged_resolution.to_csv(output_file, sep=';', index=False)
    return merged_resolution
    
def addElDensity(merged_resolution):
    file2_path = 'cyclohexane_params-sm_analysis_output.csv'
    output_file = 'CH_RMSD_Res_ED_s.csv'
    trash_file = 'trash_ED_CH-s.csv'
    
    file2 = pd.read_csv(file2_path, delimiter=';')
    new_headers = ['Entry', 'Atoms in ring']
    file2.columns = new_headers
    
    # Extract Ring ID, Ligand ID, and Coverage from Entry in file2 using regex
    file2[['Ring ID', 'Ligand ID', 'Coverage']] = file2['Entry'].str.split(',', expand=True)

    # Merge the two DataFrames on the common column "Entry ID and Ligand ID"
    merged_coverage = pd.merge(merged_resolution, file2[['Ring ID','Coverage']],
                     how='left', on=['Ring ID'])

    # Check for entries in file1 that do not have a match in file2
    unmatched_entries = merged_resolution[~merged_resolution['Ring ID'].isin(file2['Ring ID'])]

    # Save unmatched entries to the trash file
    unmatched_entries.to_csv(trash_file, sep=';', index=False)

    #Output
    merged_coverage.to_csv(output_file, sep=';', index=False)
    return merged_coverage

def Summary(ring_type, merged_coverage):
    df = merged_coverage

    # Create a folder for the output files
    output_folder = output_name
    os.makedirs(output_folder, exist_ok=True)

    if ring_type == 'cyclopentane':
        x = 5
        name_1 = 'stats_resolution_2_or_less_coverage_5.txt'
        name_2 = 'Resolution_2_or_less_Coverage_5'
        name_3 = 'Stat_Resol_2_or_less_Cover_5'
        output_name = 'Cyclopentane'
    elif ring_type == 'cyclohexane':
        x = 6
        name_1 = 'stats_resolution_2_or_less_coverage_6.txt'
        name_2 = 'Resolution_2_or_less_Coverage_6'
        name_3 = 'Stat_Resol_2_or_less_Cover_6'
        output_name = 'Cyclohexane'
    elif ring_type == 'benzene':
        x = 6
        name_1 = 'stats_resolution_2_or_less_coverage_6.txt'
        name_2 = 'Resolution_2_or_less_Coverage_6'
        name_3 = 'Stat_Resol_2_or_less_Cover_6'
        output_name = 'Benzene'
    else:
        raise ValueError("Invalid ring type. Please specify cyclopentane, cyclohexane or benzene.")

    # Create a CSV file for rows where Resolution (Å) is equal or less than 2
    output_file_path_1 = os.path.join(output_folder, 'resolution_2_or_less.csv')
    df_resolution_2_or_less = df[df['Resolution (Å)'] <= 2]
    df_resolution_2_or_less.to_csv(output_file_path_1, sep=';', index=False)

    # Create a CSV file for rows where Resolution (Å) is equal or less than 2 and Coverage is 5 or 6
    output_file_path_2 = os.path.join(output_folder, 'resolution_2_or_less_coverage_6.csv')
    df_resolution_2_or_less_coverage_5 = df[(df['Resolution (Å)'] <= 2) & (df['Coverage'] == x)] 
    df_resolution_2_or_less_coverage_5.to_csv(output_file_path_2, sep=';', index=False)

    # Create text files with statistics for different values in the Conformation column
    stats_input_file = df['Conformation'].value_counts()
    total_values_input_file = stats_input_file.sum()
    percentage_input_file = (stats_input_file / total_values_input_file) * 100
    stats_input_file_df = pd.DataFrame({
        'Conformation': stats_input_file.index,
        'Occurrences': stats_input_file.values,
        'Percentage': percentage_input_file.values
    })
    stats_input_file_df.loc[len(stats_input_file_df)] = ['Total', total_values_input_file, 100]
    stats_input_file_df.to_csv(os.path.join(output_folder, 'stats_resolution_coverage.txt'), sep=' ', index=False, float_format='%.2f')

    stats_resolution_2_or_less = df_resolution_2_or_less['Conformation'].value_counts()
    total_values_resolution_2_or_less = stats_resolution_2_or_less.sum()
    percentage_resolution_2_or_less = (stats_resolution_2_or_less / total_values_resolution_2_or_less) * 100
    stats_resolution_2_or_less_df = pd.DataFrame({
        'Conformation': stats_resolution_2_or_less.index, 
        'Occurrences': stats_resolution_2_or_less.values,
        'Percentage': percentage_resolution_2_or_less.values
    })
    stats_resolution_2_or_less_df.loc[len(stats_resolution_2_or_less_df)] = ['Total', total_values_resolution_2_or_less, 100]
    stats_resolution_2_or_less_df.to_csv(os.path.join(output_folder, 'stats_resolution_2_or_less.txt'), sep=' ', index=False, float_format='%.2f')

    stats_resolution_2_or_less_coverage_5 = df_resolution_2_or_less_coverage_5['Conformation'].value_counts()
    total_values_resolution_2_or_less_coverage_5 = stats_resolution_2_or_less_coverage_5.sum()
    percentage_resolution_2_or_less_coverage_5 = (stats_resolution_2_or_less_coverage_5 / total_values_resolution_2_or_less_coverage_5) * 100
    stats_resolution_2_or_less_coverage_5_df = pd.DataFrame({
        'Conformation': stats_resolution_2_or_less_coverage_5.index, 
        'Occurrences': stats_resolution_2_or_less_coverage_5.values,
        'Percentage': percentage_resolution_2_or_less_coverage_5.values
    })
    stats_resolution_2_or_less_coverage_5_df.loc[len(stats_resolution_2_or_less_coverage_5_df)] = ['Total', total_values_resolution_2_or_less_coverage_5, 100]
    stats_resolution_2_or_less_coverage_5_df.to_csv(os.path.join(output_folder, name_1), sep=' ', index=False, float_format='%.2f')

    # Create an Excel file with multiple sheets
    excel_file_path = os.path.join(output_folder, 'result_summary.xlsx') 
    with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Summary', index=False)
        df_resolution_2_or_less.to_excel(writer, sheet_name='Resolution_2_or_less', index=False)
        df_resolution_2_or_less_coverage_5.to_excel(writer, sheet_name=name_2, index=False)
        stats_input_file_df.to_excel(writer, sheet_name='Stat_summary', index=False)
        stats_resolution_2_or_less_df.to_excel(writer, sheet_name='Stat_Resol_2_or_less', index=False)
        stats_resolution_2_or_less_coverage_5_df.to_excel(writer, sheet_name=name_3, index=False)
    return excel_file_path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate RMSD statistics for CP or CH rings.')
    parser.add_argument('-r', '--ring', choices=['cyclopentane', 'cyclohexane', 'benzene'], help='Specify the ring type (cyclopentane, cyclohexane or benzene)')
    
    args = parser.parse_args()
    
    # Call statistic_RMSD with the specified ring type
    rmsd_result = statistic_RMSD(args.ring)
    
    # Call addResolution
    resolution_result = addResolution(rmsd_result)
    
    # Call addElDensity
    coverage_result = addElDensity(resolution_result)

    # Call Summary
    excel_file_path = Summary(args.ring, coverage_result)
