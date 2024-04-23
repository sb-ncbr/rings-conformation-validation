import pandas as pd
import argparse
import os


def statistic_RMSD(ring_type, base_dir):
    filename = "result_rmsd_chart.csv"
    file_path = os.path.join(base_dir, f'{ring_type}', "output", filename)
    #output_file = f'{ring_type}_minimum.csv'
    RMSD_data = pd.read_csv(file_path, delimiter=';', header=0)
    RMSD_data = RMSD_data.dropna(axis=1)

    values_columns = RMSD_data.columns[2:]
    RMSD_data['MinValue'] = RMSD_data[values_columns].min(axis=1)
    RMSD_data['Conformation'] = RMSD_data[values_columns].idxmin(axis=1)
    #RMSD_data.to_csv(output_file, sep=';', index=False)

    return RMSD_data


def addResolution(data_dir, ring_type, RMSD_data):
    file_information_path = os.path.join(data_dir, 'PDB_INFORMATION.csv')
    #output_file = f'{ring_type}_RMSD_resolution.csv'
    trash_file = f'trash_resolution_{ring_type}.csv'
    file2 = pd.read_csv(file_information_path, sep=';')

    RMSD_data['Entry ID_x'] = RMSD_data['Ring_ID'].str.extract(r'_(\w+)_\d+')
    RMSD_data['Entry ID'] = RMSD_data['Entry ID_x'].str.upper()

    merged_resolution = pd.merge(RMSD_data,
                                 file2[['Entry ID', 'Experimental Method', 'Deposition Date', 'Resolution (Å)']],
                                 how='left', left_on='Entry ID', right_on='Entry ID')

    unmatched_entries = RMSD_data[~RMSD_data['Entry ID'].isin(file2['Entry ID'])]
    unmatched_entries.to_csv(trash_file, sep=';', index=False)
    merged_resolution.drop(['Entry ID_x'], axis=1, inplace=True)

    xray_merged_resolution = pd.DataFrame(columns=merged_resolution.columns)
    for index, row in merged_resolution.iterrows():
        methods = row['Experimental Method']
        if isinstance(methods, str) and 'X-RAY DIFFRACTION' in methods.split(','):
            xray_merged_resolution = pd.concat([xray_merged_resolution, row.to_frame().transpose()], ignore_index=True)
    #xray_merged_resolution.to_csv(output_file, sep=';', index=False)

    return xray_merged_resolution


def addElDensity(ring_type, xray_merged_resolution):
    el_dens_result_dir = "electron_density_coverage_analysis"
    file2_path = os.path.join(el_dens_result_dir, "output", f'{ring_type}_params__analysis_output.csv')
    #output_file = f'{ring_type}_RMSD_Res_ED.csv'
    trash_file = f'trash_ED_{ring_type}.csv'

    file2 = pd.read_csv(file2_path, delimiter=';')
    new_headers = ['Entry', 'Atoms in ring']
    file2.columns = new_headers

    file2[['Ring_ID', 'Ligand_name', 'Coverage']] = file2['Entry'].str.split(',', expand=True)
    merged_coverage = pd.merge(xray_merged_resolution, file2[['Ring_ID', 'Coverage']],
                               how='left', on=['Ring_ID'])

    unmatched_entries = xray_merged_resolution[~xray_merged_resolution['Ring_ID'].isin(file2['Ring_ID'])]
    unmatched_entries.to_csv(trash_file, sep=';', index=False)
    #merged_coverage.to_csv(output_file, sep=';', index=False)

    return merged_coverage


def Summary(base_output_dir, ring_type, merged_coverage):
    df = merged_coverage
    output_folder = os.path.join(base_output_dir, ring_type, 'final_results')
    os.makedirs(output_folder, exist_ok=True)

    if ring_type == 'cyclopentane':
        x = '5'
        Sheet_1 = 'stats_resolution_2_or_less_coverage_5'
        Sheet_2 = 'Resolution_2_or_less_Coverage_5'
        Sheet_3 = 'Stat_Resol_2_or_less_Cover_5'
    elif ring_type == 'cyclohexane' or ring_type == 'benzene':
        x = '6'
        Sheet_1 = 'stats_resolution_2_or_less_coverage_6'
        Sheet_2 = 'Resolution_2_or_less_Coverage_6'
        Sheet_3 = 'Stat_Resol_2_or_less_Cover_6'

    def check_resolution(resolution):
        if isinstance(resolution, float):
            return resolution <= 2
        values = resolution.split(',')
        for value in values:
            try:
                if float(value.strip()) <= 2:
                    return True
            except ValueError:
                pass
        return False

    # Create a CSV file for rows where Resolution (Å) is equal or less than 2
    output_file_path_1 = os.path.join(output_folder, 'resolution_2_or_less.csv')
    df_resolution_2_or_less = df[df['Resolution (Å)'].apply(check_resolution)]
    df_resolution_2_or_less.to_csv(output_file_path_1, sep=';', index=False)

    # Create a CSV file for rows where Resolution (Å) is equal or less than 2 and Coverage is 5 or 6
    output_file_path_2 = os.path.join(output_folder, f'resolution_2_or_less_coverage_{x}.csv')
    df_resolution_2_or_less_coverage_x = df_resolution_2_or_less[df_resolution_2_or_less['Coverage'] == x]
    df_resolution_2_or_less_coverage_x.to_csv(output_file_path_2, sep=';', index=False)

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
    stats_input_file_df.to_csv(os.path.join(output_folder, 'stats_resolution_coverage.txt'), sep=' ', index=False,
                               float_format='%.2f')

    stats_resolution_2_or_less = df_resolution_2_or_less['Conformation'].value_counts()
    total_values_resolution_2_or_less = stats_resolution_2_or_less.sum()
    percentage_resolution_2_or_less = (stats_resolution_2_or_less / total_values_resolution_2_or_less) * 100
    stats_resolution_2_or_less_df = pd.DataFrame({
        'Conformation': stats_resolution_2_or_less.index,
        'Occurrences': stats_resolution_2_or_less.values,
        'Percentage': percentage_resolution_2_or_less.values
    })
    stats_resolution_2_or_less_df.loc[len(stats_resolution_2_or_less_df)] = ['Total', total_values_resolution_2_or_less,
                                                                             100]
    stats_resolution_2_or_less_df.to_csv(os.path.join(output_folder, 'stats_resolution_2_or_less.txt'), sep=' ',
                                         index=False, float_format='%.2f')

    stats_resolution_2_or_less_coverage_x = df_resolution_2_or_less_coverage_x['Conformation'].value_counts()
    total_values_resolution_2_or_less_coverage_x = stats_resolution_2_or_less_coverage_x.sum()
    percentage_resolution_2_or_less_coverage_x = (
                                                             stats_resolution_2_or_less_coverage_x / total_values_resolution_2_or_less_coverage_x) * 100
    stats_resolution_2_or_less_coverage_x_df = pd.DataFrame({
        'Conformation': stats_resolution_2_or_less_coverage_x.index,
        'Occurrences': stats_resolution_2_or_less_coverage_x.values,
        'Percentage': percentage_resolution_2_or_less_coverage_x.values
    })
    stats_resolution_2_or_less_coverage_x_df.loc[len(stats_resolution_2_or_less_coverage_x_df)] = ['Total',
                                                                                                   total_values_resolution_2_or_less_coverage_x,
                                                                                                   100]
    stats_resolution_2_or_less_coverage_x_df.to_csv(os.path.join(output_folder, Sheet_1), sep=' ', index=False,
                                                    float_format='%.2f')

    # Create an Excel file with multiple sheets
    excel_file_path = os.path.join(output_folder, 'result_summary.xlsx')
    with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Summary', index=False)
        df_resolution_2_or_less.to_excel(writer, sheet_name='Resolution_2_or_less', index=False)
        df_resolution_2_or_less_coverage_x.to_excel(writer, sheet_name=Sheet_2, index=False)
        stats_input_file_df.to_excel(writer, sheet_name='Stat_summary', index=False)
        stats_resolution_2_or_less_df.to_excel(writer, sheet_name='Stat_Resol_2_or_less', index=False)
        stats_resolution_2_or_less_coverage_x_df.to_excel(writer, sheet_name=Sheet_3, index=False)
    return excel_file_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate RMSD statistics for CP or CH rings.')
    parser.add_argument('-r', '--ring', choices=['cyclopentane', 'cyclohexane', 'benzene'], required=True,
                        help='Specify the ring type (cyclopentane, cyclohexane or benzene)')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="Path to the USER's main output directory. This script requires data in that directory "
                             "for analysis. That dir is the same as in th previous steps of the workflow.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to the directory with input data (local pdb, ccp4 files, etc.)')
    args = parser.parse_args()

    base_dir = os.path.join(args.output, "validation_data")

    # Call statistic_RMSD with the specified ring type
    rmsd_result = statistic_RMSD(args.ring, base_dir)

    # Call addResolution
    resolution_result = addResolution(args.input, args.ring, rmsd_result)

    # Call addElDensity
    coverage_result = addElDensity(args.ring, resolution_result)

    # Call Summary
    excel_file_path = Summary(base_dir, args.ring, coverage_result)
