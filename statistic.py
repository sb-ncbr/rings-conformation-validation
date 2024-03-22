import os
import pandas as pd

#for cyclohexane and benzene, replace 5 with 6

# Read the input CSV file
input_file_path = 'CH_RMSD_Res_ED_s.csv'
df = pd.read_csv(input_file_path, delimiter=';')

# Create a folder for the output files
output_folder = 'Cyclohexane_s'
os.makedirs(output_folder, exist_ok=True)

# Create a CSV file for rows where Resolution (Å) is equal or less than 2
output_file_path_1 = os.path.join(output_folder, 'resolution_2_or_less.csv')
df_resolution_2_or_less = df[df['Resolution'] <= 2]
df_resolution_2_or_less.to_csv(output_file_path_1, sep=';', index=False)

# Create a CSV file for rows where Resolution (Å) is equal or less than 2 and Coverage is 5 or 6
output_file_path_2 = os.path.join(output_folder, 'resolution_2_or_less_coverage_6.csv')
#df_resolution_2_or_less_coverage_5 = df[(df['Resolution'] <= 2) & (df['Coverage'] == 5)] #cyclopentane
df_resolution_2_or_less_coverage_5 = df[(df['Resolution'] <= 2) & (df['Coverage'] == 6)] #yclohexane
df_resolution_2_or_less_coverage_5.to_csv(output_file_path_2, sep=';', index=False)

# Create text files with statistics for different values in the Conformation column
# Statistics for the input file
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
stats_resolution_2_or_less_coverage_5_df.to_csv(os.path.join(output_folder, 'stats_resolution_2_or_less_coverage_6.txt'), sep=' ', index=False, float_format='%.2f')


# Create an Excel file with multiple sheets
excel_file_path = os.path.join(output_folder, 'result_summary.xlsx') 
with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
    df.to_excel(writer, sheet_name='Summary', index=False)
    df_resolution_2_or_less.to_excel(writer, sheet_name='Resolution_2_or_less', index=False)
    df_resolution_2_or_less_coverage_5.to_excel(writer, sheet_name='Resolution_2_or_less_Coverage_6', index=False)
    stats_input_file_df.to_excel(writer, sheet_name='Stat_summary', index=False)
    stats_resolution_2_or_less_df.to_excel(writer, sheet_name='Stat_Resol_2_or_less', index=False)
    stats_resolution_2_or_less_coverage_5_df.to_excel(writer, sheet_name='Stat_Resol_2_or_less_Cover_6', index=False)