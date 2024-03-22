import pandas as pd

def addResolution():
    file1_path = 'CP_minimum.csv'
    file2_path = 'CP_PDB_INFORMATION.csv'
    output_file = 'CP_RMSD_resolution.csv'
    trash_file = 'trash_resolution_CP.csv'

    # Read input CSV files into DataFrames
    file1 = pd.read_csv(file1_path, delimiter=';')
    file2 = pd.read_csv(file2_path,  sep=',', quoting=3, quotechar='"', error_bad_lines=False, warn_bad_lines=True)
    print(file2)
    #Extract Entry ID from Ring ID in file1 using regex
    file1['Entry ID_x'] = file1['Ring ID'].str.extract(r'_(\w+)_\d+')
    file1['Entry ID'] = file1['Entry ID_x'].str.upper()
    

    # Function to strip quoting from values
    def strip_quoting(value):
        if isinstance(value, str):
            return value.strip('"')
        return value
    # Apply the function to each element in the DataFrames
    file2 = file2.applymap(strip_quoting)
    file2['Resolution (Å)'] = file2['Resolution (Å)'].str.extract('(\d+\.\d+)').astype(float)
    print(file2)

    
    # Merge the two DataFrames on the common column "Entry ID"
    merged_df = pd.merge(file1, file2[['Entry ID','Experimental Method', 'Deposition Date', 'Resolution (Å)']],
                     how='left', left_on='Entry ID', right_on='Entry ID')

    # Check for entries in file1 that do not have a match in file2
    unmatched_entries = file1[~file1['Entry ID'].isin(file2['Entry ID'])]

    # Save unmatched entries to the trash file
    unmatched_entries.to_csv(trash_file, sep=';', index=False)

    # Delete the 'Entry ID_x' column
    merged_df.drop(['Entry ID_x'], axis=1, inplace=True)
    
    #Output
    merged_df.to_csv(output_file, sep=';', index=False)
    return merged_df
    
    
print(addResolution())