import pandas as pd

def addResolution():
    file1_path = 'CP_RMSD_resolution.csv'
    file2_path = 'ED_cyclopentane_params-s.csv'
    output_file = 'CP_RMSD_Res_ED_s.csv'
    trash_file = 'trash_ED_CP-s.csv'
    
    # Read input CSV files into DataFrames
    file1 = pd.read_csv(file1_path, delimiter=';')
    file2 = pd.read_csv(file2_path, delimiter=';')
    
    new_headers = ['Entry', 'Atoms in ring'] #cyclopentane
    file2 .columns = new_headers

    #Extract Entry ID from Ring ID in file1, Entry ID and Ligand ID and coverage from Entry in file2 using regex
    file1['Entry ID'] = file1['Ring ID'].str.split('_', n=1).str[1]
    
    # Extract Entry ID, Ligand ID, and Coverage from Entry in file2 using regex
    file2[['Entry ID', 'Ligand ID', 'Coverage']] = file2['Entry'].str.split(',', expand=True)
    

    # Merge the two DataFrames on the common column "Entry ID and Ligand ID"
    merged_df = pd.merge(file1, file2[['Entry ID','Ligand ID', 'Coverage']],
                     how='left', left_on=['Entry ID', 'Ligand ID'], right_on=['Entry ID', 'Ligand ID'])

    # Check for entries in file1 that do not have a match in file2
    unmatched_entries = file1[~file1['Entry ID'].isin(file2['Entry ID'])]

    # Save unmatched entries to the trash file
    unmatched_entries.to_csv(trash_file, sep=';', index=False)

    # Delete the 'Entry ID_x' column and duplicit
    #merged_df = merged_df[~merged_df['Entry ID'].isin(unmatched_entries['Entry ID'])]
    merged_df.drop(['Entry ID'], axis=1, inplace=True)
    merged_df.drop_duplicates(subset=['Ring ID'], inplace=True)

    #Output
    merged_df.to_csv(output_file, sep=';', index=False)
    return merged_df
    
print(addResolution())
