import pandas as pd
from tqdm.auto import tqdm
import re
from natsort import natsorted, ns

# Load the files into dataframes and display the head of each.
file_paths = ['MCMM_MM3_totalfre.csv', 'MCLMCS_MM3_totalfre.csv', 'Diamond_RDKit_MM3_totalfre.csv', 'LMCS_MM3_totalfre.csv', 'Diamond_MacroModel_MM3_totalfre.csv', 'SPMC_MM3_totalfre.csv']
dataframes = {}

for file in tqdm(file_paths, desc='Loading files'):
    df = pd.read_csv(file)
    dataframes[file] = df
    print(f'Head of {file}:')
    print(df.head())
    print('\n')



# Define a function to extract the 'Conformers_within_3kcal' column from each dataframe
# and rename it to match the file name without '_MM3_totalfre.csv'
def extract_and_rename(df, file_name):
    new_column_name = re.sub('_MM3_totalfre.csv', '', file_name)
    return df[['Molecule', 'Conformers_within_3kcal']].rename(columns={'Conformers_within_3kcal': new_column_name})

# Apply the function to each dataframe and merge them on 'Molecule'
merged_df = pd.DataFrame()
for file, df in tqdm(dataframes.items(), desc='Processing files'):
    if 'Molecule' in df.columns and 'Conformers_within_3kcal' in df.columns:
        temp_df = extract_and_rename(df, file)
        if merged_df.empty:
            merged_df = temp_df
        else:
            merged_df = pd.merge(merged_df, temp_df, on='Molecule', how='outer')

# Sort the 'Molecule' column naturally
merged_df['Molecule'] = pd.Categorical(merged_df['Molecule'], ordered=True, categories=natsorted(merged_df['Molecule'].unique(), alg=ns.IGNORECASE))
merged_df.sort_values('Molecule', inplace=True)

# Save the merged dataframe to a new CSV file
merged_csv_file = 'merged_conformers_3kcaldata.csv'
merged_df.to_csv(merged_csv_file, index=False)

# Display the head of the merged dataframe
print(merged_df.head())

