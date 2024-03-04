import os
import pandas as pd

def count_files_in_directory(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith(extension)]

def read_results():
    all_results = []
    for dir_name in os.listdir('.'):
        if os.path.isdir(dir_name):
            directory = dir_name
            extension = "H_fretest.log"
            files = count_files_in_directory(directory, extension)

            for file_name in files:
                with open(os.path.join(directory, file_name), "r") as f:
                    row = {
                        "Molecule": dir_name,
                        "Conf_num": file_name.split("_")[1].split("H")[0],
                        "True_minimum": 0,  # Initialize to 0
                        "Multiple_minimum": 0,  # Initialize to 0
                        "Saddle_point": 0  # Initialize to 0
                    }
                    for line in f:
                        if "Total Energy =" in line:
                            energy = float(line.split("=")[1].split("kcal")[0].strip())
                            row["Total_Energy(kcal/mol)"] = energy
                        elif "True minimum" in line:
                            row["True_minimum"] = 1
                        elif "Multiple minimum" in line:
                            row["Multiple_minimum"] = 1
                        elif "Saddle point" in line:
                            row["Saddle_point"] = 1
                    all_results.append(row)

    return pd.DataFrame(all_results)



def final_results(df, filename):
    # Calculate lowest energy for each subdirectory
    df['Lowest_energy'] = df.groupby('Molecule')['Total_Energy(kcal/mol)'].transform('min')

    # Count conformers within energy ranges
    for kcal in [1, 3, 5]:
        df[f'Conformers_within_{kcal}kcal'] = df.apply(lambda x: x['Total_Energy(kcal/mol)'] <= x['Lowest_energy'] + kcal, axis=1)

    # Aggregate results
    agg_dict = {
        'Lowest_energy': 'min',
        'Total_Energy(kcal/mol)': 'mean',
        'True_minimum': 'sum',
        'Multiple_minimum': 'sum',
        'Saddle_point': 'sum',
        'Conformers_within_1kcal': 'sum',
        'Conformers_within_3kcal': 'sum',
        'Conformers_within_5kcal': 'sum'
    }
    agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}

    aggregated_df = df.groupby('Molecule').agg(agg_dict).reset_index()

    # Add the 'Minimum' column as the sum of 'True_minimum' and 'Multiple_minimum'
    aggregated_df['Minimum'] = aggregated_df['True_minimum'] + aggregated_df['Multiple_minimum']

    aggregated_df.to_csv(filename, index=False)

# Execute functions for fretest.log
df_fretest = read_results()
final_results(df_fretest, "Diamond_OPLS4_MacroModel_totalfre.csv")

