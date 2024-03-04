import pandas as pd

# Dictionary to hold dataframes
dfs = {}

# List of method file names
method_files = [
    'MCMM_MM3_totalfre.csv',
    'MCLMCS_MM3_totalfre.csv',
    'LMCS_MM3_totalfre.csv',
    'Diamond_MM3_RDKit_totalfre.csv',
    'SPMC_MM3_totalfre.csv',
    'Diamond_MM3_MacroModel_totalfre.csv'
]

# Load the dataframes
for file in method_files:
    df = pd.read_csv(file)
    method_name = file.replace('_totalfre.csv', '')
    df['Method'] = method_name
    dfs[method_name] = df[['Molecule', 'Lowest_Total_Energy', 'Method']]

# Combine all data into one dataframe
combined_data = pd.concat(dfs.values()).reset_index(drop=True)

# Find the method with the lowest energy for each molecule
combined_data['Rank'] = combined_data.groupby('Molecule')['Lowest_Total_Energy'].rank(method='min')
combined_data_lowest = combined_data[combined_data['Rank'] == 1]

# Count the number of times each method finds the lowest energy
lowest_energy_counts = combined_data_lowest['Method'].value_counts()

# Calculate descriptive statistics for each method
stats = combined_data.groupby('Method')['Lowest_Total_Energy'].describe()

# Add this count to the statistics table
stats['Lowest_Energy_Count'] = stats.index.map(lowest_energy_counts)

# Save the final table as a CSV file
#stats.to_csv('final_statistics.csv')

print(stats)

import scipy.stats as stats

# Perform ANOVA test
anova_result = stats.f_oneway(
    combined_data[combined_data['Method'] == 'MCMM_MM3']['Lowest_Total_Energy'],
    combined_data[combined_data['Method'] == 'MCLMCS_MM3']['Lowest_Total_Energy'],
    combined_data[combined_data['Method'] == 'LMCS_MM3']['Lowest_Total_Energy'],
    combined_data[combined_data['Method'] == 'Diamond_MM3_RDKit']['Lowest_Total_Energy'],
    combined_data[combined_data['Method'] == 'SPMC_MM3']['Lowest_Total_Energy'],
    combined_data[combined_data['Method'] == 'Diamond_MM3_MacroModel']['Lowest_Total_Energy']
)

# Display the ANOVA test result
print(anova_result)


'''
import matplotlib.pyplot as plt
import seaborn as sns

# Custom order for the methods
method_order = [
    'Diamond_MM3_RDKit',
    'Diamond_MM3_MacroModel',
    'MCMM_MM3',
    'MCLMCS_MM3',
    'LMCS_MM3',
    'SPMC_MM3'
]

# Set the aesthetic style of the plots
sns.set_style('whitegrid')

# Create a boxplot to compare the distribution of 'Lowest_Total_Energy' for each method
plt.figure(figsize=(10, 6))
boxplot = sns.boxplot(x='Method', y='Lowest_Total_Energy', data=combined_data, order=method_order)
boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=45)
plt.title('Boxplot of Lowest Total Energy by Method')
plt.xlabel('Method')
plt.ylabel('Lowest Total Energy (kcal/mol)')
plt.tight_layout()

# Show the plot
plt.show()

'''













'''














## plot scatter



# Check the unique 'Source' values in the data
combined_df['Source'].unique()

# Correct the color palette
#palette = {'Diamond_MM3_MacroModel': 'red', 'Diamond_MM3_RDKit': 'blue', 'LMCS_MM3': 'lightgray', 'SPMC_MM3': 'lightgray', 'MCMM_MM3': 'lightgray', 'MCLMCS_MM3': 'lightgray'}
palette = {'Diamond_MM3_MacroModel': 'red', 'Diamond_MM3_RDKit': 'blue', 'LMCS_MM3': 'yellow', 'SPMC_MM3': 'green', 'MCMM_MM3': 'brown', 'MCLMCS_MM3': 'purple'}


# Step 1: Split the combined dataframe into two
# Diamond series dataframe
diamond_df = combined_df[combined_df['Source'].str.contains('Diamond_MM3')]

# Other sources dataframe
other_df = combined_df[~combined_df['Source'].str.contains('Diamond_MM3')]

# Step 2: Plot the data for the other sources first
plt.figure(figsize=(20, 10))
sns.scatterplot(data=other_df, x='Molecule', y='Lowest_Total_Energy', hue='Source', palette=palette)

# Step 3: Plot the Diamond series data on top
sns.scatterplot(data=diamond_df, x='Molecule', y='Lowest_Total_Energy', hue='Source', palette=palette)

# Rotate the x-axis labels for better visibility
plt.xticks(rotation=90)

# Set the title and labels
plt.title('Lowest_Total_Energy for each Molecule by Source')
plt.xlabel('Molecule')
plt.ylabel('Lowest_Total_Energy (kcal/mol)')

# Show the plot
plt.show()
'''
