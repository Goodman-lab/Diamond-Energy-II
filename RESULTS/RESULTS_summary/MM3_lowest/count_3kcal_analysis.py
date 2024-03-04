import pandas as pd
from tqdm.auto import tqdm
import seaborn as sns

# Load the CSV file into a DataFrame
file_path = 'merged_conformers_5kcaldata.csv'
df = pd.read_csv(file_path, encoding='ascii')

# Display the head of the DataFrame
df_head = df.head()
print(df_head)

# Plotting distributions
import matplotlib.pyplot as plt

# Progress bar for plotting
for column in tqdm(df.columns[1:], desc='Plotting distributions'):
    plt.figure(figsize=(10, 6))
    df[column].hist(bins=20)
    plt.title('Distribution of ' + column)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.show()


# Set the aesthetic style of the plots
sns.set_style('whitegrid')

# Plotting the number of conformers for each method
plt.figure(figsize=(14, 8))

# To compare the methods, need to melt the dataframe to have a single column for method and one for the number of conformers
melted_df = df.melt(id_vars='Molecule', var_name='Method', value_name='Number of Conformers')

# Create a boxplot to compare the number of conformers identified by each method
sns.boxplot(x='Method', y='Number of Conformers', data=melted_df)
plt.title('Comparison of the Number of Conformers Identified by Each Method')
plt.xticks(rotation=45)
plt.show()

# Calculate the mean number of conformers for each method
mean_conformers = melted_df.groupby('Method')['Number of Conformers'].mean().sort_values(ascending=False)
print(mean_conformers)

