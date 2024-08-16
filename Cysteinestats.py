import pandas as pd
import numpy as np

# Define the path to the Excel file containing cysteine counts
excel_path = '/content/cysteine_counts.xlsx'

# Read the Excel file into a DataFrame
df = pd.read_excel(excel_path)

# Ensure columns are read correctly
print("Columns in DataFrame:", df.columns)

# Prepare the list of cysteine counts by expanding the counts
all_cysteines = []
for index, row in df.iterrows():
    cysteine_count = int(row['Cysteine Residues'])  # Ensure integer type
    number_of_proteins = int(row['Number of Proteins'])
    all_cysteines.extend([cysteine_count] * number_of_proteins)

# Convert the list to a numpy array for statistical calculations
all_cysteines = np.array(all_cysteines)

# Calculate mean and median
mean_cysteines = np.mean(all_cysteines)
median_cysteines = np.median(all_cysteines)

# Output the results
print(f"Mean number of cysteine residues: {mean_cysteines:.2f}")
print(f"Median number of cysteine residues: {median_cysteines}")
