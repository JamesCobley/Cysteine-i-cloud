# scriptcys1.py
"""
Script to count the number of cysteine residues in protein sequences from a gzipped FASTA file.

Dependencies:
- Biopython
- Matplotlib (optional for plotting)
-Pandas

To install the dependencies, run:
pip install biopython
pip install matplotlib
pip install pandas
"""

import gzip
import pandas as pd

# Define the file path to the gzipped FASTA file
file_path = '/content/uniprotkb_human_AND_model_organism_9606_2024_05_08.fasta.gz'

# Dictionary to count the number of proteins by the number of cysteine residues
cysteine_counts = {}

# Open and read the gzipped FASTA file in text mode
with gzip.open(file_path, 'rt') as file:
    protein_sequence = ''
    for line in file:
        if line.startswith('>'):
            if protein_sequence:
                # Count cysteine residues in the current protein sequence
                count = protein_sequence.count('C')
                if count in cysteine_counts:
                    cysteine_counts[count] += 1
                else:
                    cysteine_counts[count] = 1
            protein_sequence = ''  # Reset the sequence for the next protein
        else:
            protein_sequence += line.strip()  # Accumulate the protein sequence

    # Handle the last protein sequence in the file
    if protein_sequence:
        count = protein_sequence.count('C')
        if count in cysteine_counts:
            cysteine_counts[count] += 1
        else:
            cysteine_counts[count] = 1

# Convert the cysteine counts dictionary to a DataFrame
df = pd.DataFrame(list(cysteine_counts.items()), columns=['Cysteine Residues', 'Number of Proteins'])

# Sort the DataFrame by the number of cysteine residues
df = df.sort_values(by='Cysteine Residues')

# Define the path to save the Excel file
excel_path = '/content/cysteine_counts.xlsx'

# Save the DataFrame to an Excel file
df.to_excel(excel_path, index=False)

# Print confirmation message
print(f"Excel file has been saved to {excel_path}")
