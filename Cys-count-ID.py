from Bio import SeqIO
import gzip
import pandas as pd

# Path to the gzipped FASTA file
filename = '/content/uniprotkb_human_AND_model_organism_9606_2024_04_21.fasta.gz'

# Dictionary to store cysteine residue classes and their corresponding accessions
cysteine_classes = {}

# Open and read the gzipped FASTA file
with gzip.open(filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Extract the accession ID from the record header
        accession = record.id.split('|')[1]  # Adjust depending on the ID format
        # Count the number of cysteine residues in the protein sequence
        cysteine_count = str(record.seq).count('C')

        # Initialize the list if this cysteine count is not yet in the dictionary
        if cysteine_count not in cysteine_classes:
            cysteine_classes[cysteine_count] = []
        
        # Append the accession ID to the corresponding cysteine count list
        cysteine_classes[cysteine_count].append(accession)

# Determine the maximum number of accessions in any cysteine class
max_accessions = max(len(accessions) for accessions in cysteine_classes.values())

# Create a DataFrame with cysteine counts as columns and accessions as rows
df = pd.DataFrame({
    count: accessions + [None] * (max_accessions - len(accessions))  # Pad shorter lists with None
    for count, accessions in sorted(cysteine_classes.items())
})

# Save the DataFrame to an Excel file
excel_path = '/content/output_file.xlsx'
df.to_excel(excel_path, index=False)

# Confirmation message
print(f"Accession list saved to {excel_path}")
