import pandas as pd
from Bio import SeqIO
import gzip

# Define the file path
fasta_file = "/content/uniprotkb_human_AND_model_organism_9606_2024_04_21.fasta.gz"

# Initialize lists to store UniProt IDs and their corresponding cysteine counts
uniprot_ids = []
cysteine_counts = []

# Open the gzipped FASTA file
with gzip.open(fasta_file, "rt") as handle:
    # Parse the FASTA file
    for record in SeqIO.parse(handle, "fasta"):
        # Extract the UniProt ID (assuming it's the first part of the description before the first space)
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        sequence = str(record.seq)
        
        # Count the number of cysteine residues (denoted by 'C')
        cysteine_count = sequence.count('C')
        
        # Store the results
        uniprot_ids.append(uniprot_id)
        cysteine_counts.append(cysteine_count)

# Create a DataFrame with the results
cysteine_df = pd.DataFrame({
    'Uniprot_ID': uniprot_ids,
    'Cysteine_count': cysteine_counts
})

# Save the DataFrame to a new Excel file
output_file_path = "/content/cysteine_counts_fasta.xlsx"
cysteine_df.to_excel(output_file_path, index=False)

# Notify the user of completion
print(f"The file has been saved as {output_file_path}")
