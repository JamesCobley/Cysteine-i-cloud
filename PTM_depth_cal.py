# Import necessary libraries
import gzip
from Bio import SeqIO
import pandas as pd
import math
from collections import Counter

# Define the input and output paths
input_fasta = "/content/uniprotkb_human_AND_model_organism_9606_2024_05_08.fasta.gz"
output_excel = "/content/PTM_space_results.xlsx"

# Define PTM counts (n-states) for each amino acid
ptm_counts = {
    "C": 50,  # Cysteine
    "K": 39,  # Lysine
    "M": 3,   # Methionine
    "Y": 8,   # Tyrosine
    "W": 4,   # Tryptophan
    "H": 3,   # Histidine
    "S": 6,   # Serine
    "T": 5,   # Threonine
    "Q": 4    # Glutamine
}

# Default PTM state for all other amino acids (2 states: modified + unmodified)
default_ptm_states = 2

# Function to calculate PTM-space for a single protein
def calculate_ptm_space(sequence):
    # Count amino acids in the sequence
    aa_counts = Counter(sequence)
    log10_ptm_space = 0  # Start in log-space to avoid overflow

    # Iterate over each amino acid and calculate its contribution
    for aa, count in aa_counts.items():
        if aa in ptm_counts:
            # Use the explicitly defined PTM count
            ptm_states = ptm_counts[aa]
        else:
            # Use the default PTM states for other amino acids
            ptm_states = default_ptm_states

        # Add the log10 contribution for this amino acid
        log10_ptm_space += count * math.log10(ptm_states)

    return log10_ptm_space

# Parse the gzipped FASTA file and calculate PTM-space
results = []  # To store UniProt entry, amino acid count, and PTM-space
total_log10_ptm_space = 0  # To store cumulative PTM-space in log-space

with gzip.open(input_fasta, "rt") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        uniprot_entry = record.id  # Get UniProt entry (record ID)
        sequence = str(record.seq)  # Get the protein sequence
        log10_ptm_space = calculate_ptm_space(sequence)  # Calculate PTM-space in log10

        # Add the protein's PTM-space to the total proteome-level space
        total_log10_ptm_space += log10_ptm_space

        # Convert PTM-space to scientific notation for readability
        ptm_space_sci = f"10^{log10_ptm_space:.3f}"

        # Append results
        results.append([uniprot_entry, len(sequence), ptm_space_sci])

# Save results to an Excel file
df = pd.DataFrame(results, columns=["UniProt entry", "Amino acid count", "PTM-space"])
df.to_excel(output_excel, index=False)

# Print the total proteome-level PTM-space in scientific notation
print(f"Total Proteome-Level PTM-Space: 10^{total_log10_ptm_space:.3f}")

print(f"Results saved to {output_excel}")
