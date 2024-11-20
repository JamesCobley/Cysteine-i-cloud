# Import necessary libraries
import gzip
from Bio import SeqIO
import pandas as pd
import math

# Define the input and output paths
input_fasta = "/content/uniprotkb_human_AND_model_organism_9606_2024_05_08.fasta.gz"
output_excel = "/content/PTM_space_results.xlsx"

# Function to calculate PTM-space using logarithms
def calculate_ptm_space_log(num_amino_acids):
    # Calculate log10(PTM-space)
    return num_amino_acids * math.log10(2)

# Parse the gzipped FASTA file
results = []  # To store UniProt entry, amino acid count, and PTM-space
log10_sum_ptm_space = None  # To store cumulative PTM-space in logarithmic form

with gzip.open(input_fasta, "rt") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        uniprot_entry = record.id  # Get UniProt entry (record ID)
        sequence = str(record.seq)  # Get the protein sequence
        amino_acid_count = len(sequence)  # Count the number of amino acids in the sequence
        log10_ptm_space = calculate_ptm_space_log(amino_acid_count)  # PTM-space in log10
        ptm_space_sci = f"10^{log10_ptm_space:.3f}"  # Represent PTM-space in scientific notation
        
        # Update log-sum of PTM-space
        if log10_sum_ptm_space is None:
            log10_sum_ptm_space = log10_ptm_space  # Initialize on first entry
        else:
            # Combine in log-space using log10(a + b) = log10(10^a + 10^b)
            if log10_ptm_space > log10_sum_ptm_space:
                log10_sum_ptm_space, log10_ptm_space = log10_ptm_space, log10_sum_ptm_space  # Ensure a >= b
            log10_sum_ptm_space += math.log10(1 + 10 ** (log10_ptm_space - log10_sum_ptm_space))
        
        results.append([uniprot_entry, amino_acid_count, ptm_space_sci])  # Append results

# Save results to an Excel file
df = pd.DataFrame(results, columns=["UniProt entry", "Amino acid integer", "PTM-space"])
df.to_excel(output_excel, index=False)

# Output the total PTM-space in scientific notation
if log10_sum_ptm_space is not None:
    total_ptm_space_sci = f"10^{log10_sum_ptm_space:.3f}"
    print(f"Total PTM-space across all entries: {total_ptm_space_sci}")

print(f"PTM-space calculations completed. Results saved to {output_excel}")
