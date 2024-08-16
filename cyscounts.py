from Bio import SeqIO
import gzip

# Define the path to the gzipped FASTA file
filename = '/content/uniprotkb_human_AND_model_organism_9606_2024_04_21.fasta.gz'  # Update this to your actual file path

# Initialize the total count of cysteine residues
total_cysteines = 0

# Open and read the gzipped FASTA file in text mode
with gzip.open(filename, "rt") as handle:
    # Process each protein record in the FASTA file
    for record in SeqIO.parse(handle, "fasta"):
        # Count cysteine residues in the current protein sequence and add to the total count
        total_cysteines += str(record.seq).count('C')

# Output the total count of cysteine residues
print(f"Total number of cysteine residues in the human proteome: {total_cysteines}")
