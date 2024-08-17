# Import necessary module from Biopython
from Bio import SeqIO

# Define the path to the FASTA file
fasta_file_path = "/content/UP000000437_7955.fasta"

# Function to count cysteine residues in a protein sequence
def count_cysteines(sequence):
    return sequence.count('C')

# Initialize variables to track the protein with the most cysteine residues
max_cysteines = 0
max_protein_id = ""

# Open and parse the FASTA file to find the protein with the most cysteine residues
with open(fasta_file_path, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq = str(record.seq)
        num_cysteines = count_cysteines(seq)

        if num_cysteines > max_cysteines:
            max_cysteines = num_cysteines
            max_protein_id = record.id

# Output the result
print(f"The protein with the most cysteine residues is {max_protein_id} with {max_cysteines} cysteines.")
