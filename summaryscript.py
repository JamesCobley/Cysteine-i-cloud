import gzip

# Define the file path
file_path = '/content/UP000008143_8364.fasta (1).gz'

# Initialize counters for cysteine residues and total amino acids
total_cysteine_count = 0
total_amino_acid_count = 0

# Open the gzipped file in text mode
with gzip.open(file_path, 'rt') as file:
    for line in file:
        if line.startswith('>'):
            continue  # Skip header lines starting with '>'
        
        # Strip whitespace and newline characters from the line
        line = line.strip()
        
        # Update the total count of amino acids and cysteine residues
        total_amino_acid_count += len(line)
        total_cysteine_count += line.count('C')

# Output the results
print(f"Total number of cysteine residues: {total_cysteine_count}")
print(f"Total number of amino acids: {total_amino_acid_count}")
