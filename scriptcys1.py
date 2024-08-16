# scriptcys1.py
"""
Script to count the number of cysteine residues in protein sequences from a gzipped FASTA file.

Dependencies:
- Biopython
- Matplotlib (optional for plotting)

To install the dependencies, run:
pip install biopython
pip install matplotlib
"""

from Bio import SeqIO
import gzip
import sys

def count_cysteines(filename):
    """
    Count cysteine residues in protein sequences from a gzipped FASTA file.

    Parameters:
    filename (str): Path to the gzipped FASTA file.

    Returns:
    dict: A dictionary with the number of cysteines as keys and counts of proteins as values.
    """
    cysteine_counts = {}

    try:
        with gzip.open(filename, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                num_cysteines = str(record.seq).count('C')
                if num_cysteines in cysteine_counts:
                    cysteine_counts[num_cysteines] += 1
                else:
                    cysteine_counts[num_cysteines] = 1
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return None

    return cysteine_counts

if __name__ == "__main__":
    # Modify this path to point to your actual file location
    filename = '/content/uniprotkb_mouse_AND_model_organism_1009_2024_04_24.fasta.gz'
    
    counts = count_cysteines(filename)

    if counts:
        # Output the results sorted by the number of cysteines
        for count in sorted(counts):
            print(f"Proteins with {count} cysteine(s): {counts[count]}")
