import pandas as pd
from Bio import SeqIO
import gzip

def extract_cysteine_counts_from_fasta(fasta_file):
    """Extract cysteine counts from a gzipped FASTA file."""
    uniprot_ids = []
    cysteine_counts = []

    # Open and parse the gzipped FASTA file
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Extract and normalize the UniProt ID
            uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
            uniprot_id = uniprot_id.split("-")[0]  # Remove version numbers or isoform info

            # Count cysteine residues (denoted by 'C')
            sequence = str(record.seq)
            cysteine_count = sequence.count('C')

            # Store the results
            uniprot_ids.append(uniprot_id)
            cysteine_counts.append(cysteine_count)

    # Create a DataFrame with the results
    return pd.DataFrame({'Uniprot_ID': uniprot_ids, 'Cysteine_count': cysteine_counts})

def match_and_assign_cysteine_counts(protein_file, cysteine_df, output_file_path):
    """Match UniProt IDs from an Excel file to cysteine counts and save the results."""
    # Load the second sheet of the Excel file containing the UniProt IDs
    proteins_df = pd.read_excel(protein_file, sheet_name=1)

    # Initialize a new column for cysteine counts
    proteins_df['Cysteine_count'] = 0

    # Iterate over each row in the proteins_df
    for index, row in proteins_df.iterrows():
        uniprot_field = row['Uniprot']
        
        if isinstance(uniprot_field, str):
            # Split the Uniprot column by semicolon
            uniprot_list = uniprot_field.split(';')
            
            # Search for the first match in the cysteine_df
            for uniprot_id in uniprot_list:
                uniprot_id = uniprot_id.split("-")[0]  # Normalize the ID
                if uniprot_id in cysteine_df['Uniprot_ID'].values:
                    cysteine_count = cysteine_df[cysteine_df['Uniprot_ID'] == uniprot_id]['Cysteine_count'].values[0]
                    proteins_df.at[index, 'Cysteine_count'] = cysteine_count
                    break  # Stop after the first match is found

    # Save the updated dataframe to an Excel file
    proteins_df.to_excel(output_file_path, index=False)
    print(f"The file has been saved as {output_file_path}")

def main():
    # Define the file paths
    fasta_file = "/content/uniprotkb_human_AND_model_organism_9606_2024_04_21.fasta.gz"
    protein_file = "/content/msb201181-sup-0002 (1).xlsx"
    output_file_path = "/content/cysteine_counts_merged_search_fixed.xlsx"

    # Extract cysteine counts from the FASTA file
    cysteine_df = extract_cysteine_counts_from_fasta(fasta_file)

    # Match and assign cysteine counts based on UniProt IDs in the Excel file
    match_and_assign_cysteine_counts(protein_file, cysteine_df, output_file_path)

if __name__ == "__main__":
    main()
