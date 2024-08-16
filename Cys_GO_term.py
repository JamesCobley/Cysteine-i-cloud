# Install the required packages
!pip install aiohttp pandas xlsxwriter nest_asyncio

import asyncio
import aiohttp
import pandas as pd
from xml.etree import ElementTree as ET
import nest_asyncio
from google.colab import files

# Allow nested event loops in Google Colab
nest_asyncio.apply()

# Load the dataset from an Excel file
file_path = '/content/cysteine_counts_by_id_human.xlsx'
data = pd.read_excel(file_path)

def create_protein_category_list(df):
    """
    Flattens the dataframe and creates a list of tuples containing (UniProt accession, Cysteine category).
    """
    protein_category_list = []
    for category in df.columns:
        for protein in df[category].dropna():
            protein_category_list.append((protein, category))
    return protein_category_list

# Create a list of proteins and their corresponding cysteine categories
protein_category_list = create_protein_category_list(data)
protein_df = pd.DataFrame(protein_category_list, columns=['UniProt', 'CysteineCategory'])

async def fetch_go_annotations(session, uniprot_id):
    """
    Asynchronously fetches GO annotations for a given UniProt accession ID.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    async with session.get(url) as response:
        if response.status == 200:
            text = await response.text()
            root = ET.fromstring(text)
            go_terms = [go_term.attrib['id'] for go_term in root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='GO']")]
            return uniprot_id, go_terms
        return uniprot_id, []

async def map_to_go_annotations_async(protein_list):
    """
    Maps UniProt accession numbers to their corresponding GO annotations asynchronously.
    """
    async with aiohttp.ClientSession() as session:
        tasks = [fetch_go_annotations(session, protein) for protein in protein_list]
        results = await asyncio.gather(*tasks)
        return dict(results)

# Map the UniProt accession numbers to GO annotations
protein_list = protein_df['UniProt'].tolist()
go_annotation_mapping = asyncio.run(map_to_go_annotations_async(protein_list))
protein_df['GO'] = protein_df['UniProt'].map(go_annotation_mapping)

# Convert lists in 'GO' column to strings for display purposes
protein_df['GO_str'] = protein_df['GO'].apply(lambda x: ', '.join(x))

# Define the list of GO terms of interest
go_terms_of_interest = [
    "GO:0003677", "GO:0003723", "GO:0003774", "GO:0003824", "GO:0003924", "GO:0005198", 
    "GO:0005215", "GO:0008092", "GO:0008289", "GO:0009975", "GO:0016209", "GO:0016491", 
    "GO:0016740", "GO:0016787", "GO:00016853", "GO:00016874", "GO:0042393", "GO:0044183", 
    "GO:0045182", "GO:0048018", "GO:0060089", "GO:0060090", "GO:0098772", "GO:0140096", 
    "GO:0140097", "GO:01400098", "GO:0140110", "GO:0140223", "GO:0140657", "GO:0098631"
]

# Define the output file path
output_path = '/content/cysteine_distribution_by_go_term_molecular.xlsx'
writer = pd.ExcelWriter(output_path, engine='xlsxwriter')

# Process each GO term of interest and save the results to separate sheets
for go_term in go_terms_of_interest:
    selected_proteins = protein_df[protein_df['GO'].apply(lambda x: go_term in x)]

    if selected_proteins.empty:
        print(f"No proteins found for GO term: {go_term}")
    else:
        # Group by cysteine category and aggregate UniProt accessions
        cysteine_distribution = selected_proteins.groupby('CysteineCategory')['UniProt'].apply(list).reset_index()
        cysteine_distribution.columns = ['CysteineResidues', 'UniProtAccessions']
        cysteine_distribution['NumberOfProteins'] = cysteine_distribution['UniProtAccessions'].apply(len)

        # Simplify the GO term for use as a sheet name
        simplified_go_term = go_term.replace("GO:", "")

        # Save the cysteine distribution to a sheet in the Excel file
        cysteine_distribution.to_excel(writer, sheet_name=simplified_go_term, index=False, columns=['CysteineResidues', 'NumberOfProteins', 'UniProtAccessions'])

# Close the Pandas Excel writer and output the Excel file
writer.close()

# Download the file
files.download(output_path)

print(f"Saved cysteine distribution to {output_path}")
