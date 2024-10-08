# Cysteine-i-cloud

## Overview

This repository provides code for analyzing cysteine redox proteoforms using Python. The analyses include calculating theoretical "i space", generating heatmaps, and defining various aspects of the redox proteoforms based on UniProt reference proteomes.

## 1. Coding Environment and Repository

All computational analyses were performed using Python in a Google Colaboratory (Colab) environment, which provides an interactive workspace for executing Python code. The primary libraries used include BioPython, Pandas, Numpy, and Scipy. The complete source code is available in our GitHub repository: [Cysteine-i-cloud](https://github.com/JamesCobley/Cysteine-i-cloud). The code is freely available under the MIT license, allowing for broad use and modification.

## 2. Calculating the Theoretical i Space of Reference Proteomes
To calculate the theoretical i space, reference proteomes were downloaded from UniProt as gzipped FASTA files. These files were uploaded onto a Colab workspace where the script scriptcys1.py was used to process the data. The steps involved:

Parsing the gzip file and extracting protein sequences.
Counting cysteine residues in each protein sequence.
Listing the number of proteins with different counts of cysteine residues.
Exporting results to an Excel file with columns for cysteine residues and protein counts.
The Excel file was processed using a “PERMUTATONIA” function to solve for the theoretical i space. For example, for a given cysteine residue count, the i space was calculated and summed to derive the theoretical value.

## Reference Proteomes
Species	Reference Proteome	Strain
E. coli	UniProt E. coli	Strain K12 (83333)
S. cerevisiae	UniProt S. cerevisiae	Strain ATCC 204508 / S288c
C. elegans	UniProt C. elegans	Bristol N2
Drosophila	UniProt Drosophila	Berkeley
D. rerio	UniProt D. rerio	n/a
X. tropicalis	UniProt X. tropicalis	n/a
M. musculus	UniProt M. musculus	C57BL/6J
H. sapiens	UniProt H. sapiens	n/a

## 3. Defining the Composition of the Human i Space
The output from scriptcys1.py was used to create visualizations displayed in Figure 2. To identify the protein with the most cysteine residues, the script Cys_R_integer_find.py was used. Supplementary data was analyzed to create a table by sorting cysteine counts.

## 4. Defining the Distribution of the Human i Space Using GO Terms
The script Cys-count-ID.py was implemented to categorize human proteins by cysteine residue counts. The resulting data file, Supplementary Data File 3, was used to sort the human i space by selected GO terms with further calculations performed manually in Excel.

## 5. Mathematically Limiting the Biological Accessible i Space
The script Cys_Expression.py was used to limit the number of downloadable cysteine redox proteoforms based on protein expression. Theoretical calculations were performed using empirical data and equation implementations to derive the biologically accessible i space value.

## 6. Factors Controlling the Number of Cysteine Redox Proteoforms
Structural Argument
An Alphafold structure of human PTEN was analyzed using Cys_SASA.py to calculate the SASA value for sulfur and alpha carbon atoms.

## Kinetic Argument
The molarity of PTEN in a HeLa cell was calculated and used to model the oxidation of PTEN. The custom script Cys_k_est.py automated the calculations, estimating the number of oxidized PTEN molecules based on reaction kinetics.

## How to use

Clone the repository: git clone https://github.com/JamesCobley/Cysteine-i-cloud.git
Install dependencies pip install -r requirements.txt

## Contact

For any questions and inquiries, please contact James Cobley (j_cobley@yahoo.com) 

The required libraries are listed in `requirements.txt`. You can install them using:
```bash
pip install -r requirements.txt

