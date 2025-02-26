# Cysteine-i-cloud

## Overview

This repository provides code for analyzing cysteine proteoforms using [Python](https://www.python.org/) as described in the [The nonlinear cysteine redox dynamics in the i-space: A proteoform-centric theory of redox regulation paper](https://www.sciencedirect.com/science/article/pii/S2213231725000369?via%3Dihub) .The analyses include calculating the theoretical "i space", generating heatmaps, and defining various aspects of the redox proteoforms based on [UniProt](https://www.uniprot.org/) reference proteomes.

## 1. Coding Environment and Repository

All computational analyses were performed using Python in a Google Colaboratory ([Colab](https://colab.google/)) environment, which provides an interactive workspace for executing Python code. The primary libraries used include BioPython, Pandas, Numpy, and Scipy. The complete source code is available in our GitHub repository: [Cysteine-i-cloud](https://github.com/JamesCobley/Cysteine-i-cloud). The code is freely available under the MIT license, allowing for broad use and modification.

## 2. Calculating the Theoretical i Space of Reference Proteomes
To calculate the theoretical i space, reference proteomes were downloaded from UniProt as gzipped FASTA files. These files were uploaded onto a Colab workspace where the script scriptcys1.py was used to process the data. 

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

The [human FASTA file](https://github.com/JamesCobley/Cysteine-i-cloud/blob/main/uniprotkb_human_AND_model_organism_9606_2024_05_08.fasta.gz) that we used has been uploaded to this repository.

To calculate the PTM-space in n-2 dimensional and n-X dimensonal phase space the PTM_space_cal.py and PTM_depth_cal.py scripts were used.

## 3. Defining the Composition of the Human i Space
The output from scriptcys1.py was used to create visualizations displayed in Figure 2. To identify the protein with the most cysteine residues, the script Cys_R_integer_find.py was used. Supplementary data was analyzed to create a table by sorting cysteine counts.

## 4. Defining the Distribution of the Human i Space Using GO Terms
The script Cys-count-ID.py was implemented to categorize human proteins by cysteine residue counts. The resulting data file, Supplementary Data File 3, was used to sort the human i space by selected GO terms with further calculations performed manually in Excel. 

## 5. Mathematically Limiting the Biological Accessible i Space
The script Cys_Expression.py was used to limit the number of downloadable cysteine redox proteoforms based on protein expression. Theoretical calculations were performed using empirical data and equation implementations to derive the biologically accessible i space value.

## 6. Factors Controlling the Number of Cysteine Redox Proteoforms
Structural Argument
An Alphafold structure of human PTEN was analyzed using Cys_SASA.py to calculate the SASA value for sulfur and alpha carbon atoms.

## 7. Kinetic Argument
The molarity of PTEN in a HeLa cell was calculated and used to model the oxidation of PTEN. The custom script Cys_k_est.py automated the calculations, estimating the number of oxidized PTEN molecules based on reaction kinetics.

## 8. Sampling the i-space

To calculate the minimum and maximum number of unique cysteine redox proteoforms the sampling_i_space_random.py and sampling_i_space.py scripts were implemented.

## 9. Monte Carlo Simulations

After elaborating PTP1B protoeforms in a matrix using the Cys_matrix.py script, Monte Carlo simulations were iterated on a virutal instance, using the [google cloud computing platform](https://cloud.google.com/?hl=en), by implementing the CysRedox_MonteCarlo_1.py script. The analysis of the resultant outputs was automated using the RedoxMonteCarlo_stats.py script.

## Note

Scripts under section 4, 6 and 7 now only apply to the [preprint](https://www.biorxiv.org/content/10.1101/2024.09.18.613618v1.abstract). 
The scripts used to generate the the figures in the manuscript are also provided where appropriate, the latest version has been [posted](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5047335). The final [paper](https://www.sciencedirect.com/science/article/pii/S2213231725000369?via%3Dihub) is available open access online.

## How to use

Clone the repository: git clone https://github.com/JamesCobley/Cysteine-i-cloud.git
Install dependencies pip install -r requirements.txt

The scripts are optimised for goolge [Colab](https://colab.google/) implementation but they can be adapted for other environments. 

## Contact

For any questions and inquiries, please contact James Cobley (j_cobley@yahoo.com) 

The required libraries are listed in `requirements.txt`. You can install them using:
```bash
pip install -r requirements.txt

