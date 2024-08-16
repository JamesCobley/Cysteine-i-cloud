import freesasa
from Bio import PDB

def calculate_cysteine_sasa(pdb_filename):
    """Calculate SASA for CA and SG atoms of cysteine residues in a PDB file."""
    
    # Load PDB file using Biopython
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_filename)

    # Initialize FreeSASA
    freesasa_structure = freesasa.Structure(pdb_filename)
    result = freesasa.calc(freesasa_structure)

    # Dictionary to store the SASA for CA and SG atoms of cysteine residues
    cysteine_sasa = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    res_id = residue.get_id()[1]
                    chain_id = chain.get_id()
                    sasa_ca = 0.0
                    sasa_sg = 0.0
                    for atom in residue:
                        atom_name = atom.get_name()
                        atom_serial = atom.get_serial_number()
                        if atom_name in ['CA', 'SG']:
                            sasa = result.atomArea(atom_serial - 1)
                            if atom_name == 'CA':
                                sasa_ca = sasa
                            elif atom_name == 'SG':
                                sasa_sg = sasa
                    cysteine_sasa[f'{chain_id}:{res_id}'] = {'CA': sasa_ca, 'SG': sasa_sg}

    return cysteine_sasa

def print_cysteine_sasa(cysteine_sasa):
    """Print the SASA for the CA and SG atoms of each cysteine residue."""
    for res, sasa in cysteine_sasa.items():
        print(f'Residue {res} - CA SASA: {sasa["CA"]}, SG SASA: {sasa["SG"]}')

def main():
    # Define the file path
    pdb_filename = '/content/AF-Q6DFD4-F1-model_v4 (1).pdb'

    # Calculate SASA for cysteine residues
    cysteine_sasa = calculate_cysteine_sasa(pdb_filename)

    # Print the SASA results
    print_cysteine_sasa(cysteine_sasa)

if __name__ == "__main__":
    main()
