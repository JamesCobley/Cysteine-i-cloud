import gzip
import pandas as pd
from Bio import SeqIO
from scipy.special import comb

# Parse the FASTA file (handles .gz)
def parse_fasta(fasta_path):
    fasta_data = {}
    with gzip.open(fasta_path, 'rt') as f:
        for record in SeqIO.parse(f, "fasta"):
            uniprot_id = record.id.split("|")[1]
            sequence = str(record.seq)
            fasta_data[uniprot_id] = sequence
    return fasta_data

# Process OxiMouse data
def process_oximouse_data(fasta_data, oximouse_path):
    oximouse_df = pd.read_excel(oximouse_path)
    results = []

    for _, row in oximouse_df.iterrows():
        uniprot_id = row.get("Uniprot ID")
        if not isinstance(uniprot_id, str) or uniprot_id not in fasta_data:
            continue

        sequence = fasta_data[uniprot_id]
        cysteine_positions = [i for i, aa in enumerate(sequence, start=1) if aa == "C"]
        R = len(cysteine_positions)

        if R == 0:
            continue

        detected_sites = [col for col in row.index[5:] if not pd.isna(row[col])]
        n_detected = len(detected_sites)
        cysteine_coverage = (n_detected / R) * 100

        detected_redox = []
        for col in detected_sites:
            try:
                value = str(row[col]).split("±")[0]
                detected_redox.append(float(value))
            except (ValueError, TypeError):
                continue

        avg_detected_redox = sum(detected_redox) / len(detected_redox) if detected_redox else 0
        total_redox_state = (n_detected * avg_detected_redox) / R

        min_k_states, min_i_states = calculate_min_states(R, total_redox_state)

        results.append({
            "Uniprot": uniprot_id,
            "R": R,
            "Cysteine coverage (%)": cysteine_coverage,
            "Redox state of Oximouse detected cysteines (%)": avg_detected_redox,
            "Redox state of the protein": total_redox_state,
            "Min K-state": min_k_states,
            "Min i-states": min_i_states,
        })

    return results

# Calculate minimum k and i states
def calculate_min_states(R, redox_state):
    min_k_states = 0
    min_i_states = 0

    for k in range(R + 1):
        percent_oxidized = (k / R) * 100
        if percent_oxidized >= redox_state:
            min_k_states = k + 1
            min_i_states = sum(int(comb(R, i)) for i in range(min_k_states))
            break

    return min_k_states, min_i_states

# Save results to Excel
def save_results_to_excel(results, output_path):
    df = pd.DataFrame(results)
    df.to_excel(output_path, index=False)
    print(f"Results saved to {output_path}")

# File paths
fasta_path = "/content/UP000000589_10090.fasta (1).gz"
oximouse_path = "/content/All_young_sites.xlsx"
output_path = "/content/oximouse_analysis_results.xlsx"

# Run the analysis
fasta_data = parse_fasta(fasta_path)
results = process_oximouse_data(fasta_data, oximouse_path)
save_results_to_excel(results, output_path)
