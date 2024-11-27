# Import required libraries
import pandas as pd
import numpy as np

# Helper functions
def generate_k_states(max_r=100, q=0.5):
    """Generates a directory of k-states for all r in the range 1 to max_r."""
    k_state_directory = {}
    for r in range(1, max_r + 1):
        total_states = 2 ** r
        k_space = r + 1  # k-space = r + 1
        k_states = [binomial_count(r, k) for k in range(k_space)]
        redox_grades = [k / r * 100 for k in range(k_space)]
        probabilities = [binomial_prob(r, k, q) for k in range(k_space)]
        k_state_directory[r] = {
            "k_states": k_states,
            "redox_grades": redox_grades,
            "probabilities": probabilities,
        }
    return k_state_directory

def binomial_count(r, k):
    """Calculates the binomial coefficient for r choose k."""
    return int(np.math.factorial(r) / (np.math.factorial(k) * np.math.factorial(r - k)))

def binomial_prob(r, k, q):
    """Calculates the binomial probability for a given k-state."""
    return binomial_count(r, k) * (q ** k) * ((1 - q) ** (r - k))

# Upload file
from google.colab import files
uploaded = files.upload()
file_path = list(uploaded.keys())[0]

# Load the dataset
df = pd.read_excel(file_path)

# Generate the k-state directory
k_state_directory = generate_k_states(max_r=100)

# Process the dataset
output_data = []

for _, row in df.iterrows():
    r = int(row['Cysteines'])
    copies = float(row['Copies per cell'])
    if r > 100 or r == 0:  # Exclude proteins outside the 1-100 cysteine range
        continue

    # Retrieve k-state data for the current r
    k_data = k_state_directory[r]
    probabilities = k_data["probabilities"]
    redox_grades = k_data["redox_grades"]
    max_proteoforms_per_k = k_data["k_states"]

    # Distribute copies across k-states
    copies_distribution = [copies * p for p in probabilities]

    # Calculate MIN proteoforms: one unique proteoform per k_state accessed
    min_proteoforms_accessed = sum(1 for c in copies_distribution if c > 0)

    # Calculate MAX proteoforms
    max_proteoforms_accessed = sum(
        min(max_p, c) for max_p, c in zip(max_proteoforms_per_k, copies_distribution)
    )

    # Ensure max does not exceed I-space
    max_proteoforms_accessed = min(max_proteoforms_accessed, 2**r)

    # Weighted mean redox state
    weighted_redox_state = sum(c * r for c, r in zip(copies_distribution, redox_grades)) / sum(copies_distribution)

    # Append data to output
    output_data.append({
        "Protein name": row["Protein name"],
        "Uniprot": row["Uniprot"],
        "Cysteines": r,
        "I space": row["I space"],
        "Equation 5": row["Equation 5"],
        "Copies per cell": copies,
        "Accessed k_range": f"{redox_grades[0]}-{redox_grades[-1]}",
        "Redox state of the molecule": weighted_redox_state,
        "MIN proteoforms": min_proteoforms_accessed,
        "MAX proteoforms": max_proteoforms_accessed,
    })

# Calculate the overall redox state of the proteome
total_copies = sum(row["Copies per cell"] for row in output_data)
total_redox_state = sum(row["Copies per cell"] * row["Redox state of the molecule"] for row in output_data) / total_copies

# Save the results to a DataFrame and output file
output_df = pd.DataFrame(output_data)
output_df["Total proteome redox state"] = total_redox_state

output_file_path = '/content/proteome_redox_analysis_results.xlsx'
output_df.to_excel(output_file_path, index=False)

# Provide the download link
print(f"Results saved to {output_file_path}")
files.download(output_file_path)
