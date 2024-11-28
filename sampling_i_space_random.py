import pandas as pd
import numpy as np
from google.colab import files

# Upload file
uploaded = files.upload()
file_path = list(uploaded.keys())[0]

# Load the dataset
df = pd.read_excel(file_path)

# Helper functions
def generate_k_states(max_r=100):
    """Generates a directory of k-states for all r in the range 1 to max_r."""
    k_state_directory = {}
    for r in range(1, max_r + 1):
        total_states = 2 ** r
        k_space = r + 1  # k-space = r + 1
        k_states = [binomial_count(r, k) for k in range(k_space)]
        redox_grades = [k / r * 100 for k in range(k_space)]
        k_state_directory[r] = {
            "k_states": k_states,
            "redox_grades": redox_grades,
        }
    return k_state_directory

def binomial_count(r, k):
    """Calculates the binomial coefficient for r choose k."""
    return int(np.math.factorial(r) / (np.math.factorial(k) * np.math.factorial(r - k)))

def assign_distribution(r, copies, profile):
    """Assigns molecules across k-states based on the selected profile."""
    k_space = r + 1
    if profile == "Random":
        probabilities = np.ones(k_space) / k_space
    elif profile == "Super Poisson (Reduced)":
        probabilities = np.exp(-np.arange(k_space))
        probabilities /= probabilities.sum()
    elif profile == "Super Poisson (Oxidized)":
        probabilities = np.exp(-np.arange(k_space)[::-1])
        probabilities /= probabilities.sum()
    elif profile == "Gaussian":
        center = k_space // 2
        probabilities = np.exp(-0.5 * ((np.arange(k_space) - center) ** 2) / (center / 2) ** 2)
        probabilities /= probabilities.sum()
    elif profile == "Polarized":
        probabilities = np.zeros(k_space)
        probabilities[0], probabilities[-1] = 0.5, 0.5
    else:
        raise ValueError("Unknown distribution profile")
    return probabilities * copies

# Generate k-state directory
max_r = 100
k_state_directory = generate_k_states(max_r=max_r)

# Define possible distribution profiles
distribution_profiles = [
    "Random",
    "Super Poisson (Reduced)",
    "Super Poisson (Oxidized)",
    "Gaussian",
    "Polarized"
]

# Sort and rank proteins by "Copies per cell" in descending order
df = df.sort_values(by="Copies per cell", ascending=False).reset_index(drop=True)

# Process the dataset
output_data = []

for index, row in df.iterrows():
    r = int(row['Cysteines'])
    copies = float(row['Copies per cell'])
    if r > max_r or r == 0:
        continue

    # Retrieve k-state data for the current r
    k_data = k_state_directory[r]
    redox_grades = k_data["redox_grades"]

    # Assign distribution profile
    if index < 100:  # Top 100 most abundant proteins
        profile = "Super Poisson (Reduced)"
    else:  # Randomly assign remaining proteins to one of the profiles
        profile = np.random.choice(distribution_profiles)

    # Distribute molecules across k-states based on the profile
    copies_distribution = assign_distribution(r, copies, profile)

    # Calculate MIN proteoforms: Count accessed k-states, capped at k + 1 and total copies
    min_proteoforms_accessed = min(sum(1 for copies_in_k in copies_distribution if copies_in_k > 0), r + 1, int(copies))

    # Calculate MAX proteoforms, capped by the total number of molecules
    max_proteoforms_accessed = sum(
        min(copies_in_k, binomial_count(r, k)) for k, copies_in_k in enumerate(copies_distribution)
    )
    max_proteoforms_accessed = min(max_proteoforms_accessed, int(copies))  # Cap at total molecules

    # Ensure MAX Proteoforms is at least as large as MIN Proteoforms
    if max_proteoforms_accessed < min_proteoforms_accessed:
        max_proteoforms_accessed = min_proteoforms_accessed

    # Weighted mean redox state
    weighted_redox_state = sum(c * r for c, r in zip(copies_distribution, redox_grades)) / sum(copies_distribution)

    # Append data to output
    output_data.append({
        "Protein name": row["Protein name"],
        "Uniprot": row["Uniprot"],
        "Cysteines": r,
        "I space": 2**r,
        "Equation 5": row["Equation 5"],
        "Copies per cell": copies,
        "Assigned profile": profile,
        "Accessed k_range": f"{redox_grades[0]}-{redox_grades[-1]}",
        "Redox state of the molecule": weighted_redox_state,
        "MIN proteoforms": min_proteoforms_accessed,
        "MAX proteoforms": max_proteoforms_accessed,
    })

# Convert to DataFrame
output_df = pd.DataFrame(output_data)

# Calculate total copies across the proteome
total_copies = output_df["Copies per cell"].sum()

# Calculate total proteome redox state
total_redox_state = (output_df["Redox state of the molecule"] * output_df["Copies per cell"]).sum() / total_copies

# Add total proteome redox state as a separate output
output_df["Total proteome redox state"] = total_redox_state

# Save results to an Excel file
output_file_path = '/content/proteome_redox_analysis_fixed.xlsx'
output_df.to_excel(output_file_path, index=False)

# Provide download link
print(f"Results saved to {output_file_path}")
files.download(output_file_path)
