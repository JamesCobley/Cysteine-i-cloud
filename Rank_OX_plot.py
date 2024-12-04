import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = "/content/proteome_redox_analysis_fixed_3000.xlsx"  # Replace with your file path
data = pd.read_excel(file_path)

# Ensure the relevant columns exist
if "Copies per cell" not in data.columns or "Redox state of the molecule" not in data.columns:
    raise ValueError("The required columns are not found in the dataset.")

# Rank molecules by "Copies per cell"
data_sorted = data.sort_values(by="Copies per cell", ascending=False)
data_sorted['Rank'] = range(1, len(data_sorted) + 1)

# Extract values for plotting
rank = data_sorted['Rank']
redox_state = data_sorted["Redox state of the molecule"]

# Plot Redox state vs. Rank
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(rank, redox_state, marker="o", linestyle="", markersize=3, alpha=0.7, color="black")
plt.title("Redox State of the Molecule vs. Rank (by Copy Number)", fontsize=14)
plt.xlabel("Rank (1 = Highest Copy Number)", fontsize=12)
plt.ylabel("Oxidation(%)", fontsize=12)
plt.grid(alpha=0.5, linestyle="--")
plt.tight_layout()

# Save and show the plot
plt.savefig("redox_state_vs_rank.png", dpi=300)
plt.show()
