import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# Data for genome size (Mb) and theoretical i-space
species = ["E. coli", "S. cerevisiae", "C. elegans", "Drosophila", 
           "D. rerio", "X. tropicalis", "M. musculus", "H. sapiens"]
genome_size = [4.6, 12.1, 100, 180, 1500, 1400, 2700, 3200]  # Genome size in Mb
i_space = [
    2.5e9, 4.95e27, 1.9e109, 6.7e190, 2.27e190, 6.79e184, 1.84e165, 3.02e169
]  # Theoretical i-space (log scale for better visualization)

# Convert i-space to log scale for better visualization
log_i_space = np.log10(i_space)

# Define colors for each species
colors = ["red", "blue", "green", "purple", "orange", "brown", "pink", "cyan"]

# Create the scatter plot with colored bubbles
plt.figure(figsize=(10, 6))
scatter = plt.scatter(genome_size, log_i_space, c=colors, edgecolor="black", s=300, alpha=0.8)


# Create a legend with species names and corresponding colors
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=10, 
                          label=f"${species[i]}$") for i in range(len(species))]
plt.legend(handles=legend_elements, title="Species", fontsize=10, title_fontsize=12)

# Customize plot
plt.title("Genome Size vs Theoretical i-Space", fontsize=14)
plt.xlabel("Genome Size (Mb)", fontsize=12)
plt.ylabel("Log10(Theoretical i-Space)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.6)

# Save as a high-resolution image
plt.savefig("genome_vs_ispace_bubbles.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()
