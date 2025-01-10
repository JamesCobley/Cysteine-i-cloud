import matplotlib.pyplot as plt
import numpy as np

# Data for genome size (Mb) and theoretical i-space
species = ["E. coli", "S. cerevisiae", "C. elegans", "Drosophila", "D. rerio", "X. tropicalis", "M. musculus", "H. sapiens"]
genome_size = [4.6, 12.1, 100, 180, 1500, 1400, 2700, 3200]  # Genome size in Mb
i_space = [
    2.5e9,
    4.95e27,
    1.9e109,
    6.7e190,
    2.27e190,
    6.79e184,
    1.84e165,
    3.02e169
]  # Theoretical i-space (log scale for better visualization)

# Convert i-space to log scale for better visualization
log_i_space = np.log10(i_space)

# Create the scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(genome_size, log_i_space, color="blue", edgecolor="black", s=100, alpha=0.8)

# Annotate each point with the species name, introducing slight offsets for overlapping points
offsets = [
    (5, -5),  # E. coli
    (5, 5),   # S. cerevisiae
    (-5, 5),  # C. elegans
    (-5, -5), # Drosophila
    (5, -10), # D. rerio
    (-5, 10), # X. tropicalis
    (5, 10),  # M. musculus
    (-10, -5) # H. sapiens
]
for i, txt in enumerate(species):
    plt.annotate(txt, (genome_size[i], log_i_space[i]), fontsize=10, 
                 xytext=offsets[i], textcoords="offset points")

# Customize plot
plt.title("Genome Size vs Theoretical i-Space", fontsize=14)
plt.xlabel("Genome Size (Mb)", fontsize=12)
plt.ylabel("Log10(Theoretical i-Space)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.6)

# Save as a high-resolution image
plt.savefig("genome_vs_ispace_fixed.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()
