import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for the grid
grid_size = 10  # 10x10 grid
protein_density = 5  # Number of protein molecules per box
reductant_density = 20  # Average number of reductants per box
oxidant_density = 20  # Average number of oxidants per box

# Step 1: Create 2D Grid with Protein, Reductants, and Oxidants
def generate_nanodomain_grid(grid_size, protein_density, reductant_density, oxidant_density):
    grid = []
    for _ in range(grid_size):
        row = []
        for _ in range(grid_size):
            proteins = protein_density
            reductants = np.random.poisson(reductant_density)
            oxidants = np.random.poisson(oxidant_density)
            row.append((proteins, reductants, oxidants))
        grid.append(row)
    return grid

# Generate the grid
nanodomain_grid = generate_nanodomain_grid(grid_size, protein_density, reductant_density, oxidant_density)

# Step 2: Compute Oxidation Levels for 3D Plot
def compute_oxidation_levels(grid):
    oxidation_levels = np.zeros((len(grid), len(grid[0])))
    for i, row in enumerate(grid):
        for j, (proteins, reductants, oxidants) in enumerate(row):
            if proteins > 0:
                # Compute the oxidation level as oxidants/(reductants + oxidants)
                oxidation_levels[i, j] = oxidants / (reductants + oxidants + 1e-6)
            else:
                oxidation_levels[i, j] = 0  # No oxidation if no proteins
    return oxidation_levels

# Calculate oxidation levels
oxidation_levels = compute_oxidation_levels(nanodomain_grid)

# Step 3: Combined Plot Function
def plot_combined_nanodomain_and_redox(grid, oxidation_levels, save_path="combined_plot.png"):
    fig = plt.figure(figsize=(16, 8), dpi=300)

    # Subplot 1: 2D Nanodomain Spray Model
    ax1 = fig.add_subplot(1, 2, 1)
    for i, row in enumerate(grid):
        for j, (proteins, reductants, oxidants) in enumerate(row):
            # Plot proteins as black dots
            ax1.scatter(
                np.random.uniform(j, j+1, proteins),
                np.random.uniform(i, i+1, proteins),
                color="black", s=10, label="Proteins" if i == 0 and j == 0 else ""
            )
            # Plot reductants as blue dots
            ax1.scatter(
                np.random.uniform(j, j+1, reductants),
                np.random.uniform(i, i+1, reductants),
                color="blue", s=5, alpha=0.6, label="Reductants" if i == 0 and j == 0 else ""
            )
            # Plot oxidants as red dots
            ax1.scatter(
                np.random.uniform(j, j+1, oxidants),
                np.random.uniform(i, i+1, oxidants),
                color="red", s=5, alpha=0.6, label="Oxidants" if i == 0 and j == 0 else ""
            )
    ax1.set_xlim(0, grid_size)
    ax1.set_ylim(0, grid_size)
    ax1.set_title("2D Nanodomain Spray Model", fontsize=14)
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.legend(loc="upper right")

    # Subplot 2: 3D Redox Profile
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    x = np.arange(oxidation_levels.shape[0])
    y = np.arange(oxidation_levels.shape[1])
    x, y = np.meshgrid(x, y)
    z = oxidation_levels.T  # Transpose to align with grid
    surf = ax2.plot_surface(x, y, z, cmap="viridis", edgecolor='k', alpha=0.8)
    ax2.set_title("3D Redox Profile", fontsize=14)
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.set_zlabel("Oxidation Level")
    fig.colorbar(surf, ax=ax2, shrink=0.5, aspect=10)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)  # Save at 300 DPI
    plt.show()

# Plot and save the combined figure
plot_combined_nanodomain_and_redox(nanodomain_grid, oxidation_levels, save_path="combined_plot.png")
