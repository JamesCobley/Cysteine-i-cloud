import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 10  # Total k-grades (from 0 to R)
total_molecules = 50000  # Total number of molecules

# Function to create custom distributions
def create_custom_distribution(type_, total_molecules, R):
    k_values = np.arange(0, R+1)
    molecules = np.zeros_like(k_values, dtype=float)
    
    if type_ == "Super-Poisson Reduced":
        # Leftward shift heavily favoring k = 0
        probs = np.exp(-0.5 * k_values)  # Exponentially decreasing
        probs /= probs.sum()
        molecules = (probs * total_molecules).astype(int)
    
    elif type_ == "Super-Poisson Oxidized":
        # Rightward shift heavily favoring k = 10
        probs = np.exp(-0.5 * (R - k_values))  # Exponentially increasing
        probs /= probs.sum()
        molecules = (probs * total_molecules).astype(int)
    
    elif type_ == "Gaussian":
        # Central peak around the middle k-value
        center = R // 2
        stddev = R / 4
        probs = np.exp(-0.5 * ((k_values - center) / stddev)**2)
        probs /= probs.sum()
        molecules = (probs * total_molecules).astype(int)
    
  
    
    elif type_ == "Polarized":
        # Bimodal distribution
        molecules[0] = 25000  # k = 0
        molecules[-1] = 25000  # k = 10
    
    else:
        raise ValueError("Unknown distribution type")
    
    return k_values, molecules

# Generate and plot all distributions in a single combined image
def plot_combined_distributions(distributions, total_molecules, R):
    fig, axes = plt.subplots(3, 2, figsize=(12, 18), dpi=300)
    axes = axes.flatten()  # Flatten the 2D grid into a 1D array for iteration

    for i, dist in enumerate(distributions):
        k_values, molecules = create_custom_distribution(dist, total_molecules, R)
        redox_percentages = k_values * (100 / R)  # Convert k to % oxidized
        mean_redox_state = np.dot(redox_percentages, molecules) / molecules.sum()

        # Plot each distribution
        ax = axes[i]
        ax.bar(k_values, molecules, color="skyblue", edgecolor="black")
        ax.set_title(f"{dist} Distribution\nAggregate Redox State: {mean_redox_state:.2f}%", fontsize=12)
        ax.set_xlabel("k-grades (% Oxidized)")
        ax.set_ylabel("Number of Molecules")
        ax.set_xticks(k_values)
        ax.grid(axis="y", linestyle="--", alpha=0.7)

    # Remove empty subplot if the number of plots is not a perfect multiple of grid size
    if len(distributions) < len(axes):
        for j in range(len(distributions), len(axes)):
            fig.delaxes(axes[j])

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig("combined_distributions.png", dpi=300)
    plt.show()

# List of distributions
distributions = ["Super-Poisson Reduced", "Super-Poisson Oxidized", "Gaussian", "Polarized"]

# Plot and combine distributions
plot_combined_distributions(distributions, total_molecules, R)
