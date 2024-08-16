# Import necessary libraries
import numpy as np

# Constants
avogadro_number = 6.022e23  # Avogadro's number (molecules/mol)
cell_volume_um3 = 2.6e3  # HeLa cell volume in cubic micrometers
cell_volume_L = cell_volume_um3 * 1e-15  # Convert volume to liters
num_pten_molecules = 3188  # Number of PTEN molecules in the cell
rate_constant = 5.717  # Rate constant in M^-1 s^-1
h2o2_concentration_nM = 10  # Hydrogen peroxide concentration in nanomoles
time_seconds = 3600  # Time in seconds (1 hour)

# Convert PTEN molecules to molarity
pten_concentration_M = num_pten_molecules / (avogadro_number * cell_volume_L)  # Molarity in M
pten_concentration_nM = pten_concentration_M * 1e9  # Convert molarity to nM

# Print PTEN concentration
print(f"PTEN concentration in the cell: {pten_concentration_nM:.2f} nM")

# Calculate the rate of oxidation
h2o2_concentration_M = h2o2_concentration_nM * 1e-9  # Convert H2O2 concentration to M
oxidation_rate = rate_constant * pten_concentration_M * h2o2_concentration_M  # Rate in s^-1

# Print oxidation rate
print(f"Oxidation rate of PTEN Cys124: {oxidation_rate:.2e} s^-1")

# Calculate the fraction of PTEN oxidized over the given time
fraction_oxidized = 1 - np.exp(-oxidation_rate * time_seconds)

# Print fraction oxidized
print(f"Fraction of PTEN oxidized over {time_seconds/3600:.1f} hours: {fraction_oxidized:.6f}")

# Calculate the number of oxidized PTEN molecules
num_oxidized_pten = num_pten_molecules * fraction_oxidized

# Print the number of oxidized PTEN molecules
print(f"Expected number of oxidized PTEN molecules: {num_oxidized_pten:.2f}")

# Summary of results
print("\nSummary:")
print(f"- PTEN concentration: {pten_concentration_nM:.2f} nM")
print(f"- Oxidation rate: {oxidation_rate:.2e} s^-1")
print(f"- Fraction oxidized over {time_seconds/3600:.1f} hours: {fraction_oxidized:.6f}")
print(f"- Number of oxidized PTEN molecules: {num_oxidized_pten:.2f}")
