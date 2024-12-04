import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Function to generate Pascal's Triangle
def generate_pascals_triangle(rows):
    triangle = [[1]]
    for i in range(1, rows):
        row = [1]
        for j in range(1, i):
            row.append(triangle[i-1][j-1] + triangle[i-1][j])
        row.append(1)
        triangle.append(row)
    return triangle

# Data for the table
data = [
    ["k", "%-oxidised", "i-states", "Structure"],
    ["0", "0", "1", "[000]"],
    ["1", "33.3", "3", "[100], [010], [001]"],
    ["2", "66.6", "3", "[110], [011], [101]"],
    ["3", "100", "1", "[111]"],
]

# Generate Pascal's Triangle
rows = 10
triangle = generate_pascals_triangle(rows)

# Create the figure and define layout
fig = plt.figure(figsize=(12, 6), dpi=300)  # High-resolution image
gs = GridSpec(1, 2, width_ratios=[2, 1])  # Two sections: figure and table

# Left panel: Pascal's Triangle
ax_triangle = fig.add_subplot(gs[0])
for i, row in enumerate(triangle):
    for j, value in enumerate(row):
        ax_triangle.text(j - i / 2, -i, str(value), ha='center', va='center', fontsize=12)

# Styling the Pascal's Triangle
ax_triangle.axis('off')
ax_triangle.set_aspect('equal')
ax_triangle.set_xlim(-5, 10)
ax_triangle.set_ylim(-10, 2)

# Right panel: Table
ax_table = fig.add_subplot(gs[1])
ax_table.axis('tight')
ax_table.axis('off')

# Create the table
table = ax_table.table(
    cellText=data,
    colLabels=None,
    cellLoc='center',
    loc='center',
    colWidths=[0.2, 0.3, 0.3, 0.8],  # Adjust column widths
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5)  # Adjust table size

# Save the combined figure
plt.tight_layout()
plt.savefig("combined_figure_and_table.png", dpi=300)
plt.show()
