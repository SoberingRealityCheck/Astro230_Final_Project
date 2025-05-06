import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from config import path

# open the table of stars, their distances, and magnitudes
target_table = Table.read(path + '../star_matching/target_table_part4.tex')

# Get the absolute magnitudes of the stars
abs_V = target_table['abs_V'].data[1:]
abs_R = target_table['abs_R'].data[1:]
abs_B = target_table['abs_B'].data[1:]

# Convert to numpy arrays
abs_V = np.array(abs_V, dtype=float)
abs_R = np.array(abs_R, dtype=float)
abs_B = np.array(abs_B, dtype=float)

# Cool dark background
plt.style.use('bmh')

# Create a scatter plot CMD
plt.scatter(abs_B - abs_V, abs_V - abs_R, label='B-V vs V-R', alpha=0.5)
#plt.scatter(abs_V - abs_R, abs_V, label='V-R', alpha=0.5)

# Add labels and title
plt.xlabel('Color Index (B-V)')
plt.ylabel('Color Index (V-R)')
plt.title('Color-Color Diagram')
plt.legend()

# Save the plot to a file
plt.savefig(path + '../Color-Color_Diagram_light.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()