import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from config import path

# open the table of stars, their distances, and magnitudes
target_table = Table.read(path + '../star_matching/target_table_part4.tex')

# Get the absolute magnitudes of the stars
abs_V = target_table['abs_V'].data
abs_R = target_table['abs_R'].data
abs_B = target_table['abs_B'].data

# Convert to numpy arrays
abs_V = np.array(abs_V, dtype=float)
abs_R = np.array(abs_R, dtype=float)
abs_B = np.array(abs_B, dtype=float)

# Cool dark background
#plt.style.use('dark_background')
plt.style.use('bmh')

distance = target_table['distance'].data  # Get the distances of the stars
distance = np.array(distance, dtype=float)  # Convert to numpy array

center_guess = np.nanmedian(distance)  # Center guess for the cluster center


B_V = abs_B - abs_V  # Color index (B-V)
V_R = abs_V - abs_R  # Color index (V-R)

sigma = target_table['sigma'].data  # Get the sigma values of the stars
sigma = np.array(sigma, dtype=float)  # Convert to numpy array
# Create a scatter plot CMD

cluster_sigma = sigma[np.where(abs(sigma) < 1)]  # Get the sigma values of stars within 1 sigma
print("Cluster Sigma Count:", len(cluster_sigma))

cluster_distance = distance[np.where(abs(sigma) < 1)]  # Get the radial distances of stars within 1 sigma
print("Cluster Distance Count:", len(cluster_distance))

cluster_color = B_V[np.where(abs(sigma) < 1)]  # Get the color index of stars within 1 sigma
print("Cluster Color Count:", len(cluster_color))

cluster_center = np.nanmedian(cluster_distance)  # Center of the cluster

radial_distance = np.abs(cluster_distance - cluster_center)  # Calculate the radial distance from the cluster center



norm = plt.Normalize(0, 2)
plt.scatter(radial_distance, cluster_color, label='B-V', alpha=0.5, c=cluster_color, cmap='twilight_shifted', norm=norm)

#plt.scatter(radial_distance, V_R, label='V-R', alpha=0.5)
#plt.scatter(abs_V - abs_R, abs_V, label='V-R', alpha=0.5)

# Add labels and title
plt.xlabel('Radial Distance from Cluster Center (parsecs)')
plt.ylabel('Color Index (B-V)')
plt.title('Radial Distance vs Color Index for Definitive($\sigma < 1$) Cluster Members')
plt.legend()

# Identify indices with actual non-NaN values
idx = np.isfinite(radial_distance) & np.isfinite(cluster_color)

# Fit a linear regression line to the data
z = np.polyfit(radial_distance[idx], cluster_color[idx], 1)
p = np.poly1d(z)
print("Linear Fit Coefficients:", z)
plt.plot(radial_distance, p(radial_distance), color='red', label='Linear Fit', alpha=0.5)
plt.legend()


# Save the plot to a file
plt.savefig(path + '../analysis/Distance-Vs-Color_light.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

