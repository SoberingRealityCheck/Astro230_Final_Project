import matplotlib.pyplot as plt 
import numpy as np
from config import path
from astropy.table import Table

#extract the x-y-z locations of the stars from the table from the spherical coordinates given
def extract_star_positions(table):
    """
    Extract star positions from the table and convert to 3D coordinates.
    
    Args:
        table (astropy.table.Table): Table containing star data.
        
    Returns:
        np.ndarray: Array of star positions in 3D coordinates.
    """
    # Extract the RA, Dec, and Parallax values from the table
    sky = table['sky_coord'].data
    parallax = table['distance'].data
    
    parallax = parallax[1:]  # Skip the first value (assumed to be a header or invalid data)
    parallax = np.array(parallax, dtype=float)
    
    sky = sky[1:]
    sky = [x.strip().split(',') for x in sky]
    sky = np.array(sky, dtype=float)
    print("ra type", type(sky))
    print("ra shape", sky.shape)
    print("ra", sky)
    
    ra = sky[:, 0]  # Right Ascension in degrees
    dec = sky[:, 1]  # Declination in degrees
    # Convert RA, Dec, and Parallax to 3D coordinates
    x = parallax * np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y = parallax * np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z = parallax * np.sin(np.radians(dec))
    
    return np.array([x, y, z]).T
    
star_table = Table.read(path + '../star_matching/target_table.tex')
star_positions = extract_star_positions(star_table)
print("Star Positions (3D):", star_positions)

print("Median Star Position:", np.median(star_positions, axis=0))

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title('3D Plot of Star Positions')
ax.scatter(star_positions[:, 0], star_positions[:, 1], star_positions[:, 2], c='b', marker='o')
ax.set_xlim3d(-139,-138)
ax.set_ylim3d(677,678)
ax.set_zlim3d(-263, -262)



plt.show()