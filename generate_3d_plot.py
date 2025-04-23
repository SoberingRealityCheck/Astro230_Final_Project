import matplotlib.pyplot as plt 
import numpy as np
from config import path
from astropy.table import Table
from scipy.stats import norm

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
    distance = table['distance'].data
    
    distance = distance[1:]  # Skip the first value (header)
    distance = np.array(distance, dtype=float)
    
    # Crop out the crazy outliers from the data. Stuff that is farther than 10000 parsecs is probably not in our cluster
    print("Distance before cropping", distance)
    
    
    for i in range(len(distance)):
        if distance[i] > 4000:
            distance[i] = np.nan
    
    # Some fun with getting a gaussian fit to our distance data to see where our cluster probably is
    mu = np.nanmedian(distance)
    std = np.nanstd(distance)
    xp = np.linspace(0, 4000, 100)
    p = norm.pdf(xp, mu, std) * 100000
    print("P", p)
    plt.style.use('dark_background')
    plt.plot(xp, p, linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.xlabel('Radial Distance (parsecs)')
    plt.ylabel('Star Density / Probability Density')
    plt.title(title)
    plt.hist(distance, bins=100)
    plt.legend(['Gaussian Fit', 'Data'])
    plt.show()
    
    # Our median distance is 739.64 parsecs and our standard deviation is 904.92 parsecs.
    # If we only keep data within 2 sigma of the median, we should get a much better picture of a cluster.
    '''
    for i in range(len(distance)):
        if distance[i] > mu + 3 * std or distance[i] < mu - 3 * std:
            distance[i] = np.nan
    '''
    
    sky = sky[1:]
    sky = [x.strip().split(',') for x in sky]
    sky = np.array(sky, dtype=float)
    print("ra type", type(sky))
    print("ra shape", sky.shape)
    print("ra", sky)
    print("Median Angle:", np.nanmedian(sky, axis=0))
    print("Median Distance:", np.nanmedian(distance))
    
    ra = sky[:, 0]  # Right Ascension in degrees
    dec = sky[:, 1]  # Declination in degrees
    
    # Convert RA, Dec, and Distance to 3D coordinates
    x = distance * np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y = distance * np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z = distance * np.sin(np.radians(dec))
    
    # Define a 'color' based on the star's sigma value from the cluster
    c = np.zeros(len(distance))
    for i in range(len(distance)):
        if distance[i] > mu + 3 * std or distance[i] < mu - 3 * std:
            c[i] = 1
        
        # Check for the nan values from our cropping
        elif np.isnan(distance[i]):
            c[i] = 1
        
        else:
            # Should range from -3 to 3
            c[i] = ((distance[i] - mu) / std)
            # normalize to 0-1
            print("type of c[i]", type(c[i]))
            print("c[i]", c[i])
            c[i] = (int(abs(c[i])) / 4)


    return np.array([x, y, z, c]).T
star_table = Table.read(path + '../star_matching/target_table.tex')
star_positions = extract_star_positions(star_table)
print("Star Positions (3D):", star_positions)

print("Median Star Position:", np.median(star_positions, axis=0))

# Dark mode plots because hell yeah
with plt.style.context('dark_background'):
    # Set the colormap to 'Spectral' for better visibility and also coolness
    plt.set_cmap('Spectral')
    plt.rcParams['grid.color'] = (0.5, 0.5, 0.5, 0)
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    
    ax.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
    ax.dist = 100000
    # Set the labels and title
    ax.set_xlabel('X-axis (parsecs)')
    ax.set_ylabel('Y-axis (parsecs)')
    ax.set_zlabel('Z-axis (parsecs)')
    ax.set_title('3D Plot of Star Positions')
    
    # Scatter plot of star positions
    ax.scatter(star_positions[:, 0], star_positions[:, 1], star_positions[:, 2], marker='o', c = star_positions[:, 3], s=1, alpha=1)
    ax.scatter(0, 0, 0, marker='o', color='yellow', s=10)  # Our Sun is at the origin
    
    # Set the limits for each axis to be, like, reasonable
    ax.set_xlim3d(-139,-138)
    ax.set_ylim3d(677,678)
    ax.set_zlim3d(-263, -262)
    
    # Set the initial viewing angle to be pretty close to our original telescope angle
    ax.view_init(elev=-200.74746374, azim=101.53983887)
    
    # Make the x, y, and z planes transparent
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')



plt.show()