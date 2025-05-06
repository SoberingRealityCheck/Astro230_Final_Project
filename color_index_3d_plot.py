import matplotlib.pyplot as plt 
import numpy as np
from config import path
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from mpl_toolkits.mplot3d import Axes3D


def cartesian_from_spherical(r, ra, dec):
    """
    Convert spherical coordinates to Cartesian coordinates.
    
    Args:
        r (array-like): Radius.
        theta (array-like): Polar angle (inclination).
        phi (array-like): Azimuthal angle (longitude).
    
    Returns:
        tuple: Cartesian coordinates (x, y, z).
    """
    x = r * np.sin(dec) * np.cos(ra)
    y = r * np.sin(dec) * np.sin(ra)
    z = r * np.cos(dec)
    return x, y, z

def extract_star_positions(table):
    """
    Extract star positions from the table and convert to Cartesian coordinates.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        np.ndarray: Array of star positions in Cartesian coordinates.
    """
    # Get the distances and angles from the table
    distance = table['distance'].data
    distance = np.array(distance, dtype=float)
    sky = table['sky_coord'].data
    sky = [x.strip().split(',') for x in sky]
    sky = np.array(sky, dtype=float)
    
    # convert the sky coord tuples to separate arrays
    ra = sky[:, 0] # Polar angle (inclination)
    dec = sky[:, 1]    # Azimuthal angle (longitude)
    
    
    # Convert spherical coordinates to Cartesian coordinates
    x, y, z = cartesian_from_spherical(distance, np.radians(ra), np.radians(dec))
    
    return np.array([x, y, z]).T

def plot_3d_scatter(x, y, z, color, title, xlabel, ylabel, zlabel):
    """
    Create a 3D scatter plot.
    
    Args:
        x (array-like): X-axis data.
        y (array-like): Y-axis data.
        z (array-like): Z-axis data.
        color (array-like): Color data for the points.
        title (str): Title of the plot.
        xlabel (str): Label for the X-axis.
        ylabel (str): Label for the Y-axis.
        zlabel (str): Label for the Z-axis.
    """ 
    # Some cool styling for the plot
    plt.style.use('dark_background')
    plt.rcParams['grid.color'] = (0.5, 0.5, 0.5, 0)
   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
    
    # Get the color range so it properly makes white stars at B-V = 1
    norm = plt.Normalize(-1, 3)
    sc = ax.scatter(x, y, z, c=color, cmap='twilight_ebushifted', alpha=.9, norm=norm)
    
    #ax.scatter(0, 0, 0, marker='o', color='yellow', s=10)  # Our Sun is at the origin
    
    # Label Stuff
    plt.colorbar(sc, label='Color Index (B-V)')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    
    # Make the x, y, and z planes transparent
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')
    
    plt.show()



if __name__ == "__main__":
    # Read the target table
    star_table = Table.read(path + '../star_matching/target_table_part4.tex')
    
    # Extract star positions in Cartesian coordinates
    star_positions = extract_star_positions(star_table)
    color = star_table['abs_B'].data -  star_table['abs_V'].data
    
    # Filter to just cluster members
    color = color[[500 < x < 1000 for x in star_table['distance'].data]]
    star_positions = star_positions[[500< x < 1000 for x in star_table['distance'].data]]
    
    print("Plotting 3D scatter plot of star positions...")
    # Plot the 3D scatter plot of star positions
    plot_3d_scatter(star_positions[:, 0], star_positions[:, 1], star_positions[:, 2],
                    color=color, 
                    title='3D Scatter Plot of Star Positions',
                    xlabel='X Coordinate (parsecs)',
                    ylabel='Y Coordinate (parsecs)',
                    zlabel='Z Coordinate (parsecs)')