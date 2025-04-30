from config import path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_superimposed_on_image(image_file, star_positions, output_file, table = None, colors=None):
    plt.style.use('dark_background')
    
    # Load the image data
    image_data = fits.getdata(image_file)
    
    # Define the logarithmic scale for the image
    log_norm = LogNorm(vmin=1, vmax=np.max(image_data))
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Display the image
    ax.imshow(image_data, cmap='gray', origin='lower', norm=log_norm)
    
    # Colorcode the stars based on some criteria if colors is not None
    if colors == 'sigma':
        # Generate colors based on sigma values
        colors = generate_sigma_colors(table)
        ax.scatter(star_positions[:, 0], star_positions[:, 1], c=colors, s=10, edgecolor='white', alpha=0.5)
        plt.legend({'really far away'},loc='upper right')
        ax.set_title('Sigma of Stars Identified in our Image')
    
    elif colors == 'distance':
        # Generate colors based on distance values
        colors = generate_distance_colors(table)
        print("Colors:", colors)
        print('type of colors[0]', type(colors[0]))
        red = np.where(colors == 'red', True, False)
        print("Red stars:", red)
        print("Red stars:", star_positions[red])
        green = np.where(colors=='green', True, False)
        blue = np.where(colors == 'blue')
        yellow = np.where(colors == 'yellow')
        
        ax.scatter(star_positions[red, 0], star_positions[red, 1], c='red', s=10, edgecolor='white', alpha=0.5)
        ax.scatter(star_positions[green, 0], star_positions[green, 1], c='green', s=10, edgecolor='white', alpha=0.5)
        ax.scatter(star_positions[blue, 0], star_positions[blue, 1], c='blue', s=10, edgecolor='white', alpha=0.5)
        ax.scatter(star_positions[yellow, 0], star_positions[yellow, 1], c='yellow', s=10, edgecolor='white', alpha=0.5)
        plt.legend(['really far away','<1500 parsecs', '<1000 parsecs', '<500 parsecs'],loc='upper right')
        ax.set_title('Distances of Stars Identified in our Image')
        
    elif colors == 'magnitude':
        # Generate colors based on magnitude values
        colors = generate_magnitude_colors(table)
        red = np.where(colors == 'red')
        print("Red stars:", red)
        print("Red stars:", star_positions[red])
        green = np.where(colors == 'green')
        blue = np.where(colors == 'blue')
        yellow = np.where(colors == 'yellow')
        
        ax.scatter(star_positions[blue, 0], star_positions[blue, 1], c='blue', s=10, marker='x', alpha=0.5, label='<10')
        ax.scatter(star_positions[green, 0], star_positions[green, 1], c='green', s=10, marker='x', alpha=0.5, label = '10-15')
        ax.scatter(star_positions[yellow, 0], star_positions[yellow, 1], c='yellow', s=10, marker='x', alpha=0.5, label='15-20')
        ax.scatter(star_positions[red, 0], star_positions[red, 1], c='red', s=10, marker='x', alpha=0.5, label = '>20')
        plt.legend(loc='upper right')
        ax.set_title('Magnitudes of Stars Identified in our Image')
    
    else:
        ax.set_title('Stars Identified in our Image')
        ax.scatter(star_positions[:, 0], star_positions[:, 1], c='red', s=10, edgecolor='white', alpha=0.5)
    
    # Set the title and labels
    
    ax.set_xlabel('X Pixel')
    ax.set_ylabel('Y Pixel')
    
    
    
    # Show the figure
    plt.show()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)

def generate_sigma_colors(table):
    """
    Generate colors based on the sigma values of the stars.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        array: numpy array of colors for each star based on their sigma values.
    """
    distance = table['distance'].data
    
    distance = distance[1:]  # Skip the first value (header)
    distance = np.array(distance, dtype=float)
    
    # Some fun with getting a gaussian fit to our distance data to see where our cluster probably is
    filtered_distance = distance[distance < 4000]  # Filter out distances greater than 4000 parsecs
    
    mu = np.nanmedian(filtered_distance)
    std = np.nanstd(filtered_distance)
    
    print("Median Distance:", mu)
    print("Standard Deviation:", std)
    
    colors = []
    
    for i in range(len(distance)):
        if distance[i] < mu - 3 * std or distance[i] > mu + 3 * std:
            colors.append('blue')
        elif distance[i] < mu - 2 * std or distance[i] > mu + 2 * std:
            colors.append('green')
        elif distance[i] < mu - std or distance[i] > mu + std:
            colors.append('yellow')
        else:
            colors.append('red')
    
    return np.array(colors)

def generate_distance_colors(table):
    """
    Generate colors based on the distance values of the stars.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        array: numpy array of colors for each star based on their distance values.
    """
    distance = table['distance'].data
    
    distance = distance[1:]  # Skip the first value (header)
    distance = np.array(distance, dtype=float)

    colors = []
    
    for value in distance:
        if value < 500:
            colors.append('blue')
        elif value < 1000:
            colors.append('green')
        elif value < 1500:
            colors.append('yellow')
        else:
            colors.append('red')
    
    return np.array(colors)

def generate_magnitude_colors(table):
    """
    Generate colors based on the magnitude values of the stars.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        array: numpy array of colors for each star based on their magnitude values.
    """
    # Assuming the table has a 'magnitude' column
    magnitudes = table['V_apparent'].data
    magnitudes = magnitudes[1:]  # Skip the first value (header)
    magnitudes = np.array(magnitudes, dtype=float)
    
    colors = []
    
    for value in magnitudes:
        if value < 10:
            colors.append('blue')
        elif value < 15:
            colors.append('green')
        elif value < 20:
            colors.append('yellow')
        else:
            colors.append('red')
    
    return np.array(colors)
    


if __name__ == "__main__":
    # Example usagew
    image_file = path + 'combined/B_20s_combined.fits'  # Replace with your image file path
    
    # Getting star positions by importing the Ra/Dec table and converting back into pixel (just to test and see)
    from astropy.table import Table
    star_table = Table.read(path + '../star_matching/target_table_part2.tex')

    from star_positions import get_pixel_from_coords
    # Get skycoords data from the star table
    sky_coords = star_table["sky_coord"].data[1:]
    
    # These are formatted as strings, so we need to convert them to a list of tuples
    sky_coords = [tuple(map(float, coord.split(','))) for coord in sky_coords]
    print("Sky Coordinates (as tuples):", sky_coords)
    
    # Convert to list of skycoord objects
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    sky_coords = [SkyCoord(ra=coord[0],dec=coord[1], unit=u.deg) for coord in sky_coords] 
    print("Sky Coordinates:", sky_coords)
    
    # Convert sky coordinates to pixel positions
    star_positions = get_pixel_from_coords(sky_coords)
    print("Star Positions (pixel):", star_positions)
    
    output_file = path + 'superimposed_stars.png'  # Output file path
    
    # Generate colors based on sigma values
    colors = generate_distance_colors(star_table)
    
    # Call the function to plot and save the image
    plot_superimposed_on_image(image_file, star_positions, output_file, table=star_table, colors='magnitude')