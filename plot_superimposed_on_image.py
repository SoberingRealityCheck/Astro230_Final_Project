from config import path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_image(image_file_V, image_file_B, image_file_R, fig, ax):
    # Load the image data
    image_V = fits.getdata(image_file_V)
    image_B = fits.getdata(image_file_B)
    image_R = fits.getdata(image_file_R)
    print("Image Data:", image_V)
    
    # Clip negative values to zero
    image_V[image_V < 0] = 0
    image_B[image_B < 0] = 0
    image_R[image_R < 0] = 0
    print("Image Data:", image_V)
    
    # Define the logarithmic scale for the image
    log_V = LogNorm(vmin=5, vmax=np.max(image_V), clip=True)
    log_B = LogNorm(vmin=5, vmax=np.max(image_B), clip=True)
    log_R = LogNorm(vmin=5, vmax=np.max(image_R), clip=True)
    
    
    #map counts in the image to a log scale from 0-255
    color_V = log_V(image_V) * 255
    color_B = log_B(image_B) * 255
    color_R = log_R(image_R) * 255
    
    # Convert to uint8 for display
    color_V = color_V.astype(np.uint8)
    color_B = color_B.astype(np.uint8)
    color_R = color_R.astype(np.uint8)
    
    print("Color V:", color_V)
    # Stack the images along the third axis
    image_data = np.array([color_R, color_V, color_B]).transpose(1, 2, 0)  # Stack the images along the third axis
    print("Image Stack:", image_data)
    
    # Display the image
    ax.imshow(image_data, origin='lower', aspect='auto')

def plot_superimposed_on_image(image_file_V, image_file_B, image_file_R, star_positions, output_file, table = None, colors=None):
    plt.style.use('dark_background')
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    
    plot_image(image_file_V, image_file_B, image_file_R, fig, ax)
    
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
        #print("Colors:", colors)
        #print('type of colors[0]', type(colors[0]))
        red = np.where(colors == 'red', True, False)
        #print("Red stars:", red)
        #print("Red stars:", star_positions[red])
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
        #print("Red stars:", red)
        #print("Red stars:", star_positions[red])
        green = np.where(colors == 'green')
        blue = np.where(colors == 'blue')
        yellow = np.where(colors == 'yellow')
        
        ax.scatter(star_positions[blue, 0], star_positions[blue, 1], c='blue', s=10, marker='x', alpha=0.75, label='<10')
        ax.scatter(star_positions[green, 0], star_positions[green, 1], c='green', s=10, marker='x', alpha=0.75, label = '10-15')
        ax.scatter(star_positions[yellow, 0], star_positions[yellow, 1], c='yellow', s=10, marker='x', alpha=0.75, label='15-20')
        ax.scatter(star_positions[red, 0], star_positions[red, 1], c='red', s=10, marker='x', alpha=0.75, label = '>20')
        plt.legend(loc='upper right')
        ax.set_title('Apparent V band Magnitudes of Stars Identified in our Image')
    elif colors == 'B-V':
        # Generate colors based on B-V values
        norm = plt.Normalize(-1, 1)
        c = generate_color_colors(table)
        sc = ax.scatter(star_positions[:,0],star_positions[:,1], c=c, cmap='twilight_shifted', alpha=.75, norm=norm, marker='x', label='B-V')
        plt.colorbar(sc, label='Color Index (B-V)', orientation='vertical')
        ax.set_title('Color Index (B-V) of Stars Identified in our Image')

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

def generate_color_colors(table):
    """
    Generate colors based on the color values of the stars.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        array: numpy array of colors for each star based on their color values.
    """
    # Assuming the table has a 'color' column
    V_apparent = table['V_apparent'].data
    B_apparent = table['B_apparent'].data

    colors = B_apparent - V_apparent  # Calculate color index (B-V)
    colors = colors[1:]  # Skip the first value (header)
    colors = np.array(colors, dtype=float)
    
    return colors


if __name__ == "__main__":
    # Example usagew
    image_file_B = path + 'combined/B_60s_combined.fits'  # Replace with your image file path
    image_file_V = path + 'combined/V_20s_combined.fits'  # Replace with your image file path
    image_file_R = path + 'combined/R_10s_combined.fits'  # Replace with your image file path
    
    # Getting star positions by importing the Ra/Dec table and converting back into pixel (just to test and see)
    from astropy.table import Table
    star_table = Table.read(path + '../star_matching/target_table_part4.tex')

    from star_positions import get_pixel_from_coords
    # Get skycoords data from the star table
    sky_coords = star_table["sky_coord"].data[1:]
    
    # These are formatted as strings, so we need to convert them to a list of tuples
    sky_coords = [tuple(map(float, coord.split(','))) for coord in sky_coords]
    #print("Sky Coordinates (as tuples):", sky_coords)
    
    # Convert to list of skycoord objects
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    sky_coords = [SkyCoord(ra=coord[0],dec=coord[1], unit=u.deg) for coord in sky_coords] 
    #print("Sky Coordinates:", sky_coords)
    
    # Convert sky coordinates to pixel positions
    star_positions = get_pixel_from_coords(sky_coords)
    #print("Star Positions (pixel):", star_positions)
    
    output_file = path + 'superimposed_stars.png'  # Output file path
    
    # Generate colors based on sigma values
    colors = generate_distance_colors(star_table)
    
    # Call the function to plot and save the image
    plot_superimposed_on_image(image_file_V, image_file_B, image_file_R, star_positions, output_file, table=star_table, colors='B-V')