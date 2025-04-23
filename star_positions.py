#get the actual x-y pixel locations of the stars in the image
# Core Imports
import numpy as np
import matplotlib.pyplot as plt

# AstroPy Imports
from astropy import table
from astropy.io import fits
from astropy.wcs import WCS
from photutils.detection import DAOStarFinder
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import fit_wcs_from_points
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

# My Imports
from get_files import fits_get
from config import path

def get_pixel_positions(file, threshold = 3, plot=False):
    """
    Get the pixel positions of stars in the image.
    
    Args:
        file (filepath): filepath to the image.
        threshold (float): Threshold for star detection.
        plot (bool): Whether to plot the images with detected stars.
    
    Returns:
        list: List of pixel positions of stars in the images.
    """
    
    # Step 1: Load the image
    image = fits.getdata(file)
    
    # Step 2: Find the stars in the images using a star detection algorithm 
    fwhm = 2.0  # Full Width at Half Maximum (FWHM) of the stars
    starfinder = DAOStarFinder(threshold=threshold, fwhm=fwhm)
    
    # Checking the code with a single image.
    star_locations = starfinder.find_stars(image)
    print("Star Locations Found:", star_locations)
    # Step 2.5(optional): Plot the images with detected stars for visual verification
    if plot:
        # Get the nice looking filename for the plot title
        filename = file.split('/')[-1]
        filename = filename.split('\\')[1]
        
        # Also define a nice logarithmic scale for our images
        from matplotlib.colors import LogNorm
        log_norm = LogNorm(vmin=10, vmax = np.max(image))
        
        plt.style.use('dark_background')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))
        ax1.imshow(image, cmap='gray', origin='lower', norm=log_norm)
        ax1.set_title(f'Original Image: \n {filename}')
        ax1.set_xlabel('X Pixel')
        ax1.set_ylabel('Y Pixel')
        
        ax2.imshow(image, cmap='gray', origin='lower', norm=log_norm)
        ax2.scatter(star_locations['xcentroid'], star_locations['ycentroid'], s=1, color='red')
        ax2.set_title(f'Detected Stars, threshold={threshold}, fwhm={fwhm}')
        ax2.set_xlabel('X Pixel')
        ax2.set_ylabel('Y Pixel')
        
        plt.show()
    
    # Step 3: Return the pixel positions of the stars
    star_positions = []
    for star in star_locations:
        star_positions.append([star['xcentroid'], star['ycentroid']])
    return np.array(star_positions)

def generate_pixel_to_world_matrix():
    '''
    Generate a pixel to world coordinate transformation matrix.
    '''
    
    # HD 49091 -- 416 x, 420 y, 06 46 03.18	-20 43 18.62, deg 6.7675500  -20.7218389
    # HD 49126 -- 461 x, 443 y, 06 46 07.20 -20 45 15.55, deg 6.7686667  -20.7543194
    # HD 49023 -- 298 x, 343 y, 06 45 35.51 -20 40 51.38, deg 6.7598639  -20.6809389 
    # HD 49317 -- 755 x, 293 y, 06:47:02.63 -20:40:32.82, deg 6.7840639  -20.6757833
    # HD 48983 -- 263 x, 505 y, 06:45:28.00 -20:50:23.00, deg 6.7577778  -20.8397222

    #found this code on stack overflow. hell yeah
    # https://stackoverflow.com/questions/63594323/astropy-wcs-transfromation-matrix
    
    stars = SkyCoord(ra=[6.7675500, 6.7686667, 6.7598639, 6.7840639, 6.7577778], 
                    dec=[-20.7218389, -20.7543194, -20.6809389, -20.6757833, -20.8397222], 
                    unit=(u.hourangle ,u.deg))
    pixels_x = np.array([416, 461, 298, 755, 263])
    pixels_y = np.array([420, 443, 343, 293, 505])
    # Create a WCS object and set the pixel coordinates and world coordinates
    pixel_to_world_wcs = fit_wcs_from_points((pixels_x, pixels_y), stars, projection='AIR')
    print("wcs", pixel_to_world_wcs)
    print("wcs_pixel_n_dim:", pixel_to_world_wcs.pixel_n_dim)
    print("wcs_world_n_dim:", pixel_to_world_wcs.world_n_dim)
    print("wcs_array_shape:", pixel_to_world_wcs.array_shape)
    
    print("WCS TEST")
    print("--------------------")
    print("star 1 pixel:", pixel_to_world_wcs.world_to_pixel(stars[0]))
    print("star 1 actual pixel:", (415.6, 419.7))
    print("")
    return pixel_to_world_wcs

def get_star_locations(star_positions, wcs):
    """
    Convert pixel locations to RA/Dec locations using a known reference point and the FOV of the sensor.

    Args:
        star_positions (Numpy array): 2d array of pixel positions of stars.
        reference_point (tuple): Reference point (RA, Dec) in degrees.
        fov (float): Field of view in degrees.

    Returns:
        list: List of RA/Dec locations of stars.
    """
    skyvals = []
    print("star_positions:", star_positions)
    print("star_positions shape:", star_positions.shape)
    print("star_positions individual value shape:", star_positions[0].shape)
    for starpos in star_positions:
        skyval_new = (pixel_to_skycoord(starpos[0], starpos[1], wcs, origin=0))
        skyvals.append(skyval_new)
    print(skyvals)  # (RA, Dec) in degrees
    return skyvals

if __name__ == "__main__":
    # Example usage
    PLOT_PIXEL = False
    PLOT_WCS = False
    
    files = fits_get(path + 'combined/')
    print("Files to be processed:", files)
    for file in files:
        filename = file.split('/')[-1]
        filename = filename.split('\\')[1]
        print("Filename:", filename)
        file_id = filename.split('.')[0]
        print("Processing file:", file)
        pixel_positions = get_pixel_positions(file, threshold=50, plot=PLOT_PIXEL)
        print("Pixel Positions of Stars:", pixel_positions)
        pixel_to_world_wcs = generate_pixel_to_world_matrix()
        skyvals = get_star_locations(pixel_positions, pixel_to_world_wcs)
        if PLOT_WCS:
            for i in range(len(skyvals)):
                plt.scatter(skyvals[i].ra, skyvals[i].dec, s=1, color='red')
            plt.xlabel('RA (degrees)')
            plt.ylabel('Dec (degrees)')
            plt.title('Star Locations in RA/Dec')
            plt.show()
            plt.xlim(103,105)
            plt.ylim(-21,-19)
        # Save the pixel positions to a numpy file
        np.save(path + f'../star_pos_identification/{file_id}_pixel_positions.npy', pixel_positions)
        np.save(path + f'../star_pos_identification/{file_id}_skyvals.npy', skyvals)