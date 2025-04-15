#get the actual x-y pixel locations of the stars in the image
# Core Imports
import numpy as np
import matplotlib.pyplot as plt

# AstroPy Imports
from astropy.io import fits
from astropy import wcs
from photutils.detection import DAOStarFinder

# My Imports
from get_files import fits_get
from config import path

def get_pixel_positions(files, threshold = 3, plot=False):
    """
    Get the pixel positions of stars in the image.
    
    Args:
        files (list): List of file paths to the images.
        threshold (float): Threshold for star detection.
        plot (bool): Whether to plot the images with detected stars.
    
    Returns:
        list: List of pixel positions of stars in the images.
    """
    
    # Step 1: Load the images
    images = [fits.getdata(file) for file in files]
    
    # Step 2: Find the stars in the images using a star detection algorithm 
    starfinder = DAOStarFinder(threshold=threshold, fwhm=2.0)
    # Assuming images[0] is gonna have the same stars as the rest of our images? need to check this
    star_locations = starfinder.find_stars(images[0])
    print("Star Locations Found:", star_locations)
    # Step 2.5(optional): Plot the images with detected stars for visual verification
    if plot:
        plt.imshow(images[0], cmap='gray', origin='lower', vmin=0, vmax=255)
        plt.scatter(star_locations['xcentroid'], star_locations['ycentroid'], s=1, color='red')
        plt.title('Detected Stars')
        plt.show()    
    
    # Step 3: Return the pixel positions of the stars
    return star_locations['xcentroid'], star_locations['ycentroid']

#convert pixel locations to ra/dec locations using some known reference point and the fov of the sensor

def get_star_locations(star_positions, reference_points, fov):
    """
    Convert pixel locations to RA/Dec locations using a known reference point and the FOV of the sensor.

    Args:
        star_positions (Numpy array): 2d array of pixel positions of stars.
        reference_point (tuple): Reference point (RA, Dec) in degrees.
        fov (float): Field of view in degrees.

    Returns:
        list: List of RA/Dec locations of stars.
    """
    
    # Step 1: Unpack the reference points (these need to be Known Standard Stars in our image that we know both the pixel x-y location and the RA/Dec location of)
    
    
    # Step 2: Calculate the initial pixel scale based on the FOV and the image size
    
    
    # Step 3: Find the average offsets of our linear model to get the most accurate pixel to WCS conversion (since there might be some inaccuracy introduced by our low resolution images)
    
    
    # Step 4: Convert pixel positions to RA/Dec using the reference points and the pixel scale
    
    
    # Step 5: Return a numpy array of RA/Dec locations of stars
    # Set the WCS information manually by setting properties of the WCS
    # object.

    #drumroll.... 
    # COPIED AND PASTED DIRECTLY FROM WCS ASTROPY DOCUMENTATION!!! I SEE NO FLAWS IN THIS PLAN 
    
    from astropy.utils.data import get_pkg_data_filename

    fn = get_pkg_data_filename(path+'reduced/B_20s_0.fits', package='astropy.wcs.tests')

    f = fits.open(fn)

    w = WCS(f[1].header)
    skyvals = []
    for starpos in star_positions:
        skyvals.append(w.pixel_to_world(starpos))
    print(skyvals)  # (RA, Dec) in degrees
    return skyvals


if __name__ == "__main__":
    # Example usage
    files = fits_get(path + 'combined/')
    pixel_positions = get_pixel_positions(files, threshold=30, plot=True)
    print("Pixel Positions of Stars:", pixel_positions)
    skyvals = get_star_locations(pixel_positions, reference_points=[(0, 0)], fov=1.0)
    
    # Save the pixel positions to a numpy file
    np.save(path + '../star_pos_identification/pixel_positions.npy', pixel_positions)