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

import gwcs

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
    
    # HD 49091           -- 416.42371 x, 420.80733 y, 06 46 03.18 -20 43 18.62, hrs 6.7675500 deg -20.7218389
    # HD 49126           -- 467.76442 x, 393.25411 y, 06 46 07.20 -20 45 15.55, hrs 6.7686667 deg -20.7543194
    # HD 49023           -- 304.31240 x, 293.98475 y, 06 45 35.51 -20 40 51.38, hrs 6.7598639 deg -20.6809389 
    # HD 49317           -- 755.22029 x, 293.52529 y, 06:47:02.63 -20:40:32.82, hrs 6.7840639 deg -20.6757833
    # HD 48983           -- 263.55565 x, 504.59604 y, 06:45:28.00 -20:50:23.00, hrs 6.7577778 deg -20.8397222
    
    # HD 48924           -- 169.31156 x, 18.231880 y, 06:45:08.95 -20:28:29.30, hrs 6.7524861 deg -20.4748056
    # CPD-20 1613        -- 487.78894 x, 472.04317 y, 06:46:11.26 -20:48:47.19, hrs 6.7697944 deg -20.8131083
    # CPD-20 1620        -- 539.10617 x, 529.76223 y, 06 46 21.27 -20 51 22.31, hrs 6.7725750 deg -20.8561972
    # CPD-20 1664        -- 941.08451 x, 588.20195 y, 06 47 39.22 -20 53 45.50, hrs 6.7942278 deg -20.8959722
    # CPD-20 1607        -- 467.97932 x, 276.33726 y, 06 46 06.97 -20 39 57.83, hrs 6.7686027 deg -20.6519361
    
    # HD 49105           -- 457.45778 x, 197.94727 y, 06 46 04.84 -20 36 24.91, hrs 6.7680111 deg -20.6069194
    # BD-20 1580         -- 921.87055 x, 186.20412 y, 06 47 34.38 -20 35 39.14, hrs 6.7928833 deg -20.5942055
    # HD 49416           -- 911.25769 x, 436.94065 y, 06 47 33.09 -20 46 58.67, hrs 6.7925249 deg -20.7829638
    # CPD-20 1659        -- 871.59616 x, 301.78246 y, 06 47 25.03 -20 40 52.93, hrs 6.7902861 deg -20.6813694
    # CPD-20 1567        -- 563.45581 x, 651.74896 y, 06 46 26.27 -20 56 52.31, hrs 6.7739638 deg -20.9478638
    
    # UCAC4 345-016898   -- 425.45777 x, 870.12045 y, 06 46 00.08 -21 06 47.28, hrs 6.7666888 deg -21.1131333
    # TYC 5961-2060-1    -- 267.63457 x, 883.33927 y, 06 45 29.55 -21 07 29.36, hrs 6.7582083 deg -21.1248222
    # Cl* NGC 2287 AR 10 -- 258.85233 x, 422.32718 y, 06 45 26.94 -20 46 39.04, hrs 6.7574833 deg -20.7775111
    # CPD-10 1623        -- 548.93354 x, 300.96545 y, 06 46 22.70 -20 41 02.10, hrs 6.7729722 deg -20.6839166
    # TYC 5961-1468-1    -- 65.694216 x, 443.14170 y, 06 44 49.75 -20 47 40.69, hrs 6.7471527 deg -20.7946361
    
    # These are our astronomical coordinates for each star
    stars = SkyCoord(
        ra=[6.7675500, 6.7686667, 6.7598639, 6.7840639, 6.7577778, 
            6.7524861, 6.7697944, 6.7725750, 6.7942278, 6.7686027,
            6.7680111, 6.7928833, 6.7925249, 6.7902861, 6.7739638,
            6.7666888, 6.7582083, 6.7574833, 6.7729722, 6.7471527], 
        dec=[-20.7218389, -20.7543194, -20.6809389, -20.6757833, -20.8397222, 
            -20.4748056, -20.8131083, -20.8561972, -20.8959722, -20.6519361,
            -20.6069194, -20.5942055, -20.7829638, -20.6813694, -20.9478638,
            -21.1131333, -21.1248222, -20.7775111, -20.6839166, -20.7946361], 
        unit=(u.hourangle ,u.deg))
    
    # These are our pixel coordinates for each star
    pixels_x = np.array([416.42371, 467.76442, 304.3124, 755.22029, 263.55565, 
                        169.31156, 487.78894, 539.10617, 941.08451, 467.97932,
                        457.45778, 921.87055, 911.25769, 871.59616, 563.45581,
                        425.45777, 267.63457, 258.85233, 548.93354, 65.694216])
    
    pixels_y = np.array([420.80733, 393.25411, 293.98475, 293.52529, 504.59604, 
                        18.231880, 472.04317, 529.76223, 588.20195, 276.33726,
                        197.94727, 186.20412, 436.94065, 301.78246, 651.74896,
                        870.12045, 883.33927, 422.32718, 300.96545, 443.14170])
    
    # Check the shapes of the arrays
    print("Stars shape:", stars.shape)
    print("pixels_x shape:", pixels_x.shape)
    print("pixels_y shape:", pixels_y.shape)
    # Create a WCS object and set the pixel coordinates and world coordinates
    pixel_to_world_wcs = fit_wcs_from_points((pixels_x, pixels_y), stars)
    print("wcs", pixel_to_world_wcs)
    print("wcs_pixel_n_dim:", pixel_to_world_wcs.pixel_n_dim)
    print("wcs_world_n_dim:", pixel_to_world_wcs.world_n_dim)
    print("wcs_array_shape:", pixel_to_world_wcs.array_shape)
    
    # Check the transform to see if running it forwards and backwards gets us close to where we started 
    print("WCS TEST")
    print("--------------------")
    print("star 1 pixel:", pixel_to_world_wcs.world_to_pixel(stars[0]))
    print("star 1 actual pixel:", pixels_x[0], pixels_y[0])
    print("")
    return pixel_to_world_wcs

def get_coords_from_pixel(star_positions, wcs):
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

def get_pixel_from_coords(skyvals):
    """
    Convert RA/Dec locations to pixel locations using a known reference point and the FOV of the sensor.

    Args:
        skyvals (list): List of RA/Dec locations of stars.
        wcs (WCS): WCS object for the image.

    Returns:
        list: List of pixel locations of stars.
    """
    wcs = generate_pixel_to_world_matrix()
    pixel_positions = []
    for skyval in skyvals:
        pixel_pos_new = (skycoord_to_pixel(skyval, wcs, origin=0))
        pixel_positions.append(pixel_pos_new)
    return np.array(pixel_positions)

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
        skyvals = get_coords_from_pixel(pixel_positions, pixel_to_world_wcs)
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