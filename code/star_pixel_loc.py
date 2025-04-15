#get the actual x-y pixel locations of the stars in the image
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from get_files import fits_get

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
    
    # Step 2: Find the stars in the images using a star detection algorithm (e.g., DAOFIND, SExtractor)
    
    
    # Step 2.5(optional): Plot the images with detected stars for visual verification
    
    
    # Step 3: Return the pixel positions of the stars
    return []