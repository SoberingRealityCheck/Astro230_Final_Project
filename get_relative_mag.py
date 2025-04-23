from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from photutils import psf
import numpy as np


def get_counts(pixel_positions, imdata):
    '''
    this function takes in an array of star positions and an array of image data and returns the photometric counts for each star.
    Most of this code is taken directly from megan-the-astronomer(my professor)'s 'whats_in_an_image' notebook because I am a lazy bastard.
    
    Args:
        pixel_positions (ndarray): List of coordinates (x, y) of stars, formatted as floats or ints.
        imdata (ndarray): image data from reduced FITS file
    
    Returns:
        counts: ndarray (ordered identically to pixel_positions) with the summed counts (background-subtracted) of each star in the image.
    '''
    # Format our pixel positions into a list of tuples so photutils is happy
    #px_pos = list(zip(pixel_positions[:,0], pixel_positions[:,1]))
    #print("Pixel Positions:", px_pos)
    #print("Pixel Positions Shape:", pixel_positions.shape)
    
    # Get FWHMs of the stars 
    fwhms = psf.fit_fwhm(imdata, xypos=pixel_positions, fit_shape=7)
    apertures = CircularAperture(pixel_positions, r=np.nanmedian(fwhms))  
    annuli = CircularAnnulus(pixel_positions, r_in=6, r_out=10)
    phot_table = aperture_photometry(pixel_positions, apertures)  
    print("Photometry Table:", phot_table)
    
    pass

def magnitudes_from_counts(counts, reference_counts, reference_magnitude):
    '''
    This function takes in the counts of stars and returns their magnitudes based on a reference star.
    Uses the classic magnitude calculation formula, the most fucked-up piece of math to ever grace the astronomical domain. 
    Don't ask me why the sun is negative. Don't ask me why it increases as stuff gets dimmer. That's just how it be. Astronomy is weird.
    
    Args:
        counts (ndarray): Array of counts for each star.
        reference_counts (float): Counts of the reference star.
        reference_magnitude (float): Magnitude of the reference star.
    
    Returns:
        magnitudes (ndarray): Array of magnitudes for each star.
    '''
    return -2.5 * np.log10(counts / reference_counts) + reference_magnitude


if __name__ == "__main__":
    # Example usage
    
    from config import path
    from get_files import fits_get
    from astropy.io import fits
    
    pixel_positions = np.load(path + '../star_pos_identification/pixel_positions.npy', allow_pickle=True)
    
    files = fits_get(path + 'combined/')
    print("Files to be processed:", files)
    
    for i, file in enumerate(files):
        image_data, header = fits.getdata(file, header=True)
        # Assuming imdata is the image data from the FITS file
        imdata = image_data
        
        counts = get_counts(pixel_positions, imdata)
        print("Counts Obtained for file", i)
        print(counts)
        
        # Save the counts to a numpy file
        np.save(path + '../star_pos_identification/counts_' + str(i) + '.npy', counts)
        
        # Get the magnitudes from the counts
        magnitudes = magnitudes_from_counts(counts)