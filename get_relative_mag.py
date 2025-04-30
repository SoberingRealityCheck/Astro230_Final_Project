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
    phot_table = aperture_photometry(imdata, apertures)  
    #print("Photometry Table:", phot_table)
    
    return phot_table

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
    from astropy.table import Table
    
    target_table = Table.read(path + '../star_matching/target_table.tex')
    # Extract the pixel positions from the target table
    pixel_positions = target_table['pixel_pos'].data[1:]
    
    # Format the pixel positions into a list of tuples so photutils is happy
    pixel_positions = [x.strip().split('..') for x in pixel_positions]
    pixel_positions = np.array(pixel_positions, dtype=float)
    
    print("Pixel Positions Loaded:", pixel_positions)
    print("Pixel Positions Shape:", pixel_positions.shape)
    # Get the list of FITS files to process
    files = fits_get(path + 'combined/')
    print("Files to be processed:", files)
    
    # Loop through each FITS file and perform photometry on each one
    for i, file in enumerate(files):
         
        imdata, header = fits.getdata(file, header=True)
        
        band = header['INSFILTE']
        exptime = int(header['EXPTIME'])
        print(f"Analyzing dataset number {i+1} from {band} band with exposure time of {exptime} seconds...")
        
        counts = get_counts(pixel_positions, imdata)
        print("Counts Obtained for file", i)
        print(counts)
        
        # Save the counts to a numpy file
        np.save(path + '../photometry/counts_' + band + '_' + str(exptime) + '.npy', counts)
        
        counts['name'] = target_table['name'][1:]
        
        print("Counts Table:", counts)
        # Save the counts table to a LaTeX file
        
        
        # Reference Magnitude:
        # HD 48924 (the reference star) has a magnitude of: 
        # 9.28 in the B band, 9.36 in the V band, unknown in the R band. Can at least use it for B and V calibration.
        
        # Going to make a nested dictionary of reference magnitudes and counts for each band and exposure time.
        # All '1' values are just placeholders for the reference counts, which I have not yet calculated. Results will thus be very inaccurate.
        # I will update this dictionary with the correct values once I have them.
        reference_dictionary = {
            'B20': {
                'reference_magnitude': 9.28,
                'reference_counts': 208880.30191664075
            },
            'B60': {
                'reference_magnitude': 9.28,
                'reference_counts': 118826.63344725754
            },
            'V20': {
                'reference_magnitude': 9.36,
                'reference_counts': 399347.81419154577
            },
            'V7': {
                'reference_magnitude': 9.36,
                'reference_counts': 192733.27095290253
            },
            'V3': {
                'reference_magnitude': 9.36,
                'reference_counts': 44683.902504858124
            },
            'R10': {
                'reference_magnitude': 9.924,
                'reference_counts': 1
            },
            'R4': {
                'reference_magnitude': 9.924,
                'reference_counts': 1
            },
        }
        
        # Loops access this dictionary with a key of the form BandExptime (e.g. B20, V7, etc.)
        key = str(band) + str(exptime)
        
        reference_magnitude = reference_dictionary[key]['reference_magnitude']
        reference_counts = reference_dictionary[key]['reference_counts']
        
        
        # Gotta format counts into just a numpy array for the function to work
        count_array = np.array([counts['aperture_sum'].data])
        count_array = count_array.astype(float)
        count_array = count_array.flatten()
        print("Count Array:", count_array)
        
        # Get the magnitudes from the counts
        magnitudes = magnitudes_from_counts(count_array, reference_counts, reference_magnitude)
        print("Magnitudes:", magnitudes)
        print("Magnitudes Shape:", magnitudes.shape)
        
        counts['AppMag'] = magnitudes
        
        np.save(path + '../photometry/magnitudes_' + band + '_' + str(exptime) + '.npy', magnitudes)
        counts.write(path + '../photometry/counts_' + band + '_' + str(exptime) + '.tex', format='latex', overwrite=True)

    # Average out the magnitudes across each band to get band magnitudes for each star
    V20 = np.load(path + '../photometry/magnitudes_V_20.npy', allow_pickle=True)
    V7 = np.load(path + '../photometry/magnitudes_V_7.npy', allow_pickle=True)
    V3 = np.load(path + '../photometry/magnitudes_V_3.npy', allow_pickle=True)
    B20 = np.load(path + '../photometry/magnitudes_B_20.npy', allow_pickle=True)
    B60 = np.load(path + '../photometry/magnitudes_B_60.npy', allow_pickle=True)
    R10 = np.load(path + '../photometry/magnitudes_R_10.npy', allow_pickle=True)
    R4 = np.load(path + '../photometry/magnitudes_R_4.npy', allow_pickle=True)
    
    v_average = (V20 + V7 + V3) / 3
    b_average = (B20 + B60) / 2
    r_average = (R10 + R4) / 2
    
    print("V Average:", v_average)
    print("V Average Shape:", v_average.shape)
    
    new_target_table = Table.read(path + '../star_matching/target_table.tex')
    # Remove the first row (header) from the target table
    new_target_table.remove_row(0)
    
    # Update our target table with the apparent magnitudes in each band
    new_target_table.add_column(v_average, name='V_apparent')
    new_target_table.add_column(b_average, name='B_apparent')
    new_target_table.add_column(r_average, name='R_apparent')
    
    print("Updated Target Table:", new_target_table)
    
    new_target_table.write(path + '../star_matching/target_table_part2.tex', format='latex', overwrite=True)