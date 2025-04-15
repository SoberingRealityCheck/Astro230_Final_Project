from get_files import fits_get
from reduce_image import reduce_image
from astropy.io import fits
import os
from code.config import crop_ranges

def reduce_everything(raw_filepath, reduced_filepath, calibration_filepath, filetypes, darks = False):
    #get our raw data files of the matching bands!
    bias_files = fits_get(calibration_filepath + "/biases")
    
    #get dark files if our name is Ryan and we actually care about thermal noise
    if darks:
        dark_files = fits_get(calibration_filepath + "/darks")
    else:
        dark_files = []
    
    for band in filetypes.keys():
        exptimes = filetypes[band]
        flat_files = fits_get(calibration_filepath + "/flats", selected_band = band)
        
        #get the files for each band and exposure time
        for exptime in exptimes:
            raw_files = fits_get(raw_filepath, band, str(exptime))
            print(f"Raw files for {band} with exposure time {exptime}: {raw_files}")
            
            for n, file in enumerate(raw_files):
                
                # Check if the file exists before processing (this should be impossible unless you delete stuff while the program is running?)
                if not os.path.isfile(file):
                    print(f"File {file} does not exist. Skipping.")
                    continue
                
                # Open the FITS file and reduce the data
                fits_data = fits.open(file)
                ndata, nheader = reduce_image(fits_data, bias_files, flat_files, dark_files, ranges = crop_ranges, debug = False, use_darks = darks, crop_data = True)
                
                # Add a new header entry to indicate that the data has been reduced
                nheader.append(('REDUCED', True, 'Data has been reduced'))
                
                #Add the original filename to the header for reference
                original_filename = os.path.basename(file)
                nheader.append(('PRE_REDU', original_filename, 'Original file name pre-reduction'))
                
                # Save the reduced data to a new FITS file
                reduced_filename = f"{band}_{exptime}s_{n}.fits"
                print(f"Saving reduced data to {reduced_filepath + reduced_filename}")
                fits.writeto(reduced_filepath + reduced_filename, ndata, header = nheader, overwrite = True)

    
    
    
    

if __name__ == "__main__":
    # the filetypes dictionary is formatted as bands as keys 
    # and a list of related exposure times as an array of values.
    from config import filetypes 
    
    # Set the paths for the raw, reduced, and calibration data. 
    # Change these paths to your own directory structure.
    from config import path

    # This is where the raw data is stored.
    inpath = path + "raw/"
    
    # This is where the reduced data will be saved.
    outpath = path + "reduced/"
    
    # Make sure this calibration directory has subdirectories for biases, flats, and darks that are named correctly.
    calpath = path + "calibration/"
    
    # This will make a bunch of new reduced FITS files in the reduced directory with new names.
    reduce_everything(inpath, outpath, calpath, filetypes)