from get_files import fits_get
from reduce_image import reduce_image
from astropy.io import fits
import os

def reduce_everything(raw_filepath, reduced_filepath, calibration_filepath, filetypes, darks = False):
    #get our raw data files of the matching bands!
    bias_files = fits_get(calibration_filepath + "/biases")
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
                fits_data = fits.open(file)
                ndata, nheader = reduce_image(fits_data, bias_files, flat_files, dark_files, ranges = None, debug = False, use_darks = darks, crop_data = False)
                
                # Save the reduced data to a new FITS file
                reduced_filename = f"{band}_{exptime}s_{n}.fits"
                print(f"Saving reduced data to {reduced_filepath + reduced_filename}")
                nheader.append(('REDUCED', True, 'Data has been reduced'))
                original_filename = os.path.basename(file)
                nheader.append(('PRE_REDU', original_filename, 'Original file name pre-reduction'))
                fits.writeto(reduced_filepath + reduced_filename, ndata, header = nheader, overwrite = True)

    
    
    
    

if __name__ == "__main__":
    filetypes = {
        "B" : [20, 60],
        "V" : [20, 7, 3],
        "R" : [10, 4],
        }
    
    path = "C:/Users/buzzs/OneDrive/Documents/Physics/Astro 230/Final_Project/data_reduction/"
    inpath = path + "raw/"
    outpath = path + "reduced/"
    calpath = path + "calibration/"
    
    reduce_everything(inpath, outpath, calpath, filetypes)