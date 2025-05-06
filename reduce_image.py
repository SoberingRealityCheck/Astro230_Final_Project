# Astropy Libraries
import astropy
from astropy.io import fits

# Base Python Libraries
import glob

# Other Libraries
import numpy as np
import matplotlib.pyplot as plt

'''
This file contains the function reduce_image, 
which exists to perform astronomical data reduction 
on a .fits image file given a collection of calibration frames.

The three types of calibration frames we use are:
    1. Bias frames
        - These are images taken with the lens cap on and a short exposure time.
        - Since they're fast, the only thing measured is the readout noise from
        the camera sensor itself.
    
    2. Flat frames
        - These are images taken with a uniform light source.
        - These can be taken by pointing the telescope at an evenly lit 
        section of the wall or something along those lines.
        - By averaging them, we can determine how the effects of the optics 
        change the way each pixel responds on average and compensate for stuff
        like changes in brightness around the edges.
    
    3. Dark frames 
        - These are images taken with the lens cap on and a long exposure time.
        - Since thermal noise is read, by averaging them we can
        determine the thermal background.
        
We take several of each kind of calibration frame before we start imaging!

This program will utilize all three kinds of frames to remove 
thermal noise, readout noise, dead pixels and optical artifacts from our image.
'''


def reduce_image(ffile, biases, flats, darks, ranges, debug=False, use_darks=True, crop_data = True):
    '''
    Applies the image reduction process to a fits file by subtracting the bias, dividing by the flat, and subtracting the dark images provided.
    
    ffile: fits file object
    
    test_biases: glob of filepaths to bias files
    test_flats: glob of filepaths to flat files
    test_darks: glob of filepaths to dark files
    
    ranges = dict containing:
        x1: first good pixel in x after trim
        x2: last pixel to include in x after trim

        y1: first good pixel in y after trim
        y2: last pixel to include in y after trim

        bx1 = first pixel of bias section
        bx2 = last pixel of bias section
    
    
    returns a tuple of reduced data, header
    '''
    
    #pull variables from ranges 
    if crop_data:
        x1 = ranges['x1']
        x2 = ranges['x2']
        y1 = ranges['y1']
        y2 = ranges['y2']
        bx1 = ranges['bx1']
        bx2 = ranges['bx2']
    
    # Get the data and header
    if crop_data:
        data = ffile[0].data[y1:y2,x1:x2]
    else:
        data = ffile[0].data

    fheader = ffile[0].header
    if debug:
        print(f"data imported. Data shape: {data.shape}")
        #plt.imshow(data, origin='lower', vmin=1500, vmax=6500)

    # Get the median of all the bias frames 
    bias_timestream = [fits.getdata(x) for x in biases]
    median_bias = np.median(bias_timestream, axis=0)
    # Subtract the median bias from the image data
    if crop_data:
        ndata_no_bias = data - median_bias[y1:y2,x1:x2]
    else:
        ndata_no_bias = data - median_bias
    
    
    if debug:
        print(f"bias subtracted. Median bias: {np.median(median_bias)}")
        print("data median: ", np.median(ndata_no_bias))
        #plt.imshow(ndata_no_bias, origin='lower', vmin=0, vmax=5000)
    del bias_timestream

    if use_darks:
        # Get the median of all the dark frames
        dark_timestream = [fits.getdata(x) for x in darks]
        mean_dark = np.mean(dark_timestream, axis=0)
        median_dark = np.median(dark_timestream, axis=0)
        stddev_dark = np.std(dark_timestream, axis=0)
        
        # Subtract the median bias from the dark
        diff = median_dark - median_bias
        
        # Subtract the dark from the image data
        if crop_data:
            ndata_no_bias = ndata_no_bias - diff[y1:y2,x1:x2]
        else:
            ndata_no_bias = ndata_no_bias - diff
        if debug:
            print(f"darks subtracted. Median dark: {np.median(median_dark)}")
            print("data median: ", np.median(ndata_no_bias))
            #plt.imshow(ndata_ff_nodark, origin='lower', vmin=1500, vmax=6500)
        del dark_timestream


    # Get the median of all the flat frames
    flat_timestream = [fits.getdata(x)for x in flats]
    mean_flat = np.mean(flat_timestream, axis=0)
    median_flat = np.median(flat_timestream, axis=0)
    
    # Subtract the median bias from the median flat
    bias_subtracted_flat = median_flat - median_bias
    
    # Normalize the flat so the pixels have a mean value of 1
    flat_mean = np.mean(bias_subtracted_flat)
    normalized_mean_flat = bias_subtracted_flat / flat_mean
    if crop_data:
        trimmed_flat = normalized_mean_flat[y1:y2,x1:x2]
    else:
        trimmed_flat = normalized_mean_flat
    
    # Divide image data by the flat
    ndata_ff = ndata_no_bias / trimmed_flat
    if debug:
        print(f"flats divided. Median flat: {np.median(median_flat)}")
        print("data median: ", np.median(ndata_ff))
        #plt.imshow(ndata_ff, origin='lower', vmin=1500, vmax=6500)
    del flat_timestream

    # Update the header
    nheader = fheader # copy the original header
    if crop_data:
        nheader['BIASSEC'] = "["+str(bx1)+":"+str(bx2)+","+str(y1)+":"+str(y2)+"]"
        nheader['TRIMSEC'] = "["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]"
    #return the reduced data and header
    return ndata_ff, nheader

if __name__ == "__main__":
    
    from config import path
    # Test the function
    test_ffile = fits.open(path + "raw/obj0016.fits")
    test_biases = glob.glob(path + 'calibration/biases/*.fits')
    test_flats = glob.glob(path + 'calibration/flats/*.fits')
    test_darks = glob.glob(path + 'calibration/darks/*.fits')

    DEBUG = True
    # Set to True to see debug information and plots
    
    DARKS = False
    # Set to True to use dark frames in the reduction process
    
    CROP_DATA = False
    # Set to True to crop the data to the specified ranges

    from config import test_ranges
    
    ndata, nheader = reduce_image(test_ffile, test_biases, test_flats, test_darks, test_ranges, debug=DEBUG, use_darks = DARKS, crop_data = CROP_DATA)
    if DEBUG:
        plt.imshow(ndata, origin='lower', vmin=1500, vmax=5500)
        plt.colorbar()
        plt.show()

    
    # If you want to save the image as a new file, uncomment the lines below
    
    #filename = "mars_1.fits"
    #fits.writeto(path+"data/reduced/"+filename, data=ndata, header=nheader, overwrite=True)

