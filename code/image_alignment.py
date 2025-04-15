# basic imports
import numpy as np
import matplotlib.pyplot as plt
import os

# spacey imports
from astropy.io import fits
from scipy.ndimage import shift

# stuff from other python documents
from config import path, filetypes
from get_files import fits_get


'''
We will be using the astroalign library. Because I am lazy and don't want to manually do this anymore. 
# let's test and see how well a pre-written alignment package handles this! lol. 

Update: astroalign sucks. Can't get it to work. Just going to match up the brightest pixel in each image, lol.
'''

def align_and_combine(filepath, output_path):
    """
    Aligns and combines a list of FITS files using astroalign.

    Parameters:
        files (list): List of file paths to the FITS files to be aligned and combined.
        output_path (str): Path where the combined FITS files will be saved.
    """
    for band_selected in filetypes.keys():
        exptimes = filetypes[band_selected]
        for exptime in exptimes:
            print("Combining Images in the", band_selected, "band with exposure time", exptime, "seconds")
            # Get the list of files for the current band and exposure time
            relevant_files = fits_get(filepath, band_selected, str(exptime))
            
            # Initialize an empty list to store the aligned images
            aligned_images = []
            
            # Align each file with the first file in the list
            reference_image, reference_header = fits.getdata(relevant_files[0], header=True)
            
            reference_image = reference_image - np.median(reference_image)  # Subtract the background
            
            reference_star_value = np.max(reference_image)
            reference_star_location = np.squeeze(np.where(reference_image == reference_star_value))
            
            #to make sure it's all working correctly, let's also track where each 'reference star' location is for each image
            reference_star_positions = []
            
            print("Reference Star location:", reference_star_location)
            
            for file in relevant_files:
                image_data, header = fits.getdata(file, header=True)
                
                #get the brightest star in the image to use as a reference point for alignment
                star_value = np.max(image_data)
                star_location = np.squeeze(np.where(image_data == star_value))
                print("Star Location:", star_location)
                
                #determine the required shift to align the images
                transform = (reference_star_location[0] - star_location[0], reference_star_location[1] - star_location[1])
                print("Transform:", transform)
                
                #subtract the estimated background from the image data
                background = np.median(image_data)
                image_data = image_data - background
                
                ## Align the image using the calculated transform
                aligned_image = shift(image_data, transform, order=1, mode='constant', cval=0.0, prefilter=False)
                aligned_images.append(aligned_image)
                reference_star_positions.append(star_location - transform)
                
            # Combine the aligned images by averaging them
            # comparing and contrasting mean vs median for image averaging
            fig, (ax1, ax2) = plt.subplots(1, 2)
            mean_images = np.mean(aligned_images, axis=0)
            im1 = ax1.imshow(mean_images, cmap='gray', origin='lower', vmin=0, vmax=100)
            median_images = np.median(aligned_images, axis=0)
            im2 = ax2.imshow(median_images, cmap ='gray', origin='lower', vmin=0, vmax=100)
            ax1.set_title('Mean Aligned Image')
            ax2.set_title('Median Aligned Image')
            fig.suptitle(f"Band: {band_selected}, Exposure Time: {exptime} seconds")
            ax1.plot(reference_star_location[1], reference_star_location[0], 'x', markersize=2, label='Reference Star')
            ax2.plot(reference_star_location[1], reference_star_location[0], 'x', markersize=2, label='Reference Star')
            ax1.scatter(reference_star_positions[0][1], reference_star_positions[0][0], c='b', s=2, label='Aligned Star')
            ax2.scatter(reference_star_positions[0][1], reference_star_positions[0][0], c='b', s=2, label='Aligned Star')
            fig.legend(loc='upper right')
            plt.show()
            
            #looks like the median gives better results.
            combined_image = median_images
            
            # Save the combined image to a new FITS file
            combined_filename = f"{band_selected}_{exptime}s_combined.fits"
            print(f"Saving combined image to {output_path + combined_filename}")
            fits.writeto(os.path.join(output_path, combined_filename), combined_image, header=reference_header, overwrite=True)

# Example usage of the align_and_combine function
if __name__ == "__main__":
    # Get the list of files to be aligned and combined
    reduced_filepath = path + "reduced/"
    output_path = path + "combined/"
    
    # Align and combine the files
    align_and_combine(reduced_filepath, output_path)