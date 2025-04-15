import os
import glob
from astropy.io import fits


def fits_get(path, selected_band = "Not Specified", selected_exptime = "Not Specified"):
    """
    Get the list of files in the specified path that match the given band and exposure time.

    Parameters:
        path (str): The directory to search for files.
        band (str) (optional): The band to match in the file names.
        exptime (int or str) (optional): The exposure time in seconds to match in the file names.

    Returns:
        list: A list of file paths that match the criteria.
    """
    # Print the initial information
    print("")
    print("=============================")
    print("Initializing FITS File Search")
    print("=============================")
    print("")
    print("Searching for files in:", path)
    print("Band Selected:", selected_band)
    print("Exposure Time Selected:", selected_exptime, "seconds")
    
    # Check if the path exists
    if not os.path.exists(path):
        raise FileNotFoundError(f"\n Hey dumbass!\n The specified path does not exist: {path} \n Check your filepath again for typos.")

    # Create a glob pattern to match files with the specified band and exposure time
    pattern = os.path.join(path, f"*.fit*")
    # Use glob to find all matching files
    all_files = glob.glob(pattern)
    
    print("")
    print(f"Found {len(all_files)} files in the directory, searching for matches...")
    print("")
    
    files = []
    for file in all_files:
        file_hdr = fits.getheader(file)
        
        if selected_band != "Not Specified":
            file_band = file_hdr['INSFILTE']
        else:
            file_band = "Not Specified"
        if selected_exptime != "Not Specified":
            file_exptime = str(int(file_hdr['EXPTIME']))
        else:
            file_exptime = "Not Specified"
        
        #this is just like the original filename
        file_number = os.path.basename(file).split('_')[0]
        
        if file_exptime == selected_exptime and file_band == selected_band:
            files.append(file)
            print(f"Matching File Found: {file_number}")
    print("")
    print("Search Complete! Found", len(files), "matching files.")
    print("")
    return files


if __name__ == "__main__":
    # Example usage from my personal computer. This will not work on your computer :)
    from config import path
    raw_path = path + "raw/"
    band = "R"
    exptime = "10"
    
    matching_files = fits_get(raw_path, selected_band = band, selected_exptime = exptime)