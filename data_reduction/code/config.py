# Configuration file for the data reduction script
# This file contains the paths and settings for the data reduction process.

# This is the local path to the directory.
path = "C:/Users/buzzs/OneDrive/Documents/Physics/Astro 230/Final_Project/data_reduction/"

# This contains info on the types of files I've captured in terms of band types and exposure times (in seconds).
# The filetypes dictionary is formatted as bands as keys
# and a list of related exposure times as an array of values.
filetypes = {
        "B" : [20, 60],
        "V" : [20, 7, 3],
        "R" : [10, 4],
        }

#McDonald 30inch telecope data ranges I use to trim off the bias strip and stuff.
crop_ranges = {
'x1' : 3,     # first good pixel in x after trim
'x2' : 1020,  # last pixel to include in x after trim
'y1' : 2,     # first good pixel in y after trim
'y2' : 1020,  # last pixel to include in y after trim
'bx1' : 1025, # first pixel of bias section
'bx2' : 1054, # last pixel of bias section}
}