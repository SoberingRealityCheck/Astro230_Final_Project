from config import path
from astropy.table import Table 
import numpy as np

star_table = Table.read(path + '../star_matching/target_table.tex')
readable_table = Table()
readable_table['name'] = star_table['name'][1:]

# Shoutout to the person on stackoverflow who provided this function to truncate floats to a certain number of decimal places
# https://stackoverflow.com/questions/8595973/truncate-to-three-decimals-in-python
import math
def truncate(number, digits) -> float:
    # Improve accuracy with floating point operations, to avoid truncate(16.4, 2) = 16.39 or truncate(-1.13, 2) = -1.12
    nbDecimals = len(str(number).split('.')[1]) 
    if nbDecimals <= digits:
        return number
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

sky = star_table['sky_coord'].data[1:]
sky = [x.strip().split(',') for x in sky]
sky = np.array(sky, dtype=float)
ra = sky[:, 0]  # Right Ascension in degrees
dec = sky[:, 1]  # Declination in degrees

distance = star_table['distance'].data[1:]
distance = np.array(distance, dtype=float)

new_ra = []
new_dec = []
new_dist = []

for i in range(len(readable_table)):
    new_ra.append(truncate(ra[i], 2))
    new_dec.append(truncate(dec[i], 2))
    new_dist.append(truncate(distance[i], 2))

print("readable_table", readable_table)
print("new_ra", len(new_ra))
readable_table['ra'] = new_ra
readable_table['dec'] = new_dec
readable_table['distance'] = new_dist
print("readable_table", readable_table)

readable_table.write(path + '../star_matching/readable_table.tex', overwrite=True)