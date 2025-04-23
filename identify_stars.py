from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import angular_separation
from config import path



def get_nearby_stars(ra, dec, radius=.1):
    """
    Get nearby stars from Simbad database based on RA and Dec.

    Args:
        ra (float): Right Ascension in degrees.
        dec (float): Declination in degrees.
        radius (float): Search radius in degrees.

    Returns:
        list: table of nearby stars with their coordinates and other info.
    """
    simbad = Simbad()
    simbad.add_votable_fields('parallax')
    result = simbad.query_region(f"{ra}d {dec}d", radius=f'{radius}d')
    
    if result is None:
        return []
    print("Simbad Search Complete. Results:", result)
    
    return result

def identify_stars(skyvals, pixel_positions, nearby_stars):
    
    """
    Identify stars in the image based on their pixel positions and sky coordinates.

    Args:
        skyvals (list): List of sky coordinates (RA, Dec) of stars.
        pixel_positions (list): List of pixel positions of stars.

    Returns:
        list: List of identified stars with their pixel and sky coordinates.
    """
    stars = Table([skyvals, pixel_positions], names=('sky_coord', 'pixel_pos'))
    
    identified_stars = []
    for i, star in enumerate(skyvals):
        # Check if the star is really close (within .001 degrees) to a known star
        print("Checking Star", i, ":", star)
        for known_star in nearby_stars:
            
            distance = star.separation(SkyCoord(ra=known_star['ra'], dec=known_star['dec'], unit=(u.deg, u.deg)))
            #print(f"Star Selected:{star}. Distance to known star {known_star} : {distance}")
            if distance < .01 * u.deg:
                # If it is, add it to the identified stars list
                print(f"Identified star: {star} matches known star: \n{known_star}")
                # Check if the known star has a parallax value
                has_parallax_value = type(known_star['plx_value']) == np.float64
                print("Has parallax value:", has_parallax_value)
                if has_parallax_value:
                    # If the known star has a parallax value, calculate the distance
                    print('known star parallax:', known_star['plx_value'], "milli arc seconds.")
                    print('parallax value type:' , type(known_star['plx_value']))
                    star_distance = 1000 / known_star['plx_value'] # Parallax to distance conversion
                    print("Star Distance:", star_distance, "parsecs")
                    identified_star = {
                    'name': known_star['main_id'],
                    'sky_coord': star,
                    'pixel_pos': pixel_positions[i],
                    'distance': star_distance,
                    }
                    identified_stars.append(identified_star)
                    print("Star Identified, added to list.")
                    break
                else:
                    # If the known star does not have a parallax value, do not add it (we don't care about it)
                    print("Known star does not have a parallax value, skipping.")
                    break

    return identified_stars


if __name__ == "__main__":
    # Example usage, this is pointed at my cluster with a reasonably large radius
    ra = 101.4991667  # Example RA in degrees
    dec = -20.7161111  # Example Dec in degrees
    radius = 0.5  # Search radius in degrees

    nearby_stars = get_nearby_stars(ra, dec, radius)
    print("Nearby Stars:", nearby_stars)

    # Assuming pixel_positions is a list of pixel coordinates of stars
    pixel_positions = np.load(path+'../star_pos_identification/pixel_positions.npy', allow_pickle=True)
    skyvals = np.load(path+'../star_pos_identification/skyvals.npy', allow_pickle=True)
    
    # Identify stars from our photo that match the nearby stars in our Simbad search
    identified_stars = identify_stars(skyvals, pixel_positions, nearby_stars)
    print("Identified Stars:", identified_stars)
    print("Number of identified stars:", len(identified_stars))
    
    star_table = Table(rows=identified_stars, names=('name','sky_coord', 'pixel_pos', 'distance'))
    print("Star Table", star_table)
    star_table.write(path+'../star_matching/target_table.tex', format='latex', overwrite=True)