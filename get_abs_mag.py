from astropy.table import Table
import numpy as np

def get_absmag(apparent_magnitude, distance):
    '''
    This function takes in the apparent magnitude and distance of a star and returns its absolute magnitude.
    
    Args:
        apparent_magnitude (float): Apparent magnitude of the star.
        distance (float): Distance to the star in parsecs.
    
    Returns:
        absolute_magnitude (float): Absolute magnitude of the star.
    '''
    return apparent_magnitude - 5 * (np.log10(distance) - 1)

def get_abs_mag_table(table):
    '''
    This function takes in a table of stars and returns a new table with the absolute magnitudes of the stars.
    
    Args:
        table (astropy.table.Table): Table containing star data.
    
    Returns:
        astropy.table.Table: New table with absolute magnitudes.
    '''
    # Create a new table to store the absolute magnitudes
    abs_mag_table = table.copy()
    
    # Remove the first row (header) from the table
    abs_mag_table.remove_row(0)
    
    # Get the apparent magnitudes and distances from the input table
    apparent_V = table['V_apparent'].data[1:]
    apparent_R = table['R_apparent'].data[1:]
    apparent_B = table['B_apparent'].data[1:]
    distances = table['distance'].data[1:]
    
    # Convert the apparent magnitudes and distances to numpy arrays
    apparent_V = np.array(apparent_V, dtype=float)
    apparent_R = np.array(apparent_R, dtype=float)
    apparent_B = np.array(apparent_B, dtype=float)
    distances = np.array(distances, dtype=float)
    
    # Calculate the absolute magnitudes
    abs_V = get_absmag(apparent_V, distances)
    abs_R = get_absmag(apparent_R, distances)
    abs_B = get_absmag(apparent_B, distances)
    
    # Add the absolute magnitudes to the new table
    abs_mag_table['abs_V'] = abs_V
    abs_mag_table['abs_R'] = abs_R
    abs_mag_table['abs_B'] = abs_B
    
    return abs_mag_table

if __name__ == "__main__":
    from config import path
    
    # Read the target table
    target_table = Table.read(path + '../star_matching/target_table_part2.tex')
    
    # Get the absolute magnitude table
    abs_mag_table = get_abs_mag_table(target_table)
    
    # Save the absolute magnitude table to a file
    abs_mag_table.write(path + '../star_matching/target_table_part3.tex', format='latex', overwrite=True)