from astropy.table import Table
import numpy as np
from scipy.stats import norm

def add_sigma_values(table):
    """
    Add sigma values to the table based on the distance of each star.
    
    Args:
        table (astropy.table.Table): Table containing star data.
        
    Returns:
        astropy.table.Table: Updated table with sigma values.
    """
    # Extract the distances from the table
    distance = table['distance'].data
    distance = np.array(distance, dtype=float)
    
    # Create a normal distribution for the distances
    center_estimate = np.nanmedian(distance)
    cluster_distance = []
    
    arbitrary_limit = 300  # Arbitrary limit to crop out the crazy outliers
    
    for i in range(len(distance)):
        if distance[i] < center_estimate + arbitrary_limit and distance[i] > center_estimate - arbitrary_limit:
            cluster_distance.append(distance[i])
    
    # Some fun with getting a gaussian fit to our distance data to see where our cluster probably is
    mu = np.nanmedian(cluster_distance)
    std = np.nanstd(cluster_distance)

    print("Mu:", mu, "Std:", std)
    
    # Create a new column for sigma values
    table['sigma'] = np.zeros(len(table), dtype=float)
    
    # Calculate sigma values based on distance
    for i in range(len(table)):
        if table['distance'][i] > 0:
            table['sigma'][i] = (table['distance'][i] - mu) / std
            print("Distance:", table['distance'][i], "Sigma:", table['sigma'][i])
        else:
            table['sigma'][i] = np.nan
    
    return table


if __name__ == "__main__":
    from config import path
    
    # Read the target table
    target_table = Table.read(path + '../star_matching/target_table_part3.tex')
    
    # Get the absolute magnitude table
    new_table = add_sigma_values(target_table)
    
    # Save the absolute magnitude table to a file
    new_table.write(path + '../star_matching/target_table_part4.tex', format='latex', overwrite=True)