#convert pixel locations to ra/dec locations using some known reference point and the fov of the sensor

def get_star_locations(star_positions, reference_points, fov):
    """
    Convert pixel locations to RA/Dec locations using a known reference point and the FOV of the sensor.

    Args:
        star_positions (list): List of pixel positions of stars.
        reference_point (tuple): Reference point (RA, Dec) in degrees.
        fov (float): Field of view in degrees.

    Returns:
        list: List of RA/Dec locations of stars.
    """
    
    # Step 1: Unpack the reference points (these need to be Known Standard Stars in our image that we know both the pixel x-y location and the RA/Dec location of)
    
    
    # Step 2: Calculate the initial pixel scale based on the FOV and the image size
    
    
    # Step 3: Find the average offsets of our linear model to get the most accurate pixel to WCS conversion (since there might be some inaccuracy introduced by our low resolution images)
    
    
    # Step 4: Convert pixel positions to RA/Dec using the reference points and the pixel scale
    
    
    # Step 5: Return a numpy array of RA/Dec locations of stars
    return []