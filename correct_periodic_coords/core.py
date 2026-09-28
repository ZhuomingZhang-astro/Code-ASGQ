import numpy as np
from numba import njit

@njit
def Correct_Periodic_Coords(pos, central_pos, radius, L_box):

    # Calculate parameters
    h = 0.68
    threshold = 30 * radius
    compensate_amp = L_box
    
    # Create result array (avoid modifying the original array)
    result = np.empty_like(pos)
    
    # Explicit loop processing (Numba optimizes loops extremely well)
    for i in range(len(pos)):
        for d in range(len(central_pos)):
            delta = pos[i, d] - central_pos[d]
            if np.abs(delta) > threshold:
                # Calculate sign and apply compensation
                sign = 1.0 if delta > 0 else -1.0
                result[i, d] = pos[i, d] - sign * compensate_amp
            else:
                result[i, d] = pos[i, d]
    
    # Apply compensation, return corrected coordinates
    return result
