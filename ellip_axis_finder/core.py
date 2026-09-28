import numpy as np
from numpy.linalg import eigh
from scipy.stats import zscore
from sklearn.decomposition import PCA
from numba import njit
from joblib import Parallel, delayed

#----------------Method 1: PCA-----------------

@njit
def Compute_Distances_Numba(particle_pos, samp_center):
    """
    Use Numba to accelerate the calculation of Euclidean distance between particle positions and the sampling center
    
    Parameters:
        particle_pos: NumPy array of shape (N, 3), 3D coordinates of the particles
        samp_center: NumPy array of shape (3,), 3D coordinates of the sampling center
    
    Returns:
        distances: NumPy array of shape (N,), distance from each particle to the sampling center
    """
    n = particle_pos.shape[0]
    distances = np.empty(n, dtype=np.float64)
    for i in range(n):
        dx = particle_pos[i, 0] - samp_center[0]
        dy = particle_pos[i, 1] - samp_center[1]
        dz = particle_pos[i, 2] - samp_center[2]
        distances[i] = np.sqrt(dx * dx + dy * dy + dz * dz)
    return distances


def Single_Sphere_Sampling(N_samp, sphere_center, radius, samp_radius, 
                          particle_pos, nearby_gal_pos, nearby_gal_rhalf, rho_N_ave, rho_N_cri_ratio):
    """
    Single-sphere sampling and calculation of ellipsoid axes (for parallel invocation)
    
    Parameters:
        N_samp: Target number of sampling points
        sphere_center: Tuple or array of shape (3,), sphere center coordinates (x, y, z)
        radius: Sphere radius
        samp_radius: Sampling sub-sphere radius (= radius * samp_radius_ratio)
        particle_pos: Array of shape (N, 3), particle coordinates
        rho_N_ave: Average particle number density
        rho_N_cri_ratio: Density threshold ratio (actual threshold = rho_N_cri_ratio * rho_N_ave)
    
    Returns:
        axes: PCA principal components (3 axis directions), shape of (3, 3)
        variances: Variance of each principal component (representing axis length), shape of (3,)
    """
    N_count = 0
    N_loop = 0
    samp_accept = []  # Store accepted sampling points (x, y, z, particle count)
    
    while N_count < N_samp:
        # Parse sphere center coordinates
        center_x, center_y, center_z = sphere_center
        
        # Generate random sampling points on the sphere (spherical coordinate system to Cartesian coordinate system)
        r = radius * np.random.random()  # Radial random value (0~radius)
        theta = np.arccos(1 - 2 * np.random.random())  # Polar angle (0~π)
        phi = 2 * np.pi * np.random.random()  # Azimuth angle (0~2π)
        
        # Convert to Cartesian coordinates
        x = center_x + r * np.sin(theta) * np.cos(phi)
        y = center_y + r * np.sin(theta) * np.sin(phi)
        z = center_z + r * np.cos(theta)
        samp_center = np.array([x, y, z])
        
        # Calculate the distance from particles to the sampling center and filter
        distances = Compute_Distances_Numba(particle_pos, samp_center)
        samp_pos = particle_pos[distances < samp_radius]  # Filter particles inside the sub-sphere
        
        # Density judgment: accept the sampling point if density exceeds the threshold
        volume = (4 / 3) * np.pi * (samp_radius **3)
        if volume > 0 and (len(samp_pos) / volume) >= rho_N_cri_ratio * rho_N_ave:
            samp_accept.append([samp_center[0], samp_center[1], samp_center[2], len(samp_pos)])
            N_count += 1
        
        N_loop = N_loop + 1
        if N_loop > 100 * N_samp:
            return None, None
    
    # Data cleaning (remove outliers)
    samp_accept = np.vstack(samp_accept)
    
    for i in range(len(nearby_gal_pos)):
        distances = Compute_Distances_Numba(samp_accept[:, 0:3], nearby_gal_pos[i, :])
        mask = distances >= 2 * nearby_gal_rhalf[i]
        samp_accept = samp_accept[mask]
        if len(samp_accept) < int(0.7 * N_samp):
            #print('The substructure is closely attached to the central galaxy. Skipping this halo!')
            return None, None
            
    z_scores = np.abs(zscore(samp_accept[:, 0:3]))  # Calculate Z-score of coordinates
    pos_clean = samp_accept[:, 0:3][(z_scores < 3).all(axis=1)]  # Keep points with Z-score < 3
    
    # PCA analysis: compute ellipsoid axis directions and variance
    pca = PCA(n_components=3)
    pca.fit(pos_clean)
    
    return pca.components_, pca.explained_variance_
                              


def Sphere_Sampling_Method(N_samp, sphere_center, radius, samp_radius_ratio, 
                           particle_pos, nearby_gal_pos, nearby_gal_rhalf, rho_N_ave, rho_N_cri_ratio, N_re=1):
    """
    Main function: multi-process parallel sphere sampling to calculate ellipsoid axis distribution
    
    Parameters:
        N_samp: Target number of sampling points per repetition
        sphere_center: Ellipsoid center coordinates (x, y, z)
        radius: Characteristic radius of all particle distribution space
        samp_radius_ratio: Sampling sub-sphere radius ratio (samp_radius = radius * this value)
        particle_pos: (N, 3) array, particle coordinates
        rho_N_ave: Average particle number density
        rho_N_cri_ratio: Density threshold ratio, accept if rho > rho_N_ave * rho_N_cri_ratio inside sampling sub-sphere
        N_re: Number of repeated samplings (executed in parallel)
    
    Returns:
        axes_collect: List, each element is a (3, 3) array storing the three ellipsoid axis directions for each repetition
        variances_collect: List, each element is a (3, N_re) array storing the PCA variances for each repetition
    """
    # Preprocess parameters
    samp_radius = radius * samp_radius_ratio
    particle_pos = np.ascontiguousarray(particle_pos)  # Ensure memory contiguity to improve Numba efficiency
    
    # Execute N_re samplings in parallel
    results = Parallel(n_jobs=N_re, verbose=0)(  # Disable verbose to avoid log interference
        delayed(Single_Sphere_Sampling)(
            N_samp=N_samp,
            sphere_center=sphere_center,
            radius=radius,
            samp_radius=samp_radius,
            particle_pos=particle_pos,
            nearby_gal_pos=nearby_gal_pos, 
            nearby_gal_rhalf=nearby_gal_rhalf, 
            rho_N_ave=rho_N_ave,
            rho_N_cri_ratio=rho_N_cri_ratio
        ) for _ in range(N_re)
    )
    
    for axes, variances in results:
        if axes is None or variances is None:
            return [None], [None]
        
    # Collect results
    axes_collect = [res[0] for res in results]
    
    # Use the first vector as a standard so that the corresponding vectors in repeated calculations form an acute angle with it
    standard = axes_collect[0]
    adjusted_axes_collect = [standard.copy()]
    
    for i in range(1, len(axes_collect)):
        current_array = axes_collect[i].copy()
        # Check each vector
        for j in range(current_array.shape[0]):
            # Calculate the dot product of the current vector and the standard vector
            dot_product = np.dot(current_array[j], standard[j])
            # If the dot product is less than 0, reverse the direction of the current vector
            if dot_product < 0:
                current_array[j] *= -1
        adjusted_axes_collect.append(current_array)
                
    variances_collect = [res[1] for res in results]
    
    return adjusted_axes_collect, variances_collect


#----------------Method 2: ITM-----------------

@njit
def Compute_Inertia_Tensor(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf): 
    for i in range(len(nearby_gal_pos)):
        distances = Compute_Distances_Numba(particle_pos[:, 0:3], nearby_gal_pos[i, :])
        mask = distances >= 2 * nearby_gal_rhalf[i]
        particle_pos = particle_pos[mask]
        particle_mass = particle_mass[mask]
    
    particle_pos = particle_pos - sphere_center
    inertia_tensor = np.zeros((3, 3), dtype=np.float64)
    
    # Iterate through each particle
    for i in range(len(particle_mass)):
        mass = particle_mass[i]
        x, y, z = particle_pos[i]
        
        # Calculate the square of the position vector
        r_squared = x**2 + y**2 + z**2
        
        # Calculate the contribution of a single particle to the inertia tensor and accumulate it
        inertia_tensor[0, 0] += mass * (r_squared - x**2)
        inertia_tensor[0, 1] += mass * (-x * y)
        inertia_tensor[0, 2] += mass * (-x * z)
        
        inertia_tensor[1, 0] += mass * (-y * x)
        inertia_tensor[1, 1] += mass * (r_squared - y**2)
        inertia_tensor[1, 2] += mass * (-y * z)
        
        inertia_tensor[2, 0] += mass * (-z * x)
        inertia_tensor[2, 1] += mass * (-z * y)
        inertia_tensor[2, 2] += mass * (r_squared - z**2)
    
    return inertia_tensor
    

def Inertia_Tensor_Method(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf):
    inertia_tensor = Compute_Inertia_Tensor(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf)
    if inertia_tensor is None:
        return [None], [None]
    moments, axes = eigh(inertia_tensor)
    return axes.T, moments
