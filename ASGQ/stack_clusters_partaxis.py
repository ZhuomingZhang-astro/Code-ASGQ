# stack_clusers.py

import numpy as np
import h5py
from sklearn.decomposition import PCA
import yt
import sys
sys.path.append("/home/zhangzhuoming/Quench_Direction") 
from ellip_axis_finder import Compute_Distances_Numba
from ellip_axis_finder import Sphere_Sampling_Method
from ellip_axis_finder import Inertia_Tensor_Method
from correct_periodic_coords import Correct_Periodic_Coords

def process_halo_data(simulation_name, h, catalog_path, snapshot_path, L, snapnum, star_num_limit=100, halo_mass_limit=1e12, include_cgal=False, part_axis='star', rho_N_cri_ratio=3, N_re=1):
    
    yt.set_log_level("warning")
    yt.set_log_level("error")
    # Iterate through each file and read particle data
    ds = yt.load(snapshot_path)
    # Get all particle data for the entire dataset
    ad = ds.all_data()
    # particle position unit: SIMBA[kpc/h], TNG100[kpc/h], EAGLE[Mpc/h], SIMBAnofb[kpc/h]
    star_pos = ad["PartType4", "Coordinates"]
    if part_axis == 'dm':
        dm_pos = ad["PartType1", "Coordinates"]
    
    file_path_catalog = catalog_path
    with h5py.File(file_path_catalog, 'r') as f:        
        # Input data
        GroupID = f['galaxy_data/GroupID'][:]
        parent_halo_index = f['galaxy_data/parent_halo_index'][:]
        central = f['galaxy_data/central'][:]
        pos = f['galaxy_data/pos'][:] # kpc
        vel = f['galaxy_data/vel'][:]
        sfr = f['galaxy_data/sfr'][()]
        slist_start = f['galaxy_data/slist_start'][()]
        slist_end = f['galaxy_data/slist_end'][()]
        slist = f['galaxy_data/lists/slist'][()]
        dicts = f['galaxy_data/dicts']
        star_mass = dicts['masses.stellar'][:]
        dm_mass = dicts['masses.dm_30kpc'][:]
        gas_mass = dicts['masses.gas'][:]
        radii_half = dicts['radii.total_half_mass'][:] # kpc
        
        halo_pos = f['halo_data/pos'][:] # kpc
        halo_dicts = f['halo_data/dicts']
        halo_radius = halo_dicts['virial_quantities.r200c'][:] # kpc
        halo_mass = halo_dicts['virial_quantities.m200c'][:]
        # dm particle
        dmlist_start = f['halo_data/dmlist_start'][()]
        dmlist_end = f['halo_data/dmlist_end'][()]
        dmlist = f['halo_data/lists/dmlist'][()]
        
    # Generate the resident array, where resident[i] provides the catalog IDs of galaxies in the halo with halo_index=i
    resident = np.empty(np.max(parent_halo_index)+1, dtype=object)

    # Initialize each element as an empty list
    for i in range(len(resident)):
        resident[i] = []

    # Write galaxies into each halo
    for i in range(len(parent_halo_index)):
        if star_mass[i] >= 10**(8):
            resident[parent_halo_index[i]].append(i)

    cluster_halo_CID = []  # Used to store the catalog indices of halos considered to be clusters

    for i in range(len(resident)):
        if len(resident[i]) >= 1:
            cluster_halo_CID.append(i)

    def central_CID_finder(index_of_XX_halo):
        N_central = 0
        galaxy_index = 0
        for i in resident[index_of_XX_halo]:
            if (central[i] == True):
                galaxy_index = i
                N_central = N_central + 1
        accept = True
        if (N_central > 1):
            accept = False
        if (N_central == 0):
            accept = False
        return galaxy_index, accept

    def central_rCID_finder(index_of_XX_halo):
        N_central = 0
        galaxy_index = 0
        for i in range(len(resident[index_of_XX_halo])):
            if (central[resident[index_of_XX_halo][i]] == True):
                galaxy_index = i
                N_central = N_central + 1
        accept = True
        if (N_central > 1):
            accept = False
        if (N_central == 0):
            accept = False
        return galaxy_index, accept


    sgaldata = [[] for _ in range(N_re+1)]
    CID_vec_var = []

    for index in range(len(cluster_halo_CID)):
        # Galaxy cluster/group index
        # Read the stellar coordinates of the central galaxy    
        gal_index = np.array(resident[cluster_halo_CID[index]])  # List of CIDs of all member galaxies in the current halo
        central_rCID, accept = central_rCID_finder(cluster_halo_CID[index])  # rCID of the central galaxy in the current halo
        if not accept:
            continue
        central_CID, accept = central_CID_finder(cluster_halo_CID[index])  # CID of the central galaxy in the current halo
        cgal_pos = pos[central_CID]
        R200c = halo_radius[cluster_halo_CID[index]]  # Distance threshold [kpc]
        M200c = halo_mass[cluster_halo_CID[index]]
        
        # Read all stellar coordinates of the central galaxy from the snapshot
        star_cgal_SID = slist[slist_start[central_CID]:slist_end[central_CID]]
        if (R200c != 0) and (M200c >= halo_mass_limit) and (len(star_cgal_SID) > star_num_limit):
            mask = (GroupID != central_CID) # & (star_mass >= 1e8)
            CID_other = GroupID[mask]
            pos_other = pos[mask, :]
            radii_half_other = radii_half[mask]
            radius_cluster = 10 * R200c
            
            for i in range(len(cgal_pos)):
                offset = np.zeros(pos_other.shape[1])
                if cgal_pos[i] < radius_cluster:
                    offset[i] = L
                    CID_other = np.concatenate([CID_other, CID_other])
                    pos_other = np.concatenate([pos_other, pos_other - offset])
                    radii_half_other = np.concatenate([radii_half_other, radii_half_other])
                if cgal_pos[i] > L - radius_cluster:
                    offset[i] = L
                    CID_other = np.concatenate([CID_other, CID_other])
                    pos_other = np.concatenate([pos_other, pos_other + offset])
                    radii_half_other = np.concatenate([radii_half_other, radii_half_other])

            # Calculate the Euclidean distance from each point in the neighboring galaxies to the center point
            distances = Compute_Distances_Numba(pos_other, cgal_pos)
            # Filter out points with distances less than the threshold
            nearby_gal_CID = CID_other[distances < radius_cluster]
            nearby_gal_pos = pos_other[distances < radius_cluster]
            nearby_gal_rhalf = radii_half_other[distances < radius_cluster]
            
            
            if part_axis == 'star':
                # star particles' positions
                star_cgal_pos_unscale = star_pos[star_cgal_SID]
                star_cgal_pos_unscale = np.column_stack((star_cgal_pos_unscale[:,0].astype(np.float64), star_cgal_pos_unscale[:,1].astype(np.float64), star_cgal_pos_unscale[:,2].astype(np.float64)))
                
                if simulation_name == 'SIMBA':
                    star_cgal_pos = star_cgal_pos_unscale / h # kpc/h -> kpc
                if simulation_name == 'TNG100':
                    star_cgal_pos = star_cgal_pos_unscale / h # kpc/h -> kpc
                if simulation_name == 'EAGLE':
                    star_cgal_pos = star_cgal_pos_unscale * 1000 / h # Mpc/h -> kpc
                if simulation_name == 'SIMBAnofb':
                    star_cgal_pos = star_cgal_pos_unscale / h # kpc/h -> kpc
                    
                radius_star = 2 * radii_half[central_CID]
                star_cgal_pos = Correct_Periodic_Coords(pos=star_cgal_pos, central_pos=cgal_pos, radius=radius_star, L_box=L)
                star_cgal_pos_samp = star_cgal_pos
                star_rho_N_ave = len(star_cgal_pos_samp) / (4 / 3 * np.pi * radius_star ** 3)
    
                star_axes_collect, star_variances_collect = Sphere_Sampling_Method(N_samp=5000, sphere_center=cgal_pos, radius=radius_star, samp_radius_ratio=1/20,
                                                           particle_pos=star_cgal_pos_samp, nearby_gal_pos=nearby_gal_pos, nearby_gal_rhalf=nearby_gal_rhalf,
                                                           rho_N_ave=star_rho_N_ave, rho_N_cri_ratio=rho_N_cri_ratio, N_re=N_re)
                axes_collect = star_axes_collect
                variances_collect = star_variances_collect
                    

            if part_axis == 'dm':
                # dark matter particles' positions
                dm_cgal_SID = dmlist[dmlist_start[cluster_halo_CID[index]]:dmlist_end[cluster_halo_CID[index]]]
                dm_cgal_pos_unscale = dm_pos[dm_cgal_SID] # Mpc / h
                dm_cgal_pos_unscale = np.column_stack((dm_cgal_pos_unscale[:,0].astype(np.float64), dm_cgal_pos_unscale[:,1].astype(np.float64), dm_cgal_pos_unscale[:,2].astype(np.float64)))
                
                if simulation_name == 'SIMBA':
                    dm_cgal_pos = dm_cgal_pos_unscale / h # kpc/h -> kpc
                if simulation_name == 'TNG100':
                    dm_cgal_pos = dm_cgal_pos_unscale / h # kpc/h -> kpc
                if simulation_name == 'EAGLE':
                    dm_cgal_pos = dm_cgal_pos_unscale * 1000 / h # Mpc/h -> kpc
                if simulation_name == 'SIMBAnofb':
                    dm_cgal_pos = dm_cgal_pos_unscale / h # kpc/h -> kpc
                
                radius_dm = R200c
                dm_cgal_pos = Correct_Periodic_Coords(pos=dm_cgal_pos, central_pos=cgal_pos, radius=radius_dm, L_box=L)
                distances = np.sqrt(np.sum((dm_cgal_pos - cgal_pos) ** 2, axis=1))
                dm_cgal_pos = dm_cgal_pos[distances < radius_dm]

                if len(dm_cgal_pos) > 100000:
                    dm_cgal_pos_samp = dm_cgal_pos[np.random.choice(len(dm_cgal_pos), 100000, replace=False)]
                if len(dm_cgal_pos) <= 100000:
                    dm_cgal_pos_samp = dm_cgal_pos
                dm_rho_N_ave = len(dm_cgal_pos_samp) / (4 / 3 * np.pi * radius_dm ** 3)
                dm_axes_collect, dm_variances_collect = Sphere_Sampling_Method(N_samp=5000, sphere_center=cgal_pos, radius=radius_dm, samp_radius_ratio=1/20,
                                                           particle_pos=dm_cgal_pos_samp, nearby_gal_pos=nearby_gal_pos, nearby_gal_rhalf=nearby_gal_rhalf,
                                                           rho_N_ave=dm_rho_N_ave, rho_N_cri_ratio=rho_N_cri_ratio, N_re=N_re)
                axes_collect = dm_axes_collect
                variances_collect = dm_variances_collect
                    
                
            if axes_collect[0] is None:
                continue
                    
            # Average over repeated calculations (axis=0 indicates along the "repetition" dimension)
            ave_vectors = np.array(axes_collect).mean(axis=0)
            ave_variances = np.array(variances_collect).mean(axis=0)
            
            CID_vec_var.append([central_CID, axes_collect, variances_collect])
            
            for N_i in range(len(sgaldata)):
                if (N_i == 0):
                    # Extract the average values of the three directions
                    vec1, vec2, vec3 = ave_vectors
                if (N_i != 0):
                    # Extract the average values of the three directions
                    vec1, vec2, vec3 = axes_collect[N_i-1]
                    
                # Rotate the coordinate system using the three axes of the ellipsoid (from longest to shortest) as the new x, y, z axes
                # Construct the rotation matrix
                # Ensure that vec1, vec2, and vec3 are unit vectors
                vec1 = vec1 / np.linalg.norm(vec1)
                vec2 = vec2 / np.linalg.norm(vec2)
                vec3 = vec3 / np.linalg.norm(vec3)
                # Construct the rotation matrix
                rotation_matrix = np.vstack([vec1, vec2, vec3]).T

                if (include_cgal == True):
                    gal_index_select = np.append(nearby_gal_CID, central_CID)
                    sgal_abspos = np.vstack((nearby_gal_pos, cgal_pos))
                if (include_cgal == False):
                    gal_index_select = nearby_gal_CID 
                    sgal_abspos = nearby_gal_pos

                    
                cgal_abspos = np.column_stack((np.repeat(cgal_pos[0], len(gal_index_select)), 
                               np.repeat(cgal_pos[1], len(gal_index_select)), 
                               np.repeat(cgal_pos[2], len(gal_index_select))))
                # Relative coordinates
                relative_pos = sgal_abspos - cgal_pos
                # Radial distance r
                r = np.linalg.norm(relative_pos, axis=1)
                mask = (r > 2 * radii_half[central_CID]) | (gal_index_select == central_CID) # Keep only satellite galaxies with r > 2 * r_half
                gal_index_select = gal_index_select[mask]
                r = r[mask]
                sgal_abspos = sgal_abspos[mask]
                cgal_abspos = cgal_abspos[mask]
                relative_pos = relative_pos[mask]
                
                sgal_absvel = vel[gal_index_select]
                cgal_absvel = np.column_stack((np.repeat(vel[central_CID][0], len(gal_index_select)), 
                               np.repeat(vel[central_CID][1], len(gal_index_select)), 
                               np.repeat(vel[central_CID][2], len(gal_index_select))))
                
                # Apply the rotation matrix
                rotated_pos = np.dot(relative_pos, rotation_matrix)

                relative_vel = vel[gal_index_select] - vel[central_CID]
                rotated_vel = np.dot(relative_vel, rotation_matrix)

                # Calculate spherical coordinates
                ra_deg = np.zeros_like(r)
                dec_deg = np.zeros_like(r)
                for i in range(len(r)): 
                    if (r[i] != 0):
                        # Right ascension RA (radians)
                        ra_i = np.arctan2(rotated_pos[i, 1], rotated_pos[i, 0])
                        # Declination DEC (radians)
                        dec_i = np.arcsin(rotated_pos[i, 2] / r[i])
                        # Convert RA and DEC from radians to degrees
                        ra_deg[i] = np.mod(np.degrees(ra_i), 360)
                        dec_deg[i] = np.degrees(dec_i)

                    if (r[i] == 0):
                        # Right ascension RA (radians)
                        ra_deg[i] = np.nan
                        # Declination DEC (radians)
                        dec_deg[i] = np.nan

                # Record the star formation rate
                sgal_sfr = sfr[gal_index_select]

                # Record the stellar mass of each satellite galaxy
                sgal_dmmass = dm_mass[gal_index_select]
                sgal_gasmass = gas_mass[gal_index_select]
                sgal_starmass = star_mass[gal_index_select]
                sgal_radiihalf = radii_half[gal_index_select]

                sgal_mass = sgal_dmmass + sgal_gasmass + sgal_starmass

                # Record the Group ID of the central galaxy in the current cluster
                central_CID_array = np.full(len(gal_index_select), central_CID)

                # Record the radius of the current cluster
                Radius = np.full(len(gal_index_select), R200c)
                Mass = np.full(len(gal_index_select), M200c)

                # Others
                progen_CID = np.full(len(gal_index_select), np.nan)
                progen_cgalCID = np.full(len(gal_index_select), np.nan)
                descen_CID_select = np.full(len(gal_index_select), np.nan)
                descen_cgalCID_select  = np.full(len(gal_index_select), np.nan)
                merge = np.full(len(gal_index_select), np.nan)


                #(x, y, z, r, RA, DEC sfr, starmass, CatalogueID, central_galaxy_CatalogueID, progen_CatalogueID, progen_central_galaxy_CatalogueID, descen_CatalogueID, descen_central_galaxy_CatalogueID , vx, vy, vz, R200, Rgal_half, sgal_abspos(x, y, z), cgal_abspos(x, y, z), sgal_absvel(x, y, z), cgal_absvel(x, y, z), gal_total_mass, Halo_M200c, whether_merge, gal_dmmass, gal_gasmass)
                sgaldata[N_i].append([rotated_pos[:, 0], rotated_pos[:, 1], rotated_pos[:, 2], r, ra_deg, dec_deg, sgal_sfr, sgal_starmass, gal_index_select, central_CID_array, progen_CID, progen_cgalCID, descen_CID_select, descen_cgalCID_select, rotated_vel[:, 0], rotated_vel[:, 1], rotated_vel[:, 2], Radius, sgal_radiihalf, sgal_abspos[:, 0], sgal_abspos[:, 1], sgal_abspos[:, 2],
    cgal_abspos[:, 0], cgal_abspos[:, 1], cgal_abspos[:, 2], sgal_absvel[:, 0], sgal_absvel[:, 1], sgal_absvel[:, 2], cgal_absvel[:, 0], cgal_absvel[:, 1], cgal_absvel[:, 2], sgal_mass, Mass, merge, sgal_dmmass, sgal_gasmass])
    
    return sgaldata, CID_vec_var
