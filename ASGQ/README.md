# ASGQ Core Pipeline

Core statistical pipeline for primary data processing, statistical analysis, and figure generation for the project.

---

## 📂 Code Scripts & Notebooks

1. [`stack_clusters_partaxis.py`](stack_clusters_partaxis.py) Note: Utility script (imported and used as a module).
   * **Description:** Computes ellipsoidal triaxiality based on "central galaxy stellar particle data" / "host halo dark matter particles within $R_\mathrm{200c}$" across $N_\mathrm{re}$ independent repeated runs. Within `stack_clusters_partaxis.py`, `galaxy_ave.npy` uses the triaxial vectors averaged over the $N_\mathrm{re} = 50$ independent repeated calculations (`ave_vectors = np.array(axes_collect).mean(axis=0)`), whereas each individual run file (`galaxy_{i}.npy`) utilizes the actual triaxial directions independently computed in that specific run.
   * **Output Data Structure:**
     * **1. Galaxy and Environment Data (`sgaldata`):** Records information for every central galaxy and all surrounding galaxies within $10 R_\mathrm{200c}$ and beyond 2 times half mass radius of central galaxy of its position across independent runs.
       
| Column / Index | Variable Name in Script | Physical Meaning & Description |
| :--- | :--- | :--- |
| **0 – 2** | `rotated_pos[:, 0:3]` | Coordinates $(x, y, z)$ relative to the central galaxy's reference frame, aligned with its major, intermediate, and minor axes (kpc) |
| **3** | `r` | Radial distance ($r$) (kpc) |
| **4** | `ra_deg` | Right Ascension (RA) relative to the central galaxy reference frame (deg) |
| **5** | `dec_deg` | Declination (DEC) relative to the central galaxy reference frame (deg) |
| **6** | `sgal_sfr` | Star Formation Rate (SFR) ($M_\odot  · \mathrm{yr}^{-1}$) |
| **7** | `sgal_starmass` | Stellar mass of the galaxy ($M_\odot$) |
| **8** | `gal_index_select` | Galaxy Catalogue ID |
| **9** | `central_CID_array` | Central galaxy Catalogue ID |
| **10** | `progen_CID` | Progenitor galaxy Catalogue ID |
| **11** | `progen_cgalCID` | Progenitor central galaxy Catalogue ID |
| **12** | `descen_CID_select` | Descendant galaxy Catalogue ID |
| **13** | `descen_cgalCID_select` | Descendant central galaxy Catalogue ID |
| **14 – 16** | `rotated_vel[:, 0:3]` | Velocity components $(v_x, v_y, v_z)$ relative to the central galaxy reference frame ($\mathrm{km · s^{-1}}$) |
| **17** | `Radius` | Host halo $R_\mathrm{200}$ radius (kpc) |
| **18** | `sgal_radiihalf` | Galaxy half-mass radius (kpc) |
| **19 – 21** | `sgal_abspos[:, 0:3]` | Coordinates of galaxy in simulation box (kpc) |
| **22 – 24** | `cgal_abspos[:, 0:3]` | Coordinates of central galaxy in simulation box (kpc) |
| **25 – 27** | `sgal_absvel[:, 0:3]` | Velocity of galaxy in simulation box ($\mathrm{km · s^{-1}}$) |
| **28 – 30** | `cgal_absvel[:, 0:3]` | Velocity of central galaxy in simulation box ($\mathrm{km · s^{-1}}$) |
| **31** | `sgal_mass` | Total mass of the galaxy ($M_\odot$) |
| **32** | `Mass` | Host halo mass within $R_\mathrm{200c}$ ($M_\odot$) |
| **33** | `merge` | Merger status flag |
| **34** | `sgal_dmmass` | Dark matter mass of the galaxy ($M_\odot$) |
| **35** | `sgal_gasmass` | Gas mass of the galaxy ($M_\odot$) |

   * **2. Halo Triaxiality and Variance Data:** Compiles halo-level attributes across the $N_{re}$ runs:
     * **Central Galaxy ID:** Unique identifier mapping each central galaxy and its host halo.
     * **Ellipsoidal Triaxial Orientation:** Calculated directional vectors of the cental galaxy's / host halo's triaxial ellipsoid.
     * **Variance:** The variance (eigenvalues) along the principal axes derived from the PCA algorithm, quantifying the shape and dispersion of the distribution.

2. [`SFGVQ.py`](SFGVQ.py) Note: Utility script (imported and used as a module).
* **Description:** Takes filtered satellite galaxy data as input and extracts satellite galaxies located in different directions within *the major and minor axis planes*, based on their coordinates relative to the central galaxy or halo reference frame (ellipsoidal triaxial coordinate system). Each directional subset consists of galaxies locating inside a cone defined by the given direction as its axis and a half-opening angle `sample_angle` (set to $45^\circ$ in our paper). The orientation array is generated via `orient = np.linspace(0, 360, int((360/orient_angle)+1))`, where each `orient` value corresponds to a specific directional angle (e.g., with an `orient_angle` of $30^\circ$, `orient = 0^\circ` and $180^\circ$ correspond to the major axis, while $90^\circ$ and $270^\circ$ correspond to the minor axis). For each orientation, it automatically aggregates and sums up the total counts across all halos for total galaxies, star-forming galaxies, green valley galaxies, quenched galaxies, and computes the corresponding quenched fraction within each conical volume.
   * **Output Data Structure:**
     * **1. Directional Statistics and Quenched Fraction Distribution Data (`N_distri`):** Saved as a `.npy` file (file path format: `.../NQfrac_{bins_num}_{snap_num:03}.npy`), with an array shape of `(6, N_orient)` (6 rows and `N_orient` columns).

| Row Index | Variable Name in Script | Data Type | Physical Meaning & Description |
| :--- | :--- | :--- | :--- |
| **0** | `Ndata` | `float / int` | Total number of galaxies within the conical volume for each direction (`SF + GV + Q`) |
| **1** | `QNdata` | `float / int` | Total Number of quenched galaxies within the conical volume for each direction |
| **2** | `GVNdata` | `float / int` | Total Number of green valley galaxies within the conical volume for each direction |
| **3** | `SFNdata` | `float / int` | Total Number of star-forming galaxies within the conical volume for each direction |
| **4** | `Qfracdata` | `float` | Quenched fraction, defined as `QNdata / Ndata` (set to `NaN` if the total number is 0) |
| **5** | `orient` | `float` | Orientation angle in degrees (e.g., generated by `np.linspace` with a given `orient_angle`, representing specific directions relative to the axes) |

3. [`classifier.ipynb`](classifier.ipynb)
   * **Description:** Calls `stack_clusters_partaxis.py` to iterate across four cosmological simulations (SIMBA, TNG100, EAGLE, and SIMBA-nofb) with $N_\mathrm{re} = 50$ independent repeated runs. For each snapshot (`snapnum`), it extracts and aggregates surrounding galaxy data from all qualified halos within $10 R_\mathrm{200c}$ of the central galaxy (applying a stellar number limit of $100$ and a halo mass limit of $1\times 10^{12} M_\odot$). 
   * **Output Data Structure:**
     * **1. Directory Structure:** Saved under `sgaldata_save_path/GalInHalo_{part_axis}_axis/{snapnum:03}/`, where `part_axis` can be `'star'` or `'dm'` denoting the stellar particle component or dark matter particle component used to define the triaxial frame.
     * **2. Central Galaxy ID and Triaxial Directions (`CID_vec.pkl`):** A pickled list storing the central galaxy catalog IDs and their ellipsoidal triaxial directions across all repetitions (`CID_vec_var`) for the processed snapshot.
     * **3. Averaged-Axis Galaxy Data (`galaxy_ave.npy`):** Saved as a `.npy` file containing galaxy data combined from all halos, rotated and aligned with the triaxial coordinates averaged over the $N_\mathrm{re} = 50$ independent runs (`N_i = 0`).
     * **4. Individual Run Galaxy Data (`galaxy_{i}.npy`):** Saved as `.npy` files containing galaxy data combined from all halos for each individual repeated run, rotated using the actual triaxial vectors of that specific run (`N_i` ranging from $1$ to $N_\mathrm{re}$, corresponding to loop indices `i` from $0$ to $N_\mathrm{re}-1$).

| File / Output Name | File Format | Path / Naming Convention | Physical Meaning & Description |
| :--- | :--- | :--- | :--- |
| **`CID_vec.pkl`** | Pickle (`.pkl`) | `.../GalInHalo_{part_axis}_axis/{snapnum:03}/CID_vec.pkl` | List of central galaxy Catalogue IDs and their ellipsoidal triaxial directions across repetitions (`CID_vec_var`) |
| **`galaxy_ave.npy`** | NumPy Array (`.npy`) | `.../GalInHalo_{part_axis}_axis/{snapnum:03}/galaxy_ave.npy` | Galaxy data combined from all halos, aligned with the triaxial coordinates averaged over $N_\mathrm{re} = 50$ independent runs (`N_i = 0`) |
| **`galaxy_{i}.npy`** | NumPy Array (`.npy`) | `.../GalInHalo_{part_axis}_axis/{snapnum:03}/galaxy_{i}.npy` | Galaxy data combined from all halos, aligned with the actual triaxial coordinates from individual run `i` (`i` from $0$ to $N_\mathrm{re}-1$) |

4. [`geometry.ipynb`](geometry.ipynb)
   * **Description:** Processes the previously generated `GalInHalo_{part_axis}` data across the four cosmological simulations (SIMBA, TNG100, EAGLE, and SIMBA-nofb). It first loads and matches the central galaxy IDs and triaxial directions from both stellar (`star`) and dark matter (`dm`) pickle files (`CID_vec.pkl`), applies a halo mass cut ($M_\mathrm{200c} \ge 1 \times 10^{12} \, M_\odot$), and iterates across the $N_\mathrm{re} = 50$ independent repetitions to compute alignment angles between the stellar and dark matter ellipsoidal principal axes.
   * **Output Data Structure:**
     * **1. Directory Structure:** Saved under `cluster_data/{sim_name}/{snapnum}/`.
     * **2. Consolidated Cluster Data (`all_data.pkl`):** A pickled list (`all_data`) storing structured records for each qualified halo. Each entry contains the central galaxy ID (`CID`), mean alignment angles (`ave_angle_i`) and their standard deviations (`std_angle_i`) between stellar and dark matter principal axes across repetitions, the average stellar triaxial directions (`ave_vec_star`) and variances (`ave_var_star`), and the average dark matter triaxial directions (`ave_vec_dm`) and variances (`ave_var_dm`).

| File / Output Name | File Format | Path / Naming Convention | Physical Meaning & Description |
| :--- | :--- | :--- | :--- |
| **`all_data.pkl`** | Pickle (`.pkl`) | `cluster_data/{sim_name}/{snapnum}/all_data.pkl` | Comprehensive structured list containing central galaxy IDs, stellar-DM axis alignment angles (mean and standard deviation), and averaged triaxial vectors/variances for both stellar and dark matter components across the $N_\mathrm{re} = 50$ runs |

*****

5. [`ASGQ_Mstar.ipynb`](ASGQ_Mstar.ipynb)
   * **Description:** Generates Figure 1 by utilizing the qualified halo data from `all_data.pkl` (produced by `geometry.ipynb`). For each cosmological simulation (SIMBA, TNG100, EAGLE, and SIMBA-nofb) and across three radial bins defined by projected distance ratios ($R/R_\mathrm{200c} \in (2r_\mathrm{half,c}/R_\mathrm{200c}, 1)$, $(1, 3)$, and $(3, 5)$), it filters satellite galaxies located along the major and minor axes of the host halo (within a sample angle of $45^\circ$). It then classifies galaxies into Star-Forming (SF) and Quenched (Q) populations based on a redshift-dependent SFR threshold, computes the quenched fraction $F_\mathrm{Q}$ for both axes, and calculates the interpolated difference curve ($\Delta F_\mathrm{Q} = F_\mathrm{Q, major} - F_\mathrm{Q, minor}$), with asymmetric error bars corresponding to the $3\sigma$ / 99.7% confidence interval) to quantify anisotropic satellite quenching across stellar mass bins.
   * **Output Data Structure:**
     * **1. Intermediate Data Directory:** Saved under `fig_data_starmass/{sim_name}/` containing compressed NumPy arrays (`.npz`) for each radial bin (e.g., `all_R_0_to_1.npz`, `all_R_1_to_3.npz`, `all_R_3_to_5.npz`).
     * **2. Figure Results Directory:** Saved under `fig_result/` containing publication-ready PDF plots for each simulation (e.g., `{sim_name}_BCG_all.pdf`).

| File / Output Name | File Format | Path / Naming Convention | Physical Meaning & Description |
| :--- | :--- | :--- | :--- |
| **`{axis_align}_R_{low}_to_{high}.npz`** | NumPy Archive (`.npz`) | `fig_data_starmass/{sim_name}/all_R_{low}_to_{high}.npz` | Contains interpolated stellar masses (`x_star`), median quenched fraction differences (`star_median`), and associated asymmetric errors (`star_err`) between major and minor axes for the given $R/R_\mathrm{200c}$ radial interval |
| **`{sim_name}_BCG_all.pdf`** | PDF Document (`.pdf`) | `fig_result/{sim_name}_BCG_all.pdf` | Final multi-radial-bin plot illustrating the quenched fraction difference ($\Delta F_\mathrm{Q}$) as a function of stellar mass ($\log_{10}(M_* / M_\odot)$) for the specified cosmological simulation |

6. [`ASGQ_Radius.ipynb`](ASGQ_Radius.ipynb)
   * **Description:** Examines radial dependencies and scale-segregation profiles.

7. [`ASGQ_Nratio.ipynb`](ASGQ_Nratio.ipynb)
   * **Description:** Computes major-to-minor axis galaxy number ratios and related environmental metrics.

8. [`mass_function.ipynb`](mass_function.ipynb)
   * **Description:** Computes halo mass functions and associated statistical distributions.

---

## 📁 Data Directories

* `data_input/` — Raw or pre-filtered observational/simulation catalogs used as inputs.
* `data_processed/` — Intermediate processed outputs, stacked profiles, or matching results.
* `outputs/` — Generated figures, tables, and numerical results for the paper.

---

## ⚙️ Usage Notes

Make sure you have activated the project's Conda environment (`workenv`) from the root directory before running these scripts or notebooks:
