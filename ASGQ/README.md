# ASGQ Core Pipeline

Core statistical pipeline for primary data processing, statistical analysis, and figure generation for the project.

---

## 📂 Code Scripts & Notebooks

1. [`stack_clusters_partaxis.py`](stack_clusters_partaxis.py) Note: Utility script (imported and used as a module).
   * **Description:** Computes ellipsoidal triaxiality based on "central galaxy stellar particle data" / "host halo dark matter particles within $R_\mathrm{200c}$" across $N_\mathrm{re}$ independent repeated runs.
   * **Output Data Structure:**
     * **1. Galaxy and Environment Data (`sgaldata`):** Records information for every central galaxy and all surrounding galaxies within $10 R_\mathrm{200c}$ of its position across independent runs.
       
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
     * **1. Directional Statistics and Quenched Fraction Distribution Data (`N_distri`):** with an array shape of `(6, N_orient)`.

| Row Index | Variable Name in Script | Data Type | Physical Meaning & Description |
| :--- | :--- | :--- | :--- |
| **0** | `Ndata` | `float / int` | Total number of galaxies within the conical volume for each direction (`SF + GV + Q`) |
| **1** | `QNdata` | `float / int` | Total Number of quenched galaxies within the conical volume for each direction |
| **2** | `GVNdata` | `float / int` | Total Number of green valley galaxies within the conical volume for each direction |
| **3** | `SFNdata` | `float / int` | Total Number of star-forming galaxies within the conical volume for each direction |
| **4** | `Qfracdata` | `float` | Quenched fraction, defined as `QNdata / Ndata` (set to `NaN` if the total number is 0) |
| **5** | `orient` | `float` | Orientation angle in degrees (e.g., generated by `np.linspace` with a given `orient_angle`, representing specific directions relative to the axes) |

4. [`classifier.ipynb`](classifier.ipynb)
   * **Description:** Jupyter notebook for galaxy classification and population filtering.

5. [`geometry.ipynb`](geometry.ipynb)
   * **Description:** Notebook for geometric analyses, alignment angles, and coordinate transformations.

6. [`ASGQ_Mstar.ipynb`](ASGQ_Mstar.ipynb)
   * **Description:** Analyzes stellar mass dependent trends and property distributions in ASGQ studies.

7. [`ASGQ_Radius.ipynb`](ASGQ_Radius.ipynb)
   * **Description:** Examines radial dependencies and scale-segregation profiles.

8. [`ASGQ_Nratio.ipynb`](ASGQ_Nratio.ipynb)
   * **Description:** Computes major-to-minor axis galaxy number ratios and related environmental metrics.

9. [`mass_function.ipynb`](mass_function.ipynb)
   * **Description:** Computes halo mass functions and associated statistical distributions.

---

## 📁 Data Directories

* `data_input/` — Raw or pre-filtered observational/simulation catalogs used as inputs.
* `data_processed/` — Intermediate processed outputs, stacked profiles, or matching results.
* `outputs/` — Generated figures, tables, and numerical results for the paper.

---

## ⚙️ Usage Notes

Make sure you have activated the project's Conda environment (`workenv`) from the root directory before running these scripts or notebooks:
