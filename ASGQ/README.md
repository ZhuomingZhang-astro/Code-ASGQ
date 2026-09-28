# ASGQ Core Pipeline

Core statistical pipeline for primary data processing, statistical analysis, and figure generation for the project.

---

## 📂 Code Scripts & Notebooks

1. [`stack_clusters_partaxis.py`](stack_clusters_partaxis.py)
   * **Description:** Computes ellipsoidal triaxiality based on "central galaxy stellar particle data" / "host halo dark matter particles within $R_\mathrm{200c}$" across $N_\mathrm{re}$ independent repeated runs.
   * **Output Data Structure:**
     * **1. Galaxy and Environment Data (`sgaldata`):** Records information for every central galaxy and all surrounding galaxies within $10 R_\mathrm{200c}$ of its position across independent runs.
       
       | Column / Index | Variable Name in Script | Physical Meaning & Description |
       | :--- | :--- | :--- |
       | **0 – 2** | `rotated_pos[:, 0:3]` | Rotated 3D position coordinates $(x, y, z)$ |
       | **3** | `r` | Radial distance ($r$) |
       | **4** | `ra_deg` | Right Ascension (RA in degrees) |
       | **5** | `dec_deg` | Declination (DEC in degrees) |
       | **6** | `sgal_sfr` | Star Formation Rate (SFR) |
       | **7** | `sgal_starmass` | Stellar mass of the galaxy |
       | **8** | `gal_index_select` | Galaxy Catalogue ID (`CatalogueID`) |
       | **9** | `central_CID_array` | Central galaxy Catalogue ID (`central_galaxy_CatalogueID`) |
       | **10** | `progen_CID` | Progenitor Catalogue ID (`progen_CatalogueID`) |
       | **11** | `progen_cgalCID` | Progenitor central galaxy Catalogue ID (`progen_central_galaxy_CatalogueID`) |
       | **12** | `descen_CID_select` | Descendant Catalogue ID (`descen_CatalogueID`) |
       | **13** | `descen_cgalCID_select` | Descendant central galaxy Catalogue ID (`descen_central_galaxy_CatalogueID`) |
       | **14 – 16** | `rotated_vel[:, 0:3]` | Rotated 3D velocity components $(vx, vy, vz)$ |
       | **17** | `Radius` | Host halo $R_{200}$ radius ($R_{200}$) |
       | **18** | `sgal_radiihalf` | Galaxy half-mass radius ($R_{gal\_half}$) |
       | **19 – 21** | `sgal_abspos[:, 0:3]` | Absolute position of the galaxy $(x, y, z)$ |
       | **22 – 24** | `cgal_abspos[:, 0:3]` | Absolute position of the central galaxy $(x, y, z)$ |
       | **25 – 27** | `sgal_absvel[:, 0:3]` | Absolute velocity of the galaxy $(vx, vy, vz)$ |
       | **28 – 30** | `cgal_absvel[:, 0:3]` | Absolute velocity of the central galaxy $(vx, vy, vz)$ |
       | **31** | `sgal_mass` | Total mass of the galaxy (`gal_total_mass`) |
       | **32** | `Mass` | Host halo mass within $R_{200c}$ (`Halo_M200c`) |
       | **33** | `merge` | Merger status flag (`whether_merge`) |
       | **34** | `sgal_dmmass` | Dark matter mass of the galaxy (`gal_dmmass`) |
       | **35** | `sgal_gasmass` | Gas mass of the galaxy (`gal_gasmass`) |

     * **2. Halo Triaxiality and Variance Data:** Compiles halo-level attributes across the $N_{re}$ runs:
       * **Central Galaxy / Host Halo ID:** Unique identifier mapping each host halo and its central galaxy.
       * **Ellipsoidal Triaxial Orientation:** Calculated directional vectors / axes of the host halo's triaxial ellipsoid.
       * **Variance:** Statistical variance associated with the triaxiality computations across the repeated runs.

2. [`SFGVQ.py`](SFGVQ.py)
   * **Description:** Core script for evaluating spatial segregation and anisotropic satellite quenching statistics.

3. [`classifier.ipynb`](classifier.ipynb)
   * **Description:** Jupyter notebook for galaxy classification and population filtering.

4. [`geometry.ipynb`](geometry.ipynb)
   * **Description:** Notebook for geometric analyses, alignment angles, and coordinate transformations.

5. [`ASGQ_Mstar.ipynb`](ASGQ_Mstar.ipynb)
   * **Description:** Analyzes stellar mass dependent trends and property distributions in ASGQ studies.

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
