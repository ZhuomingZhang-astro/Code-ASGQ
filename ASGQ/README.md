# ASGQ Core Pipeline

Core statistical pipeline for primary data processing, statistical analysis, and figure generation for the project.

---

## 📂 Code Scripts

* [`stack_clusters_partaxis.py`](stack_clusters_partaxis.py) — Computes ellipsoidal triaxiality based on central galaxy stellar particle data and host halo dark matter particles within $R_{200c}$ across independent repeated runs.

---

## 📊 Output Data Structure

The pipeline computes ellipsoidal triaxiality across $N_{re}$ independent repeated runs, organizing results into two main output categories:

### 1. Galaxy and Environment Data
For each of the $N_{re}$ independent repeated runs, this dataset records information for every central galaxy and all surrounding galaxies within $10 R_{200c}$ of its position. 

The data fields (columns) stored in `sgaldata[N_i]` correspond to the following physical parameters:

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

---

### 2. Halo Triaxiality and Variance Data
For each of the $N_{re}$ independent repeated runs, this dataset compiles halo-level attributes:
* **Central Galaxy / Host Halo ID**: Unique identifier mapping each host halo and its central galaxy.
* **Ellipsoidal Triaxial Orientation**: The calculated directional vectors / axes of the host halo's triaxial ellipsoid.
* **Variance**: Statistical variance associated with the triaxiality computations across the repeated runs.

* [`SFGVQ.py`](SFGVQ.py) — Computes halo mass functions and associated statistical distributions.
* [`classifier.ipynb`](classifier.ipynb) — Computes halo mass functions and associated statistical distributions.
* [`geometry.ipynb`](geometry.ipynb) — Computes halo mass functions and associated statistical distributions.
* [`ASGQ_Mstar.ipynb`](ASGQ_Mstar.ipynb) — Computes halo mass functions and associated statistical distributions.
* [`ASGQ_Radius.ipynb`](ASGQ_Radius.ipynb) — Computes halo mass functions and associated statistical distributions.
* [`ASGQ_Nratio.ipynb`](ASGQ_Nratio.ipynb) — Computes halo mass functions and associated statistical distributions.
* [`mass_function.ipynb`](mass_function.ipynb) — Computes halo mass functions and associated statistical distributions.


## 📁 Data Directories

* `data_input/` — *(Optional example)* Contains raw or pre-filtered observational/simulation catalogs used as inputs.
* `data_processed/` — Stores intermediate processed outputs, stacked profiles, or matching results.
* `outputs/` — Generated figures, tables, and numerical results for the paper.

---

## ⚙️ Usage Notes

Make sure you have activated the project's Conda environment (`workenv`) from the root directory before running these scripts or notebooks:
