# Code-ASGQ (Anisotropic Satellite Galaxy Quenching)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)

This repository contains the core analysis scripts and numerical tools used in our research investigating **Anisotropic Satellite Galaxy Quenching (ASGQ)**, Paper title: **Spatial segregation as the origin of anisotropic satellite quenching in haloes and beyond**

---

## 📂 Repository Structure

* `caesar/` — Processing snapshot/particle data from simulation for halo and galaxy catalogs.
* `pygadgetreader/` — Utility modules for reading snapshot and particle data from cosmological hydrodynamic simulations.
* `correct_periodic_coords/` — Routines for handling periodic boundary conditions and coordinate transformations.
* `ellip_axis_finder/` — Tools to determine the morphological and spatial axes of central galaxies and host halos.
* `ASGQ/` — Core statistical pipeline for measuring anisotropic satellite quenching and spatial segregation.
* `environment.yml` — Conda environment specification file containing all required dependencies.

---

## 🚀 Installation & Setup

To set up the environment and run the analysis tools, clone this repository and create the Conda environment using the provided YAML file:

```bash
# Clone the repository
git clone [https://github.com/ZhuomingZhang-astro/Code-ASGQ.git](https://github.com/ZhuomingZhang-astro/Code-ASGQ.git)
cd Code-ASGQ

# Create and activate the conda environment
conda env create -f environment.yml
conda activate workenv
