"""
ellip_axis_finder: Sphere sampling and ellipsoid axis calculation toolkit
Computes the direction and length of ellipsoid axes for particle distributions via [spherical random sampling and PCA analysis] or [moment of inertia principal axis calculation]
"""

from .core import (
    Compute_Distances_Numba,
    Single_Sphere_Sampling,
    Sphere_Sampling_Method,
    Compute_Inertia_Tensor,
    Inertia_Tensor_Method
)

__all__ = [
    "Compute_Distances_Numba",
    "Single_Sphere_Sampling",
    "Sphere_Sampling_Method",
    "Compute_Inertia_Tensor",
    "Inertia_Tensor_Method"
]

__version__ = "0.5.7"  # Package version number
