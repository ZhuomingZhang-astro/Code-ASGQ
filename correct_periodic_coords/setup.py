from setuptools import setup, find_packages

setup(
    name="correct_periodic_coords",  # Package name
    version="0.1.1",  # Version number
    packages=find_packages(),  # Automatically discover all modules in the package
    install_requires=[  # Dependent third-party libraries
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "numba>=0.55.0",
        "joblib>=1.1.0"
    ],
    author="Zhuoming Zhang",  # Author name (can be modified)
    description="A tool for finding ellipsoid axes via sphere sampling and PCA or ITM",  # Package description
    url="",  # Project URL (optional)
)
