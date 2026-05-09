"""
Setup script for the SEM package with DOLFINx support and integrated pdb2xyz and gen_dist scripts.
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = """
    SEM (Steric Exclusion Model) calculation package for nanopore analytics.
    
    This package provides tools for calculating vertical movement of analytes through nanopores
    using modern finite element methods with DOLFINx and MPI parallelism support for multiple pore geometries.
    
    Features:
    - DOLFINx-based finite element solver with MPI parallelism
    - Multiple pore geometries (cylindrical, double cone, biological, binary files)
    - van der Waals radii support for accurate atomic modeling
    - Integrated PDB to XYZ conversion and distance field generation tools
    - High-performance parallel computing capabilities
    
    INSTALLATION:
    
    1. Create and activate conda environment:
       conda env create -f environment.yml
       conda activate sem-env
    
    2. Install this package in development mode:
       pip install -e .
    
    For manual installation:
    1. Install DOLFINx via conda-forge:
       conda install conda-forge::fenics-dolfinx
    
    2. Install other dependencies:
       conda install conda-forge::mpi4py conda-forge::petsc4py
       
    3. Install this package:
       pip install -e .
    
    Note: DOLFINx must be installed via conda-forge as it requires complex compiled dependencies.
    """

extras_require = {
    "viz": [
        "matplotlib>=3.5.0",
    ],
    "structure": [
        "pdbfixer>=1.8.0",
        "openmm>=8.0.0",
        "pdb2pqr>=3.0.0",
    ],
    "fem": [
        "mpi4py>=3.1.0",
        "petsc4py>=3.18.0",
        # Install DOLFINx separately (typically via conda-forge::fenics-dolfinx)
    ],
    "dev": [
        "pytest>=6.0",
        "pytest-cov>=2.10",
        "pytest-mpi>=0.4",
        "black>=22.0",
        "flake8>=4.0",
        "mypy>=0.910",
        "pre-commit>=2.15.0",
    ],
    "docs": [
        "sphinx>=4.0",
        "sphinx-rtd-theme>=1.0",
        "sphinx-autodoc-typehints>=1.12",
        "nbsphinx>=0.8.0",
    ],
    "examples": [
        "jupyter>=1.0.0",
        "ipywidgets>=7.6.0",
    ],
}

extras_require["all"] = sorted(
    {
        dependency
        for key, dependencies in extras_require.items()
        if key != "all"
        for dependency in dependencies
    }
)

setup(
    name="nanopore-sem",
    version="0.1.0",
    author="SEM Development Team",
    author_email="pinhao2@illinois.edu",
    description="Steric Exclusion Model calculations for nanopore analytics with DOLFINx and MPI support",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    python_requires=">=3.9",
    
    # Core dependencies that can be installed via pip
    install_requires=[
        # Numerical core
        "numpy>=1.23.0",
        "scipy>=1.9.0",
        "numba>=0.56.0",
        # Molecular structure handling
        "MDAnalysis>=2.3.0",
        "biopython>=1.81",
        "gridDataFormats>=1.0.1",
    ],
    
    extras_require=extras_require,
    
    entry_points={
        "console_scripts": [
            # `sem` works for both serial and parallel runs — mpi4py inside
            # sem.cli:main auto-detects the MPI rank, so users invoke
            # `sem run config.json` directly or `mpirun -n N sem run config.json`
            # for parallel. No separate `sem-mpi` entry is needed.
            "sem=sem.cli:main",
            "pdb2xyz=sem.scripts.pdb2xyz:main",
            "gen_dist=sem.scripts.gen_dist:main"
        ],
    },
    
    include_package_data=True,
    package_data={
        "sem": [
            "*.json", 
            "examples/*",
            "examples/**/*",
            "scripts/*.py",
            "data/*",
            "tests/data/*",
        ],
    },
    
    keywords="nanopore electrochemistry finite-element simulation molecular-dynamics pdb xyz distance-field dolfinx mpi parallel-computing",
    
    project_urls={
        "Bug Reports": "https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit/issues",
        "Source": "https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit",
        "Documentation": "https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit#readme",
        "DOLFINx Documentation": "https://docs.fenicsproject.org/dolfinx/",
    },
    
    zip_safe=False,
    platforms=["Linux", "MacOS"],
)
