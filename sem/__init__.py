"""
SEM (Steric Exclusion Model) calculation package.

This package provides tools for calculating vertical movement of analytes through nanopores
using finite element methods. It supports multiple pore geometries including:

- Cylindrical pores (with optional corner rounding)
- Double cone (hourglass) pores  
- Biological pores from PDB structures
- Binary file-based pore structures

Key Features:
- Van der Waals radii assignment for accurate atomic volumes
- Multiple pore geometry support
- FEniCS-based finite element calculations
- JSON configuration management
- Command line interface

Usage:
    from sem import VerticalMovementSEM
    from sem.config import load_config
    
    # Load configuration
    config = load_config("config.json")
    
    # Create SEM instance
    sem = VerticalMovementSEM(**config_params)
    
    # Run simulation
    results = sem.run()

Command Line Usage:
    python -m sem config.json run           # Run simulation
    python -m sem config.json preview_only  # Preview only
"""

import logging

try:
    from mpi4py import MPI
    _MPI_RANK = MPI.COMM_WORLD.Get_rank()
except Exception:
    _MPI_RANK = 0

# Configure logging for the package
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s: %(levelname)s: %(message)s"
)

_root_logger = logging.getLogger()


class _RootRankFilter(logging.Filter):
    """Allow log records from MPI rank 0; always keep errors."""

    def filter(self, record: logging.LogRecord) -> bool:
        if _MPI_RANK == 0:
            return True
        return record.levelno >= logging.ERROR


_root_logger.addFilter(_RootRankFilter())

if _MPI_RANK != 0:
    _root_logger.setLevel(logging.ERROR)
    for _handler in _root_logger.handlers:
        _handler.setLevel(logging.ERROR)

# Import main classes and functions
from .vertical_movement_sem import VerticalMovementSEM, AnalyteOverlapError
from .van_der_waals import VanDerWaalsRadii
from .pore_geometry import PoreGeometry
from .conductivity_models import SimpleConductivityModel
from .config import load_config, validate_config, print_config_summary, create_example_config
from .rotation import RotationSpec, rotate_pdb_to_grid_center, parse_angle_file, random_uniform_rotations
from .cli import main, create_sem_from_config

# Package metadata
__version__ = "1.0.0"
__author__ = "Pinhao Gu"
__email__ = "pinhao2@illinois.edu"
__description__ = "Steric Exclusion Model calculations for nanopore analytics"

# Export main classes and functions
__all__ = [
    'VerticalMovementSEM',
    'AnalyteOverlapError',
    'VanDerWaalsRadii', 
    'PoreGeometry',
    'SimpleConductivityModel',
    'load_config',
    'validate_config',
    'print_config_summary',
    'create_example_config',
    'RotationSpec',
    'rotate_pdb_to_grid_center',
    'parse_angle_file',
    'random_uniform_rotations',
    'main'
]
