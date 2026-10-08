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

# Package metadata
__version__ = "0.1.0"
__author__ = "Pinhao Gu"
__email__ = "pinhao2@illinois.edu"
__description__ = "Steric Exclusion Model calculations for nanopore analytics"

# Public API — see _LAZY below for the (name -> module) routing.
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
    'main',
    'create_sem_from_config',
]

# PEP 562 lazy attribute loading. Why:
#   - `VerticalMovementSEM`, `PoreGeometry`, `cli.main` transitively import
#     `dolfinx`, which only ships via conda-forge. Eager imports here mean
#     `import sem` fails for any pip-only environment (e.g. our CI smoke
#     test, downstream tools that just want VanDerWaalsRadii or load_config).
#   - Lazy `__getattr__` defers every public symbol until first access. The
#     pip-only smoke test gets a successful `import sem` + `__version__`
#     readout; real users running `from sem import VerticalMovementSEM`
#     trigger the dolfinx import exactly when the symbol is requested, with
#     a clean ImportError if dolfinx isn't installed.
#   - This is the modern Python idiom for optional heavy dependencies
#     (numpy, sklearn, dask all use the same pattern).
_LAZY = {
    'VerticalMovementSEM':       ('.vertical_movement_sem', 'VerticalMovementSEM'),
    'AnalyteOverlapError':       ('.vertical_movement_sem', 'AnalyteOverlapError'),
    'VanDerWaalsRadii':          ('.van_der_waals',         'VanDerWaalsRadii'),
    'PoreGeometry':              ('.pore_geometry',         'PoreGeometry'),
    'SimpleConductivityModel':   ('.conductivity_models',   'SimpleConductivityModel'),
    'load_config':               ('.config',                'load_config'),
    'validate_config':           ('.config',                'validate_config'),
    'print_config_summary':      ('.config',                'print_config_summary'),
    'create_example_config':     ('.config',                'create_example_config'),
    'RotationSpec':              ('.rotation',              'RotationSpec'),
    'rotate_pdb_to_grid_center': ('.rotation',              'rotate_pdb_to_grid_center'),
    'parse_angle_file':          ('.rotation',              'parse_angle_file'),
    'random_uniform_rotations':  ('.rotation',              'random_uniform_rotations'),
    'main':                      ('.cli',                   'main'),
    'create_sem_from_config':    ('.cli',                   'create_sem_from_config'),
}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        module_path, attr = _LAZY[name]
        mod = importlib.import_module(module_path, package=__package__)
        value = getattr(mod, attr)
        # Cache on the module so subsequent accesses skip __getattr__.
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    # Make tab-completion / dir(sem) report the public names.
    return sorted(set(globals()) | set(_LAZY))
