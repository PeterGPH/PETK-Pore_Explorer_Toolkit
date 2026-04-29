"""
Utility scripts for SEM package.
Makes pdb2xyz.py and gen_dist.py discoverable by the main package.
"""

from pathlib import Path

# Make script paths available to the package
SCRIPT_DIR = Path(__file__).parent
PDB2XYZ_SCRIPT = SCRIPT_DIR / "pdb2xyz.py"
GEN_DIST_SCRIPT = SCRIPT_DIR / "gen_dist.py"

# Verify scripts exist
def check_scripts():
    """Check if all required scripts are present"""
    missing = []
    if not PDB2XYZ_SCRIPT.exists():
        missing.append("pdb2xyz.py")
    if not GEN_DIST_SCRIPT.exists():
        missing.append("gen_dist.py")
    
    if missing:
        import warnings
        warnings.warn(f"Missing scripts in {SCRIPT_DIR}: {missing}")
    
    return len(missing) == 0

# Check on import
_scripts_available = check_scripts()

__all__ = [
    'SCRIPT_DIR', 
    'PDB2XYZ_SCRIPT', 
    'GEN_DIST_SCRIPT', 
    'BIN_COMPARE_SCRIPT'
]