"""
Tiny dolfinx-free provenance helper shared by `sem.cli` (the open-pore
result file/console header) and `sem.scripts.derive_geometry` (the
`*_derivation.json` provenance object).

Kept in its own module -- rather than folded into `sem.cli` -- so
`sem.scripts.derive_geometry`, which deliberately imports only numpy,
scipy, `sem.grid_io` and `sem.geometry_profiles` (see that module's
docstring), does not have to pull in `sem.cli`'s heavier import surface
(argparse, `sem.config`, `sem.rotation`, ...) just to stamp a commit hash.
"""

import subprocess
from pathlib import Path

__all__ = ["git_commit"]


def git_commit(cwd=None, timeout=2.0):
    """
    Short git commit SHA for provenance stamping, or ``"unknown"`` if it
    cannot be determined. Never raises: any failure (git not installed, not
    a git checkout, timeout, permissions, unexpected output, ...) is
    swallowed and mapped to ``"unknown"`` so it is always safe to call from
    an output-writing code path.

    `cwd` defaults to the directory containing this file (the `sem`
    package directory); git walks upward from there to find the repository
    root, so this works whether `sem` is installed editable from a git
    checkout or not (in which case it falls back to ``"unknown"``).
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=cwd or Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except Exception:
        return "unknown"
    if result.returncode != 0:
        return "unknown"
    sha = result.stdout.strip()
    return sha if sha else "unknown"
