"""
Import-surface tests for the `sem` package.

`sem/__init__.py` uses a PEP 562 lazy `_LAZY` map + `__getattr__` so that
merely `import sem` (and importing the numpy-only submodules `sem.grid_io`
and `sem.pore_geometry`, plus `sem.config`) never touches dolfinx. Only
accessing a dolfinx-backed attribute (e.g. `sem.VerticalMovementSEM`)
should trigger the heavy import. These tests run the surface in a
subprocess with `dolfinx` forced unavailable so they verify the real
behaviour a pip-only / conda-less environment would see, regardless of
whether dolfinx happens to be installed in the environment running the
test suite.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_import_surface_without_dolfinx():
    """sem, sem.grid_io, sem.pore_geometry, sem.config must import cleanly
    even when dolfinx cannot be imported at all."""
    code = (
        "import sys\n"
        "sys.modules['dolfinx'] = None\n"
        "import sem\n"
        "import sem.grid_io\n"
        "import sem.pore_geometry\n"
        "import sem.config\n"
        "print('IMPORT_SURFACE_OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        "dolfinx-free import surface failed:\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "IMPORT_SURFACE_OK" in result.stdout


def test_accessing_dolfinx_backed_attribute_fails_cleanly_without_dolfinx():
    """Lazy attributes that need dolfinx should raise (not silently no-op)
    once dolfinx is genuinely unavailable, proving the surface is lazy
    rather than accidentally eager."""
    code = (
        "import sys\n"
        "sys.modules['dolfinx'] = None\n"
        "import sem\n"
        "sem.VerticalMovementSEM\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode != 0
    assert "dolfinx" in result.stderr.lower()


def _version_from_setup_py():
    """Ask setup.py for its computed version via the standard distutils
    `--version` global option, rather than importing it (setup.py is not
    import-safe as a module) or regexing a literal that no longer exists
    (setup.py derives its version from sem/__init__.py at run time)."""
    result = subprocess.run(
        [sys.executable, "setup.py", "--version"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"setup.py --version failed:\n{result.stderr}"
    lines = [line.strip() for line in result.stdout.strip().splitlines() if line.strip()]
    assert lines, "setup.py --version produced no output"
    return lines[-1]


def _version_from_citation_cff():
    text = (REPO_ROOT / "CITATION.cff").read_text()
    match = re.search(r'^version:\s*["\']?([^"\'\n]+)["\']?\s*$', text, re.MULTILINE)
    assert match, "could not find a 'version:' line in CITATION.cff"
    return match.group(1).strip()


def test_version_matches_setup_py_and_citation_cff():
    import sem

    assert sem.__version__ == _version_from_setup_py()
    assert sem.__version__ == _version_from_citation_cff()


def test_version_is_0_2_0():
    import sem

    assert sem.__version__ == "0.2.0"
