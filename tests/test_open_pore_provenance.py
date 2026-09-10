"""
Unit tests for the provenance helpers used by `sem.cli`: `_git_commit`
(imported into `sem.cli` as an alias of `sem.provenance.git_commit`) and
`_open_pore_result_lines` (defined in `sem.cli`). Both are pure functions
(no dolfinx import, no real `VerticalMovementSEM` instance required), so
they can — and must — be exercised without dolfinx: `_open_pore_result_lines`
is fed a lightweight stand-in object exposing the same public attributes
`VerticalMovementSEM.calculate_open_pore_current`'s caller reads, instead
of a real (dolfinx-backed) SEM instance.
"""

import subprocess
import types

from sem.cli import _git_commit, _open_pore_result_lines
from sem.grid_io import RAMP_MAX, RAMP_MIN
import sem


class _FakeSEM(types.SimpleNamespace):
    """Minimal stand-in exposing exactly the attributes
    `_open_pore_result_lines` reads for a given `pore_type`."""


def _cylindrical_sem(**overrides):
    defaults = dict(
        pore_type="cylindrical",
        pore_radius=100.0,
        corner_radius=0.0,
        chamfer_depth=None,
        distance_metric="euclidean",
        membrane_thickness=200.0,
        voltage=0.5,
        bulk_conductivity=1.660843,
        grid_resolution=1.0,
    )
    defaults.update(overrides)
    return _FakeSEM(**defaults)


# ---------------------------------------------------------------------------
# _git_commit
# ---------------------------------------------------------------------------
def test_git_commit_returns_short_sha_in_this_repo():
    sha = _git_commit()
    assert sha != "unknown"
    assert 4 <= len(sha) <= 40
    int(sha, 16)  # must be valid hex


def test_git_commit_never_raises_on_subprocess_error(monkeypatch):
    def _boom(*args, **kwargs):
        raise OSError("git not found")

    monkeypatch.setattr(subprocess, "run", _boom)
    assert _git_commit() == "unknown"


def test_git_commit_never_raises_on_timeout(monkeypatch):
    def _timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="git", timeout=2)

    monkeypatch.setattr(subprocess, "run", _timeout)
    assert _git_commit() == "unknown"


def test_git_commit_returns_unknown_on_nonzero_exit(monkeypatch):
    class _Result:
        returncode = 128
        stdout = ""

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: _Result())
    assert _git_commit() == "unknown"


# ---------------------------------------------------------------------------
# _open_pore_result_lines
# ---------------------------------------------------------------------------
def test_header_contains_sem_version_git_commit_and_ramp():
    fake = _cylindrical_sem()
    lines = _open_pore_result_lines(fake, 4.5e-1)

    assert f"# sem_version: {sem.__version__}" in lines
    assert any(line.startswith("# git_commit: ") for line in lines)
    assert f"# ramp: {RAMP_MIN} {RAMP_MAX}" in lines


def test_header_preserves_existing_fields():
    fake = _cylindrical_sem()
    lines = _open_pore_result_lines(fake, 4.5e-1)

    assert "# Open pore current calculation results" in lines
    assert "# Pore type: cylindrical" in lines
    assert "# Pore radius: 100.0 Angstrom" in lines
    assert "# Distance metric: euclidean" in lines
    assert "# Membrane thickness: 200.0 Angstrom" in lines
    # Last line is the data row (no leading '#'), matching the original
    # single-precision-tolerant format.
    assert lines[-1] == f"{4.5e-1:.6e}"
    assert lines[-2] == "# Open_pore_current(nA)"


def test_header_includes_corner_and_chamfer_when_set():
    fake = _cylindrical_sem(corner_radius=20.0, chamfer_depth=15.0)
    lines = _open_pore_result_lines(fake, 1.0)

    assert "# Corner radius: 20.0 Angstrom" in lines
    assert "# Chamfer depth: 15.0 Angstrom" in lines


def test_header_all_lines_are_comments_except_the_last():
    fake = _cylindrical_sem()
    lines = _open_pore_result_lines(fake, 1.0)

    assert all(line.startswith("#") for line in lines[:-1])
    assert not lines[-1].startswith("#")
