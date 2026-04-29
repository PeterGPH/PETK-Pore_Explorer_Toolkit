"""Validate that the mesh-based overlap path correctly handles ConicalPore.

What this test demonstrates
---------------------------
The default overlap-protection branch in
``VerticalMovementSEM._update_conductivity_with_analyte`` (vertical_movement_sem.py
lines ~1206-1225) inspects the *base* conductivity grid that
``ConicalPore.get_conductivity_interpolator()`` produces. The check at
line 1211 is::

    membrane_mask = base_vals <= self.membrane_conductivity

so for the check to behave correctly the conductivity grid must:

1. Sit at ~``bulk_conductivity`` everywhere outside the membrane and
   inside the local pore opening.
2. Drop sharply across the wall and bottom out at ~``membrane_conductivity``
   inside the solid.
3. Use the *asymmetric* local pore radius
   ``local_pore_radius(z) = bottom_radius
       + (top_radius - bottom_radius) * (z + half) / (2 * half)``
   so the same radial distance can be inside the pore at +half and
   inside the wall at -half.

This script verifies all three numerically.

A subtlety worth knowing
------------------------
``_conductivity_from_distance`` runs the distance-to-wall through
``condfrac`` (sem/utils.py:15), which is a linear ramp from 1.3 Å
(fraction=0) to 4.1 Å (fraction=1) with a *floor of 1e-7*, not 0. So a
mesh point sitting exactly on or just inside the wall reports::

    σ_wall_floor = membrane_conductivity + 1e-7 * (bulk - membrane_conductivity)

i.e. about ``1e-7 * bulk_conductivity``. This is **not** a conical
issue — it's the same for the cylindrical and double-cone classes —
but it means the strict ``σ ≤ membrane_conductivity`` test only fires
when the user's ``membrane_conductivity`` is chosen to sit above that
floor. The conductivity grid itself is correct; the threshold is what
gets sensitive.

Run directly::

    python tests/test_conical_mesh_overlap.py

or under pytest::

    pytest tests/test_conical_mesh_overlap.py -v
"""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]


def _stub_fem_modules() -> None:
    """Make ``import sem.pore_geometry`` work without dolfinx/mpi4py."""
    for name in (
        "dolfinx",
        "dolfinx.fem",
        "dolfinx.io",
        "dolfinx.mesh",
        "mpi4py",
        "mpi4py.MPI",
        "petsc4py",
        "petsc4py.PETSc",
        "ufl",
    ):
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["mpi4py.MPI"].COMM_WORLD = types.SimpleNamespace(
        rank=0, size=1, Get_rank=lambda: 0
    )
    sys.modules["mpi4py"].MPI = sys.modules["mpi4py.MPI"]


def _load_pore_geometry():
    """Load sem/pore_geometry.py without triggering sem/__init__.py."""
    _stub_fem_modules()
    pkg = types.ModuleType("sem")
    pkg.__path__ = [str(REPO_ROOT / "sem")]
    sys.modules["sem"] = pkg

    def _load(name: str):
        path = REPO_ROOT / "sem" / f"{name}.py"
        spec = importlib.util.spec_from_file_location(f"sem.{name}", path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[f"sem.{name}"] = mod
        spec.loader.exec_module(mod)
        return mod

    _load("utils")
    return _load("pore_geometry")


# ---------------------------------------------------------------------
# Test parameters — chosen so asymmetry is unambiguous
# ---------------------------------------------------------------------
TOP_RADIUS = 30.0      # wide mouth at z = +HALF
BOTTOM_RADIUS = 10.0   # narrow mouth at z = -HALF
HALF = 20.0            # membrane half-thickness
BULK_SIGMA = 10.5      # S/m
MEMBRANE_SIGMA = 1e-7  # S/m

GRID_HALF_EXTENT = 60.0
GRID_N = 121  # ~1 Å spacing — enough so condfrac's transition is resolved

# How wide is the soft-transition zone in σ-space?
# condfrac maps distance ∈ [1.3, 4.1] linearly to [0, 1] of
# (bulk - membrane). So σ = bulk only for ``distance >= 4.1`` Å, and
# σ ≈ floor for ``distance <= 1.3`` Å.
PORE_BAND_DEPTH = 4.5   # Å inside the pore where σ should be ≈ bulk
WALL_BAND_DEPTH = 4.5   # Å inside the wall where σ should be ≈ floor

# Practical thresholds derived from condfrac
SIGMA_BULK_THRESHOLD = 0.95 * BULK_SIGMA      # "≈ bulk"
SIGMA_FLOOR_THRESHOLD = 1e-3 * BULK_SIGMA     # "≈ floor"


def _local_radius_analytical(z: float) -> float:
    """Closed-form local pore radius at axial coordinate z (frustum)."""
    t = np.clip((z + HALF) / (2.0 * HALF), 0.0, 1.0)
    return BOTTOM_RADIUS + (TOP_RADIUS - BOTTOM_RADIUS) * t


def _build_pore(pg):
    """Construct the ConicalPore on a uniform 3D grid."""
    ax = np.linspace(-GRID_HALF_EXTENT, GRID_HALF_EXTENT, GRID_N)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    return pg.PoreGeometry.create_pore(
        "conical",
        X=X, Y=Y, Z=Z,
        top_radius=TOP_RADIUS,
        bottom_radius=BOTTOM_RADIUS,
        membrane_half_thickness=HALF,
        bulk_conductivity=BULK_SIGMA,
        membrane_conductivity=MEMBRANE_SIGMA,
    )


def _signed_distance_to_wall(xyz: np.ndarray) -> float:
    """Positive inside the pore, negative inside the wall, big positive in bulk.

    Bulk (|z| > HALF) returns a large positive sentinel because we never want
    to claim a bulk point is "near the wall".
    """
    x, y, z = xyz
    if abs(z) > HALF:
        return 1e3  # bulk
    r = np.hypot(x, y)
    return _local_radius_analytical(z) - r


def _classify_expected(xyz: np.ndarray) -> str:
    """Ground-truth label using the soft transition width of condfrac."""
    sd = _signed_distance_to_wall(xyz)
    if sd >= PORE_BAND_DEPTH:
        return "bulk_or_pore"
    if sd <= -WALL_BAND_DEPTH:
        return "wall_floor"
    return "transition"


def _classify_observed(interp, xyz: np.ndarray) -> tuple[float, str]:
    """Sample the conductivity grid the same way the SEM solver does."""
    sigma = float(interp([xyz])[0])
    if sigma >= SIGMA_BULK_THRESHOLD:
        return sigma, "bulk_or_pore"
    if sigma <= SIGMA_FLOOR_THRESHOLD:
        return sigma, "wall_floor"
    return sigma, "transition"


# ---------------------------------------------------------------------
# Test cases — picked to exercise every regime
# ---------------------------------------------------------------------
TEST_POSITIONS = [
    # Bulk above and below the membrane (σ must equal bulk)
    ("bulk above membrane",                np.array([ 0.0, 0.0,  40.0]),  "bulk_or_pore"),
    ("bulk below membrane",                np.array([ 0.0, 0.0, -40.0]),  "bulk_or_pore"),
    ("on-axis at z=0",                     np.array([ 0.0, 0.0,   0.0]),  "bulk_or_pore"),
    ("on-axis at z=+15",                   np.array([ 0.0, 0.0,  15.0]),  "bulk_or_pore"),
    ("on-axis at z=-15",                   np.array([ 0.0, 0.0, -15.0]),  "bulk_or_pore"),

    # Asymmetry pair: same r=20, opposite z. Wide end → in pore; narrow → in wall.
    ("ASYM r=20  near +half  (inside wide mouth)",  np.array([20.0, 0.0,  18.0]), "bulk_or_pore"),
    ("ASYM r=20  near -half  (deep inside narrow wall)", np.array([20.0, 0.0, -18.0]), "wall_floor"),

    # Deep inside each mouth (well within the local opening)
    ("deep inside top mouth   (r=22, z=+18)", np.array([22.0, 0.0,  18.0]), "bulk_or_pore"),
    ("deep inside bottom mouth (r=4,  z=-18)", np.array([ 4.0, 0.0, -18.0]), "bulk_or_pore"),

    # Deep inside the wall, at each end
    ("deep wall  (r=40, z=+18)",  np.array([40.0, 0.0,  18.0]), "wall_floor"),
    ("deep wall  (r=20, z=-18)",  np.array([20.0, 0.0, -18.0]), "wall_floor"),

    # Slanted-wall midplane: at z=0 the local radius is (top+bottom)/2 = 20.
    ("midplane r=14 (deep inside)",  np.array([14.0, 0.0, 0.0]), "bulk_or_pore"),
    ("midplane r=26 (deep inside wall)", np.array([26.0, 0.0, 0.0]), "wall_floor"),
]


def main() -> int:
    pg = _load_pore_geometry()
    pore = _build_pore(pg)
    interp = pore.get_conductivity_interpolator()

    print("=" * 96)
    print(" Conical pore: validity of mesh-based overlap detection")
    print("=" * 96)
    print(
        f"Geometry:   top_radius={TOP_RADIUS} Å   "
        f"bottom_radius={BOTTOM_RADIUS} Å   half_thickness={HALF} Å"
    )
    print(
        f"Materials:  bulk σ = {BULK_SIGMA:.3f} S/m   "
        f"membrane σ = {MEMBRANE_SIGMA:.0e} S/m   "
        f"(condfrac floor ≈ {1e-7 * BULK_SIGMA:.2e} S/m)"
    )
    print(
        f"Grid:       {GRID_N}^3 over ±{GRID_HALF_EXTENT} Å, "
        f"spacing ≈ {2 * GRID_HALF_EXTENT / (GRID_N - 1):.2f} Å"
    )
    print(
        f"Bands:      bulk_or_pore if σ ≥ {SIGMA_BULK_THRESHOLD:.3f} S/m "
        f"(= 0.95 × bulk),"
    )
    print(
        f"            wall_floor   if σ ≤ {SIGMA_FLOOR_THRESHOLD:.3e} S/m "
        f"(= 1e-3 × bulk),"
    )
    print(
        f"            otherwise: condfrac transition zone "
        f"({PORE_BAND_DEPTH} Å either side of wall)"
    )
    print()
    print(
        f"{'#':>3}  {'case':<46}  {'pos (Å)':<22}  "
        f"{'σ@atom':>11}  {'expected':<14}  {'observed':<14}  result"
    )
    print("-" * 130)

    failures: list[str] = []
    for i, (label, xyz, expected) in enumerate(TEST_POSITIONS, start=1):
        sigma, observed = _classify_observed(interp, xyz)
        analytical = _classify_expected(xyz)

        # Both the analytical ground truth and the empirical readback
        # must agree. If either disagrees the mesh-based check is
        # operating on a wrong conductivity grid.
        ok = (analytical == expected) and (observed == expected)
        verdict = "PASS" if ok else "FAIL"
        if not ok:
            failures.append(
                f"#{i} {label}: expected={expected} "
                f"analytical={analytical} observed={observed} σ={sigma:.3e}"
            )

        pos_str = f"({xyz[0]:+.1f},{xyz[1]:+.1f},{xyz[2]:+.1f})"
        print(
            f"{i:>3}  {label:<46}  {pos_str:<22}  "
            f"{sigma:>11.3e}  {analytical:<14}  {observed:<14}  {verdict}"
        )

    print("-" * 130)

    # Asymmetry probe: same r=20, opposite z. The wide end should be at
    # bulk; the narrow end should be at the wall floor.
    plus = float(interp([[20.0, 0.0, +18.0]])[0])
    minus = float(interp([[20.0, 0.0, -18.0]])[0])
    ratio = plus / max(minus, 1e-30)
    print(
        f"\nAsymmetry probe at r=20:  σ(z=+18) = {plus:.3e} S/m,  "
        f"σ(z=-18) = {minus:.3e} S/m  →  ratio = {ratio:.3e}"
    )
    print(
        "  Symmetric DoubleConePore would give ratio=1 here. "
        "ratio ≫ 1 confirms the conical asymmetry is faithfully "
        "encoded in the conductivity grid the overlap check inspects."
    )

    # Slope probe along z at fixed r: σ should rise as z goes -half → +half.
    zs = np.linspace(-HALF + 1, HALF - 1, 9)
    print(
        f"\nσ(r=15, z) sweep across membrane thickness "
        f"(should be monotone increasing from wall→pore as z grows):"
    )
    last = -np.inf
    monotone = True
    for z in zs:
        sigma = float(interp([[15.0, 0.0, z]])[0])
        marker = " " if sigma >= last else " ← non-monotone!"
        if sigma < last:
            monotone = False
        print(f"    z = {z:+5.1f} Å    σ = {sigma:.3e} S/m{marker}")
        last = sigma
    if not monotone:
        failures.append("σ(r=15, z) is not monotone non-decreasing.")
    else:
        print("    monotone non-decreasing — consistent with frustum geometry.")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        print(f"\n{len(failures)}/{len(TEST_POSITIONS) + 1} checks failed.")
        return 1

    print(f"\nAll {len(TEST_POSITIONS)} classification cases + asymmetry + monotonicity PASS.")
    print(
        "The conductivity grid that vertical_movement_sem.py:1206-1225 reads is "
        "correctly built from the asymmetric conical geometry, so the mesh-based "
        "overlap check operates on faithful inputs for ConicalPore."
    )
    return 0


# ---------------------------------------------------------------------
# pytest entry points
# ---------------------------------------------------------------------
try:  # pragma: no cover
    import pytest  # type: ignore

    @pytest.fixture(scope="module")
    def conical_interp():
        pg = _load_pore_geometry()
        pore = _build_pore(pg)
        return pore.get_conductivity_interpolator()

    @pytest.mark.parametrize("label,xyz,expected", TEST_POSITIONS)
    def test_classification(conical_interp, label, xyz, expected):
        assert _classify_expected(xyz) == expected, (
            f"Analytical model disagrees with expectation for {label}"
        )
        _, observed = _classify_observed(conical_interp, xyz)
        assert observed == expected, (
            f"ConicalPore conductivity grid says {observed!r} for {label}, "
            f"expected {expected!r}"
        )

    def test_asymmetry(conical_interp):
        plus = float(conical_interp([[20.0, 0.0, +18.0]])[0])
        minus = float(conical_interp([[20.0, 0.0, -18.0]])[0])
        # Wide-end σ ≥ near-bulk; narrow-end σ ≤ floor band.
        assert plus >= SIGMA_BULK_THRESHOLD
        assert minus <= SIGMA_FLOOR_THRESHOLD
        assert plus / max(minus, 1e-30) > 1e6

    def test_monotonicity_along_z(conical_interp):
        zs = np.linspace(-HALF + 1, HALF - 1, 9)
        sigmas = [float(conical_interp([[15.0, 0.0, z]])[0]) for z in zs]
        # At r=15 we cross from inside-the-wall (z near -half) to
        # inside-the-pore (z near +half) → σ should be non-decreasing.
        diffs = np.diff(sigmas)
        assert np.all(diffs >= -1e-12), f"σ(z) not monotone: σ list = {sigmas}"
except ImportError:  # pragma: no cover
    pass


if __name__ == "__main__":
    raise SystemExit(main())
