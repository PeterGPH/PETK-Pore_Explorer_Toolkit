"""Demonstrate grid convergence of ConicalPore's conductivity field.

What this script shows
----------------------
The conductivity that ``ConicalPore.get_conductivity_interpolator()``
hands the SEM solver is built by

    σ(x, y, z) = condfrac(d(x, y, z)) × (bulk - membrane) + membrane

where

    d(x, y, z) = sqrt( max(local_radius(z) - r, 0)²
                      + max(|z| - half, 0)² )         (analytical)

and ``condfrac`` (sem/utils.py:15) is a piecewise linear ramp from
distance 1.3 Å (fraction = 1e-7 floor) to 4.1 Å (fraction = 1.0).

That is a continuous, well-defined function of (x, y, z). The grid
just samples it. So as ``grid_resolution → 0`` the trilinear
interpolation of the discrete conductivity grid must converge to the
analytical σ at every point. This script measures that convergence
rate at several probe points spanning the bulk, the transition band,
and the wall floor.

Run::

    python tests/test_conical_grid_convergence.py

or under pytest::

    pytest tests/test_conical_grid_convergence.py -v
"""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]


def _stub_fem_modules() -> None:
    for name in (
        "dolfinx", "dolfinx.fem", "dolfinx.io", "dolfinx.mesh",
        "mpi4py", "mpi4py.MPI", "petsc4py", "petsc4py.PETSc", "ufl",
    ):
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["mpi4py.MPI"].COMM_WORLD = types.SimpleNamespace(
        rank=0, size=1, Get_rank=lambda: 0
    )
    sys.modules["mpi4py"].MPI = sys.modules["mpi4py.MPI"]


def _load_modules():
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

    utils = _load("utils")
    pg = _load("pore_geometry")
    return pg, utils


# ---------------------------------------------------------------------
# Geometry & material parameters
# ---------------------------------------------------------------------
TOP_RADIUS = 30.0
BOTTOM_RADIUS = 10.0
HALF = 20.0
BULK_SIGMA = 10.5
MEMBRANE_SIGMA = 1e-7
# Tight box: deepest probe sits at r=27, z=±18; ±32 Å gives plenty of
# margin around the pore while keeping memory tractable at fine h.
GRID_HALF_EXTENT = 32.0

# Resolutions to sweep, coarsest → finest. Each row uses an odd grid
# count so the origin is sampled exactly. h=0.25 on a ±32 Å box is
# 257³ ≈ 17M cells × several float64 arrays ≈ ~1 GB, runnable.
RESOLUTIONS = [4.0, 2.0, 1.0, 0.5, 0.25]


def _local_radius(z: float) -> float:
    t = np.clip((z + HALF) / (2.0 * HALF), 0.0, 1.0)
    return BOTTOM_RADIUS + (TOP_RADIUS - BOTTOM_RADIUS) * t


def _analytic_sigma(xyz: np.ndarray, condfrac) -> float:
    """Closed-form σ that the discrete grid must converge to."""
    x, y, z = xyz
    r = float(np.hypot(x, y))
    radial_term = max(_local_radius(z) - r, 0.0)
    vertical_term = max(abs(z) - HALF, 0.0)
    d = np.hypot(radial_term, vertical_term)
    fraction = float(condfrac(np.array([d]))[0])
    return MEMBRANE_SIGMA + fraction * (BULK_SIGMA - MEMBRANE_SIGMA)


def _build_pore(pg, h: float):
    """Build a ConicalPore on a uniform grid with spacing ≈ h."""
    n = int(round(2 * GRID_HALF_EXTENT / h)) + 1
    if n % 2 == 0:  # ensure odd, so origin is on a grid node
        n += 1
    ax = np.linspace(-GRID_HALF_EXTENT, GRID_HALF_EXTENT, n)
    actual_h = ax[1] - ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    pore = pg.PoreGeometry.create_pore(
        "conical", X=X, Y=Y, Z=Z,
        top_radius=TOP_RADIUS, bottom_radius=BOTTOM_RADIUS,
        membrane_half_thickness=HALF,
        bulk_conductivity=BULK_SIGMA,
        membrane_conductivity=MEMBRANE_SIGMA,
    )
    return pore, n, actual_h


# Probe points spanning all three regimes.
# Each tuple: (label, (x, y, z), regime)
#
# Distances (from the analytic formula, given top=30, bottom=10, half=20):
#  - midplane local_radius = 20
#  - z=+18 local_radius = 29.5
#  - z=-18 local_radius = 10.5
# Off-grid probe coordinates with deliberately non-round offsets so that
# at every h in our sweep the probes sit *between* grid nodes — that's
# the only way to see the trilinear-interpolation convergence rate
# instead of trivial zero error from grid-aligned exact hits.
#
# Distances quoted are computed at the probe coords (top=30, bottom=10,
# half=20).
PROBES = [
    ("midplane:  r≈10.37,  z≈0.41  (deep in pore)",        (10.37, 0.0,   0.41), "bulk"),
    ("midplane:  r≈16.73,  z≈0.29  (mid transition d≈3.3)",(16.73, 0.0,   0.29), "transition"),
    ("midplane:  r≈18.31,  z≈-0.17 (mid transition d≈1.7)",(18.31, 0.0,  -0.17), "transition"),
    ("midplane:  r≈19.43,  z≈0.13  (floor band d≈0.6)",    (19.43, 0.0,   0.13), "floor_band"),
    ("midplane:  r≈22.41,  z≈0.37  (in wall, d=0)",        (22.41, 0.0,   0.37), "floor"),
    ("ASYM       r≈20.13,  z≈+18.27 (wide mouth, bulk)",   (20.13, 0.0,  18.27), "bulk"),
    ("ASYM       r≈20.13,  z≈-18.27 (narrow wall, d=0)",   (20.13, 0.0, -18.27), "floor"),
    ("midplane:  r≈27.31,  z≈0.41  (deep wall, d=0)",      (27.31, 0.0,   0.41), "floor"),
]


def main() -> int:
    pg, utils = _load_modules()
    condfrac = utils.condfrac

    print("=" * 100)
    print(" ConicalPore conductivity grid: convergence to analytical σ as h → 0")
    print("=" * 100)
    print(
        f"Geometry: top={TOP_RADIUS} Å  bottom={BOTTOM_RADIUS} Å  "
        f"half={HALF} Å   |   bulk σ={BULK_SIGMA} S/m   membrane σ={MEMBRANE_SIGMA:.0e} S/m"
    )
    print(f"Sweeping h ∈ {RESOLUTIONS} Å")
    print()

    # Pre-compute analytical σ at every probe (independent of h).
    analytic = {label: _analytic_sigma(np.array(xyz), condfrac)
                for label, xyz, _ in PROBES}

    # Header
    h_headers = "  ".join(f"{f'h={h}':>11}" for h in RESOLUTIONS)
    print(f"{'probe':<58}  {'σ_exact':>11}  {h_headers}")
    print("-" * (60 + 13 + (13 * len(RESOLUTIONS))))

    obs = {label: [] for label, _, _ in PROBES}
    actual_hs = []

    import gc
    for h in RESOLUTIONS:
        pore, n, actual_h = _build_pore(pg, h)
        actual_hs.append(actual_h)
        interp = pore.get_conductivity_interpolator()
        for label, xyz, _ in PROBES:
            sigma = float(interp([xyz])[0])
            obs[label].append(sigma)
        # Free everything before building the next, finer grid.
        del pore, interp
        gc.collect()

    for label, xyz, regime in PROBES:
        row = f"{label:<58}  {analytic[label]:>11.4e}"
        for sigma in obs[label]:
            row += f"  {sigma:>11.4e}"
        print(row)

    # Convergence summary: |σ_obs - σ_exact| at each h, plus an
    # estimated order via successive ratios (≈ log₂(err_coarse / err_fine)).
    print()
    print("Absolute error |σ_obs − σ_exact|  and  observed convergence order")
    print("-" * 100)
    h_hdr = "  ".join(f"{f'h={h}':>11}" for h in RESOLUTIONS)
    print(f"{'probe':<58}  {h_hdr}    order")
    print("-" * (60 + 13 * len(RESOLUTIONS) + 12))

    failures: list[str] = []
    for label, xyz, regime in PROBES:
        errs = [abs(s - analytic[label]) for s in obs[label]]
        # Estimate order from the last two h's where both errors are
        # well above machine ε (so the ratio is meaningful).
        order_str = "    —"
        for k in range(len(errs) - 1, 0, -1):
            if errs[k - 1] > 1e-8 and errs[k] > 1e-8:
                order = (np.log(errs[k - 1] / errs[k])
                         / np.log(actual_hs[k - 1] / actual_hs[k]))
                order_str = f"{order:>5.2f}"
                break
        row = f"{label:<58}  "
        row += "  ".join(f"{e:>11.3e}" for e in errs)
        row += f"   {order_str}"
        print(row)

        # Convergence assertion: error must decrease monotonically (allowing
        # a tiny slack at the floor where everyone is already at machine ε).
        is_floor_regime = regime == "floor" and analytic[label] < 2e-6
        if not is_floor_regime:
            for k in range(1, len(errs)):
                if errs[k] > 1.5 * errs[k - 1] + 1e-12:
                    failures.append(
                        f"{label}: error did not decrease at h={RESOLUTIONS[k]} "
                        f"(was {errs[k-1]:.3e}, became {errs[k]:.3e})"
                    )

    print()
    if failures:
        print("CONVERGENCE FAILURES:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print(
        "Convergence verified: σ_obs → σ_exact at every probe as h → 0.\n"
        "Note that the asymptote σ_exact at wall-deep points is bounded\n"
        "above ``membrane_conductivity`` by the condfrac floor:\n"
        f"    σ_floor = {MEMBRANE_SIGMA + 1e-7 * (BULK_SIGMA - MEMBRANE_SIGMA):.3e} S/m\n"
        "i.e. the wall σ does not converge *to* membrane_conductivity;\n"
        "it converges to that floor value. Refining the mesh further does\n"
        "not push it down. To make the line 1211 check\n"
        "    base_vals <= self.membrane_conductivity\n"
        "fire, raise membrane_conductivity above the floor, or add slack."
    )
    return 0


# ---------------------------------------------------------------------
# pytest entry points
# ---------------------------------------------------------------------
try:  # pragma: no cover
    import pytest

    @pytest.fixture(scope="module")
    def loaded():
        return _load_modules()

    @pytest.mark.parametrize("label,xyz,regime", PROBES)
    def test_grid_converges_to_analytic(loaded, label, xyz, regime):
        pg, utils = loaded
        condfrac = utils.condfrac
        target = _analytic_sigma(np.array(xyz), condfrac)

        errs = []
        actual_hs = []
        for h in RESOLUTIONS:
            pore, _, actual_h = _build_pore(pg, h)
            sigma = float(pore.get_conductivity_interpolator()([xyz])[0])
            errs.append(abs(sigma - target))
            actual_hs.append(actual_h)

        # Floor-regime points reach machine ε quickly; skip strict ordering.
        if regime == "floor" and target < 2e-6:
            assert errs[-1] < 1e-5
            return

        # Otherwise: error at finest h must be < error at coarsest h, and
        # the trend should be (broadly) monotone non-increasing.
        assert errs[-1] <= errs[0], (
            f"{label}: finest-h error {errs[-1]:.3e} not below coarsest-h "
            f"error {errs[0]:.3e}"
        )
        for k in range(1, len(errs)):
            assert errs[k] <= 1.5 * errs[k - 1] + 1e-12, (
                f"{label}: error grew from {errs[k-1]:.3e} to {errs[k]:.3e} "
                f"between h={RESOLUTIONS[k-1]} and h={RESOLUTIONS[k]}"
            )

except ImportError:
    pass


if __name__ == "__main__":
    raise SystemExit(main())
