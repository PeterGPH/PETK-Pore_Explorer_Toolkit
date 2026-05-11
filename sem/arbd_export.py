"""
ARBD-compatible grid export.

Exports the FEM-solved electrostatic potential phi (in volts) and derived
per-ion potential / steric grids as OpenDX (.dx) volumes that ARBD reads
natively (via Grid pair-potentials in the BD config file).

Conventions
-----------
* phi from dolfinx is in *volts* (V): voltage BCs are mV-converted-to-V upstream.
* Per-ion potential energy is U_ion(r) = z_ion * e * phi(r), then converted to
  kcal/mol via Faraday's constant:
      U[kcal/mol] = z_ion * 23.060548867 * phi[V]
  (F = 96485.33212 C/mol; 1 kcal = 4184 J)
* The steric grid is a hard-wall potential: ``wall_height`` (kcal/mol) wherever
  the membrane material lives (distance-to-membrane <= 0), 0 elsewhere.

Configuration block (under config["output"])
-------------------------------------------
::

    "arbd_export": {
        "ions":          [["POT", 1], ["CLA", -1]],
        "stride":        0,         # 0 = every step; N>0 = every N-th step
        "wall_height":   100.0,     # kcal/mol, steric exclusion height
        "temperature_K": 295.0,     # metadata only (not used in scaling)
        "include_phi":   true,      # also dump raw phi.dx (volts)
        "include_steric": true      # also dump the steric grid
    }

Each enabled SEM run produces, per Z position (and one "openpore" baseline):

* ``{prefix}_{pore_type}_z{+/-####.#A}_open_pore_phi.dx``  -- phi in V
* ``{prefix}_{pore_type}_z{+/-####.#A}_{ION}.dx``          -- per-ion U in kcal/mol
* ``{prefix}_{pore_type}_z{+/-####.#A}_steric.dx``         -- steric U in kcal/mol
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


# Faraday constant / (1 kcal in J) = 96485.33212 / 4184 = 23.06054886...
# Multiplies phi[V] * valence to give U_ion in kcal/mol.
_V_TO_KCALMOL_PER_E = 23.060548867


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _z_tag(z_pos: Optional[float]) -> str:
    """Filename tag for a Z position. ``None`` -> ``"openpore"``."""
    if z_pos is None:
        return "openpore"
    sign = "+" if z_pos >= 0 else "-"
    return f"z{sign}{abs(z_pos):07.1f}A"


def _build_grid_axes(sem_instance, custom_resolution=None):
    """Return (x_range, y_range, z_range, delta, origin, shape) for a regular grid
    that covers the dolfinx mesh domain.
    """
    if not (hasattr(sem_instance, "domain_min")
            and hasattr(sem_instance, "domain_max")
            and hasattr(sem_instance, "num_cells")):
        raise RuntimeError(
            "ARBD export needs sem_instance.domain_min/max/num_cells to be set. "
            "These are populated when the dolfinx mesh is built."
        )

    if custom_resolution is not None:
        num_cells = [
            int(round((sem_instance.domain_max[i] - sem_instance.domain_min[i])
                      / custom_resolution))
            for i in range(3)
        ]
        delta = [float(custom_resolution)] * 3
    else:
        num_cells = list(sem_instance.num_cells)
        delta = [
            (sem_instance.domain_max[i] - sem_instance.domain_min[i]) / num_cells[i]
            for i in range(3)
        ]
    num_points = [n + 1 for n in num_cells]

    x_range = np.linspace(sem_instance.domain_min[0], sem_instance.domain_max[0], num_points[0])
    y_range = np.linspace(sem_instance.domain_min[1], sem_instance.domain_max[1], num_points[1])
    z_range = np.linspace(sem_instance.domain_min[2], sem_instance.domain_max[2], num_points[2])

    origin = (float(sem_instance.domain_min[0]),
              float(sem_instance.domain_min[1]),
              float(sem_instance.domain_min[2]))
    shape = (num_points[0], num_points[1], num_points[2])
    return x_range, y_range, z_range, delta, origin, shape


def _sample_dolfinx_function_on_grid(uh, x_range, y_range, z_range, fill_value=0.0):
    """Evaluate a dolfinx Function ``uh`` on the regular grid x*y*z.

    Points outside the mesh receive ``fill_value`` (default 0.0 V).
    Returns a numpy array with shape ``(nx, ny, nz)``.
    """
    # Import locally so the module is still importable when dolfinx is absent
    # (e.g. for the lazy public API + compare-only consumers).
    from dolfinx import geometry as dgeom  # type: ignore

    nx, ny, nz = len(x_range), len(y_range), len(z_range)
    X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]).astype(np.float64)

    mesh = uh.function_space.mesh
    bb_tree = dgeom.bb_tree(mesh, mesh.topology.dim)
    cell_candidates = dgeom.compute_collisions_points(bb_tree, points)
    colliding_cells = dgeom.compute_colliding_cells(mesh, cell_candidates, points)

    # Pick the first colliding cell per point; mark points with no cell as missing.
    cells = np.full(len(points), -1, dtype=np.int32)
    for i in range(len(points)):
        link = colliding_cells.links(i)
        if len(link) > 0:
            cells[i] = link[0]

    valid_mask = cells >= 0
    out = np.full(len(points), fill_value, dtype=np.float64)
    if valid_mask.any():
        valid_points = points[valid_mask]
        valid_cells = cells[valid_mask]
        values = uh.eval(valid_points, valid_cells).reshape(-1)
        out[valid_mask] = values
    n_missing = int((~valid_mask).sum())
    if n_missing:
        logger.debug("ARBD export: %d/%d grid points fell outside the mesh (filled with %.3g).",
                     n_missing, len(points), fill_value)

    return out.reshape((nx, ny, nz))


def _write_dx(filename, grid_data, delta, origin):
    """Write a 3D numpy array as an OpenDX file via the ``gridData`` package."""
    try:
        from gridData import Grid  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "ARBD export requires the 'gridData' package: pip install gridData"
        ) from exc

    g = Grid(grid=grid_data, delta=delta, origin=origin)
    g.export(filename)
    logger.info("Wrote DX grid: %s  shape=%s  delta=%s",
                filename, grid_data.shape, tuple(float(d) for d in delta))


# --------------------------------------------------------------------------- #
# Public exports
# --------------------------------------------------------------------------- #


def export_phi_dx(sem_instance, uh, filename, custom_resolution=None):
    """Sample phi (volts) on the mesh-spanning regular grid and dump to DX."""
    x_range, y_range, z_range, delta, origin, _ = _build_grid_axes(
        sem_instance, custom_resolution
    )
    phi_grid_V = _sample_dolfinx_function_on_grid(uh, x_range, y_range, z_range, fill_value=0.0)
    _write_dx(filename, phi_grid_V, delta, origin)
    if sem_instance.rank == 0:
        logger.info("  phi range: [%.4g, %.4g] V", float(phi_grid_V.min()), float(phi_grid_V.max()))
    return phi_grid_V, delta, origin


def export_ion_potential_dx(phi_grid_V, valence, filename, delta, origin):
    """U_ion(r) = z * F * phi(r) -> kcal/mol, written as DX."""
    u_kcalmol = valence * _V_TO_KCALMOL_PER_E * phi_grid_V
    _write_dx(filename, u_kcalmol, delta, origin)
    return u_kcalmol


def export_steric_dx(sem_instance, wall_height, filename, custom_resolution=None):
    """Hard-wall steric potential: ``wall_height`` kcal/mol inside the membrane,
    0 elsewhere.

    Uses the pore object's distance-to-membrane field when available, otherwise
    falls back to a (membrane_half_thickness, R > pore_radius) test for cylindrical.
    """
    x_range, y_range, z_range, delta, origin, _ = _build_grid_axes(
        sem_instance, custom_resolution
    )
    X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
    grid_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    inside_membrane = _membrane_mask(sem_instance, grid_points).reshape(X.shape)
    steric = np.where(inside_membrane, float(wall_height), 0.0).astype(np.float64)
    _write_dx(filename, steric, delta, origin)
    return steric


def _membrane_mask(sem_instance, grid_points):
    """Boolean mask: True where the membrane material is.

    Uses ``base_dist_interp`` when present (treats distance<=0 as inside wall);
    otherwise builds an analytic mask from pore geometry parameters.
    """
    if getattr(sem_instance, "base_dist_interp", None) is not None:
        d = sem_instance.base_dist_interp(grid_points)
        d = np.asarray(d).ravel()
        # If the interpolator returns NaN for outside-mesh points, treat as bulk (no wall).
        d = np.where(np.isnan(d), 1.0, d)
        return d <= 0.0

    pore_type = getattr(sem_instance, "pore_type", "cylindrical")
    half_t = float(sem_instance.membrane_thickness) / 2.0
    if half_t <= 0:
        return np.zeros(len(grid_points), dtype=bool)

    R = np.sqrt(grid_points[:, 0] ** 2 + grid_points[:, 1] ** 2)
    abs_z = np.abs(grid_points[:, 2])
    in_slab = abs_z <= half_t

    if pore_type == "cylindrical":
        return in_slab & (R > float(sem_instance.pore_radius))

    if pore_type == "double_cone":
        z_frac = np.clip(abs_z / half_t, 0.0, 1.0)
        local_r = sem_instance.pore_radius + (sem_instance.outer_radius - sem_instance.pore_radius) * z_frac
        return in_slab & (R > local_r)

    if pore_type == "conical":
        signed_z = grid_points[:, 2]
        t = np.clip((signed_z + half_t) / (2.0 * half_t), 0.0, 1.0)
        local_r = sem_instance.bottom_radius + (sem_instance.top_radius - sem_instance.bottom_radius) * t
        return in_slab & (R > local_r)

    # Fallback: no steric region known.
    return np.zeros(len(grid_points), dtype=bool)


# --------------------------------------------------------------------------- #
# Orchestrator -- called once per (z-position or openpore)
# --------------------------------------------------------------------------- #


def export_arbd_grids(sem_instance, uh, z_pos, arbd_cfg):
    """Write phi + per-ion + steric DX grids for one solve step.

    Parameters
    ----------
    sem_instance : VerticalMovementSEM
        Provides ``output_prefix``, ``pore_type``, mesh domain, distance interpolator.
    uh : dolfinx.fem.Function
        Solved electrostatic potential (volts).
    z_pos : float | None
        Analyte z-position in angstroms (used for filename tag). Pass ``None``
        for the open-pore baseline.
    arbd_cfg : dict
        The ``output.arbd_export`` config block.
    """
    if not arbd_cfg:
        return None
    if sem_instance.rank != 0:
        # Only rank 0 writes files. Other ranks still participate in uh.eval via the
        # caller, but we keep DX writes serial.
        return None

    tag = _z_tag(z_pos)
    prefix = f"{sem_instance.output_prefix}_{sem_instance.pore_type}_{tag}"
    custom_res = arbd_cfg.get("resolution", None)

    include_phi = bool(arbd_cfg.get("include_phi", True))
    include_steric = bool(arbd_cfg.get("include_steric", True))
    ions = arbd_cfg.get("ions", []) or []

    phi_path = f"{prefix}_open_pore_phi.dx"
    logger.info("ARBD export @ %s -> %s", tag, prefix)

    # phi.dx is needed either way (to derive per-ion grids).
    phi_grid_V, delta, origin = export_phi_dx(sem_instance, uh, phi_path, custom_res)

    if not include_phi:
        try:
            Path(phi_path).unlink()
        except OSError:
            pass

    for ion in ions:
        try:
            name, valence = ion[0], float(ion[1])
        except (TypeError, ValueError, IndexError):
            logger.warning("Skipping malformed ARBD ion spec: %r", ion)
            continue
        ion_path = f"{prefix}_{name}.dx"
        export_ion_potential_dx(phi_grid_V, valence, ion_path, delta, origin)

    if include_steric:
        wall_h = float(arbd_cfg.get("wall_height", 100.0))
        steric_path = f"{prefix}_steric.dx"
        export_steric_dx(sem_instance, wall_h, steric_path, custom_res)

    return prefix


def should_export_at_step(step_index, arbd_cfg):
    """Honor the ``stride`` setting: 0 = every step, N>0 = every N-th step."""
    if not arbd_cfg:
        return False
    stride = int(arbd_cfg.get("stride", 0) or 0)
    if stride <= 0:
        return True
    return (step_index % stride) == 0
