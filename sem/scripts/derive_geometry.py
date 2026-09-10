#!/usr/bin/env python3
"""
Derive a parametric `pore_geometry` block (or a radius-vs-z `profile`
table) from an all-atom (`gen_dist`) conductivity/distance map.

The all-atom path measures a true 3-D distance-to-solid field on a grid and
converts it to a conductivity fraction with the `condfrac` ramp
(`sem.grid_io`). The parametric path builds an analytic wall (`PoreProfile`
in `sem.geometry_profiles`) and applies the same ramp. This module inverts
that process: given an all-atom `.bin` map, it recovers the ramp-weighted
equivalent radius at every z-slice (`R_map`), fits it to the requested
parametric pore type (or hands back the raw radius-vs-z table for the
`profile` pore type), and refines the fit against the true Euclidean
distance metric with `scipy.optimize.least_squares`.

Only numpy, scipy, `sem.grid_io` and `sem.geometry_profiles` are imported
(no dolfinx, no `sem.utils`), so this module -- and the standalone CLI it
exposes -- works in pip-only environments that never install the
DOLFINx/MPI stack.

Library entry point: `derive(bin_path, pore_type, **opts) -> DerivedGeometry`.
CLI: `python -m sem.scripts.derive_geometry BIN --pore-type ... `.
"""

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import least_squares

# sem.grid_io / sem.geometry_profiles are numpy-only siblings of this
# module. This is normally imported as part of the `sem` package
# (`sem.scripts.derive_geometry`), but fall back to an absolute import so
# `python sem/scripts/derive_geometry.py` also works standalone (matching
# the convention in sem/scripts/gen_dist.py).
try:
    from ..grid_io import (
        RAMP_MAX,
        RAMP_MIN,
        area_coefficients,
        bin_axes,
        condfrac,
        radius_from_area,
        readbinGrid,
    )
except ImportError:  # pragma: no cover - relative import fallback
    from sem.grid_io import (
        RAMP_MAX,
        RAMP_MIN,
        area_coefficients,
        bin_axes,
        condfrac,
        radius_from_area,
        readbinGrid,
    )

try:
    from ..geometry_profiles import PoreProfile
except ImportError:  # pragma: no cover - relative import fallback
    from sem.geometry_profiles import PoreProfile

# `sem.provenance` (git_commit) and `sem.__version__` are dolfinx-free, so
# importing them here does not compromise this module's own dolfinx-free
# guarantee (see the module docstring / test_module_source_has_no_dolfinx_or_utils_import).
try:
    from .. import __version__ as _SEM_VERSION
    from ..provenance import git_commit as _git_commit
except ImportError:  # pragma: no cover - relative import fallback
    from sem import __version__ as _SEM_VERSION
    from sem.provenance import git_commit as _git_commit


PARAMETRIC_PORE_TYPES = ("cylindrical", "double_cone", "conical")
PORE_TYPES = PARAMETRIC_PORE_TYPES + ("profile",)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------
@dataclass
class SliceProfile:
    """
    Per-z-slice quantities measured from an all-atom fraction field.
    `R_map`, `area_param`, `r_fit` and `used_in_fit` are filled in later
    (once the parametric fit / profile table is known) so this can double
    as the row source for `P_profile.csv`.
    """

    z: np.ndarray
    area_bin: np.ndarray
    area_pore: np.ndarray
    f_far: np.ndarray
    r_eq: np.ndarray
    spacing: np.ndarray
    far_radius: float
    R_map: Optional[np.ndarray] = None
    area_param: Optional[np.ndarray] = None
    r_fit: Optional[np.ndarray] = None
    used_in_fit: Optional[np.ndarray] = None


@dataclass
class DerivedGeometry:
    """Library-level result of `derive()`."""

    pore_geometry: dict
    profile_table: tuple
    slices: SliceProfile
    derivation: dict


# ---------------------------------------------------------------------------
# Loading and slicing
# ---------------------------------------------------------------------------
def load_fraction_field(path, units="distance", bulk_conductivity=None):
    """
    Load a `.bin` grid and convert it to a conductivity fraction field
    `f in [0, 1]` (0 = solid, 1 = bulk).

    `units="distance"`: the bin holds a distance-to-solid-surface field;
    `f = condfrac(distance)`.
    `units="conductivity"`: the bin holds conductivity values `sigma`;
    `f = clip((sigma - sigma_min) / (sigma_bulk - sigma_min), 0, 1)` where
    `sigma_bulk = bulk_conductivity` if given else `max(sigma)`, and
    `sigma_min = min(sigma)`.

    Returns `(f3d, axes, spacing, meta)`: `f3d` shaped `(nx, ny, nz)`,
    `axes = (x, y, z)` 1-D coordinate arrays, `spacing` the `(hx, hy, hz)`
    grid spacing, and `meta` a dict carrying the raw grid metadata plus
    the resolved `units`/`bulk_conductivity` (and `sigma_min` for
    conductivity units).
    """
    val3d, _dims, shape, metadata = readbinGrid(path, return_metadata=True)
    val3d = val3d.astype(np.float64)
    axes = bin_axes(metadata, shape)
    spacing = np.asarray(metadata["spacing"], dtype=float)

    units = (units or "distance").lower()
    meta = dict(metadata)
    meta["units"] = units

    if units == "distance":
        f3d = condfrac(val3d)
        meta["bulk_conductivity"] = bulk_conductivity
    elif units == "conductivity":
        sigma_min = float(val3d.min())
        sigma_bulk = float(bulk_conductivity) if bulk_conductivity is not None else float(val3d.max())
        denom = sigma_bulk - sigma_min
        if denom <= 0:
            raise ValueError(
                f"Cannot derive a fraction field: bulk conductivity ({sigma_bulk}) "
                f"must exceed the minimum value in the bin ({sigma_min})"
            )
        f3d = np.clip((val3d - sigma_min) / denom, 0.0, 1.0)
        meta["bulk_conductivity"] = sigma_bulk
        meta["sigma_min"] = sigma_min
    else:
        raise ValueError(f"Unknown units {units!r}; expected 'distance' or 'conductivity'")

    return f3d, axes, spacing, meta


def slice_profile(f3d, axes, spacing, far_radius=None):
    """
    Reduce a 3-D fraction field to per-z-slice quantities: bin area
    `A_bin(z) = h**2 * sum(f)`, far-field fraction `f_far(z)` (mean of `f`
    over `R > far_radius`, default `0.9 * half box width`), pore area
    `A_pore(z) = h**2 * sum((f - f_far) / (1 - f_far))` (the far-field
    baseline subtracted off, for face-band flare), and the naive
    equivalent radius `r_eq(z) = sqrt(A_pore(z) / pi)`.
    """
    x, y, z = axes
    hx, hy, hz = float(spacing[0]), float(spacing[1]), float(spacing[2])
    pix_area = hx * hy

    X, Y = np.meshgrid(x, y, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2)

    half_width_x = (float(x.max()) - float(x.min())) / 2.0
    half_width_y = (float(y.max()) - float(y.min())) / 2.0
    half_box_width = min(half_width_x, half_width_y)
    if far_radius is None:
        far_radius = 0.9 * half_box_width

    far_mask = R > far_radius
    n_far = int(far_mask.sum())

    nz = f3d.shape[2]
    area_bin = np.empty(nz)
    area_pore = np.empty(nz)
    f_far = np.empty(nz)

    for k in range(nz):
        f_slice = f3d[:, :, k]
        area_bin[k] = pix_area * float(f_slice.sum())
        f_far_k = float(f_slice[far_mask].mean()) if n_far > 0 else 1.0
        f_far[k] = f_far_k
        if f_far_k < 0.999:
            area_pore[k] = pix_area * float(np.sum((f_slice - f_far_k) / (1.0 - f_far_k)))
        else:
            area_pore[k] = area_bin[k]

    r_eq = np.sqrt(np.maximum(area_pore, 0.0) / np.pi)

    return SliceProfile(
        z=np.asarray(z, dtype=float),
        area_bin=area_bin,
        area_pore=area_pore,
        f_far=f_far,
        r_eq=r_eq,
        spacing=np.array([hx, hy, hz]),
        far_radius=float(far_radius),
    )


def _crossings(z, values, level):
    """z positions (linearly interpolated) where `values` crosses `level`."""
    diffs = np.asarray(values, dtype=float) - level
    signs = np.sign(diffs)
    idx = np.where(np.diff(signs) != 0)[0]
    result = []
    for i in idx:
        z0, z1 = z[i], z[i + 1]
        v0, v1 = diffs[i], diffs[i + 1]
        if v1 == v0:
            continue
        t = -v0 / (v1 - v0)
        result.append(float(z0 + t * (z1 - z0)))
    return result


def membrane_extent(prof):
    """
    Estimate the membrane thickness two ways: `L_integral` (exact for a
    flat face under the condfrac ramp) from integrating `1 - f_far(z)`
    over the whole z range, and `L_crossing` from the distance between the
    two z where `f_far` crosses 0.5. Also returns `z_center`, the midpoint
    of those two crossings.
    """
    hz = float(prof.spacing[2])
    L_integral = float(np.sum(1.0 - prof.f_far) * hz - (RAMP_MIN + RAMP_MAX))

    crossings = _crossings(prof.z, prof.f_far, 0.5)
    if len(crossings) < 2:
        raise ValueError(
            "Could not find two f_far=0.5 crossings along z; the bin must "
            "contain both membrane (f_far < 0.5) and bulk (f_far > 0.5) slices"
        )
    z_lo, z_hi = crossings[0], crossings[-1]
    L_crossing = float((z_hi - z_lo) - (RAMP_MIN + RAMP_MAX))
    z_center = float((z_lo + z_hi) / 2.0)
    return L_integral, L_crossing, z_center


# ---------------------------------------------------------------------------
# Ramp-weighted radius mapping
# ---------------------------------------------------------------------------
def map_radius_from_area(area, cos_theta=1.0):
    """
    Exact slope-aware inverse of the ramp-weighted equivalent area (see
    `sem.grid_io.radius_from_area`), with the input area clamped to the
    quadratic's valid domain so noisy/near-zero areas never produce a
    negative discriminant.
    """
    c1, c0 = area_coefficients(cos_theta)
    floor_area = np.pi * (c0 - c1 * c1 / 4.0)
    safe_area = np.maximum(np.asarray(area, dtype=float), floor_area + 1e-9)
    return radius_from_area(safe_area, cos_theta)


def _moving_average3(a):
    a = np.asarray(a, dtype=float)
    if a.size < 3:
        return a.copy()
    out = a.copy()
    out[1:-1] = (a[:-2] + a[1:-1] + a[2:]) / 3.0
    return out


def _compute_R_map(area, z, iterations=2, window=3):
    """
    `R_map(z) = map_radius_from_area(A_pore(z), cos_theta)` where
    `cos_theta = 1 / sqrt(1 + (dr/dz)**2)` is derived from the local slope
    of a smoothed `R_map(z)`, iterated `iterations` times starting from
    `cos_theta = 1`. For a plain cylinder this converges immediately
    (`cos_theta = 1` exactly, zero slope).
    """
    area = np.asarray(area, dtype=float)
    z = np.asarray(z, dtype=float)
    cos_theta = np.ones_like(area)
    R_map = map_radius_from_area(area, cos_theta)
    for _ in range(iterations):
        R_smooth = _moving_average3(R_map) if window == 3 else R_map
        slope = np.gradient(R_smooth, z) if z.size >= 2 else np.zeros_like(R_map)
        cos_theta = 1.0 / np.sqrt(1.0 + slope ** 2)
        R_map = map_radius_from_area(area, cos_theta)
    return R_map


def profile_from_slices(prof, L, z_center, *, smooth=3):
    """
    The `profile` pore-type output: `(z, r)` restricted to
    `|z - z_center| <= L / 2`, with the end points set exactly to `+-L/2`
    by linear interpolation and `z` expressed relative to `z_center` (so
    it feeds `PoreProfile.from_table` directly). `r` is `R_map(z)` after
    the slope iteration.
    """
    half_L = L / 2.0
    z_lo = z_center - half_L
    z_hi = z_center + half_L

    # Compute R_map over the *whole* z axis (not just the restricted
    # window) so np.interp has real bracketing samples on both sides of
    # +-L/2 -- restricting first would silently clamp (rather than
    # interpolate) the end points whenever +-L/2 falls between grid
    # points, biasing the mouth radius low.
    R_map_full = _compute_R_map(prof.area_pore, prof.z, window=smooth)

    mask = (prof.z >= z_lo) & (prof.z <= z_hi)
    if mask.sum() < 2:
        raise ValueError(
            "Not enough membrane slices within the fitted thickness to build a profile table"
        )
    z_sel = prof.z[mask]
    R_sel = R_map_full[mask]

    r_lo = float(np.interp(z_lo, prof.z, R_map_full))
    r_hi = float(np.interp(z_hi, prof.z, R_map_full))

    # The end points are set to exactly +-half_L (not `z_lo/z_hi - z_center`,
    # which -- since z_lo/z_hi were themselves built as `z_center -+ half_L`
    # -- can drift from +-half_L by a floating-point rounding ULP or two
    # once z_center != 0). `half_L = L / 2.0` here is bit-identical to
    # `membrane_thickness / 2.0` wherever the caller reports `L` as
    # `membrane_thickness`, so the emitted table's span matches the
    # emitted pore_geometry block's half_thickness by construction rather
    # than by chance -- this is what PoreProfile.from_table's +-h span
    # check (and hence validate_config) relies on.
    interior = (z_sel > z_lo) & (z_sel < z_hi)
    z_out = np.concatenate(([-half_L], z_sel[interior] - z_center, [half_L]))
    r_out = np.concatenate(([r_lo], R_sel[interior], [r_hi]))
    return z_out, r_out


# ---------------------------------------------------------------------------
# Parametric fits
# ---------------------------------------------------------------------------
def fit_cylindrical(prof, L, z_center, *, face_exclude=5.0, fit_corner=False):
    """`R = mean(R_map)` over the core (membrane slices, `face_exclude`
    Angstrom excluded at each face). With `fit_corner`, also estimate
    `corner_radius`/`chamfer_depth` from the face-band flare (both are
    only initial guesses -- `refine` does the precision work)."""
    half_L = L / 2.0
    dz = np.abs(prof.z - z_center)
    membrane = prof.f_far < 0.5
    core = membrane & (dz <= half_L - face_exclude)
    if core.sum() < 2:
        raise ValueError("Not enough core membrane slices to fit a cylinder (check face_exclude)")

    R_map_core = _compute_R_map(prof.area_pore[core], prof.z[core])
    R = float(np.mean(R_map_core))

    warnings_out = []
    core_slope = 0.0
    if prof.z[core].size >= 2:
        core_slope = float(np.polyfit(prof.z[core] - z_center, R_map_core, 1)[0])
        if abs(core_slope) > 0.02:
            warnings_out.append(
                f"cylinder core slope {core_slope:.4f} exceeds 0.02; "
                "wall may not be a plain cylinder"
            )

    result = {
        "pore_radius": R,
        "core_mask": core,
        "core_slope": core_slope,
        "warnings": warnings_out,
    }

    if fit_corner:
        face_band = membrane & (dz > half_L - face_exclude) & (dz <= half_L)
        if face_band.sum() >= 1:
            R_map_face = map_radius_from_area(prof.area_pore[face_band])
            corner_radius = float(np.mean(R_map_face)) - R

            depth = half_L - dz[face_band]
            order = np.argsort(depth)
            depth_sorted = depth[order]
            R_face_sorted = R_map_face[order]
            exceeds = R_face_sorted > (R + 1.0)
            if not exceeds.any():
                chamfer_depth = 0.0
            elif exceeds.all():
                chamfer_depth = float(depth_sorted[-1])
            else:
                first_false = int(np.argmax(~exceeds))
                if first_false == 0:
                    chamfer_depth = 0.0
                else:
                    d0, d1 = depth_sorted[first_false - 1], depth_sorted[first_false]
                    r0, r1 = R_face_sorted[first_false - 1], R_face_sorted[first_false]
                    t = 0.0 if r1 == r0 else (R + 1.0 - r0) / (r1 - r0)
                    chamfer_depth = float(d0 + t * (d1 - d0))
        else:
            corner_radius = 0.0
            chamfer_depth = float(face_exclude)
        result["corner_radius"] = max(corner_radius, 1e-3)
        result["chamfer_depth"] = max(chamfer_depth, 1e-3)
        result["face_band_mask"] = face_band

    return result


def _two_pass_linear_fit(area, x):
    """Fit `R_map(x)` (x = signed or absolute z offset) linearly, derive a
    single wall-slope `cos_theta` from the fitted slope, re-invert the
    area with that `cos_theta`, and refit once more."""
    if len(x) < 2:
        raise ValueError("Not enough slices to fit a linear wall (check face_exclude/apex_exclude)")
    R = map_radius_from_area(area, 1.0)
    slope, _intercept = np.polyfit(x, R, 1)
    cos_theta = 1.0 / np.sqrt(1.0 + slope ** 2)
    R2 = map_radius_from_area(area, cos_theta)
    slope2, intercept2 = np.polyfit(x, R2, 1)
    return float(slope2), float(intercept2), float(cos_theta)


def fit_double_cone(prof, L, z_center, *, face_exclude=5.0, apex_exclude=5.0):
    """Two-pass slope fit of `R_map` vs `|z - z_center|` on the sloped
    wall (apex and faces excluded): slope -> `cos_theta` -> re-invert ->
    refit once. `inner = intercept`, `outer = inner + slope * L/2`. Also
    reports the bottom/top asymmetry from independent per-side fits."""
    half_L = L / 2.0
    dz = np.abs(prof.z - z_center)
    membrane = prof.f_far < 0.5
    core = membrane & (dz <= half_L - face_exclude) & (dz >= apex_exclude)
    if core.sum() < 2:
        raise ValueError(
            "Not enough core membrane slices to fit double_cone (check face_exclude/apex_exclude)"
        )

    slope, intercept, cos_theta = _two_pass_linear_fit(prof.area_pore[core], dz[core])
    inner = intercept
    outer = inner + slope * half_L

    result = {
        "pore_radius": max(inner, 1e-3),
        "outer_radius": max(outer, 1e-3),
        "cos_theta": cos_theta,
        "slope": slope,
        "core_mask": core,
    }

    signed = prof.z - z_center
    bottom = core & (signed < 0)
    top = core & (signed > 0)
    if bottom.sum() >= 2 and top.sum() >= 2:
        sb, ib, _ = _two_pass_linear_fit(prof.area_pore[bottom], dz[bottom])
        st, it, _ = _two_pass_linear_fit(prof.area_pore[top], dz[top])
        outer_b = ib + sb * half_L
        outer_t = it + st * half_L
        result["asymmetry"] = {"inner": ib - it, "outer": outer_b - outer_t}

    return result


def fit_conical(prof, L, z_center, *, face_exclude=5.0):
    """Two-pass slope fit of `R_map` vs signed `z - z_center` (asymmetric,
    unlike double_cone): `bottom = R(-L/2)`, `top = R(+L/2)`."""
    half_L = L / 2.0
    signed = prof.z - z_center
    membrane = prof.f_far < 0.5
    core = membrane & (np.abs(signed) <= half_L - face_exclude)
    if core.sum() < 2:
        raise ValueError("Not enough core membrane slices to fit conical (check face_exclude)")

    slope, intercept, cos_theta = _two_pass_linear_fit(prof.area_pore[core], signed[core])
    bottom = intercept + slope * (-half_L)
    top = intercept + slope * half_L

    return {
        "bottom_radius": max(bottom, 1e-3),
        "top_radius": max(top, 1e-3),
        "cos_theta": cos_theta,
        "slope": slope,
        "core_mask": core,
    }


# ---------------------------------------------------------------------------
# Parametric-model area profile and its comparison to the bin
# ---------------------------------------------------------------------------
def parametric_area_profile(profile, axes, spacing, metric="euclidean"):
    """
    `A_param(z) = h**2 * sum_xy condfrac(distance)` for a `PoreProfile`,
    evaluated one z-slice at a time (bounded memory for large grids).
    `axes = (x, y, z)`: `x`/`y` should already be shifted so the profile's
    (0, 0) axis lines up with the true pore center, and `z` should already
    be expressed relative to the profile's own z=0 (its membrane midplane).
    """
    x, y, z = axes
    pix_area = float(spacing[0]) * float(spacing[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    R_flat = np.sqrt(X ** 2 + Y ** 2).ravel()

    z = np.atleast_1d(np.asarray(z, dtype=float))
    A = np.empty(z.size)
    for i, zi in enumerate(z):
        Zq = np.full(R_flat.shape, zi, dtype=float)
        d = profile.distance(R_flat, Zq, metric=metric)
        A[i] = pix_area * float(condfrac(d).sum())
    return A


def predicted_deviation_pct(A_bin, A_param, spacing):
    """
    `100 * (sum(h / A_bin) - sum(h / A_param)) / sum(h / A_param)`: the
    fractional difference between the 1-D axial resistance integral of the
    measured bin area and of the fitted parametric model, over whatever
    slices the caller selects (typically the membrane slices). This 1-D
    integral predicts the full FEM current deviation to about 0.5%.
    """
    hz = float(spacing[2]) if hasattr(spacing, "__len__") else float(spacing)
    A_bin = np.asarray(A_bin, dtype=float)
    A_param = np.asarray(A_param, dtype=float)
    R_bin = np.sum(hz / A_bin)
    R_param = np.sum(hz / A_param)
    return float(100.0 * (R_bin / R_param - 1.0))


# ---------------------------------------------------------------------------
# Refinement
# ---------------------------------------------------------------------------
def _param_names_for(pore_type, fit_corner):
    if pore_type == "cylindrical":
        return ["pore_radius"] + (["corner_radius", "chamfer_depth"] if fit_corner else [])
    if pore_type == "double_cone":
        return ["pore_radius", "outer_radius"]
    if pore_type == "conical":
        return ["bottom_radius", "top_radius"]
    raise ValueError(f"refine() does not support pore_type={pore_type!r}")


def _refine_params(
    pore_type,
    params,
    L,
    prof,
    axes_xy,
    offset,
    z_center,
    *,
    refine_thickness=False,
    fit_corner=False,
    metric="euclidean",
):
    """
    `scipy.optimize.least_squares` on `(A_param - A_bin) / A_bin` over
    slices with `|z - z_center| <= L/2 + RAMP_MAX + 1` (wide enough to
    cover a face-band chamfer/corner), refining the fitted radii (and `L`
    itself when `refine_thickness`). All parameters are bounded > 0.
    """
    x, y = axes_xy
    x0, y0 = offset
    x_shift = x - x0
    y_shift = y - y0

    names = _param_names_for(pore_type, fit_corner)
    x0_vec = [float(params[n]) for n in names]
    if refine_thickness:
        names_all = names + ["membrane_thickness"]
        x0_vec = x0_vec + [float(L)]
    else:
        names_all = names

    half_window = L / 2.0 + RAMP_MAX + 1.0
    mask = np.abs(prof.z - z_center) <= half_window
    z_sel = prof.z[mask] - z_center
    A_bin_sel = prof.area_bin[mask]

    static_kwargs = {}
    if pore_type == "cylindrical" and not fit_corner:
        static_kwargs["corner_radius"] = params.get("corner_radius", 0.0)
        static_kwargs["chamfer_depth"] = params.get("chamfer_depth", None)

    def build(vec):
        kwargs = dict(zip(names, vec[: len(names)]))
        kwargs.update(static_kwargs)
        Lcur = vec[len(names)] if refine_thickness else L
        kwargs["membrane_thickness"] = Lcur
        return PoreProfile.from_params(pore_type, **kwargs)

    def residual(vec):
        profile_obj = build(vec)
        A_param = parametric_area_profile(
            profile_obj, (x_shift, y_shift, z_sel), prof.spacing, metric=metric
        )
        return (A_param - A_bin_sel) / A_bin_sel

    result = least_squares(residual, x0_vec, bounds=(1e-6, np.inf))
    refined = {name: float(val) for name, val in zip(names_all, result.x)}
    residual_rms = float(np.sqrt(np.mean(result.fun ** 2)))
    return refined, residual_rms


# Public alias matching the task-brief function list (kept distinct from
# the private implementation so `derive()`'s `refine: bool` keyword
# argument can shadow the name without losing access to the function).
refine = _refine_params


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _axis_symmetry_warnings(axes):
    warnings_out = []
    for name, axis in zip(("x", "y", "z"), axes):
        lo, hi = float(axis.min()), float(axis.max())
        scale = max(abs(lo), abs(hi), 1.0)
        if abs(lo + hi) > 1e-2 * scale:
            warnings_out.append(
                f"bin axis {name} is not symmetric about 0 (range [{lo:.3f}, {hi:.3f}]); "
                "BinFilePore recentres this grid when loaded for simulation"
            )
    return warnings_out


def _axis_offset(f3d, x, y, core_mask):
    if core_mask.sum() < 1:
        return 0.0, 0.0
    X, Y = np.meshgrid(x, y, indexing="ij")
    f_block = f3d[:, :, core_mask]
    weight = float(f_block.sum())
    if weight <= 0:
        return 0.0, 0.0
    x0 = float((X[:, :, None] * f_block).sum() / weight)
    y0 = float((Y[:, :, None] * f_block).sum() / weight)
    return x0, y0


def derive(
    bin_path,
    pore_type,
    *,
    units="distance",
    bulk_conductivity=None,
    face_exclude=5.0,
    apex_exclude=5.0,
    fit_corner=False,
    refine=True,
    refine_thickness=False,
):
    """
    Derive a `pore_geometry` block (or a `profile` radius-vs-z table) from
    an all-atom `.bin` conductivity/distance map. See the module docstring
    for the overall approach.
    """
    pore_type = pore_type.lower()
    if pore_type not in PORE_TYPES:
        raise ValueError(f"Unknown pore_type {pore_type!r}; expected one of {PORE_TYPES}")

    f3d, axes, spacing, meta = load_fraction_field(bin_path, units=units, bulk_conductivity=bulk_conductivity)
    x, y, z = axes

    warnings_list = _axis_symmetry_warnings(axes)

    prof = slice_profile(f3d, axes, spacing)
    L_integral, L_crossing, z_center = membrane_extent(prof)
    L = L_integral

    membrane_mask = prof.f_far < 0.5
    # Restricted to the *core* (face_exclude Å excluded at each face): the
    # face itself is a solid surface too, so f_far necessarily ramps
    # between 0 and 1 over a band of order RAMP_MIN..RAMP_MAX Å straddling
    # each true membrane face regardless of box size -- checking the full
    # f_far < 0.5 mask would flag that ordinary ramp overhang as a false
    # "far field not solid" positive on every pore. (The notes' literal
    # "f_far < 0.999" read against the full membrane mask would fire on
    # every membrane slice; f_far > 0.001 on the face-excluded core is the
    # intended, non-vacuous reading of that check.)
    core_for_offset = membrane_mask & (np.abs(prof.z - z_center) <= L / 2.0 - face_exclude)
    if np.any(prof.f_far[core_for_offset] > 0.001):
        max_f_far = float(np.max(prof.f_far[core_for_offset]))
        warnings_list.append(
            f"far field not solid for at least one core membrane slice (max f_far={max_f_far:.4f}); "
            "box may be too small relative to the pore"
        )

    x0, y0 = _axis_offset(f3d, x, y, core_for_offset)
    offset_mag = float(np.hypot(x0, y0))
    if offset_mag > 0.5:
        warnings_list.append(
            f"pore axis offset from box center is {offset_mag:.3f} A (x0={x0:.3f}, y0={y0:.3f}); "
            "consider re-centring the source structure"
        )

    params = None
    refined_params = None
    residual_rms = None
    final_L = L
    fit_diagnostics = {}

    if pore_type in PARAMETRIC_PORE_TYPES:
        if pore_type == "cylindrical":
            fit_result = fit_cylindrical(prof, L, z_center, face_exclude=face_exclude, fit_corner=fit_corner)
            warnings_list.extend(fit_result.get("warnings", []))
            params = {"pore_radius": fit_result["pore_radius"]}
            if fit_corner:
                params["corner_radius"] = fit_result["corner_radius"]
                params["chamfer_depth"] = fit_result["chamfer_depth"]
            else:
                params["corner_radius"] = 0.0
                params["chamfer_depth"] = None
        elif pore_type == "double_cone":
            fit_result = fit_double_cone(prof, L, z_center, face_exclude=face_exclude, apex_exclude=apex_exclude)
            params = {"pore_radius": fit_result["pore_radius"], "outer_radius": fit_result["outer_radius"]}
            if "asymmetry" in fit_result:
                # Per-side (bottom vs top) fit diagnostic, not a fitted
                # parameter: report it alongside the residuals rather than
                # feeding it into params/refine.
                fit_diagnostics["asymmetry"] = fit_result["asymmetry"]
        else:  # conical
            fit_result = fit_conical(prof, L, z_center, face_exclude=face_exclude)
            params = {"bottom_radius": fit_result["bottom_radius"], "top_radius": fit_result["top_radius"]}

        # Every parametric params dict above uses exactly the keyword names
        # PoreProfile.from_params (and the pore_geometry block itself) want,
        # so the refine/build step below is shared across all three types.
        final_params = dict(params)
        if refine:
            refined_params, residual_rms = _refine_params(
                pore_type, params, L, prof, (x, y), (x0, y0), z_center,
                refine_thickness=refine_thickness, fit_corner=fit_corner,
            )
            final_params.update({k: v for k, v in refined_params.items() if k != "membrane_thickness"})
            if refine_thickness:
                final_L = refined_params["membrane_thickness"]

        pore_geometry = {
            "pore_type": pore_type,
            **final_params,
            "membrane_thickness": final_L,
            "distance_metric": "euclidean",
        }
        profile_obj = PoreProfile.from_params(pore_type, membrane_thickness=final_L, **final_params)
        fit_result_core = fit_result

    else:  # profile
        z_tab, r_tab = profile_from_slices(prof, L, z_center)
        pore_geometry = {
            "pore_type": "profile",
            "profile_path": None,
            "membrane_thickness": final_L,
            "distance_metric": "euclidean",
        }
        profile_obj = PoreProfile.from_table(z_tab, r_tab, half_thickness=final_L / 2.0)
        fit_result_core = None

    # Always produce the raw radius-vs-z reconstruction: written as
    # *_profile_pore.csv regardless of pore_type, and it *is* the fit for
    # pore_type == "profile" (recomputed above with the final L already,
    # so just reuse it there instead of doing the work twice).
    if pore_type == "profile":
        pass  # z_tab, r_tab already computed against final_L above
    else:
        z_tab, r_tab = profile_from_slices(prof, final_L, z_center)

    # Full-axis diagnostics for P_profile.csv.
    prof.R_map = _compute_R_map(prof.area_pore, prof.z)
    prof.area_param = parametric_area_profile(
        profile_obj, (x - x0, y - y0, prof.z - z_center), prof.spacing
    )
    prof.r_fit = profile_obj.local_radius(prof.z - z_center)

    if pore_type == "profile":
        used = np.abs(prof.z - z_center) <= final_L / 2.0
    elif fit_result_core is not None and "core_mask" in fit_result_core:
        used = fit_result_core["core_mask"].copy()
        if "face_band_mask" in fit_result_core:
            used = used | fit_result_core["face_band_mask"]
    else:
        used = np.zeros_like(prof.z, dtype=bool)
    prof.used_in_fit = used

    deviation_pct = predicted_deviation_pct(
        prof.area_bin[membrane_mask], prof.area_param[membrane_mask], prof.spacing
    )

    derivation = {
        "inputs": {
            "bin_path": str(bin_path),
            "pore_type": pore_type,
            "units": meta.get("units"),
            "bulk_conductivity": meta.get("bulk_conductivity"),
        },
        "ramp": {"min": RAMP_MIN, "max": RAMP_MAX},
        "membrane_thickness_integral": L_integral,
        "membrane_thickness_crossing": L_crossing,
        "membrane_thickness_used": final_L,
        "z_center": z_center,
        "axis_offset": [x0, y0],
        "far_radius": prof.far_radius,
        "fit_parameters": params,
        "fit_diagnostics": fit_diagnostics,
        "refined_parameters": refined_params,
        "residual_rms": residual_rms,
        "predicted_deviation_pct": deviation_pct,
        "warnings": warnings_list,
        "provenance": {
            "sem_version": _SEM_VERSION,
            "git_commit": _git_commit(),
            "ramp": [RAMP_MIN, RAMP_MAX],
        },
    }
    if "sigma_min" in meta:
        derivation["inputs"]["sigma_min"] = meta["sigma_min"]

    return DerivedGeometry(
        pore_geometry=pore_geometry,
        profile_table=(z_tab, r_tab),
        slices=prof,
        derivation=derivation,
    )


# ---------------------------------------------------------------------------
# JSON / CSV output helpers
# ---------------------------------------------------------------------------
def _jsonable(obj):
    """Recursively convert numpy scalars/arrays (and dataclass-free dicts
    that may contain them) into plain-Python JSON-serialisable values."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _write_profile_csv(path, prof):
    header = "z,area_bin,area_pore,area_param,f_far,r_eq,R_map,r_fit,used_in_fit\n"
    with open(path, "w") as fh:
        fh.write(header)
        for i in range(prof.z.size):
            fh.write(
                "{:.6f},{:.6f},{:.6f},{:.6f},{:.8f},{:.6f},{:.6f},{:.6f},{}\n".format(
                    prof.z[i],
                    prof.area_bin[i],
                    prof.area_pore[i],
                    prof.area_param[i],
                    prof.f_far[i],
                    prof.r_eq[i],
                    prof.R_map[i],
                    prof.r_fit[i],
                    bool(prof.used_in_fit[i]),
                )
            )


def _write_profile_pore_csv(path, z_tab, r_tab):
    # Full precision (not the "%.6f" used for the diagnostic P_profile.csv):
    # this table's z=+-half_thickness end points must round-trip back to
    # the exact float used for pore_geometry["membrane_thickness"] / 2, or
    # re-loading it via PoreProfile.from_table/validate_config can reject
    # the tool's own output over a rounding-level shortfall (see
    # profile_from_slices for the matching in-memory endpoint pinning).
    # "%.17g" is always enough significant digits to round-trip a double.
    with open(path, "w") as fh:
        fh.write("z,r\n")
        for zi, ri in zip(z_tab, r_tab):
            fh.write("{:.17g},{:.17g}\n".format(float(zi), float(ri)))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="sem derive_geometry",
        description=(
            "Fit a parametric pore_geometry block (or a radius-vs-z profile "
            "table) to an all-atom (gen_dist) conductivity/distance map."
        ),
    )
    parser.add_argument("bin_path", metavar="BIN", help="Path to the all-atom .bin grid file")
    parser.add_argument(
        "--pore-type",
        required=True,
        choices=PORE_TYPES,
        help="Parametric pore type to fit, or 'profile' for a radius-vs-z table",
    )
    parser.add_argument(
        "--units", choices=("distance", "conductivity"), default="distance",
        help="What the bin stores (default: distance)",
    )
    parser.add_argument(
        "--bulk-conductivity", type=float, default=None,
        help="Bulk conductivity S/m (conductivity units only; default: max value in the bin)",
    )
    parser.add_argument("--face-exclude", type=float, default=5.0, help="Å excluded at each face (default: 5.0)")
    parser.add_argument(
        "--apex-exclude", type=float, default=5.0,
        help="Å excluded around the apex/center for double_cone (default: 5.0)",
    )
    parser.add_argument(
        "--fit-corner", action="store_true",
        help="Also fit a corner_radius/chamfer_depth for a cylindrical pore",
    )
    parser.add_argument("--no-refine", action="store_true", help="Skip the least_squares refinement step")
    parser.add_argument(
        "--refine-thickness", action="store_true",
        help="Also refine membrane_thickness during least_squares refinement",
    )
    parser.add_argument(
        "--output-prefix", default="derived",
        help="Prefix for output files (default: derived)",
    )
    parser.add_argument("--template", default=None, help="Existing config JSON to merge the derived geometry into")
    parser.add_argument(
        "--write-config", default=None,
        help="Where to write the merged config (requires --template)",
    )
    return parser


def main(argv=None):
    """
    CLI entry point. Returns a process exit code (0 on success).
    """
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    if bool(args.template) != bool(args.write_config):
        print("Error: --template and --write-config must be given together", file=sys.stderr)
        return 2

    try:
        result = derive(
            args.bin_path,
            args.pore_type,
            units=args.units,
            bulk_conductivity=args.bulk_conductivity,
            face_exclude=args.face_exclude,
            apex_exclude=args.apex_exclude,
            fit_corner=args.fit_corner,
            refine=not args.no_refine,
            refine_thickness=args.refine_thickness,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    prefix = args.output_prefix
    geom_path = f"{prefix}_pore_geometry.json"
    profile_csv_path = f"{prefix}_profile.csv"
    profile_pore_csv_path = f"{prefix}_profile_pore.csv"
    derivation_path = f"{prefix}_derivation.json"

    pore_geometry = dict(result.pore_geometry)
    if pore_geometry.get("pore_type") == "profile":
        pore_geometry["profile_path"] = profile_pore_csv_path

    Path(geom_path).parent.mkdir(parents=True, exist_ok=True)
    with open(geom_path, "w") as fh:
        json.dump(_jsonable(pore_geometry), fh, indent=2)

    _write_profile_csv(profile_csv_path, result.slices)

    z_tab, r_tab = result.profile_table
    _write_profile_pore_csv(profile_pore_csv_path, z_tab, r_tab)

    with open(derivation_path, "w") as fh:
        json.dump(_jsonable(result.derivation), fh, indent=2)

    if args.template:
        with open(args.template) as fh:
            cfg = json.load(fh)
        cfg["pore_geometry"] = _jsonable(pore_geometry)
        with open(args.write_config, "w") as fh:
            json.dump(cfg, fh, indent=2)

    print(f"Derived {args.pore_type} geometry from {args.bin_path}")
    for key, value in pore_geometry.items():
        print(f"  {key}: {value}")
    print(f"  predicted_deviation_pct: {result.derivation['predicted_deviation_pct']:.4f}")
    for warning in result.derivation["warnings"]:
        print(f"WARNING: {warning}")
    print(f"Wrote {geom_path}, {profile_csv_path}, {profile_pore_csv_path}, {derivation_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
