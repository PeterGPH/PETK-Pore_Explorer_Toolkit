"""
Axisymmetric reference solver for the open-pore conductance of a nanopore.

Solves the steady conduction problem  div(sigma grad phi) = 0  in cylindrical
coordinates (r, z) for a membrane of half-thickness ``h`` pierced by a pore
whose wall is a surface of revolution r = a(z), with bilinear (Q1) finite
elements on a graded tensor-product grid. Only numpy/scipy are required, so
it runs without the DOLFINx stack and is intended as an independent,
grid-converged reference for the 3D SEM pipeline and for the analytical
open-pore models in :mod:`sem.open_pore_theory`.

Conventions
-----------
* Lengths in Angstrom, bulk conductivity sigma_0 = 1. Resistances are
  returned in units of ``rho / Angstrom``; multiply by ``1e10 / sigma`` (sigma
  in S/m) for Ohms.
* ``far_field="infinite"`` models unbounded reservoirs with a Robin condition
  that is exact for the monopole far field phi = phi_inf - C/s (s measured
  from the pore mouth), leaving an O((a/S)^2) truncation error.
* ``far_field="box"`` reproduces a finite simulation box: phi fixed on the
  top/bottom faces at |z| = z_box/2, insulating side wall at r = r_box.
* Symmetric pores (``profile_symmetric=True``) are solved on z >= 0 with
  phi = 0 on the midplane; asymmetric ones on the full domain.

The wall distance entering the SEM ramp can be either the Euclidean distance
to the membrane (``distance="euclidean"``) or the horizontal distance
a(z) - r used by the tapered pores of :mod:`sem.pore_geometry`
(``distance="horizontal"``). For a straight cylinder the two coincide.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse import csgraph

# SEM conductivity ramp (sem.utils.condfrac)
R_MIN = 1.3
R_CUT = 4.1


def ramp_fraction(D, r_min: float = R_MIN, r_cut: float = R_CUT):
    """Linear SEM ramp: 0 for D <= r_min, 1 for D >= r_cut (a step if equal)."""
    if r_cut < r_min or r_min < 0:
        raise ValueError("ramp needs 0 <= r_min <= r_cut")
    D = np.asarray(D, float)
    if r_cut == r_min:
        return (D >= r_cut).astype(float)
    return np.clip((D - r_min) / (r_cut - r_min), 0.0, 1.0)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
@dataclass
class AxisymmetricPore:
    """Membrane {|z| <= h, r >= a(z)} with a piecewise-linear wall.

    ``wall_z`` must run from -h to +h (increasing) and ``wall_r`` gives the
    pore radius at those heights.
    """

    wall_z: np.ndarray
    wall_r: np.ndarray

    def __post_init__(self):
        self.wall_z = np.asarray(self.wall_z, float)
        self.wall_r = np.asarray(self.wall_r, float)
        if self.wall_z.ndim != 1 or self.wall_z.shape != self.wall_r.shape:
            raise ValueError("wall_z and wall_r must be 1D arrays of equal length")
        if np.any(np.diff(self.wall_z) <= 0):
            raise ValueError("wall_z must be strictly increasing")
        if not np.isclose(self.wall_z[0], -self.wall_z[-1]):
            raise ValueError("wall_z must span [-h, h]")
        if np.any(self.wall_r <= 0):
            raise ValueError("wall radii must be positive")

    # -- constructors ------------------------------------------------------
    @classmethod
    def cylinder(cls, radius: float, thickness: float) -> "AxisymmetricPore":
        h = 0.5 * thickness
        return cls([-h, h], [radius, radius])

    @classmethod
    def double_cone(cls, inner_radius: float, outer_radius: float,
                    thickness: float) -> "AxisymmetricPore":
        h = 0.5 * thickness
        return cls([-h, 0.0, h], [outer_radius, inner_radius, outer_radius])

    @classmethod
    def cone(cls, bottom_radius: float, top_radius: float,
             thickness: float) -> "AxisymmetricPore":
        h = 0.5 * thickness
        return cls([-h, h], [bottom_radius, top_radius])

    # -- properties --------------------------------------------------------
    @property
    def h(self) -> float:
        return float(self.wall_z[-1])

    @property
    def symmetric(self) -> bool:
        return bool(np.allclose(self.wall_r, self.wall_r[::-1])
                    and np.allclose(self.wall_z, -self.wall_z[::-1]))

    def radius_at(self, z):
        return np.interp(np.clip(z, -self.h, self.h), self.wall_z, self.wall_r)

    def is_solid(self, r, z):
        return (np.abs(z) <= self.h) & (r >= self.radius_at(z))

    def distance(self, r, z, mode: str = "euclidean"):
        """Distance from (r, z) to the solid membrane (0 inside it)."""
        r = np.asarray(r, float)
        z = np.asarray(z, float)
        if mode == "horizontal":
            rad = np.maximum(self.radius_at(z) - r, 0.0)
            ver = np.maximum(np.abs(z) - self.h, 0.0)
            return np.sqrt(rad ** 2 + ver ** 2)
        if mode != "euclidean":
            raise ValueError("mode must be 'euclidean' or 'horizontal'")
        h = self.h
        # boundary of the solid: top face ray, wall polyline, bottom face ray
        d2 = np.full(r.shape, np.inf)
        # face rays: points (x, +-h) with x >= a(+-h)
        for zf, af in ((h, self.wall_r[-1]), (-h, self.wall_r[0])):
            dx = np.maximum(af - r, 0.0)
            d2 = np.minimum(d2, dx ** 2 + (z - zf) ** 2)
        # wall segments
        for k in range(len(self.wall_z) - 1):
            p0 = np.array([self.wall_r[k], self.wall_z[k]])
            p1 = np.array([self.wall_r[k + 1], self.wall_z[k + 1]])
            v = p1 - p0
            t = ((r - p0[0]) * v[0] + (z - p0[1]) * v[1]) / (v @ v)
            t = np.clip(t, 0.0, 1.0)
            d2 = np.minimum(d2, (r - p0[0] - t * v[0]) ** 2 + (z - p0[1] - t * v[1]) ** 2)
        dist = np.sqrt(d2)
        return np.where(self.is_solid(r, z), 0.0, dist)


# --------------------------------------------------------------------------
# Grids
# --------------------------------------------------------------------------
def _graded(length: float, h0: float, hmax: float, growth: float) -> np.ndarray:
    """Offsets 0..length with spacing h0 growing geometrically to hmax."""
    if not (h0 > 0 and growth > 1.0):
        raise ValueError("grid needs h_fine > 0 and growth > 1")
    pts = [0.0]
    step = h0
    while pts[-1] + step < length - 1e-12:
        pts.append(pts[-1] + step)
        step = min(step * growth, hmax)
    if len(pts) > 1 and length - pts[-1] < 0.3 * step:
        pts[-1] = length
    else:
        pts.append(length)
    return np.asarray(pts)


def graded_axis(lo: float, hi: float, fine_intervals: Sequence[tuple],
                h_fine: float, h_max: float, growth: float = 1.06) -> np.ndarray:
    """1D node set on [lo, hi]: uniform spacing ``h_fine`` on each interval in
    ``fine_intervals`` (clipped to [lo, hi], merged) and geometric grading
    (factor ``growth``, capped at ``h_max``) in between."""
    if any(a > b for a, b in fine_intervals):
        raise ValueError("fine intervals must satisfy a <= b")
    iv = sorted((max(lo, a), min(hi, b)) for a, b in fine_intervals if b > lo and a < hi)
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    nodes = []
    cursor = lo
    for a, b in merged:
        if a > cursor:  # gap [cursor, a]: grade from both ends toward the middle
            gap = a - cursor
            if cursor == lo and lo not in [m[0] for m in merged]:
                # grade away from the fine block only (coarse at the outer end)
                seg = a - _graded(gap, h_fine, h_max, growth)[::-1]
            else:
                half = 0.5 * gap
                left = cursor + _graded(half, h_fine, h_max, growth)
                right = a - _graded(half, h_fine, h_max, growth)[::-1]
                seg = np.concatenate([left, right])
            nodes.append(seg)
        n = max(1, int(np.ceil((b - a) / h_fine - 1e-9)))
        nodes.append(np.linspace(a, b, n + 1))
        cursor = b
    if hi > cursor:
        nodes.append(cursor + _graded(hi - cursor, h_fine, h_max, growth))
    x = np.unique(np.round(np.concatenate(nodes), 9))
    return x


def _with_nodes(x: np.ndarray, extra) -> np.ndarray:
    """Insert the geometric breakpoints ``extra`` into the sorted node set
    ``x``, replacing any node closer than 1e-6 of the local spacing so that
    no degenerate (near-zero-width) element is created."""
    x = np.asarray(x, float)
    for e in np.atleast_1d(np.asarray(extra, float)):
        if e < x[0] or e > x[-1]:
            continue
        i = np.searchsorted(x, e)
        h = np.min(np.diff(x[max(i - 1, 0):i + 2])) if len(x) > 1 else 1.0
        near = np.abs(x - e) <= 1e-6 * h
        x = np.sort(np.concatenate([x[~near], [e]]))
    return x


# --------------------------------------------------------------------------
# Solver
# --------------------------------------------------------------------------
_GAUSS = {
    2: (np.array([-1.0, 1.0]) / np.sqrt(3.0), np.array([1.0, 1.0])),
    4: (np.array([-0.8611363115940526, -0.3399810435848563,
                  0.3399810435848563, 0.8611363115940526]),
        np.array([0.3478548451374538, 0.6521451548625461,
                  0.6521451548625461, 0.3478548451374538])),
}


def _assemble(rn, zn, sigma_fn, nq=4):
    nr, nz = len(rn), len(zn)
    R0, Z0 = np.meshgrid(rn[:-1], zn[:-1], indexing="ij")
    hr = np.diff(rn)[:, None] * np.ones((1, nz - 1))
    hz = np.ones((nr - 1, 1)) * np.diff(zn)[None, :]
    xg, wg = _GAUSS[nq]
    Ke = np.zeros(R0.shape + (4, 4))
    for gx, wx in zip(xg, wg):
        x = 0.5 * (gx + 1.0)
        for gy, wy in zip(xg, wg):
            y = 0.5 * (gy + 1.0)
            r = R0 + x * hr
            z = Z0 + y * hz
            w = 0.25 * wx * wy * hr * hz * r * sigma_fn(r, z)
            dx = np.array([-(1 - y), (1 - y), -y, y])
            dy = np.array([-(1 - x), -x, (1 - x), x])
            for p in range(4):
                for q in range(p, 4):
                    val = w * (dx[p] * dx[q] / hr ** 2 + dy[p] * dy[q] / hz ** 2)
                    Ke[..., p, q] += val
                    if q != p:
                        Ke[..., q, p] += val
    I, J = np.meshgrid(np.arange(nr - 1), np.arange(nz - 1), indexing="ij")
    loc = np.stack([I * nz + J, (I + 1) * nz + J, I * nz + J + 1, (I + 1) * nz + J + 1], -1)
    rows = np.broadcast_to(loc[..., :, None], Ke.shape).ravel()
    cols = np.broadcast_to(loc[..., None, :], Ke.shape).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(nr * nz, nr * nz))
    return K


def _robin_edges(boundary_edges, centre_z, phi_inf, N, sigma_fn):
    """Robin far-field term on a list of boundary edges.
    Each edge: (node0, node1, (r0, z0), (r1, z1)). Returns (Kb, rhs)."""
    rows, cols, vals = [], [], []
    rhs = np.zeros(N)
    xg, wg = _GAUSS[2]
    for n0, n1, p0, p1 in boundary_edges:
        p0 = np.asarray(p0, float)
        p1 = np.asarray(p1, float)
        L = np.hypot(*(p1 - p0))
        normal = np.array([p0[1] - p1[1], p1[0] - p0[0]]) / L
        mid = 0.5 * (p0 + p1)
        if normal @ (mid - np.array([0.0, centre_z])) < 0:
            normal = -normal
        for gx, wx in zip(xg, wg):
            t = 0.5 * (gx + 1.0)
            p = p0 + t * (p1 - p0)
            s = p - np.array([0.0, centre_z])
            coef = float(sigma_fn(np.array([p[0]]), np.array([p[1]]))[0]) * (normal @ s) / (s @ s)
            w = 0.5 * wx * L * p[0] * coef
            Nv = (1.0 - t, t)
            nn = (n0, n1)
            for a in range(2):
                rhs[nn[a]] += w * Nv[a] * phi_inf
                for b in range(2):
                    rows.append(nn[a])
                    cols.append(nn[b])
                    vals.append(w * Nv[a] * Nv[b])
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), rhs


def _components(K):
    return csgraph.connected_components(K != 0, directed=False)[1]


def _solve(K, rhs, fixed, fixed_val, anchors):
    """Solve K phi = rhs with phi[fixed] = fixed_val. Nodes with zero
    conductivity and conducting islands that touch no anchor (Dirichlet or
    Robin node) carry no current and are left at phi = 0."""
    N = K.shape[0]
    diag = K.diagonal()
    free_mask = diag > 0
    labels = _components(K)
    free_mask &= np.isin(labels, labels[anchors])
    free_mask[fixed] = False
    free = np.flatnonzero(free_mask)
    phi = np.zeros(N)
    phi[fixed] = fixed_val
    b = rhs[free] - K[free][:, fixed] @ fixed_val if len(fixed) else rhs[free]
    phi[free] = spla.spsolve(K[free][:, free].tocsc(), b)
    return phi


@dataclass
class SolveResult:
    resistance: float          # full-pore resistance, rho / Angstrom
    n_nodes: int
    rn: np.ndarray
    zn: np.ndarray
    phi: Optional[np.ndarray] = None

    @property
    def conductance(self) -> float:
        return 1.0 / self.resistance


def _solve_pore(
    pore: AxisymmetricPore,
    wall: str = "ramp",
    *,
    r_min: float = R_MIN,
    r_cut: float = R_CUT,
    step_distance: Optional[float] = None,
    distance: str = "euclidean",
    h_fine: float = 0.1,
    band: float = 8.0,
    far_field: str = "infinite",
    far_factor: float = 200.0,
    box_radius: Optional[float] = None,
    box_height: Optional[float] = None,
    growth: float = 1.03,
    keep_field: bool = False,
    sigma_fn: Optional[Callable] = None,
    use_symmetry: bool = True,
) -> SolveResult:
    """Single solve on one grid; see :func:`pore_resistance`."""
    h = pore.h
    a_lo, a_hi = float(pore.wall_r.min()), float(pore.wall_r.max())
    symmetric = pore.symmetric and use_symmetry
    kinks = ()
    if sigma_fn is None and wall == "ramp":
        kinks = (r_min, r_cut)
    elif sigma_fn is None and wall == "step":
        kinks = (step_distance,)

    if sigma_fn is None:
        if wall == "sharp":
            def sigma_fn(r, z):
                return np.where(pore.is_solid(r, z), 0.0, 1.0)
        elif wall == "ramp":
            def sigma_fn(r, z):
                return ramp_fraction(pore.distance(r, z, distance), r_min, r_cut)
        elif wall == "step":
            if step_distance is None:
                raise ValueError("wall='step' needs step_distance")

            def sigma_fn(r, z):
                return ((pore.distance(r, z, distance) >= step_distance)
                        & ~pore.is_solid(r, z)).astype(float)
        else:
            raise ValueError(f"unknown wall {wall!r}")

    if far_field == "infinite":
        far = far_factor * max(a_hi, h, 1.0)
        r_out, z_out = far, h + far
    elif far_field == "box":
        if box_radius is None or box_height is None:
            raise ValueError("box far field needs box_radius and box_height")
        r_out, z_out = box_radius, 0.5 * box_height
        if r_out <= a_hi or z_out <= h:
            raise ValueError("box does not contain the pore")
    else:
        raise ValueError("far_field must be 'infinite' or 'box'")

    band_r = min(band, 0.9 * a_lo)
    band_z = min(band, 0.9 * h) if h > 0 else band
    rn = graded_axis(0.0, r_out, [(a_lo - band_r, a_hi + band)], h_fine,
                     h_max=max(r_out / 8.0, h_fine), growth=growth)
    rn = _with_nodes(rn, pore.wall_r)
    # z: fine along the wall's z-extent if it is tapered, else near the faces
    if np.allclose(pore.wall_r, pore.wall_r[0]):
        zf = [(h - band_z, h + band)]
        if not symmetric:
            zf.append((-h - band, -h + band_z))
    else:
        zf = [(-h - band, h + band)]
    z_lo = 0.0 if symmetric else -z_out
    zn = graded_axis(z_lo, z_out, zf + ([(-band_z, band_z)] if not symmetric else []),
                     h_fine, h_max=max(z_out / 8.0, h_fine), growth=growth)
    zn = _with_nodes(zn, pore.wall_z[pore.wall_z >= z_lo])
    if symmetric:
        zn = _with_nodes(zn, [0.0])
    # put the kinks of the wall ramp on grid lines where the wall is straight
    # (faces always; the bore of a cylinder); otherwise Gauss quadrature of the
    # kinks adds a grid-alignment-dependent error of ~1e-6-1e-5
    if kinks:
        zn = _with_nodes(zn, [s * (h + k) for k in kinks for s in (1, -1)])
        if a_lo == a_hi:
            rn = _with_nodes(rn, [a_lo - k for k in kinks if a_lo - k > 0])

    nr, nz = len(rn), len(zn)
    N = nr * nz
    K = _assemble(rn, zn, sigma_fn)

    def idx(i, j):
        return i * nz + j

    rhs = np.zeros(N)
    fixed = []
    fixed_val = []
    if far_field == "infinite":
        edges_top = [(idx(i, nz - 1), idx(i + 1, nz - 1), (rn[i], zn[-1]), (rn[i + 1], zn[-1]))
                     for i in range(nr - 1)]
        edges_top += [(idx(nr - 1, j), idx(nr - 1, j + 1), (rn[-1], zn[j]), (rn[-1], zn[j + 1]))
                      for j in range(nz - 1) if zn[j] >= h - 1e-9]
        Kb, rb = _robin_edges(edges_top, h, 1.0, N, sigma_fn)
        K = K + Kb
        rhs += rb
        robin_top_rows = np.unique([e[0] for e in edges_top] + [e[1] for e in edges_top])
        terminal_a = robin_top_rows
        if not symmetric:
            edges_bot = [(idx(i, 0), idx(i + 1, 0), (rn[i], zn[0]), (rn[i + 1], zn[0]))
                         for i in range(nr - 1)]
            edges_bot += [(idx(nr - 1, j), idx(nr - 1, j + 1), (rn[-1], zn[j]), (rn[-1], zn[j + 1]))
                          for j in range(nz - 1) if zn[j + 1] <= -h + 1e-9]
            Kb2, rb2 = _robin_edges(edges_bot, -h, -1.0, N, sigma_fn)
            K = K + Kb2
            rhs += rb2
            terminal_b = np.unique([e[0] for e in edges_bot] + [e[1] for e in edges_bot])
    else:
        top = [idx(i, nz - 1) for i in range(nr)]
        fixed += top
        fixed_val += [1.0] * len(top)
        terminal_a = np.asarray(top)
        if not symmetric:
            bot = [idx(i, 0) for i in range(nr)]
            fixed += bot
            fixed_val += [-1.0] * len(bot)
            terminal_b = np.asarray(bot)

    if symmetric:
        a_mid = float(pore.radius_at(0.0))
        mid = [idx(i, 0) for i in range(nr) if rn[i] <= a_mid + 1e-9]
        fixed += mid
        fixed_val += [0.0] * len(mid)
        terminal_b = np.asarray(mid)
    fixed = np.asarray(fixed, int)
    fixed_val = np.asarray(fixed_val, float)

    # a pore closed by the membrane (or by the zero-conductivity layer) carries
    # no current: report R = inf instead of a round-off value
    labels = _components(K)
    if not np.intersect1d(labels[terminal_a], labels[terminal_b]).size:
        return SolveResult(np.inf, N, rn, zn, None)
    anchors = np.unique(np.concatenate([fixed, terminal_a, terminal_b]))
    phi = _solve(K, rhs, fixed, fixed_val, anchors)

    # current through the domain
    if symmetric:
        mid = fixed[fixed_val == 0.0]
        current = abs(2 * np.pi * (K[mid] @ phi - rhs[mid]).sum())
        resistance = 2.0 * 1.0 / current          # half-pore drop is 1
    else:
        if far_field == "infinite":
            # consistent Robin flux on the top boundary (phi_inf = +1)
            current = abs(2 * np.pi * (rb[robin_top_rows] - (Kb @ phi)[robin_top_rows]).sum())
        else:
            top_nodes = fixed[fixed_val == 1.0]
            current = abs(2 * np.pi * (K[top_nodes] @ phi - rhs[top_nodes]).sum())
        resistance = 2.0 / current                # total drop is 2
    return SolveResult(resistance, N, rn, zn, phi.reshape(nr, nz) if keep_field else None)


def pore_resistance(pore: AxisymmetricPore, wall: str = "ramp", *,
                    extrapolate: bool = False, **kwargs) -> SolveResult:
    """Open-pore resistance of ``pore`` (rho / Angstrom; full pore).

    wall: "sharp" (sigma = 1 in the fluid, 0 in the membrane), "ramp" (SEM
    linear ramp from r_min to r_cut in the wall distance), or "step" (sigma
    = 1 where the wall distance >= step_distance, i.e. a sharp wall on the
    Minkowski offset of the membrane). ``sigma_fn(r, z)`` overrides all.

    The grid is uniform (``h_fine``) within ``band`` of the wall and the
    faces and graded outward. A tapered wall is fine-gridded over its whole
    z-extent, so for long cones start from h_fine = 0.2-0.5 A; h_fine = 0.1 A
    and growth = 1.03 give the ramp conductance of a cylinder to ~2e-5
    (the error is set by the graded bulk, ~ (growth - 1)^2, not by h_fine).
    ``extrapolate=True`` adds a solve with h_fine/2 and (growth - 1)/2 and
    Richardson-extrapolates assuming second order, which brings smooth (ramp)
    maps to ~1e-6. Sharp corners converge as h^(4/3); for sharp cylinders use
    :mod:`sem.cylinder_mode_matching` instead.
    """
    res = _solve_pore(pore, wall, **kwargs)
    if not extrapolate:
        return res
    fine = dict(kwargs)
    fine["h_fine"] = 0.5 * kwargs.get("h_fine", 0.1)
    fine["growth"] = 1.0 + 0.5 * (kwargs.get("growth", 1.03) - 1.0)
    res2 = _solve_pore(pore, wall, **fine)
    R = res2.resistance + (res2.resistance - res.resistance) / 3.0
    return SolveResult(R, res.n_nodes + res2.n_nodes, res2.rn, res2.zn, res2.phi)
