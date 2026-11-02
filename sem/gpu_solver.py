"""
Batched GPU backend for SEM rotation scans.

Solves -div(sigma grad phi) = 0 on the configured box with phi = voltage on the top face,
phi = 0 on the bottom face and insulating side faces, for many analyte orientations at once.

Discretisations (`method`):
- "fem" (default): the same P1 finite elements as VerticalMovementSEM. DOLFINx create_box splits
  every cube into the 6 Kuhn tetrahedra around its main diagonal, with one conductivity per
  tetrahedron (DG0, sampled at the centroid). On that mesh the P1 stiffness matrix couples a
  vertex only to its 6 axis neighbours: each tetrahedron adds -sigma_t * (area / length) / 6 to
  each of its 3 axis edges and nothing to its diagonal edges. The GPU assembles exactly that
  7-point vertex stencil, so results agree with the FEM path to solver tolerance. The current
  is PETK's top-face flux (sigma * d(phi)/dz over the top boundary facets).
- "fv": cell-centred finite volumes on the box cells, harmonic-mean face conductances, current
  from dissipated power. Faster per orientation, but a different discretisation of the model.

Solver: conjugate gradients preconditioned with a multigrid V-cycle (Galerkin coarse operators
on blocks of 2 unknowns per axis, overcorrected coarse-grid correction, weighted Jacobi
smoothing, dense coarsest solve), started from the open-pore solution and run for a fixed
number of iterations.

Conductivity follows the FEM path exactly: base (pore) conductivity is the pore's own
interpolator at the sample points, and the analyte modification reproduces
VerticalMovementSEM.calculate_analyte_conductivity_modification (nearest atom centre within
`cutoff`, minus that atom's radius, then the conductivity model, then the minimum with the
base). Units are PETK's: Angstrom, volts, conductivity in nS/A, current in nA.
"""

import itertools
import logging
import math
import time
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn.functional as F

logger = logging.getLogger(__name__)

OMEGA = 0.8   # weighted-Jacobi damping
NU = 2        # pre- and post-smoothing sweeps
# Scale of the coarse-grid correction. Piecewise-constant Galerkin coarse operators are too
# stiff by about 2, so the correction is scaled up (overcorrection, Braess 1995).
OVERCORRECTION = 1.8
DENSE_MAX = 1000  # largest coarsest level solved with a dense inverse

# Kuhn tetrahedra of a cube, as the order of axis steps from the cube origin to the opposite
# corner (DOLFINx create_box). The centroid sits at 3/4, 1/2, 1/4 of a cell along those axes.
KUHN = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))


def kuhn_centroid(perm):
    off = [0.0, 0.0, 0.0]
    for frac, axis in zip((0.75, 0.5, 0.25), perm):
        off[axis] = frac
    return tuple(off)


# ---------------------------------------------------------------- grid
@dataclass
class BoxGrid:
    """The cells of the DOLFINx box mesh of VerticalMovementSEM.setup_dolfinx."""
    n: tuple   # cells per axis
    h: tuple   # spacing per axis (A)
    lo: tuple  # lower corner (A)

    @classmethod
    def from_sem(cls, sem):
        bd = sem.box_dimensions
        ext = [bd["x"][1] - bd["x"][0] + sem.xy_margin,
               bd["y"][1] - bd["y"][0] + sem.xy_margin,
               bd["z"][1] - bd["z"][0]]
        n = tuple(max(int(np.ceil(e / sem.grid_resolution)), 1) for e in ext)
        h = tuple(e / k for e, k in zip(ext, n))
        lo = tuple(-e / 2.0 for e in ext)
        return cls(n, h, lo)

    def points(self, offset=(0.5, 0.5, 0.5)):
        """(N, 3) points at `offset` (in cell units) inside every cell, C order (x slowest)."""
        axes = [self.lo[a] + (np.arange(self.n[a]) + offset[a]) * self.h[a] for a in range(3)]
        X, Y, Z = np.meshgrid(*axes, indexing="ij")
        return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    def centres(self):
        return self.points()


# ---------------------------------------------------------------- operator
def _harm(a, b):
    return 2.0 * a * b / (a + b)


def _pad(t, axis, before, after):
    """Zero-pad tensor t along dimension `axis` (>= 1; dimension 0 is the batch)."""
    pads = [0] * (2 * (t.dim() - axis))
    pads[-2], pads[-1] = before, after
    return F.pad(t, pads)


def _diag(L):
    D = (_pad(L["gx"], 1, 1, 0) + _pad(L["gx"], 1, 0, 1)
         + _pad(L["gy"], 2, 1, 0) + _pad(L["gy"], 2, 0, 1)
         + _pad(L["gz"], 3, 1, 0) + _pad(L["gz"], 3, 0, 1))
    nz = D.shape[3]
    D = D + _pad(L["gt"][..., None], 3, nz - 1, 0) + _pad(L["gb"][..., None], 3, 0, nz - 1)
    L["D"] = D
    L["Dinv"] = 1.0 / D


def fv_level(sig, h):
    """Finite volumes: face conductances from cell conductivity sig (B, nx, ny, nz)."""
    hx, hy, hz = h
    L = {
        "gx": (hy * hz / hx) * _harm(sig[:, 1:], sig[:, :-1]),
        "gy": (hx * hz / hy) * _harm(sig[:, :, 1:], sig[:, :, :-1]),
        "gz": (hx * hy / hz) * _harm(sig[..., 1:], sig[..., :-1]),
        "gt": (2.0 * hx * hy / hz) * sig[..., -1],
        "gb": (2.0 * hx * hy / hz) * sig[..., 0],
    }
    _diag(L)
    return L


def _kuhn_edges(sig6, h, a):
    """P1 coupling of the vertex edges along axis a: (area/length)/6 * sum of the 6 tetrahedra
    around each edge. sig6 is (B, 6, nx, ny, nz); the result has n_a edges along a and n+1
    vertices along the other two axes."""
    others = [d for d in range(3) if d != a]
    w = h[others[0]] * h[others[1]] / h[a] / 6.0
    C = 0.0
    for r in itertools.product((0, 1), repeat=2):
        # the edge sits at relative position r (along the other two axes) in a cube; the
        # tetrahedra containing it step along exactly the axes where r = 1 before stepping along a
        before = {others[i] for i in range(2) if r[i]}
        tets = [t for t, p in enumerate(KUHN) if set(p[:p.index(a)]) == before]
        S = sig6[:, tets].sum(1)
        for i, d in enumerate(others):   # r = 0: edge vertex owns cube v; r = 1: cube v - 1
            S = _pad(S, d + 1, r[i], 1 - r[i])
        C = C + S
    return w * C


def fem_level(sig6, h):
    """P1 finite elements on the Kuhn mesh: 7-point stencil on the vertex unknowns (all vertices
    except the Dirichlet top and bottom layers), plus the top-face flux weights."""
    Cx, Cy, Cz = (_kuhn_edges(sig6, h, a) for a in range(3))
    L = {
        "gx": Cx[..., 1:-1],
        "gy": Cy[..., 1:-1],
        "gz": Cz[..., 1:-1],
        "gb": Cz[..., 0],
        "gt": Cz[..., -1],
    }
    _diag(L)
    hx, hy, hz = h
    zfirst = [t for t, p in enumerate(KUHN) if p[0] == 2]  # tetrahedra with a facet on the top face
    L["wt"] = sig6[:, zfirst, :, :, -1].sum(1) * (hx * hy / (2.0 * hz))
    return L


def apply(L, u):
    gx, gy, gz = L["gx"], L["gy"], L["gz"]
    return (L["D"] * u
            - _pad(gx * u[:, 1:], 1, 0, 1) - _pad(gx * u[:, :-1], 1, 1, 0)
            - _pad(gy * u[:, :, 1:], 2, 0, 1) - _pad(gy * u[:, :, :-1], 2, 1, 0)
            - _pad(gz * u[..., 1:], 3, 0, 1) - _pad(gz * u[..., :-1], 3, 1, 0))


def _block_sum(t, axis, f):
    """Sum consecutive blocks of f entries along `axis`; a short last block is zero-padded."""
    if f == 1:
        return t
    n = t.shape[axis]
    m = -(-n // f)
    if m * f != n:
        t = _pad(t, axis, 0, m * f - n)
    return t.unflatten(axis, (m, f)).sum(axis + 1)


def _interface(gl, axis, f):
    """Fine links that cross the faces between consecutive blocks of f unknowns along `axis`."""
    if f == 1:
        return gl
    g = gl.movedim(axis, -1)
    n = g.shape[-1] + 1
    return g[..., f - 1:n - 1:f].movedim(-1, axis)


def _factors(n):
    return tuple(2 if k > 5 else 1 for k in n)


def coarsen(L, f):
    """Galerkin coarse operator P^T A P for piecewise-constant P on blocks of f unknowns.

    A coarse link is the sum of the fine links crossing the face between two blocks, so a
    conducting channel never closes on a coarser level, however it bends inside a block.
    """
    fx, fy, fz = f
    C = {
        "gx": _interface(_block_sum(_block_sum(L["gx"], 2, fy), 3, fz), 1, fx),
        "gy": _interface(_block_sum(_block_sum(L["gy"], 1, fx), 3, fz), 2, fy),
        "gz": _interface(_block_sum(_block_sum(L["gz"], 1, fx), 2, fy), 3, fz),
        "gt": _block_sum(_block_sum(L["gt"], 1, fx), 2, fy),
        "gb": _block_sum(_block_sum(L["gb"], 1, fx), 2, fy),
    }
    _diag(C)
    return C


def _restrict(r, f):
    for a, k in zip((1, 2, 3), f):
        r = _block_sum(r, a, k)
    return r


def _prolong(e, f, shape):
    for a, k in zip((1, 2, 3), f):
        if k > 1:
            e = e.repeat_interleave(k, a).narrow(a, 0, shape[a])
    return e


def _dense_inverse(L):
    """Inverse of the coarsest-level matrix, assembled directly from its conductances (float64)."""
    B, nx, ny, nz = L["D"].shape
    k = nx * ny * nz
    dev, d = L["D"].device, torch.float64
    idx = torch.arange(k, device=dev).reshape(nx, ny, nz)
    A = torch.zeros(B, k, k, device=dev, dtype=d)
    A[:, idx.ravel(), idx.ravel()] = L["D"].reshape(B, k).to(d)
    for g, a, b in ((L["gx"], idx[:-1], idx[1:]), (L["gy"], idx[:, :-1], idx[:, 1:]),
                    (L["gz"], idx[..., :-1], idx[..., 1:])):
        a, b, v = a.ravel(), b.ravel(), -g.reshape(B, -1).to(d)
        A[:, a, b] = v
        A[:, b, a] = v
    L["Ainv"] = torch.linalg.inv(A).to(L["D"].dtype)


def setup_hierarchy(fine):
    levels, factors = [fine], []
    while True:
        f = _factors(levels[-1]["D"].shape[1:])
        if f == (1, 1, 1):
            break
        factors.append(f)
        levels.append(coarsen(levels[-1], f))
    k = int(np.prod(levels[-1]["D"].shape[1:]))
    if k > DENSE_MAX:
        raise ValueError(f"coarsest level has {k} unknowns (> {DENSE_MAX})")
    _dense_inverse(levels[-1])
    return levels, factors


def vcycle(levels, factors, l, r):
    L = levels[l]
    if l == len(levels) - 1:
        B = r.shape[0]
        return (L["Ainv"] @ r.reshape(B, -1, 1)).reshape(r.shape)
    x = OMEGA * L["Dinv"] * r
    for _ in range(NU - 1):
        x = x + OMEGA * L["Dinv"] * (r - apply(L, x))
    rc = _restrict(r - apply(L, x), factors[l])
    x = x + OVERCORRECTION * _prolong(vcycle(levels, factors, l + 1, rc), factors[l], r.shape)
    for _ in range(NU):
        x = x + OMEGA * L["Dinv"] * (r - apply(L, x))
    return x


def rhs(L, voltage):
    nz = L["D"].shape[3]
    return _pad((L["gt"] * voltage)[..., None], 3, nz - 1, 0)


def power_current(L, u, voltage):
    """Current (nA) from dissipated power; accumulated in float64."""
    d = torch.float64
    u = u.to(d)
    gx, gy, gz = L["gx"].to(d), L["gy"].to(d), L["gz"].to(d)
    e = (gx * (u[:, 1:] - u[:, :-1]) ** 2).sum((1, 2, 3))
    e = e + (gy * (u[:, :, 1:] - u[:, :, :-1]) ** 2).sum((1, 2, 3))
    e = e + (gz * (u[..., 1:] - u[..., :-1]) ** 2).sum((1, 2, 3))
    e = e + (L["gt"].to(d) * (voltage - u[..., -1]) ** 2).sum((1, 2))
    e = e + (L["gb"].to(d) * u[..., 0] ** 2).sum((1, 2))
    return e / voltage


def flux_current(L, u, voltage):
    """PETK's current for the FEM level: sigma * d(phi)/dz integrated over the top facets."""
    wt = L["wt"].to(torch.float64)
    nx, ny = wt.shape[1:]
    return (wt * (voltage - u[:, :nx, :ny, -1].to(torch.float64))).sum((1, 2))


def pcg(levels, factors, b, x, iters, current=None, track=False):
    """Fixed-iteration multigrid-preconditioned CG; with track, records current(x) per iteration."""
    L, red = levels[0], (1, 2, 3)
    r = b - apply(L, x)
    z = vcycle(levels, factors, 0, r)
    p = z.clone()
    rz = (r * z).sum(red, keepdim=True)
    hist = [current(x)] if track else None
    tiny = torch.finfo(x.dtype).tiny
    for _ in range(iters):
        q = apply(L, p)
        alpha = rz / ((p * q).sum(red, keepdim=True) + tiny)
        x = x + alpha * p
        r = r - alpha * q
        z = vcycle(levels, factors, 0, r)
        rz_new = (r * z).sum(red, keepdim=True)
        p = z + (rz_new / (rz + tiny)) * p
        rz = rz_new
        if track:
            hist.append(current(x))
    return x, hist


# ---------------------------------------------------------------- analyte conductivity
class AnalyteSigma:
    """Exact GPU copy of calculate_analyte_conductivity_modification, evaluated at one point per
    grid cell (at `offset` in cell units).

    With `membrane_conductivity` set, `batch(..., with_overlap=True)` also flags orientations in
    which an atom's hard core (surface distance <= model.min_distance) reaches a membrane point,
    the check prevent_analyte_overlap makes in the FEM path.
    """

    def __init__(self, grid, base_sigma, radii, cutoff, model, device, offset=(0.5, 0.5, 0.5),
                 membrane_conductivity=None):
        self.grid = grid
        self.device = device
        self.base = torch.as_tensor(base_sigma, dtype=torch.float64, device=device)
        self.radii = torch.as_tensor(np.asarray(radii, dtype=np.float64), device=device)
        self.cutoff = float(cutoff)
        self.min_distance = float(model.min_distance)
        self.model_cutoff = float(model.cutoff)
        self.bulk = float(model.bulk_conductivity)
        self.membrane = None if membrane_conductivity is None else float(membrane_conductivity)
        self.n = torch.tensor(grid.n, device=device)
        self.h = torch.tensor(grid.h, dtype=torch.float64, device=device)
        self.lo = torch.tensor(grid.lo, dtype=torch.float64, device=device)
        self.off = torch.tensor(offset, dtype=torch.float64, device=device)
        k = [int(math.ceil(self.cutoff / hh)) + 1 for hh in grid.h]
        offs = torch.stack(torch.meshgrid(*(torch.arange(-kk + 1, kk + 1, device=device) for kk in k),
                                          indexing="ij"), -1).reshape(-1, 3)
        self.offsets = offs  # candidate cells around each atom

    def model(self, s):
        frac = ((s - self.min_distance) / (self.model_cutoff - self.min_distance)).clamp(0, 1)
        return torch.where(s < self.model_cutoff, 1e-7 * self.bulk + frac * self.bulk,
                           torch.full_like(s, self.bulk))

    def __call__(self, positions):
        """positions: (M, 3) float64 array of atom centres -> flat conductivity (float64)."""
        return self.batch([positions])[0]

    def batch(self, positions, with_overlap=False):
        """positions: list of (M, 3) arrays, one per orientation -> (B, ncell) conductivity
        (and a (B,) bool tensor of hard-core membrane overlaps when with_overlap)."""
        p = torch.as_tensor(np.stack(positions), dtype=torch.float64, device=self.device)  # (B, M, 3)
        B, M = p.shape[:2]
        ncell = int(np.prod(self.grid.n))
        base_idx = torch.floor((p - self.lo) / self.h - self.off).long()   # sample point just below
        cells = base_idx[:, :, None, :] + self.offsets                      # (B, M, K, 3)
        inside = ((cells >= 0) & (cells < self.n)).all(-1)
        point = self.lo + (cells.double() + self.off) * self.h
        d = torch.sqrt(((point - p[:, :, None, :]) ** 2).sum(-1))
        keep = inside & (d < self.cutoff)
        atom = torch.arange(M, device=self.device)[None, :, None].expand_as(d)[keep]
        orient = torch.arange(B, device=self.device)[:, None, None].expand_as(d)[keep]
        d = d[keep]
        c = cells[keep]
        flat = orient * ncell + (c[:, 0] * self.grid.n[1] + c[:, 1]) * self.grid.n[2] + c[:, 2]
        dmin = torch.full((B * ncell,), float("inf"), dtype=torch.float64, device=self.device)
        dmin = dmin.scatter_reduce(0, flat, d, reduce="amin")
        nearest = d == dmin[flat]
        r = torch.full((B * ncell,), float("inf"), dtype=torch.float64, device=self.device)
        r = r.scatter_reduce(0, flat[nearest], self.radii[atom[nearest]], reduce="amin")
        hit = torch.isfinite(dmin)
        s = (dmin[hit] - r[hit]).clamp_min(0.0)
        sig = self.base.repeat(B)
        base_hit = sig[hit]
        sig[hit] = torch.minimum(base_hit, self.model(s))
        sig = sig.reshape(B, ncell)
        if not with_overlap:
            return sig
        overlap = torch.zeros(B, dtype=torch.bool, device=self.device)
        if self.membrane is not None:
            bad = (s <= self.min_distance) & (base_hit <= self.membrane)
            which = torch.nonzero(hit).squeeze(1)[bad] // ncell
            overlap[which] = True
        return sig, overlap


# ---------------------------------------------------------------- potential on a grid
class GridPotential:
    """P1 potential on the vertices of the Kuhn box mesh, for the ARBD export.

    `sample_on_grid` evaluates the P1 function exactly as dolfinx Function.eval does: the Kuhn
    tetrahedron containing a point is the one whose step order is the order of the point's
    fractional coordinates within its cube, largest first.
    """

    def __init__(self, grid, phi):
        self.grid = grid
        self.phi = np.asarray(phi, dtype=np.float64)       # (nx+1, ny+1, nz+1) volts

    def sample_on_grid(self, x_range, y_range, z_range, fill_value=0.0):
        X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
        P = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        n, h, lo = np.array(self.grid.n), np.array(self.grid.h), np.array(self.grid.lo)
        rel = (P - lo) / h
        tol = 1e-9
        inside = ((rel >= -tol) & (rel <= n + tol)).all(1)
        cube = np.clip(np.floor(rel), 0, n - 1).astype(int)
        frac = np.clip(rel - cube, 0.0, 1.0)
        order = np.argsort(-frac, axis=1, kind="stable")
        rows = np.arange(len(P))
        v = cube.copy()
        out = self.phi[v[:, 0], v[:, 1], v[:, 2]].copy()
        for k in range(3):
            a = order[:, k]
            nxt = v.copy()
            nxt[rows, a] += 1
            out += frac[rows, a] * (self.phi[nxt[:, 0], nxt[:, 1], nxt[:, 2]]
                                    - self.phi[v[:, 0], v[:, 1], v[:, 2]])
            v = nxt
        out[~inside] = fill_value
        return out.reshape(X.shape)


# ---------------------------------------------------------------- scan driver
class _NoMeshSEM:
    """VerticalMovementSEM constructed without its DOLFINx mesh (analyte, pore and model only)."""

    @staticmethod
    def create(config, prepare_analyte=True):
        from . import vertical_movement_sem as vms
        from .cli import create_sem_from_config

        original = vms.VerticalMovementSEM.setup_dolfinx
        vms.VerticalMovementSEM.setup_dolfinx = lambda self: None
        try:
            sem, _ = create_sem_from_config(config, prepare_analyte=prepare_analyte)
        finally:
            vms.VerticalMovementSEM.setup_dolfinx = original
        return sem


class GPUSolver:
    """Batched GPU solver for one configuration (one analyte, one pore, one grid).

    It solves any list of (rotation matrix, z) cases; `run_scan` reproduces
    VerticalMovementSEM.run (the `run` mode) and `run` the fixed-z rotation scan.
    """

    def __init__(self, config, device="cuda", dtype=None, iters=None, open_iters=100,
                 compile=False, method="fem", current=None, prepare_analyte=True):
        if method not in ("fem", "fv"):
            raise ValueError("method must be 'fem' or 'fv'")
        current = current or ("flux" if method == "fem" else "power")
        if current not in ("flux", "power") or (method == "fv" and current == "flux"):
            raise ValueError("current must be 'flux' (fem only) or 'power'")
        # The flux reads phi next to the top electrode, where it is within ~0.1 % of the voltage,
        # so its error is ~1e4 times the potential error there: it needs a float64 solve (14
        # iterations: <= 2e-4 blockade points on the test set). The power current is quadratic
        # in the error and reaches <= 1e-3 points in float32 within 8 iterations.
        if dtype is None:
            dtype = torch.float64 if current == "flux" else torch.float32
        if iters is None:
            iters = 14 if current == "flux" else 8
        t0 = time.time()
        self.method = method
        self.current = current
        self.compiled = compile
        self._iterate_fn = torch.compile(self._iterate, dynamic=False) if compile else self._iterate
        self.sem = sem = _NoMeshSEM.create(config, prepare_analyte)
        if prepare_analyte and getattr(sem, "use_charges", False):
            raise NotImplementedError("GPU backend does not support use_charges")
        if sem.mesh_engine != "dolfinx":
            raise NotImplementedError("GPU backend supports the uniform dolfinx box mesh only")
        self.device = torch.device(device)
        self.dtype = dtype
        self.iters = iters
        self.grid = BoxGrid.from_sem(sem)
        self.voltage = float(sem.voltage)
        # the domain attributes setup_dolfinx would set (the ARBD export reads them)
        sem.domain_min = list(self.grid.lo)
        sem.domain_max = [lo + n * h for lo, n, h in zip(self.grid.lo, self.grid.n, self.grid.h)]
        sem.num_cells = list(self.grid.n)
        offsets = [kuhn_centroid(p) for p in KUHN] if method == "fem" else [(0.5, 0.5, 0.5)]
        self.sample_offsets = offsets
        self.base_sigma = []
        self.sigma = []
        has_analyte = prepare_analyte and len(sem.moving_radii) > 0
        membrane = sem.membrane_conductivity if sem.prevent_analyte_overlap else None
        for off in offsets:
            base = np.asarray(sem.base_cond_interp(self.grid.points(off)), dtype=np.float64)
            base[np.isnan(base)] = sem.bulk_conductivity          # as utils.loadFunc
            self.base_sigma.append(base)
            if has_analyte:
                self.sigma.append(AnalyteSigma(self.grid, base, sem.moving_radii, sem.cutoff,
                                               sem.conductivity_model, self.device, offset=off,
                                               membrane_conductivity=membrane))
        self.z = float(self.z_positions()[-1])  # run() reports the last z position
        self.timing = {"prepare_s": time.time() - t0}
        t1 = time.time()
        # open-pore solution in float64: every case is then solved for its change from phi0
        base = torch.as_tensor(np.stack(self.base_sigma), device=self.device)
        L0 = self._fine(base.reshape(1, len(offsets), *self.grid.n))
        levels, factors = setup_hierarchy(L0)
        b = rhs(L0, self.voltage)
        self.phi0, _ = pcg(levels, factors, b, torch.zeros_like(b), open_iters)
        self.open_current = float(self._current(L0, self.phi0)[0])
        self._sync()
        self.timing["open_pore_s"] = time.time() - t1

    # ------------------------------------------------------------ helpers
    def _sync(self):
        if self.device.type == "cuda":
            torch.cuda.synchronize()

    def z_positions(self):
        """The z positions VerticalMovementSEM.run visits."""
        sem = self.sem
        num_steps = int(abs(sem.z_end - sem.z_start) / sem.z_step) + 1
        direction = np.sign(sem.z_end - sem.z_start)
        return sem.z_start + np.arange(num_steps) * (direction * sem.z_step)

    def _fine(self, sig):
        """Finest level, in float64, from conductivity (B, k, nx, ny, nz), k = 6 (fem) or 1 (fv)."""
        sig = sig.to(torch.float64)
        return fem_level(sig, self.grid.h) if self.method == "fem" else fv_level(sig[:, 0], self.grid.h)

    def _current(self, L, x):
        if self.current == "flux":
            return flux_current(L, x, self.voltage)
        return power_current(L, x, self.voltage)

    def _prepare(self, sig):
        """float64 operator and the residual of phi0; the solve hierarchy in self.dtype.

        Solving A psi = b - A phi0 for psi = phi - phi0 keeps the float32 solve away from the
        cancellation in (voltage - phi) next to the top electrode, where phi is within ~0.1 % of
        the voltage; the current is then evaluated in float64 from phi0 + psi.
        """
        L = self._fine(sig)
        r0 = rhs(L, self.voltage) - apply(L, self.phi0.expand(L["D"].shape))
        Ls = {k: v.to(self.dtype) for k, v in L.items()}
        levels, factors = setup_hierarchy(Ls)
        return L, levels, factors, r0

    def vertex_potential(self, phi):
        """Full P1 vertex potential (nx+1, ny+1, nz+1), Dirichlet layers included (fem only)."""
        if self.method != "fem":
            raise NotImplementedError("vertex potentials exist for --gpu-method fem only")
        phi = phi.to(torch.float64)
        full = _pad(_pad(phi, 3, 1, 0), 3, 0, 1)
        full[..., -1] = self.voltage
        return full.cpu().numpy()

    def open_pore_potential(self):
        return GridPotential(self.grid, self.vertex_potential(self.phi0)[0])

    # ------------------------------------------------------------ cases
    def positions(self, rotation_matrix=None, z=None):
        """Atom centres for one case, computed exactly as VerticalMovementSEM.run does: the
        rotation about the centre of mass, then the centre of mass moved to z."""
        z = self.z if z is None else float(z)
        if rotation_matrix is not None:
            self.sem.set_analyte_rotation_matrix(rotation_matrix)
        try:
            pos = self.sem._get_rotated_analyte_positions()
        finally:
            if rotation_matrix is not None:
                self.sem.reset_analyte_rotation()
        return pos + np.array([0.0, 0.0, z - self.sem.moving_com[2]])

    def _radius_overlap(self, pos):
        """PETK's geometric overlap check (prevent_analyte_overlap with use_radius_overlap_check)."""
        from .vertical_movement_sem import AnalyteOverlapError
        sem = self.sem
        if not (sem.prevent_analyte_overlap and sem.use_radius_overlap_check):
            return None
        try:
            sem._assert_radius_overlap(pos, sem.moving_radii)
        except AnalyteOverlapError as exc:
            return str(exc)
        return None

    def conductivity(self, cases, chunk=16, with_overlap=False):
        """(B, k, nx, ny, nz) conductivity at the k sample points per cell (float64) for cases
        given as rotation matrices or (rotation matrix, z) pairs; with_overlap also returns a list
        with an overlap message (or None) per case."""
        if not self.sigma:
            raise RuntimeError("this solver was built without an analyte (prepare_analyte=False)")
        cases = [c if isinstance(c, tuple) else (c, None) for c in cases]
        out, notes = [], []
        for s in range(0, len(cases), chunk):
            pos = [self.positions(R, z) for R, z in cases[s:s + chunk]]
            parts = [sg.batch(pos, with_overlap=True) for sg in self.sigma]
            out.append(torch.stack([p[0] for p in parts], 1))
            if with_overlap:
                hard = torch.stack([p[1] for p in parts]).any(0).cpu().numpy()
                for P, h in zip(pos, hard):
                    note = self._radius_overlap(P)
                    if note is None and h:
                        note = "analyte hard core overlaps the membrane"
                    notes.append(note)
        sig = torch.cat(out).reshape(len(cases), len(self.sigma), *self.grid.n)
        return (sig, notes) if with_overlap else sig

    def _iterate(self, levels, factors, r):
        psi, _ = pcg(levels, factors, r, torch.zeros_like(r), self.iters)
        return psi

    def solve(self, sig, iters=None, track=False, return_phi=False):
        """CG in self.dtype for psi = phi - phi0; the current is evaluated in float64 from phi0 + psi.
        Returns (current, history) or, with return_phi, (current, history, phi)."""
        L, levels, factors, r0 = self._prepare(sig)
        r = r0.to(self.dtype)
        phi = lambda u: self.phi0 + u.to(torch.float64)
        if iters is None and not track:
            psi, hist = self._iterate_fn(levels, factors, r), None
        else:
            psi, hist = pcg(levels, factors, r, torch.zeros_like(r), self.iters if iters is None else iters,
                            lambda u: self._current(L, phi(u)), track)
        full = phi(psi)
        current = self._current(L, full)
        return (current, hist, full) if return_phi else (current, hist)

    def solve_cases(self, cases, batch_size=64, keep_phi=()):
        """Solve (rotation matrix, z) cases in batches.

        Returns a dict with the current per case (NaN where prevent_analyte_overlap rejects the
        case), the overlap message per case, per-case conductivity and solve times (the batch
        time shared evenly), and {index: GridPotential} for the cases in keep_phi.
        """
        keep_phi = set(keep_phi)
        n_all = len(cases)
        currents = np.full(n_all, np.nan)
        notes = [None] * n_all
        t_sig = np.zeros(n_all)
        t_solve = np.zeros(n_all)
        phis = {}
        for s in range(0, n_all, batch_size):
            idx = list(range(s, min(s + batch_size, n_all)))
            batch = [cases[i] for i in idx]
            if self.compiled and len(batch) < batch_size:  # keep one shape for the compiled solve
                batch = batch + [batch[-1]] * (batch_size - len(batch))
            t0 = time.time()
            sig, note = self.conductivity(batch, with_overlap=True)
            self._sync()
            t1 = time.time()
            want = any(i in keep_phi for i in idx)
            out = self.solve(sig, return_phi=want)
            cur = out[0].cpu().numpy()
            self._sync()
            t2 = time.time()
            for j, i in enumerate(idx):
                notes[i] = note[j]
                currents[i] = np.nan if note[j] else cur[j]
                t_sig[i] = (t1 - t0) / len(idx)
                t_solve[i] = (t2 - t1) / len(idx)
                if want and i in keep_phi and not note[j]:
                    phis[i] = GridPotential(self.grid, self.vertex_potential(out[2][j:j + 1])[0])
        return {"currents": currents, "overlap": notes, "sigma_s": t_sig, "solve_s": t_solve, "phi": phis}

    def run(self, rotation_matrices, batch_size=64):
        """Currents (nA) at the last z position for each rotation matrix, in order."""
        t0 = time.time()
        res = self.solve_cases([(R, None) for R in rotation_matrices], batch_size)
        self.timing.update(sigma_s=float(res["sigma_s"].sum()), solve_s=float(res["solve_s"].sum()),
                           n=len(rotation_matrices), run_s=time.time() - t0)
        return res["currents"]

    # ------------------------------------------------------------ the `run` mode
    def run_scan(self, batch_size=64):
        """VerticalMovementSEM.run on the GPU: the z scan of the analyte in its PDB orientation.

        Writes the same files (<prefix>_results.txt, and the per-position, timing-analysis and
        timing-summary files when verbose_output is on), runs the ARBD export when configured,
        and returns the same results dict.
        """
        from .arbd_export import export_arbd_grids, should_export_at_step
        sem = self.sem
        t_start = time.time()
        z = self.z_positions()
        arbd = getattr(sem, "arbd_export_config", None)
        if arbd and self.method != "fem":
            raise NotImplementedError("ARBD export with the GPU backend needs --gpu-method fem")
        if arbd:
            export_arbd_grids(sem, self.open_pore_potential(), None, arbd)
        keep = [i for i in range(len(z)) if arbd and should_export_at_step(i, arbd)]
        t_sim = time.time()
        res = self.solve_cases([(None, zz) for zz in z], batch_size, keep_phi=keep)
        for i, note in enumerate(res["overlap"]):
            if note:
                logger.warning("Skipping position %d/%d (Z=%.2f A) due to analyte overlap: %s",
                               i + 1, len(z), z[i], note)
        for i in keep:
            if i in res["phi"]:
                export_arbd_grids(sem, res["phi"][i], float(z[i]), arbd)
        simulation_time = time.time() - t_sim
        return write_run_outputs(sem, z, res, self.open_current, self.timing["open_pore_s"],
                                 simulation_time, time.time() - t_start + self.timing["prepare_s"])


GPURotationScan = GPUSolver  # earlier name


def write_run_outputs(sem, z_positions, res, open_current, open_pore_time, simulation_time, total_time):
    """The output files and results dict of VerticalMovementSEM.run, from a GPU z scan."""
    currents = np.asarray(res["currents"], dtype=float)
    ok = np.isfinite(currents)
    cond_t = np.where(ok, res["sigma_s"], np.nan)
    solve_t = np.where(ok, res["solve_s"], np.nan)
    mesh_t = np.where(ok, 0.0, np.nan)
    pos_t = cond_t + solve_t
    normalized = currents / open_current
    blockages = (1 - normalized) * 100
    n = len(z_positions)

    def fmt(seconds):
        if seconds < 60:
            return f"{seconds:.1f}s"
        if seconds < 3600:
            return f"{seconds / 60:.1f}m"
        return f"{seconds / 3600:.1f}h"

    avg_pos = float(np.nanmean(pos_t)) if ok.any() else float("nan")
    per_hour = 3600 / avg_pos if ok.any() and avg_pos > 0 else 0.0
    prefix = sem.output_prefix
    if sem.verbose_output:
        for i in np.flatnonzero(ok):
            with open(f"{prefix}_position_{i:04d}.dat", "w") as f:
                f.write(f"# Position {i + 1}/{n}\n")
                f.write(f"# Z_position: {z_positions[i]} Å\n")
                f.write(f"# Current: {currents[i]:.6e} nA\n")
                f.write(f"# Blockage: {blockages[i]:.2f}%\n")
                f.write(f"# Mesh_time: {mesh_t[i]:.3f} s\n")
                f.write(f"# Conductivity_time: {cond_t[i]:.3f} s\n")
                f.write(f"# Solver_time: {solve_t[i]:.3f} s\n")
                f.write(f"# Total_time: {pos_t[i]:.3f} s\n")
                f.write(f"# Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"{z_positions[i]:.1f} {currents[i]:.6e} {blockages[i]:.2f} {mesh_t[i]:.3f} "
                        f"{cond_t[i]:.3f} {solve_t[i]:.3f} {pos_t[i]:.3f}\n")
    header = "\n".join(f"# {line}" for line in [
        "Z_position(Å) Current(nA) Normalized_Current Blockage(%) Mesh_time(s) Conductivity_time(s) "
        "Solver_time(s) Total_time(s)",
        f"Pore_type: {sem.pore_type}",
        f"Open_pore_current: {open_current:.6e} A",
        f"Total_simulation_time: {fmt(total_time)}",
        f"Average_time_per_position: {avg_pos:.3f}s",
        f"Throughput: {per_hour:.1f} positions/hour",
        f"Use_VdW_radii: {sem.use_vdw_radii}",
        "Backend: gpu",
    ])
    np.savetxt(f"{prefix}_results.txt",
               np.column_stack([z_positions, currents, normalized, blockages, mesh_t, cond_t, solve_t, pos_t]),
               header=header, fmt=["%.1f", "%.6e", "%.6f", "%.2f", "%.3f", "%.3f", "%.3f", "%.3f"])
    if sem.verbose_output:
        timing_header = "\n".join(f"# {line}" for line in [
            "COMPREHENSIVE TIMING ANALYSIS (GPU backend: per-position times are batch times shared evenly)",
            f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
            f"Total positions: {n}",
            f"Total simulation time: {fmt(total_time)}",
            f"Open pore calculation: {fmt(open_pore_time)}",
            f"Main simulation: {fmt(simulation_time)}",
            f"Average per position: {avg_pos:.3f}s",
            f"Throughput: {per_hour:.1f} positions/hour",
            "",
            "Position Z_position(Å) Mesh_time(s) Conductivity_time(s) Solver_time(s) Total_time(s)",
        ])
        np.savetxt(f"{prefix}_timing_analysis.txt",
                   np.column_stack([np.arange(1, n + 1), z_positions, mesh_t, cond_t, solve_t, pos_t]),
                   header=timing_header, fmt=["%d", "%.1f", "%.6f", "%.6f", "%.6f", "%.6f"])
        with open(f"{prefix}_timing_summary.txt", "w") as f:
            f.write("SEM Simulation Timing Summary (GPU backend)\n")
            f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"{'=' * 50}\n\n")
            f.write(f"Total Time: {fmt(total_time)}\n")
            f.write(f"Positions: {n}\n")
            f.write(f"Throughput: {per_hour:.1f} positions/hour\n")
            f.write(f"Open Pore: {fmt(open_pore_time)}\n")
            f.write(f"Simulation: {fmt(simulation_time)}\n\n")
            f.write("Per Position (average):\n")
            f.write(f"  Conductivity: {np.nanmean(cond_t):.4f}s\n")
            f.write(f"  GPU Solver:   {np.nanmean(solve_t):.4f}s\n")
            f.write(f"  Total:        {avg_pos:.4f}s\n")
    if ok.any():
        k = int(np.nanargmax(blockages))
        logger.info("Simulation complete! Maximum blockage: %.1f%% at Z = %.1f A", blockages[k], z_positions[k])
    logger.info("Results saved to: %s_results.txt", prefix)
    return {
        "z_positions": z_positions,
        "currents": currents,
        "normalized_currents": normalized,
        "blockages": blockages,
        "open_current": open_current,
        "pore_type": sem.pore_type,
        "pore_radius": sem.pore_radius,
        "outer_radius": sem.outer_radius if sem.pore_type == "double_cone" else None,
        "corner_radius": sem.corner_radius if sem.pore_type == "cylindrical" else None,
        "biological_pore_pdb": sem.biological_pore_pdb if sem.pore_type == "biological" else None,
        "bin_file_path": sem.bin_file_path if sem.pore_type == "bin_file" else None,
        "use_vdw_radii": sem.use_vdw_radii,
        "timing": {
            "total_time": total_time,
            "open_pore_time": open_pore_time,
            "simulation_time": simulation_time,
            "position_times": pos_t,
            "mesh_times": mesh_t,
            "conductivity_times": cond_t,
            "solver_times": solve_t,
            "avg_position_time": avg_pos,
            "positions_per_hour": per_hour,
        },
    }
