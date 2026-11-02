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
    grid cell (at `offset` in cell units)."""

    def __init__(self, grid, base_sigma, radii, cutoff, model, device, offset=(0.5, 0.5, 0.5)):
        self.grid = grid
        self.device = device
        self.base = torch.as_tensor(base_sigma, dtype=torch.float64, device=device)
        self.radii = torch.as_tensor(np.asarray(radii, dtype=np.float64), device=device)
        self.cutoff = float(cutoff)
        self.min_distance = float(model.min_distance)
        self.model_cutoff = float(model.cutoff)
        self.bulk = float(model.bulk_conductivity)
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

    def batch(self, positions):
        """positions: list of (M, 3) arrays, one per orientation -> (B, ncell) conductivity."""
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
        sig[hit] = torch.minimum(sig[hit], self.model(s))
        return sig.reshape(B, ncell)


# ---------------------------------------------------------------- scan driver
class _NoMeshSEM:
    """VerticalMovementSEM constructed without its DOLFINx mesh (analyte, pore and model only)."""

    @staticmethod
    def create(config):
        from . import vertical_movement_sem as vms
        from .cli import create_sem_from_config

        original = vms.VerticalMovementSEM.setup_dolfinx
        vms.VerticalMovementSEM.setup_dolfinx = lambda self: None
        try:
            sem, _ = create_sem_from_config(config, prepare_analyte=True)
        finally:
            vms.VerticalMovementSEM.setup_dolfinx = original
        return sem


class GPURotationScan:
    """Batched rotation scan for one protein and one pore configuration."""

    def __init__(self, config, device="cuda", dtype=torch.float32, iters=12, open_iters=100,
                 compile=False, method="fem", current=None, refine=None):
        if method not in ("fem", "fv"):
            raise ValueError("method must be 'fem' or 'fv'")
        current = current or ("flux" if method == "fem" else "power")
        if current not in ("flux", "power") or (method == "fv" and current == "flux"):
            raise ValueError("current must be 'flux' (fem only) or 'power'")
        t0 = time.time()
        self.method = method
        self.current = current
        # refinement passes: the flux reads phi next to the top electrode with ~1e4 amplification,
        # so it needs two float32 passes; the power current is quadratic in the error and needs one
        self.refine = refine or (2 if current == "flux" and dtype != torch.float64 else 1)
        self.iters = iters
        self._inner = self._steps(iters)[0]
        self.compiled = compile
        self._iterate_fn = torch.compile(self._iterate, dynamic=False) if compile else self._iterate
        self.sem = sem = _NoMeshSEM.create(config)
        for flag in ("prevent_analyte_overlap", "use_charges"):
            if getattr(sem, flag, False):
                raise NotImplementedError(f"GPU backend does not support {flag}")
        if sem.mesh_engine != "dolfinx":
            raise NotImplementedError("GPU backend supports the uniform dolfinx box mesh only")
        self.device = torch.device(device)
        self.dtype = dtype
        self.iters = iters
        self.grid = BoxGrid.from_sem(sem)
        self.voltage = float(sem.voltage)
        offsets = [kuhn_centroid(p) for p in KUHN] if method == "fem" else [(0.5, 0.5, 0.5)]
        self.sample_offsets = offsets
        self.base_sigma = []
        self.sigma = []
        for off in offsets:
            base = np.asarray(sem.base_cond_interp(self.grid.points(off)), dtype=np.float64)
            base[np.isnan(base)] = sem.bulk_conductivity          # as utils.loadFunc
            self.base_sigma.append(base)
            self.sigma.append(AnalyteSigma(self.grid, base, sem.moving_radii, sem.cutoff,
                                           sem.conductivity_model, self.device, offset=off))
        num_steps = int(abs(sem.z_end - sem.z_start) / sem.z_step) + 1
        direction = np.sign(sem.z_end - sem.z_start)
        self.z = float(sem.z_start + (num_steps - 1) * direction * sem.z_step)  # run() keeps the last z
        self.timing = {"prepare_s": time.time() - t0}
        t1 = time.time()
        # open-pore solution in float64: every orientation is then solved for its change from phi0
        base = torch.as_tensor(np.stack(self.base_sigma), device=self.device)
        L0 = self._fine(base.reshape(1, len(offsets), *self.grid.n))
        levels, factors = setup_hierarchy(L0)
        b = rhs(L0, self.voltage)
        self.phi0, _ = pcg(levels, factors, b, torch.zeros_like(b), open_iters)
        self.open_current = float(self._current(L0, self.phi0)[0])
        self._sync()
        self.timing["open_pore_s"] = time.time() - t1

    def _sync(self):
        if self.device.type == "cuda":
            torch.cuda.synchronize()

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

    def positions(self, rotation_matrix):
        """Atom centres for one orientation, computed exactly as VerticalMovementSEM.run does."""
        self.sem.set_analyte_rotation_matrix(rotation_matrix)
        pos = self.sem._get_rotated_analyte_positions()
        self.sem.reset_analyte_rotation()
        return pos + np.array([0.0, 0.0, self.z - self.sem.moving_com[2]])

    def conductivity(self, rotation_matrices, chunk=16):
        """(B, k, nx, ny, nz) conductivity at the k sample points per cell (float64)."""
        out = []
        for s in range(0, len(rotation_matrices), chunk):
            pos = [self.positions(R) for R in rotation_matrices[s:s + chunk]]
            out.append(torch.stack([sg.batch(pos) for sg in self.sigma], 1))
        return torch.cat(out).reshape(len(rotation_matrices), len(self.sigma), *self.grid.n)

    def _iterate(self, levels, factors, r):
        d, _ = pcg(levels, factors, r, torch.zeros_like(r), self._inner)
        return d

    def _steps(self, iters):
        """Iterations per refinement pass (the total is split as evenly as possible)."""
        n = max(1, min(self.refine, iters))
        return [iters // n + (1 if k < iters % n else 0) for k in range(n)]

    def solve(self, sig, iters=None, track=False):
        """Mixed-precision iterative refinement: each pass runs CG in self.dtype on the float64
        residual of the previous one, so the float32 rounding floor shrinks with the residual."""
        L, levels, factors, r0 = self._prepare(sig)
        psi = torch.zeros_like(r0, dtype=torch.float64)
        r, hist = r0, ([] if track else None)
        steps = self._steps(self.iters if iters is None else iters)
        for k, n in enumerate(steps):
            rs = r.to(self.dtype)
            if not track and n == self._inner:
                d = self._iterate_fn(levels, factors, rs)
            else:
                base = psi
                d, h = pcg(levels, factors, rs, torch.zeros_like(rs), n,
                           lambda u: self._current(L, self.phi0 + base + u.to(torch.float64)), track)
                if track:
                    hist.extend(h if k == 0 else h[1:])
            psi = psi + d.to(torch.float64)
            if k < len(steps) - 1:
                r = r0 - apply(L, psi)
        return self._current(L, self.phi0 + psi), hist

    def run(self, rotation_matrices, batch_size=64):
        """Currents (nA) for each rotation matrix, in order."""
        out, t_sig, t_solve = [], 0.0, 0.0
        for s in range(0, len(rotation_matrices), batch_size):
            t0 = time.time()
            mats = rotation_matrices[s:s + batch_size]
            n = len(mats)
            if self.compiled and n < batch_size:  # keep one shape so the compiled solve is reused
                mats = list(mats) + [mats[-1]] * (batch_size - n)
            sig = self.conductivity(mats)
            self._sync()
            t1 = time.time()
            current, _ = self.solve(sig)
            out.append(current[:n].cpu().numpy())
            t2 = time.time()
            t_sig += t1 - t0
            t_solve += t2 - t1
        self.timing.update(sigma_s=t_sig, solve_s=t_solve, n=len(rotation_matrices))
        return np.concatenate(out) if out else np.zeros(0)
