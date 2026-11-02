"""
Batched GPU backend for SEM rotation scans.

Solves -div(sigma grad phi) = 0 on the configured box with phi = voltage on the top face,
phi = 0 on the bottom face and insulating side faces, for many analyte orientations at once.

- Discretisation: cell-centred finite volumes on the box grid (one cell per DOLFINx box cell),
  harmonic-mean face conductances.
- Solver: conjugate gradients preconditioned with a geometric multigrid V-cycle
  (series-parallel coarsening by 2 or 3 per axis, weighted Jacobi smoothing, dense coarsest
  solve), started from the open-pore solution and run for a fixed number of iterations.
- Current: dissipated power divided by the voltage, whose error is quadratic in the
  solution error.

Conductivity follows the FEM path exactly: base (pore) conductivity is the pore's own
interpolator evaluated at the cell centres, and the analyte modification reproduces
VerticalMovementSEM.calculate_analyte_conductivity_modification (nearest atom centre within
`cutoff`, minus that atom's radius, then the conductivity model, then the minimum with the
base). Units are PETK's: Angstrom, volts, conductivity in nS/A, current in nA.
"""

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
DENSE_MAX = 1000  # largest coarsest level solved with a dense inverse


# ---------------------------------------------------------------- grid
@dataclass
class BoxGrid:
    """Cell-centred grid matching the DOLFINx box mesh of VerticalMovementSEM.setup_dolfinx."""
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

    def axis_centres(self, a):
        return self.lo[a] + (np.arange(self.n[a]) + 0.5) * self.h[a]

    def centres(self):
        """(N, 3) cell centres in C order (x slowest)."""
        X, Y, Z = np.meshgrid(*(self.axis_centres(a) for a in range(3)), indexing="ij")
        return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)


# ---------------------------------------------------------------- operator
def _harm(a, b):
    return 2.0 * a * b / (a + b)


def _pad(t, axis, before, after):
    """Zero-pad a (B, nx, ny, nz) tensor along spatial axis 1..3."""
    pads = [0, 0, 0, 0, 0, 0]
    k = 2 * (3 - axis)
    pads[k], pads[k + 1] = before, after
    return F.pad(t, pads)


def _diag(L):
    D = (_pad(L["gx"], 1, 1, 0) + _pad(L["gx"], 1, 0, 1)
         + _pad(L["gy"], 2, 1, 0) + _pad(L["gy"], 2, 0, 1)
         + _pad(L["gz"], 3, 1, 0) + _pad(L["gz"], 3, 0, 1))
    nz = D.shape[3]
    D = D + _pad(L["gt"][..., None], 3, nz - 1, 0) + _pad(L["gb"][..., None], 3, 0, nz - 1)
    L["D"] = D
    L["Dinv"] = 1.0 / D


def fine_level(sig, h):
    """Face conductances of the finest level from cell conductivity sig (B, nx, ny, nz)."""
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


def apply(L, u):
    gx, gy, gz = L["gx"], L["gy"], L["gz"]
    return (L["D"] * u
            - _pad(gx * u[:, 1:], 1, 0, 1) - _pad(gx * u[:, :-1], 1, 1, 0)
            - _pad(gy * u[:, :, 1:], 2, 0, 1) - _pad(gy * u[:, :, :-1], 2, 1, 0)
            - _pad(gz * u[..., 1:], 3, 0, 1) - _pad(gz * u[..., :-1], 3, 1, 0))


def _block_sum(t, axis, f):
    if f == 1:
        return t
    n = t.shape[axis]
    return t.unflatten(axis, (n // f, f)).sum(axis + 1)


def _series(gl, axis, f):
    """Coarse link conductances along `axis` from fine links gl (length n-1 on that axis)."""
    if f == 1:
        return gl
    g = gl.movedim(axis, -1)
    n = g.shape[-1] + 1
    if f == 2:   # coarse centre on the face between the pair: half + full + half link
        R = 0.5 / g[..., 0:n - 3:2] + 1.0 / g[..., 1:n - 2:2] + 0.5 / g[..., 2:n - 1:2]
    else:        # f == 3, coarse centre on the middle fine cell: three full links
        R = 1.0 / g[..., 1:n - 4:3] + 1.0 / g[..., 2:n - 3:3] + 1.0 / g[..., 3:n - 2:3]
    return (1.0 / R).movedim(-1, axis)


def _boundary(gz_lat, gb_lat, f, top):
    """Coarse boundary conductance on the top or bottom face."""
    if f == 1:
        return gb_lat
    link = gz_lat[..., -1] if top else gz_lat[..., 0]
    w = 0.5 if f == 2 else 1.0
    return 1.0 / (w / link + 1.0 / gb_lat)


def _factors(n):
    out = []
    for k in n:
        if k <= 6:
            out.append(1)
        elif k % 2 == 0:
            out.append(2)
        elif k % 3 == 0:
            out.append(3)
        else:
            out.append(1)
    return tuple(out)


def coarsen(L, f):
    fx, fy, fz = f
    gx = _block_sum(_block_sum(L["gx"], 2, fy), 3, fz)
    gy = _block_sum(_block_sum(L["gy"], 1, fx), 3, fz)
    gz = _block_sum(_block_sum(L["gz"], 1, fx), 2, fy)
    gt = _block_sum(_block_sum(L["gt"], 1, fx), 2, fy)
    gb = _block_sum(_block_sum(L["gb"], 1, fx), 2, fy)
    C = {
        "gx": _series(gx, 1, fx),
        "gy": _series(gy, 2, fy),
        "gz": _series(gz, 3, fz),
        "gt": _boundary(gz, gt, fz, top=True),
        "gb": _boundary(gz, gb, fz, top=False),
    }
    _diag(C)
    return C


def _restrict(r, f):
    for a, k in zip((1, 2, 3), f):
        r = _block_sum(r, a, k)
    return r


def _prolong(e, f):
    for a, k in zip((1, 2, 3), f):
        if k > 1:
            e = e.repeat_interleave(k, a)
    return e


def _dense_inverse(L):
    B, nx, ny, nz = L["D"].shape
    k = nx * ny * nz
    eye = torch.eye(k, device=L["D"].device, dtype=L["D"].dtype)
    cols = [apply(L, eye[j].expand(B, k).reshape(B, nx, ny, nz)).reshape(B, k) for j in range(k)]
    A = torch.stack(cols, dim=2).cpu().double()
    L["Ainv"] = torch.linalg.inv(A).to(L["D"].device, L["D"].dtype)


def setup_hierarchy(sig, h):
    levels, factors = [fine_level(sig, h)], []
    while True:
        n = levels[-1]["D"].shape[1:]
        f = _factors(n)
        if f == (1, 1, 1):
            break
        factors.append(f)
        levels.append(coarsen(levels[-1], f))
    k = int(np.prod(levels[-1]["D"].shape[1:]))
    if k > DENSE_MAX:
        raise ValueError(f"grid {tuple(sig.shape[1:])} coarsens only to {k} cells; "
                         f"choose cell counts with factors 2 and 3")
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
    x = x + _prolong(vcycle(levels, factors, l + 1, rc), factors[l])
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


def pcg(levels, factors, b, x, iters, voltage=None, track=False):
    """Fixed-iteration multigrid-preconditioned CG; optionally records the current per iteration."""
    L, red = levels[0], (1, 2, 3)
    r = b - apply(L, x)
    z = vcycle(levels, factors, 0, r)
    p = z.clone()
    rz = (r * z).sum(red, keepdim=True)
    hist = [power_current(L, x, voltage)] if track else None
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
            hist.append(power_current(L, x, voltage))
    return x, hist


# ---------------------------------------------------------------- analyte conductivity
class AnalyteSigma:
    """Exact GPU copy of calculate_analyte_conductivity_modification on a BoxGrid."""

    def __init__(self, grid, base_sigma, radii, cutoff, model, device):
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
        k = [int(math.ceil(self.cutoff / hh)) + 1 for hh in grid.h]
        offs = torch.stack(torch.meshgrid(*(torch.arange(-kk + 1, kk + 1, device=device) for kk in k),
                                          indexing="ij"), -1).reshape(-1, 3)
        self.offsets = offs  # candidate cells around each atom

    def model(self, s):
        frac = ((s - self.min_distance) / (self.model_cutoff - self.min_distance)).clamp(0, 1)
        return torch.where(s < self.model_cutoff, 1e-7 * self.bulk + frac * self.bulk,
                           torch.full_like(s, self.bulk))

    def __call__(self, positions):
        """positions: (M, 3) float64 array of atom centres -> flat cell conductivity (float64)."""
        p = torch.as_tensor(positions, dtype=torch.float64, device=self.device)
        base_idx = torch.floor((p - self.lo) / self.h - 0.5).long()        # cell centre just below
        cells = base_idx[:, None, :] + self.offsets[None]                   # (M, K, 3)
        inside = ((cells >= 0) & (cells < self.n)).all(-1)
        centre = self.lo + (cells.double() + 0.5) * self.h
        d = torch.sqrt(((centre - p[:, None, :]) ** 2).sum(-1))
        keep = inside & (d < self.cutoff)
        atom = torch.arange(len(p), device=self.device)[:, None].expand_as(d)[keep]
        d = d[keep]
        c = cells[keep]
        flat = (c[:, 0] * self.grid.n[1] + c[:, 1]) * self.grid.n[2] + c[:, 2]
        ncell = int(np.prod(self.grid.n))
        dmin = torch.full((ncell,), float("inf"), dtype=torch.float64, device=self.device)
        dmin = dmin.scatter_reduce(0, flat, d, reduce="amin")
        nearest = d == dmin[flat]
        r = torch.full((ncell,), float("inf"), dtype=torch.float64, device=self.device)
        r = r.scatter_reduce(0, flat[nearest], self.radii[atom[nearest]], reduce="amin")
        hit = torch.isfinite(dmin)
        s = (dmin[hit] - r[hit]).clamp_min(0.0)
        sig = self.base.clone()
        sig[hit] = torch.minimum(sig[hit], self.model(s))
        return sig


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

    def __init__(self, config, device="cuda", dtype=torch.float32, iters=6, open_iters=60):
        t0 = time.time()
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
        centres = self.grid.centres()
        base = np.asarray(sem.base_cond_interp(centres), dtype=np.float64)
        base[np.isnan(base)] = sem.bulk_conductivity          # as utils.loadFunc
        self.base_sigma = base
        self.sigma = AnalyteSigma(self.grid, base, sem.moving_radii, sem.cutoff,
                                  sem.conductivity_model, self.device)
        num_steps = int(abs(sem.z_end - sem.z_start) / sem.z_step) + 1
        direction = np.sign(sem.z_end - sem.z_start)
        self.z = float(sem.z_start + (num_steps - 1) * direction * sem.z_step)  # run() keeps the last z
        self.timing = {"prepare_s": time.time() - t0}
        t1 = time.time()
        sig0 = torch.as_tensor(base, device=self.device).reshape(1, *self.grid.n).to(dtype)
        levels, factors = setup_hierarchy(sig0, self.grid.h)
        b = rhs(levels[0], self.voltage)
        self.phi0, _ = pcg(levels, factors, b, torch.zeros_like(b), open_iters)
        self.open_current = float(power_current(levels[0], self.phi0, self.voltage)[0])
        self._sync()
        self.timing["open_pore_s"] = time.time() - t1

    def _sync(self):
        if self.device.type == "cuda":
            torch.cuda.synchronize()

    def positions(self, rotation_matrix):
        """Atom centres for one orientation, computed exactly as VerticalMovementSEM.run does."""
        self.sem.set_analyte_rotation_matrix(rotation_matrix)
        pos = self.sem._get_rotated_analyte_positions()
        self.sem.reset_analyte_rotation()
        return pos + np.array([0.0, 0.0, self.z - self.sem.moving_com[2]])

    def conductivity(self, rotation_matrices):
        sig = [self.sigma(self.positions(R)) for R in rotation_matrices]
        return torch.stack(sig).reshape(len(sig), *self.grid.n)

    def solve(self, sig, iters=None, track=False):
        levels, factors = setup_hierarchy(sig.to(self.dtype), self.grid.h)
        b = rhs(levels[0], self.voltage)
        x0 = self.phi0.expand_as(b).clone()
        x, hist = pcg(levels, factors, b, x0, self.iters if iters is None else iters,
                      self.voltage, track)
        return power_current(levels[0], x, self.voltage), hist

    def run(self, rotation_matrices, batch_size=64):
        """Currents (nA) for each rotation matrix, in order."""
        out, t_sig, t_solve = [], 0.0, 0.0
        for s in range(0, len(rotation_matrices), batch_size):
            t0 = time.time()
            sig = self.conductivity(rotation_matrices[s:s + batch_size])
            self._sync()
            t1 = time.time()
            current, _ = self.solve(sig)
            out.append(current.cpu().numpy())
            t2 = time.time()
            t_sig += t1 - t0
            t_solve += t2 - t1
        self.timing.update(sigma_s=t_sig, solve_s=t_solve, n=len(rotation_matrices))
        return np.concatenate(out) if out else np.zeros(0)
