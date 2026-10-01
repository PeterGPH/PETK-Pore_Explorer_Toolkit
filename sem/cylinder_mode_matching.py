"""
Exact open-pore resistance of a SHARP-walled cylindrical pore by Ritz mode
matching (reference for :mod:`sem.open_pore_theory`).

Geometry: a cylinder of radius a and length L through an insulating membrane
between two unbounded half-space reservoirs (conductivity 1). Writing

    R = rho (L + 2 kappa(L/a) a) / (pi a^2),

this module computes the end-correction coefficient kappa(L/a) and the
amplitude k of the r^(2/3) potential singularity at the rim,

    phi ~ K rho^(2/3) cos(2 theta / 3),   K = k E a^(1/3),   E = I / (sigma pi a^2),

(rho: distance from the rim edge, theta measured from the membrane face).

Method: the unknown is the axial flux density w(r) through the pore mouth,
expanded as w = sum_m c_m (1 - r^2/a^2)^(nu_m) with nu_m = (2m - 1)/3, which
carries the exact r^(-1/3) edge behaviour of the 270-degree fluid corner and
its sub-leading powers. By Thomson's principle the dissipation at fixed total
current is minimised by the true flux, so the Ritz value is an upper bound on
the end resistance that converges from above. The half-space dissipation has
a closed (Weber-Schafheitlin) form; the tube side is a sum over the Neumann
modes J0(j_{1,n} r/a) of the cylinder, with tanh(j L / 2a) for the finite
half-length (the midplane is equipotential by antisymmetry) and a Hurwitz-
zeta tail correction. Converges to ~1e-9 in kappa for L/a >~ 0.05 with the
default 12 basis functions; for thinner membranes k converges slowly (it
approaches the thin-plate law 0.5607 (L/a)^(-1/6)).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import gammaln, j0, j1, jn_zeros, jv, zeta

#: Flanged-tube end correction per end, units of a (kappa(L/a -> inf)).
KAPPA_INF_EXACT = 0.8216742


def _j1_zeros(n: int) -> np.ndarray:
    """First n positive zeros of J1 (McMahon start, Newton polish)."""
    k = np.arange(1, n + 1, dtype=float)
    beta = (k + 0.25) * np.pi
    x = beta - 3.0 / (8.0 * beta) + 36.0 / (1536.0 * beta ** 3)
    m = min(n, 200)
    x[:m] = jn_zeros(1, m)
    for _ in range(4):
        f = j1(x)
        x = x - f / (j0(x) - f / x)
    return x


def _hankel_basis(nu: float, k):
    """Order-0 Hankel transform of (1 - r^2)^nu on r < 1."""
    k = np.asarray(k, float)
    c = np.exp(nu * np.log(2.0) + gammaln(nu + 1.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = c * jv(nu + 1.0, k) / k ** (nu + 1.0)
    return np.where(k == 0, 1.0 / (2.0 * (nu + 1.0)), out)


def _half_space_entry(nu1: float, nu2: float) -> float:
    """2 pi int_0^inf W_nu1(k) W_nu2(k) dk in closed form."""
    mu, mup = nu1 + 1.0, nu2 + 1.0
    lam = mu + mup
    log_ws = (gammaln(lam) + 0.5 * np.log(np.pi) - lam * np.log(2.0)
              - gammaln(mup + 0.5) - gammaln(mu + mup + 0.5) - gammaln(mu + 0.5))
    pref = (nu1 + nu2) * np.log(2.0) + gammaln(nu1 + 1.0) + gammaln(nu2 + 1.0)
    return 2.0 * np.pi * np.exp(pref + log_ws)


class _Ritz:
    def __init__(self, n_basis: int = 12, n_modes: int = 200000):
        self.nus = np.array([(2 * m - 1) / 3.0 for m in range(n_basis)])
        self.jz = _j1_zeros(n_modes)
        self.j0sq = j0(self.jz) ** 2
        self.H = np.array([[_half_space_entry(p, q) for q in self.nus] for p in self.nus])
        self.g = 1.0 / (2.0 * (self.nus + 1.0))      # int_0^1 (1 - r^2)^nu r dr
        self.W = np.array([_hankel_basis(nu, self.jz) for nu in self.nus])

    def tube(self, half_length: float) -> np.ndarray:
        n = len(self.nus)
        th = np.tanh(self.jz * half_length) if np.isfinite(half_length) else 1.0
        base = 4.0 * np.pi * th / (self.j0sq * self.jz)
        T = np.zeros((n, n))
        nn = np.arange(len(self.jz) - 63, len(self.jz) + 1, dtype=float)
        for i in range(n):
            for j in range(i, n):
                terms = base * self.W[i] * self.W[j]
                q = self.nus[i] + self.nus[j] + 3.0           # terms ~ C (n + 1/4)^-q
                C = np.mean(terms[-64:] * (nn + 0.25) ** q)
                T[i, j] = T[j, i] = terms.sum() + C * zeta(q, len(self.jz) + 1.25)
        return T

    def solve(self, lam: float):
        mat = self.H + (self.tube(0.5 * lam) if lam > 0 else 0.0)
        x = np.linalg.solve(mat, self.g)
        gx = self.g @ x
        kappa = 1.0 / (4.0 * np.pi * gx)              # = pi * R_end (rho = a = 1)
        # w ~ c_0 (2 rho / a)^(-1/3) near the rim and dphi/dz = (sqrt(3)/3) K rho^(-1/3)
        k_rim = np.sqrt(3.0) * 2.0 ** (-1.0 / 3.0) * abs(x[0]) / (2.0 * gx)
        return kappa, k_rim


@lru_cache(maxsize=4)
def _ritz(n_basis: int, n_modes: int) -> _Ritz:
    return _Ritz(n_basis, n_modes)


def sharp_cylinder_coefficients(L_over_a, n_basis: int = 12, n_modes: int = 200000):
    """(kappa, k_rim) of a sharp-walled cylinder for each L/a (array-like)."""
    rz = _ritz(n_basis, n_modes)
    lam = np.atleast_1d(np.asarray(L_over_a, float))
    out = np.array([rz.solve(x) for x in lam])
    return out[:, 0], out[:, 1]


def sharp_cylinder_resistance_exact(d, L, sigma: float = 1.0, **kwargs):
    """Exact (Ritz) resistance in Ohm of a sharp cylinder; d, L in Angstrom."""
    a = 0.5 * float(d)
    kappa, _ = sharp_cylinder_coefficients(float(L) / a, **kwargs)
    return (float(L) + 2.0 * kappa[0] * a) / (np.pi * a ** 2) * 1e10 / sigma
