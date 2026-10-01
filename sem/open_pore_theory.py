"""
Closed-form open-pore conductance of cylindrical solid-state pores, and the
mapping between the SEM soft-wall conductivity map and sharp-wall theory.

Why this module exists
----------------------
SEM does not place a sharp insulating wall at the nominal pore surface: the
conductivity ramps linearly from 0 at ``r_min`` to bulk at ``r_cut`` in the
Euclidean distance D to the membrane (:func:`sem.utils.condfrac`). Comparing
an SEM conductance with a sharp-wall formula therefore needs (i) the
sharp-wall geometry that is electrically equivalent to the ramped map and
(ii) a sharp-wall formula that is itself accurate. The usual choice,

    R = rho * (4 L / (pi d^2) + 1 / d)            (Hall access + channel)

evaluated at d - 2t, L + 2t with t = (r_min + r_cut)/2 = 2.7 A, gets (i)
almost right but not (ii), and is 2-3 % off for d ~ L.

What the reference solver (:mod:`sem.axisym_conductance`) shows
------------------------------------------------------------------
1. Channel: inside the bore the field is exactly axial, so the equal-area
   (sigma-weighted) radius ``a_e`` is exact for the channel resistance.
   Its curvature term is tiny (a - a_e = 2.685 A at a = 25 A, -> 2.7 A).
2. Faces / access region ("dr" part): current runs tangentially along the
   membrane faces, so the face ramp acts like the face moved out by the same
   deficit length t (L -> L + 2t is exact for a sharp-cornered offset by
   translation invariance). The Euclidean distance, however, rounds the rim
   with radius t, which removes a small end length

       dl = c(lambda) * t**(4/3) * a_e**(-1/3)      per end,

   the t**(4/3) scaling being that of rounding a 270-degree corner, whose
   sharp solution is singular as rho**(2/3). c -> 0.53 for L >~ a.
3. Theory: Hall's access term rho/(4a) per side is the thin-membrane
   (equipotential disc) value. For a membrane of finite thickness the exact
   end correction per end is kappa(L/a) * a, rising from pi/4 = 0.7854 at
   L/a -> 0 to the flanged-tube value 0.82155 already at L/a ~ 1, so Hall
   underestimates the access resistance of practical pores by 4.6 %. This,
   not the soft wall, is what makes the "required" offset drift with d.

Putting the three together,

    R_SEM = rho * (L_e + 2 kappa(L_e/a_e) a_e) / (pi a_e^2),
    a_e = equal-area radius,  L_e = L + 2t - 2 dl,

reproduces the converged SEM ramp conductance to 0.03 % (max) over
d = 2-50 nm, L = 1-40 nm (0.16 % at d = 1.5 nm).

Units: lengths in Angstrom, conductivity in S/m, resistance in Ohm,
conductance in S.
"""
from __future__ import annotations

import numpy as np

# SEM ramp defaults (sem.utils.condfrac)
R_MIN = 1.3
R_CUT = 4.1

#: End correction per end, in units of the pore radius, of a long cylinder in
#: a thick insulating plane (flanged tube). Grid-converged with the
#: axisymmetric reference solver: 0.82155 +- 0.00002; cf. 0.8216 in the
#: acoustics literature (Norris & Sheng 1989).
KAPPA_INF = 0.82155
#: Thin-membrane limit (equipotential disc, Hall 1975): pi/4.
KAPPA_THIN = np.pi / 4.0
#: First zero of J1: decay rate (in L/a) of the interaction between the two
#: pore ends through the first non-trivial mode of an insulating cylinder.
J1_ZERO = 3.8317059702075125

# Fitted shape of kappa(lambda) between the two exact limits (max abs error
# 6e-5 against the grid-converged table, lambda = L/a in [0.02, 20]).
_KAPPA_W = 0.49634
_KAPPA_C = 4.19081

# Rim (rounded-corner) end-length deficit, fitted to the converged SEM ramp
# conductance over d = 2-50 nm, L = 1-40 nm (r_min = 1.3 A, r_cut = 4.1 A).
_RIM_C = 0.5289
_RIM_B = 0.1631
_RIM_S = 3.1853
# Dependence of the rim coefficient on the ramp shape, through the relative
# ramp width x = (r_cut - r_min) / t (x = 0 is a sharp step at D = t). Fitted
# to thick-membrane solves (a = 100 A, L = 400 A) for 0 <= r_min/r_cut <= 0.9;
# no linear term, since t is the first moment of the ramp deficit.
_RIM_SHAPE = (0.45515, 0.08743, -0.01614, 0.00346)

_A_PER_M = 1e10


def deficit_length(r_min: float = R_MIN, r_cut: float = R_CUT) -> float:
    """Conductance deficit of the ramp on a flat wall, int_0^inf (1 - f) dD.

    For the linear ramp this is the ramp midpoint, (r_min + r_cut) / 2.
    """
    return 0.5 * (r_min + r_cut)


def equal_area_radius(a, r_min: float = R_MIN, r_cut: float = R_CUT):
    """Radius of the sharp pore with the same sigma-weighted cross-section.

    pi a_e^2 = int_0^a f(a - r) 2 pi r dr, valid for a >= r_cut. This is the
    exact channel-equivalent radius (the bore field is purely axial).
    """
    a = np.asarray(a, float)
    if np.any(a < r_cut):
        raise ValueError("equal_area_radius requires a >= r_cut")
    dr = r_cut - r_min
    c0 = r_cut ** 2 - dr * r_min - (2.0 / 3.0) * dr ** 2
    return np.sqrt(a ** 2 - (r_cut + r_min) * a + c0)


def end_correction_coefficient(lam):
    """kappa(L/a): end correction per end, in units of a, of a sharp-walled
    cylinder of radius a in an insulating membrane of thickness L.

    R_sharp = rho (L + 2 kappa a) / (pi a^2). kappa(0) = pi/4 (Hall's thin
    aperture, exact), kappa(inf) = 0.82155 (flanged tube). Accurate to 6e-5.
    """
    lam = np.asarray(lam, float)
    if np.any(lam < 0):
        raise ValueError("L/a must be non-negative")
    shape = np.exp(-J1_ZERO * lam) * (_KAPPA_W + (1.0 - _KAPPA_W)
                                      * np.exp(-_KAPPA_C * lam ** (2.0 / 3.0)))
    return KAPPA_INF - (KAPPA_INF - KAPPA_THIN) * shape


def hall_resistance(d, L, sigma: float = 1.0):
    """Kowalczyk et al. (2011) Eq. 3: rho (4L/(pi d^2) + 1/d). d, L in A."""
    d = np.asarray(d, float)
    L = np.asarray(L, float)
    return (4.0 * L / (np.pi * d ** 2) + 1.0 / d) * _A_PER_M / sigma


def sharp_cylinder_resistance(d, L, sigma: float = 1.0):
    """Sharp-walled cylinder with the exact finite-thickness access term.

    Same as :func:`hall_resistance` with 1/d replaced by 4 kappa(2L/d)/(pi d);
    within 1e-4 of the exact Laplace solution for any L/d.
    """
    d = np.asarray(d, float)
    L = np.asarray(L, float)
    a = 0.5 * d
    kappa = end_correction_coefficient(L / a)
    return (L + 2.0 * kappa * a) / (np.pi * a ** 2) * _A_PER_M / sigma


def _rim_shape_factor(r_min: float, r_cut: float) -> float:
    t = deficit_length(r_min, r_cut)
    x = (r_cut - r_min) / t
    c0, c2, c3, c4 = _RIM_SHAPE
    return c0 + c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4


def rim_length_deficit(a_e, L_e0, r_min: float = R_MIN, r_cut: float = R_CUT):
    """End-length deficit per end from the rim rounding of the SEM map.

    The Euclidean-distance ramp makes the effective wall (roughly) the
    Minkowski offset of the membrane, whose rim is rounded with radius t;
    relative to a sharp-cornered offset this shortens each end by
    c(lambda) t^(4/3) a_e^(-1/3), lambda = L_e0 / a_e, L_e0 = L + 2t.
    The t^(4/3) a^(-1/3) scaling is that of rounding a 270-degree corner;
    c depends on the ramp shape (0.53 for the default ramp in a thick
    membrane, 0.46 for a sharp step) and grows for thin membranes.
    """
    a_e = np.asarray(a_e, float)
    lam = np.asarray(L_e0, float) / a_e
    t = deficit_length(r_min, r_cut)
    shape = _rim_shape_factor(r_min, r_cut) / _rim_shape_factor(R_MIN, R_CUT)
    c = _RIM_C * shape * (1.0 + _RIM_B * np.exp(-_RIM_S * lam))
    return c * t ** (4.0 / 3.0) * a_e ** (-1.0 / 3.0)


def effective_cylinder(d, L, r_min: float = R_MIN, r_cut: float = R_CUT):
    """Sharp-walled (d_e, L_e) that is electrically equivalent to an SEM
    cylinder of nominal diameter d and thickness L (all in A).

    d_e = 2 a_e (equal-area), L_e = L + 2t - 2 dl (face offset minus rim).
    Evaluate :func:`sharp_cylinder_resistance` at (d_e, L_e) -- not Hall.
    """
    a = 0.5 * np.asarray(d, float)
    t = deficit_length(r_min, r_cut)
    a_e = equal_area_radius(a, r_min, r_cut)
    L_e0 = np.asarray(L, float) + 2.0 * t
    L_e = L_e0 - 2.0 * rim_length_deficit(a_e, L_e0, r_min, r_cut)
    return 2.0 * a_e, L_e


def sem_cylinder_resistance(d, L, sigma: float = 1.0,
                            r_min: float = R_MIN, r_cut: float = R_CUT,
                            model: str = "full"):
    """Predicted open-pore resistance (Ohm) of an SEM cylindrical pore with
    nominal diameter d and membrane thickness L (A), infinite reservoirs.

    model="full":  sharp_cylinder_resistance(*effective_cylinder(d, L)) --
                   0.03 % of the converged ramp solution (d >= 2 nm).
    model="paper": hall_resistance(d - 2t, L + 2t) -- the area-conserving
                   offset with Hall's access term (2-3 % high in G at d ~ L).
    """
    if model == "full":
        d_e, L_e = effective_cylinder(d, L, r_min, r_cut)
        return sharp_cylinder_resistance(d_e, L_e, sigma)
    if model == "paper":
        t = deficit_length(r_min, r_cut)
        return hall_resistance(np.asarray(d, float) - 2 * t,
                               np.asarray(L, float) + 2 * t, sigma)
    raise ValueError("model must be 'full' or 'paper'")


def sem_cylinder_conductance(d, L, sigma: float = 1.0, **kwargs):
    """Conductance (S) counterpart of :func:`sem_cylinder_resistance`."""
    return 1.0 / sem_cylinder_resistance(d, L, sigma, **kwargs)


def sem_cylinder_for_sharp_pore(d_sharp, L_sharp, r_min: float = R_MIN,
                                r_cut: float = R_CUT, tol: float = 1e-10):
    """Inverse map: nominal SEM (d, L) whose ramped map conducts exactly like
    a sharp-walled cylinder (d_sharp, L_sharp), e.g. an experimental pore.

    Solved by fixed-point iteration on effective_cylinder (converges in a
    few steps since the offsets depend only weakly on d and L).
    """
    d_sharp = float(d_sharp)
    L_sharp = float(L_sharp)
    d, L = d_sharp + 2 * deficit_length(r_min, r_cut), L_sharp
    for _ in range(100):
        d_e, L_e = effective_cylinder(d, L, r_min, r_cut)
        d_new, L_new = d + (d_sharp - d_e), L + (L_sharp - L_e)
        if abs(d_new - d) < tol and abs(L_new - L) < tol:
            return float(d_new), float(L_new)
        d, L = d_new, L_new
    raise RuntimeError("inverse mapping did not converge")
