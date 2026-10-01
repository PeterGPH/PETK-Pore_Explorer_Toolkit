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
almost right but not (ii): it overestimates G by up to 3 % (2.0 % at
d = 25 nm, L = 20 nm).

The mapping (checked against :mod:`sem.axisym_conductance` and the exact
sharp-cylinder solution of :mod:`sem.cylinder_mode_matching`)
------------------------------------------------------------------------
1. Bore: away from the ends the field is axial, so the conducting layers add
   in parallel and the equal-area radius ``a_e`` is exact for the channel:
   a_e^2 = (a - t)^2 + (r_cut - r_min)^2 / 12 (ramp mean and variance).
2. Faces / access region (the radial "dr" part): current runs tangentially
   along the faces, where a ramp acts like the face moved out by the same
   deficit length t, with no correction at second order. For a sharp-cornered
   offset this is exactly L -> L + 2t.
3. Rim: the Euclidean distance rounds the rim, so the equivalent wall is the
   rounded (Minkowski) offset. Relative to the sharp-cornered offset this
   removes an end length per end

       dl = C (k(lambda)/k_inf)^2 t^(4/3) a_e^(-1/3),   lambda = (L + 2t)/a_e,

   from the r^(2/3) singularity of the 270-degree fluid corner: C is a local
   constant of the ramp shape (0.535 for the default ramp, 0.454 for a sharp
   step at D = t), k(lambda) the rim singularity strength of the sharp
   cylinder (k_inf = 0.7549; grows as lambda^(-1/6) for thin membranes).
4. Theory: Hall's access term rho/(4a) per side is the thin-plate limit. In a
   membrane of finite thickness the end correction per end is kappa(L/a) a,
   rising from pi/4 = 0.7854 at L/a -> 0 to the flanged-tube value 0.82167
   already at L/a ~ 1, so Hall underestimates the access resistance of
   practical pores by 4.6 %. This, not the soft wall, is what makes the
   "required" offset drift with d.

Putting them together,

    R_SEM = rho * (L + 2t + 2 kappa(lambda) a_e - 2 dl) / (pi a_e^2),

with no parameter fitted to SEM solutions; it reproduces the converged SEM
ramp conductance over d = 2-50 nm, L = 1-40 nm (see
validation/soft_wall_mapping/README.md for the error budget).

Units: lengths in Angstrom, conductivity in S/m, resistance in Ohm,
conductance in S.
"""
from __future__ import annotations

import numpy as np

# SEM ramp defaults (sem.utils.condfrac)
R_MIN = 1.3
R_CUT = 4.1

#: End correction per end, in units of the pore radius, of a long cylinder in
#: a thick insulating plane (flanged tube), from sem.cylinder_mode_matching
#: (variational, converged to 1e-9); 0.8216 in Norris & Sheng (1989).
KAPPA_INF = 0.8216742
#: Thin-membrane limit (equipotential aperture, Hall 1975): pi/4.
KAPPA_THIN = np.pi / 4.0
#: First zero of J1: decay rate (in L/a) of the interaction between the two
#: pore ends through the first non-trivial mode of an insulating cylinder.
J1_ZERO = 3.8317059702075125
#: Rim singularity strength k = K / (E a^(1/3)) of a long sharp cylinder,
#: phi ~ K rho^(2/3) cos(2 theta/3), E = I/(sigma pi a^2).
K_RIM_INF = 0.75485

# kappa(lambda) between the exact limits: kappa_inf - (kappa_inf - pi/4)
# e^(-j lambda) [1 + e^(-c lambda^(2/3))] / 2; one constant, fitted to the
# mode-matching table; max abs error 5e-5 for lambda in [0.02, 20].
_KAPPA_C = 4.313

# k(lambda)/k_inf from sem.cylinder_mode_matching (14 basis functions;
# converged to <= 0.2 %). Below the table k follows the thin-plate law
# k ~ lambda^(-1/6); above it k = k_inf.
_K_RIM_LAM = (0.05, 0.06297, 0.0793, 0.09988, 0.12578, 0.15841, 0.1995, 0.25125,
              0.31643, 0.3985, 0.50188, 0.63206, 0.79602, 1.0025, 1.26254, 1.59005,
              2.0025, 2.52194, 3.17612, 4.0)
_K_RIM_RATIO = (1.27370, 1.23366, 1.19618, 1.16225, 1.13199, 1.10508, 1.08131, 1.06071,
                1.04336, 1.02931, 1.01848, 1.01065, 1.00546, 1.00241, 1.00088, 1.00025,
                1.00005, 1.00001, 1.00000, 1.00000)

# Thick-membrane rim coefficient C as a function of the relative ramp width
# x = (r_cut - r_min)/t (x = 0: sharp step at D = t), from axisymmetric
# solves at a = 100 A, L = 400 A (sem.scripts.soft_wall_mapping rimshape);
# C(x = 1.037) = 0.5349 for the default ramp, cf. -2 pi beta k_inf^2 = 0.5349
# with the local corner constant beta = -0.1494.
# No linear term: t is the first moment of the ramp deficit.
_RIM_SHAPE = (0.45511, 0.08670, -0.01542, 0.00326)

_A_PER_M = 1e10


def deficit_length(r_min: float = R_MIN, r_cut: float = R_CUT) -> float:
    """Conductance deficit of the ramp on a flat wall, int_0^inf (1 - f) dD.

    For the linear ramp this is the ramp midpoint, (r_min + r_cut) / 2.
    """
    if r_cut < r_min or r_min < 0:
        raise ValueError("ramp needs 0 <= r_min <= r_cut")
    return 0.5 * (r_min + r_cut)


def equal_area_radius(a, r_min: float = R_MIN, r_cut: float = R_CUT):
    """Radius of the sharp pore with the same sigma-weighted cross-section.

    pi a_e^2 = int_0^a f(a - r) 2 pi r dr = pi [(a - t)^2 + (r_cut - r_min)^2/12],
    valid for a >= r_cut. This is the exact channel-equivalent radius.
    """
    a = np.asarray(a, float)
    if np.any(a < r_cut):
        raise ValueError("equal_area_radius requires a >= r_cut")
    t = deficit_length(r_min, r_cut)
    return np.sqrt((a - t) ** 2 + (r_cut - r_min) ** 2 / 12.0)


def end_correction_coefficient(lam):
    """kappa(L/a): end correction per end, in units of a, of a sharp-walled
    cylinder of radius a in an insulating membrane of thickness L.

    R_sharp = rho (L + 2 kappa a) / (pi a^2). kappa(0) = pi/4 (Hall's thin
    aperture, exact), kappa(inf) = 0.82167 (flanged tube). Accurate to 5e-5.
    """
    lam = np.asarray(lam, float)
    if np.any(lam < -1e-9):
        raise ValueError("L/a must be non-negative")
    lam = np.maximum(lam, 0.0)
    shape = 0.5 * np.exp(-J1_ZERO * lam) * (1.0 + np.exp(-_KAPPA_C * lam ** (2.0 / 3.0)))
    return KAPPA_INF - (KAPPA_INF - KAPPA_THIN) * shape


def rim_strength(lam):
    """Rim singularity strength k(L/a) of a sharp cylinder (k_inf = 0.7549)."""
    lam = np.asarray(lam, float)
    lo, hi = _K_RIM_LAM[0], _K_RIM_LAM[-1]
    ratio = np.interp(np.log(np.clip(lam, lo, hi)), np.log(_K_RIM_LAM), _K_RIM_RATIO)
    thin = _K_RIM_RATIO[0] * (np.maximum(lam, 1e-12) / lo) ** (-1.0 / 6.0)
    return K_RIM_INF * np.where(lam < lo, thin, ratio)


def rim_shape_coefficient(r_min: float = R_MIN, r_cut: float = R_CUT) -> float:
    """Thick-membrane rim coefficient C of a linear ramp (0.535 by default)."""
    x = (r_cut - r_min) / deficit_length(r_min, r_cut)
    c0, c2, c3, c4 = _RIM_SHAPE
    return c0 + c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4


def rim_coefficient(lam, r_min: float = R_MIN, r_cut: float = R_CUT):
    """c(lambda) = C (k(lambda)/k_inf)^2 in dl = c t^(4/3) a_e^(-1/3)."""
    return rim_shape_coefficient(r_min, r_cut) * (rim_strength(lam) / K_RIM_INF) ** 2


def rim_length_deficit(a_e, lam, r_min: float = R_MIN, r_cut: float = R_CUT):
    """End-length deficit per end (A) from the rounded rim of the SEM map,
    relative to the sharp-cornered offset cylinder (a_e, L + 2t);
    lam = (L + 2t)/a_e."""
    t = deficit_length(r_min, r_cut)
    return rim_coefficient(lam, r_min, r_cut) * t ** (4.0 / 3.0) * np.asarray(a_e, float) ** (-1.0 / 3.0)


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


def sem_cylinder_resistance(d, L, sigma: float = 1.0,
                            r_min: float = R_MIN, r_cut: float = R_CUT,
                            model: str = "full"):
    """Predicted open-pore resistance (Ohm) of an SEM cylindrical pore with
    nominal diameter d and membrane thickness L (A), infinite reservoirs.

    model="full":  rho (L + 2t + 2 kappa(lambda) a_e - 2 dl) / (pi a_e^2).
    model="paper": hall_resistance(d - 2t, L + 2t) -- the area-conserving
                   offset with Hall's access term (up to 3 % high in G).
    """
    d = np.asarray(d, float)
    L = np.asarray(L, float)
    t = deficit_length(r_min, r_cut)
    if model == "paper":
        return hall_resistance(d - 2 * t, L + 2 * t, sigma)
    if model != "full":
        raise ValueError("model must be 'full' or 'paper'")
    a_e = equal_area_radius(0.5 * d, r_min, r_cut)
    lam = (L + 2 * t) / a_e
    length = (L + 2 * t + 2 * end_correction_coefficient(lam) * a_e
              - 2 * rim_length_deficit(a_e, lam, r_min, r_cut))
    return length / (np.pi * a_e ** 2) * _A_PER_M / sigma


def sem_cylinder_conductance(d, L, sigma: float = 1.0, **kwargs):
    """Conductance (S) counterpart of :func:`sem_cylinder_resistance`."""
    return 1.0 / sem_cylinder_resistance(d, L, sigma, **kwargs)


def effective_cylinder(d, L, r_min: float = R_MIN, r_cut: float = R_CUT):
    """Sharp-walled (d_e, L_e), in A, that conducts exactly like the SEM map
    of a cylinder with nominal diameter d and thickness L.

    d_e = 2 a_e (equal-area); L_e solves sharp_cylinder_resistance(d_e, L_e)
    = sem_cylinder_resistance(d, L), i.e. L_e ~ L + 2t - 2 dl.
    """
    d = np.asarray(d, float)
    L = np.asarray(L, float)
    a_e = equal_area_radius(0.5 * d, r_min, r_cut)
    target = sem_cylinder_resistance(d, L, 1.0, r_min, r_cut) / _A_PER_M * np.pi * a_e ** 2
    L_e = np.asarray(L + 2 * deficit_length(r_min, r_cut), float)
    for _ in range(50):   # Newton on L_e + 2 kappa(L_e/a_e) a_e = target
        lam = L_e / a_e
        f = L_e + 2 * end_correction_coefficient(lam) * a_e - target
        h = 1e-6 * np.maximum(lam, 1e-3)
        dk = (end_correction_coefficient(lam + h) - end_correction_coefficient(np.maximum(lam - h, 0)))
        df = 1 + 2 * dk / (lam + h - np.maximum(lam - h, 0))
        step = f / df
        L_e = np.maximum(L_e - step, 1e-9)
        if np.all(np.abs(step) < 1e-10):
            break
    return 2 * a_e, L_e


def sem_cylinder_for_sharp_pore(d_sharp, L_sharp, r_min: float = R_MIN,
                                r_cut: float = R_CUT):
    """Inverse map: nominal SEM (d, L) whose ramped map conducts like a
    sharp-walled cylinder (d_sharp, L_sharp), e.g. an experimental pore, with
    matching channel area (d_e = d_sharp) and resistance (L_e = L_sharp).

    Raises ValueError when no SEM pore can do so: d_sharp below the smallest
    equal-area diameter the ramp allows, or a sharp membrane thinner than
    the ramped faces themselves (L_sharp < L_e(d, 0) ~ 2t - 2 dl, ~4 A).
    """
    from scipy.optimize import brentq

    d_sharp = float(d_sharp)
    L_sharp = float(L_sharp)
    t = deficit_length(r_min, r_cut)
    d_min = 2 * float(equal_area_radius(r_cut, r_min, r_cut))
    if d_sharp <= d_min:
        raise ValueError(f"d_sharp must exceed {d_min:.2f} A for this ramp")
    # diameter: invert a_e^2 = (a - t)^2 + (r_cut - r_min)^2 / 12
    d = 2 * (t + np.sqrt((0.5 * d_sharp) ** 2 - (r_cut - r_min) ** 2 / 12.0))
    # thickness: L_e(d, L) increases monotonically with L
    L_min = float(effective_cylinder(d, 0.0, r_min, r_cut)[1])
    if L_sharp < L_min:
        raise ValueError(f"a sharp membrane thinner than {L_min:.2f} A cannot be "
                         f"represented by the ramped map at d = {d:.1f} A")
    L = brentq(lambda x: float(effective_cylinder(d, x, r_min, r_cut)[1]) - L_sharp,
               0.0, L_sharp + 2 * t, xtol=1e-12)
    return float(d), float(L)
