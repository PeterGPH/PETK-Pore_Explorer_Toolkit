"""Soft-wall -> sharp-wall mapping (sem.open_pore_theory), its exact
sharp-cylinder reference (sem.cylinder_mode_matching) and the axisymmetric
solver (sem.axisym_conductance). Grids are coarse to keep the suite fast; the
converged comparison lives in validation/soft_wall_mapping/."""
import math

import numpy as np
import pytest

from sem import open_pore_theory as opt
from sem.axisym_conductance import AxisymmetricPore, pore_resistance, ramp_fraction
from sem.cylinder_mode_matching import KAPPA_INF_EXACT, sharp_cylinder_coefficients


# ------------------------------------------------------------ exact sharp cylinder
def test_mode_matching_limits():
    kappa, k_rim = sharp_cylinder_coefficients([0.0, 20.0])
    assert kappa[0] == pytest.approx(math.pi / 4, abs=5e-5)     # thin aperture (Rayleigh/Hall)
    assert kappa[1] == pytest.approx(KAPPA_INF_EXACT, abs=1e-6)  # flanged tube
    assert kappa[1] == pytest.approx(0.8216, abs=1e-4)          # Norris & Sheng (1989)
    assert k_rim[1] == pytest.approx(opt.K_RIM_INF, rel=1e-4)


def test_closed_forms_match_mode_matching():
    lam = np.array([0.03, 0.1, 0.3, 0.56, 1.0, 1.7, 5.0])
    kappa, k_rim = sharp_cylinder_coefficients(lam)
    np.testing.assert_allclose(opt.end_correction_coefficient(lam), kappa, atol=6e-5)
    sel = lam >= 0.05
    np.testing.assert_allclose(opt.rim_strength(lam[sel]), k_rim[sel], rtol=3e-3)


def test_kappa_limits_and_monotone():
    assert opt.end_correction_coefficient(0.0) == pytest.approx(math.pi / 4, abs=1e-12)
    assert opt.end_correction_coefficient(50.0) == pytest.approx(opt.KAPPA_INF, abs=1e-12)
    lam = np.linspace(0, 5, 200)
    assert np.all(np.diff(opt.end_correction_coefficient(lam)) >= 0)


def test_sharp_cylinder_reduces_to_hall_for_thin_membrane():
    # Hall's 1/d access term is exact for an aperture in a thin plate...
    d = 100.0
    assert opt.sharp_cylinder_resistance(d, 0.0) == pytest.approx(opt.hall_resistance(d, 0.0), rel=1e-12)
    # ...and 4.6 % short for a thick one.
    acc_exact = opt.sharp_cylinder_resistance(d, 1000.0) - opt.hall_resistance(d, 1000.0) + 1e10 / d
    assert acc_exact / (1e10 / d) == pytest.approx(4 * opt.KAPPA_INF / math.pi, rel=1e-9)
    assert 4 * opt.KAPPA_INF / math.pi == pytest.approx(1.0462, abs=1e-4)


# ------------------------------------------------------------------ ramp mapping
@pytest.mark.parametrize("a", [5.0, 25.0, 125.0])
def test_equal_area_radius_matches_quadrature(a):
    edges = np.linspace(0.0, a, 200001)
    r = 0.5 * (edges[1:] + edges[:-1])          # midpoint rule
    area = np.sum(ramp_fraction(a - r) * 2 * np.pi * r) * (edges[1] - edges[0])
    assert math.pi * opt.equal_area_radius(a) ** 2 == pytest.approx(area, rel=1e-6)


def test_effective_cylinder_is_consistent():
    for d, L in [(50.0, 10.0), (100.0, 200.0), (500.0, 20.0)]:
        d_e, L_e = opt.effective_cylinder(d, L)
        assert opt.sharp_cylinder_resistance(d_e, L_e) == pytest.approx(
            opt.sem_cylinder_resistance(d, L), rel=1e-10)
        assert d - d_e == pytest.approx(2 * opt.deficit_length(), abs=0.05)
        assert 0 < L_e - L < 2 * opt.deficit_length()        # the rim gives back part of 2t


def test_inverse_map_roundtrip():
    d, L = opt.sem_cylinder_for_sharp_pore(100.0, 200.0)
    d_e, L_e = opt.effective_cylinder(d, L)
    assert (float(d_e), float(L_e)) == pytest.approx((100.0, 200.0), abs=1e-8)
    assert d > 100.0 and L < 200.0


def test_rim_coefficient_theory_values():
    # -2 pi beta k_inf^2 with beta = -0.1494 (local corner problem) for the default ramp
    assert opt.rim_shape_coefficient() == pytest.approx(2 * math.pi * 0.1494 * opt.K_RIM_INF ** 2, rel=0.01)
    assert opt.rim_coefficient(4.0) == pytest.approx(opt.rim_shape_coefficient(), rel=1e-4)
    assert opt.rim_coefficient(0.1) > opt.rim_coefficient(1.0)      # thin plates: stronger rim


# ----------------------------------------------------------------------- solver
def test_long_sharp_cylinder_end_correction():
    a, L = 10.0, 200.0
    R = pore_resistance(AxisymmetricPore.cylinder(a, L), "sharp", h_fine=0.1, band=4.0).resistance
    kappa = (R * math.pi * a ** 2 - L) / (2 * a)
    assert kappa == pytest.approx(opt.KAPPA_INF, abs=1e-3)   # corner: converges as h^(4/3)


def test_non_representable_thickness_regression():
    # L = 0.56 * 100 = 56.00000000000001 used to create a degenerate element
    a, lam = 100.0, 0.56
    L = lam * a
    R = pore_resistance(AxisymmetricPore.cylinder(a, L), "sharp", h_fine=0.2, band=8.0).resistance
    kappa, _ = sharp_cylinder_coefficients(lam)
    assert (R * math.pi * a ** 2 - L) / (2 * a) == pytest.approx(kappa[0], abs=5e-4)


@pytest.mark.parametrize("d_nm, L_nm", [(5, 20), (25, 20), (10, 2)])
def test_full_mapping_matches_ramp_solution(d_nm, L_nm):
    R_ramp = pore_resistance(AxisymmetricPore.cylinder(5 * d_nm, 10 * L_nm), "ramp",
                             h_fine=0.2).resistance * 1e10
    R_full = opt.sem_cylinder_resistance(10 * d_nm, 10 * L_nm, model="full")
    R_paper = opt.sem_cylinder_resistance(10 * d_nm, 10 * L_nm, model="paper")
    assert R_full / R_ramp == pytest.approx(1.0, abs=4e-4)
    if d_nm == 25:   # the Hall offset overestimates G by ~2 % at d ~ L
        assert R_ramp / R_paper - 1 == pytest.approx(0.020, abs=0.002)


def test_symmetric_and_full_domain_agree():
    pore = AxisymmetricPore.cylinder(25.0, 100.0)
    R_half = pore_resistance(pore, "ramp", h_fine=0.2).resistance
    R_full = pore_resistance(pore, "ramp", h_fine=0.2, use_symmetry=False).resistance
    assert R_full == pytest.approx(R_half, rel=1e-4)


def test_distance_modes_coincide_for_cylinder():
    pore = AxisymmetricPore.cylinder(30.0, 60.0)
    r, z = np.meshgrid(np.linspace(0, 60, 61), np.linspace(-60, 60, 121), indexing="ij")
    np.testing.assert_allclose(pore.distance(r, z, "euclidean"), pore.distance(r, z, "horizontal"),
                               atol=1e-12)


def test_finite_box_approaches_infinite_reservoirs():
    pore = AxisymmetricPore.cylinder(50.0, 200.0)
    R_inf = pore_resistance(pore, "ramp", h_fine=0.2).resistance
    errs = [abs(pore_resistance(pore, "ramp", h_fine=0.2, far_field="box", box_radius=s,
                                box_height=2 * s + 200).resistance / R_inf - 1)
            for s in (400.0, 4000.0)]
    assert errs[1] < errs[0] and errs[1] < 1e-3


# ------------------------------------------------------------- edge cases
def test_inverse_map_rejects_unrepresentable_pores():
    with pytest.raises(ValueError, match="thinner"):
        opt.sem_cylinder_for_sharp_pore(500.0, 2.0)        # graphene-thin sharp membrane
    with pytest.raises(ValueError, match="d_sharp"):
        opt.sem_cylinder_for_sharp_pore(2.0, 100.0)


def test_closed_pore_has_infinite_resistance():
    barrel = AxisymmetricPore([-10.0, 0.0, 10.0], [1.0, 10.0, 1.0])   # neck narrower than r_min
    for sym in (True, False):
        assert math.isinf(pore_resistance(barrel, "ramp", h_fine=0.2, use_symmetry=sym).resistance)


def test_step_wall_limits_and_bulk_conductivity_scaling():
    pore = AxisymmetricPore.cylinder(10.0, 20.0)
    R_sharp = pore_resistance(pore, "sharp", h_fine=0.2).resistance
    assert pore_resistance(pore, "step", step_distance=0.0, h_fine=0.2).resistance == pytest.approx(R_sharp, rel=1e-12)
    R1 = pore_resistance(pore, "sharp", h_fine=0.2).resistance
    R2 = pore_resistance(pore, sigma_fn=lambda r, z: 10.5 * (~pore.is_solid(r, z)), h_fine=0.2).resistance
    assert R2 * 10.5 == pytest.approx(R1, rel=1e-9)        # far-field Robin term scales with sigma


def test_grid_parameters_are_validated():
    pore = AxisymmetricPore.cylinder(10.0, 20.0)
    for bad in (dict(growth=0.9), dict(h_fine=0.0)):
        with pytest.raises(ValueError):
            pore_resistance(pore, "ramp", **bad)
    with pytest.raises(ValueError):
        ramp_fraction(1.0, r_min=4.1, r_cut=1.3)
