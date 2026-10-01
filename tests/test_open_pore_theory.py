"""Soft-wall -> sharp-wall mapping (sem.open_pore_theory) and its reference
solver (sem.axisym_conductance). Grids are coarse to keep the suite fast; the
converged comparison lives in validation/soft_wall_mapping/."""
import math

import numpy as np
import pytest

from sem import open_pore_theory as opt
from sem.axisym_conductance import AxisymmetricPore, pore_resistance, ramp_fraction


def test_kappa_limits():
    assert opt.end_correction_coefficient(0.0) == pytest.approx(math.pi / 4, abs=1e-12)
    assert opt.end_correction_coefficient(50.0) == pytest.approx(opt.KAPPA_INF, abs=1e-12)
    lam = np.linspace(0, 5, 200)
    assert np.all(np.diff(opt.end_correction_coefficient(lam)) >= 0)


def test_sharp_cylinder_reduces_to_hall_for_thin_membrane():
    # Hall's 1/d access term is exact for an aperture in a thin plate.
    d = 100.0
    assert opt.sharp_cylinder_resistance(d, 0.0) == pytest.approx(opt.hall_resistance(d, 0.0), rel=1e-12)
    # ...and 4.6 % short for a thick one.
    acc_exact = opt.sharp_cylinder_resistance(d, 1000.0) - opt.hall_resistance(d, 1000.0) + 1e10 / d
    assert acc_exact / (1e10 / d) == pytest.approx(opt.KAPPA_INF / (math.pi / 4), rel=1e-9)


@pytest.mark.parametrize("a", [5.0, 25.0, 125.0])
def test_equal_area_radius_matches_quadrature(a):
    edges = np.linspace(0.0, a, 200001)
    r = 0.5 * (edges[1:] + edges[:-1])          # midpoint rule
    area = np.sum(ramp_fraction(a - r) * 2 * np.pi * r) * (edges[1] - edges[0])
    assert math.pi * opt.equal_area_radius(a) ** 2 == pytest.approx(area, rel=1e-6)


def test_inverse_map_roundtrip():
    d, L = opt.sem_cylinder_for_sharp_pore(100.0, 200.0)
    d_e, L_e = opt.effective_cylinder(d, L)
    assert (float(d_e), float(L_e)) == pytest.approx((100.0, 200.0), abs=1e-8)
    assert d > 100.0 and L < 200.0


def test_long_sharp_cylinder_end_correction():
    a, L = 10.0, 200.0
    R = pore_resistance(AxisymmetricPore.cylinder(a, L), "sharp", h_fine=0.1, band=4.0).resistance
    kappa = (R * math.pi * a ** 2 - L) / (2 * a)
    assert kappa == pytest.approx(opt.KAPPA_INF, abs=2e-3)   # converges as h^(4/3)


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
