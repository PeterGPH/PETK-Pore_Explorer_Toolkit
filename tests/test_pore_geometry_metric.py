"""
Behavioural tests for the `distance_metric` wire-up in `sem.pore_geometry`:
`CylindricalPore`/`DoubleConePore`/`ConicalPore` built on `PoreProfile`, and
the new `ProfilePore`. Imports only `sem.pore_geometry` (numpy + scipy, no
dolfinx), matching Task 0's guarantee that this module is importable
without dolfinx.
"""

import numpy as np
import pytest

from sem.pore_geometry import CylindricalPore, DoubleConePore, ProfilePore


def _grid():
    x = np.linspace(-90.0, 90.0, 13)
    y = np.linspace(-90.0, 90.0, 13)
    z = np.linspace(-60.0, 60.0, 13)
    return np.meshgrid(x, y, z, indexing="ij")


def _grid_points(X, Y, Z):
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])


def _query_conductivities(interp, X, Y, Z):
    return interp(_grid_points(X, Y, Z))


def test_cylindrical_pore_corner_zero_legacy_vs_euclidean_identical():
    X, Y, Z = _grid()
    common = dict(
        X=X, Y=Y, Z=Z, pore_radius=50.0, membrane_half_thickness=50.0,
        corner_radius=0.0, bulk_conductivity=10.5, membrane_conductivity=0.0001,
    )
    pore_legacy = CylindricalPore(distance_metric="legacy", **common)
    pore_euclid = CylindricalPore(distance_metric="euclidean", **common)

    vals_legacy = _query_conductivities(pore_legacy.get_conductivity_interpolator(), X, Y, Z)
    vals_euclid = _query_conductivities(pore_euclid.get_conductivity_interpolator(), X, Y, Z)

    np.testing.assert_array_equal(vals_legacy, vals_euclid)
    assert pore_euclid.get_distance_interpolator() is None


def test_cylindrical_pore_corner_positive_legacy_vs_euclidean_differ():
    X, Y, Z = _grid()
    common = dict(
        X=X, Y=Y, Z=Z, pore_radius=50.0, membrane_half_thickness=50.0,
        corner_radius=20.0, chamfer_depth=20.0,
        bulk_conductivity=10.5, membrane_conductivity=0.0001,
    )
    pore_legacy = CylindricalPore(distance_metric="legacy", **common)
    pore_euclid = CylindricalPore(distance_metric="euclidean", **common)

    vals_legacy = _query_conductivities(pore_legacy.get_conductivity_interpolator(), X, Y, Z)
    vals_euclid = _query_conductivities(pore_euclid.get_conductivity_interpolator(), X, Y, Z)

    assert not np.array_equal(vals_legacy, vals_euclid)


def test_double_cone_pore_legacy_matches_old_formula_inlined():
    """Regression: distance_metric='legacy' must exactly reproduce the
    pre-Euclidean DoubleConePore.get_conductivity_interpolator formula."""
    X, Y, Z = _grid()
    inner_radius, outer_radius, half_thickness = 25.0, 50.0, 50.0
    bulk_conductivity, membrane_conductivity = 10.5, 0.0001

    pore = DoubleConePore(
        X, Y, Z, inner_radius=inner_radius, outer_radius=outer_radius,
        membrane_half_thickness=half_thickness,
        bulk_conductivity=bulk_conductivity, membrane_conductivity=membrane_conductivity,
        distance_metric="legacy",
    )
    vals = _query_conductivities(pore.get_conductivity_interpolator(), X, Y, Z)

    # Old formula, inlined verbatim (pre-Task-1 DoubleConePore body).
    R = np.sqrt(X ** 2 + Y ** 2)
    abs_z = np.abs(Z)
    z_fraction = np.clip(abs_z / half_thickness, 0.0, 1.0)
    local_pore_radius = inner_radius + (outer_radius - inner_radius) * z_fraction
    radial_term = np.maximum(local_pore_radius - R, 0.0)
    vertical_term = np.maximum(abs_z - half_thickness, 0.0)
    distance_map = np.sqrt(radial_term ** 2 + vertical_term ** 2)
    fraction = np.clip((distance_map - 1.3) / (4.1 - 1.3), 1e-7, 1.0)
    expected_grid = membrane_conductivity + fraction * (bulk_conductivity - membrane_conductivity)
    expected = expected_grid.ravel()

    np.testing.assert_allclose(vals, expected, atol=1e-10)


def test_profile_pore_from_3_vertex_table_matches_double_cone_euclidean():
    X, Y, Z = _grid()
    inner_radius, outer_radius, half_thickness = 25.0, 50.0, 50.0

    dc_pore = DoubleConePore(
        X, Y, Z, inner_radius=inner_radius, outer_radius=outer_radius,
        membrane_half_thickness=half_thickness,
        distance_metric="euclidean",
    )
    z_table = np.array([-half_thickness, 0.0, half_thickness])
    r_table = np.array([outer_radius, inner_radius, outer_radius])
    profile_pore = ProfilePore(
        X, Y, Z, profile_table=(z_table, r_table),
        membrane_half_thickness=half_thickness,
        distance_metric="euclidean",
    )

    vals_dc = _query_conductivities(dc_pore.get_conductivity_interpolator(), X, Y, Z)
    vals_profile = _query_conductivities(profile_pore.get_conductivity_interpolator(), X, Y, Z)
    np.testing.assert_array_equal(vals_dc, vals_profile)
    assert profile_pore.get_distance_interpolator() is None
