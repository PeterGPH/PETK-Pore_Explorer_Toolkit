"""
Behavioural tests for `sem.geometry_profiles`, the numpy-only module
holding the axisymmetric wall-distance geometry (`PoreProfile`) shared by
every parametric pore type. None of this module may import scipy or
dolfinx, and it must not import from `sem.utils`.
"""

import re
from pathlib import Path

import numpy as np
import pytest

import sem.geometry_profiles
from sem.geometry_profiles import (
    PoreProfile,
    euclidean_distance_to_profile,
    legacy_distance,
)


def test_geometry_profiles_module_source_has_no_scipy_or_dolfinx_import():
    """Static guard: sem.geometry_profiles must stay numpy-only."""
    source = Path(sem.geometry_profiles.__file__).read_text()
    assert not re.search(r"^\s*(import|from)\s+(scipy|dolfinx)\b", source, re.MULTILINE)
    assert not re.search(r"^\s*(import|from)\s+.*\butils\b", source, re.MULTILINE)


# ---------------------------------------------------------------------------
# Cylinder: euclidean must be bit-identical to legacy (vertical wall).
# ---------------------------------------------------------------------------
def test_cylinder_corner_zero_euclidean_matches_legacy_exactly():
    profile = PoreProfile.cylindrical(pore_radius=50.0, half_thickness=100.0)
    assert profile.is_vertical_cylinder()

    rng = np.random.default_rng(0)
    R = rng.uniform(0, 150, size=2000)
    z = rng.uniform(-150, 150, size=2000)

    d_euclid = profile.distance(R, z, metric="euclidean")
    d_legacy = profile.distance(R, z, metric="legacy")
    assert np.array_equal(d_euclid, d_legacy)


# ---------------------------------------------------------------------------
# Double cone: sloped-wall interior distance vs. legacy's radial-only gap.
# ---------------------------------------------------------------------------
def test_double_cone_interior_distance_matches_geometric_formula():
    inner, outer, h = 25.0, 50.0, 50.0
    profile = PoreProfile.double_cone(inner_radius=inner, outer_radius=outer, half_thickness=h)

    z = 25.0  # strictly between the apex (z=0) and rim (z=h): on the upper wall segment
    r_at_z = inner + (outer - inner) * (z / h)  # 37.5
    delta = 5.0
    R = r_at_z - delta  # inside the bore, near the sloped wall

    tan_theta = (outer - inner) / h  # 0.5, slope dr/dz of the wall
    cos_theta = 1.0 / np.sqrt(1.0 + tan_theta ** 2)

    d_euclid = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    d_legacy = profile.distance(np.array([R]), np.array([z]), metric="legacy")[0]

    assert d_legacy == pytest.approx(r_at_z - R, abs=1e-9)
    assert d_euclid == pytest.approx(delta * cos_theta, abs=1e-9)
    assert d_euclid < d_legacy


# ---------------------------------------------------------------------------
# Conical: asymmetric signed-z interpolation (not abs(z) like double_cone).
# ---------------------------------------------------------------------------
def test_conical_signed_z_local_radius():
    bottom, top, h = 30.0, 60.0, 50.0
    profile = PoreProfile.conical(bottom_radius=bottom, top_radius=top, half_thickness=h)

    assert profile.local_radius(np.array([-h])) == pytest.approx(bottom)
    assert profile.local_radius(np.array([h])) == pytest.approx(top)
    assert profile.local_radius(np.array([0.0])) == pytest.approx((bottom + top) / 2.0)

    # Legacy formula regression: t=(z+h)/(2h), r = bottom + (top-bottom)*t
    z = np.array([-40.0, -10.0, 0.0, 20.0, 45.0])
    t = np.clip((z + h) / (2.0 * h), 0.0, 1.0)
    expected_local_radius = bottom + (top - bottom) * t
    np.testing.assert_allclose(profile.local_radius(z), expected_local_radius)

    R = np.full_like(z, 5.0)  # well inside the bore everywhere
    expected_legacy = legacy_distance(R, np.abs(z), expected_local_radius, h)
    np.testing.assert_allclose(profile.distance(R, z, metric="legacy"), expected_legacy)


# ---------------------------------------------------------------------------
# Mouth cases (cylinder): above the face, inside vs. outside the bore, and
# inside the solid (including inside a chamfer).
# ---------------------------------------------------------------------------
def test_mouth_above_face_inside_bore():
    r0, h = 50.0, 100.0
    profile = PoreProfile.cylindrical(pore_radius=r0, half_thickness=h)
    R, z = 20.0, 110.0  # above the top face, radially inside the bore
    expected = np.sqrt((r0 - R) ** 2 + (z - h) ** 2)
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == pytest.approx(expected)


def test_mouth_above_face_outside_bore():
    r0, h = 50.0, 100.0
    profile = PoreProfile.cylindrical(pore_radius=r0, half_thickness=h)
    R, z = 80.0, 110.0  # above the top face, radially beyond the pore mouth
    expected = z - h
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == pytest.approx(expected)


def test_inside_solid_is_zero_distance():
    r0, h = 50.0, 100.0
    profile = PoreProfile.cylindrical(pore_radius=r0, half_thickness=h)
    R, z = 80.0, 50.0  # within the slab, beyond the pore radius: inside the solid
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == 0.0


def test_inside_chamfer_solid_is_zero_distance():
    profile = PoreProfile.cylindrical(
        pore_radius=50.0, half_thickness=100.0, corner_radius=20.0, chamfer_depth=20.0
    )
    # z=95 is inside the chamfer zone (80 < |z| <= 100); local_radius there is 65.
    assert profile.local_radius(np.array([95.0]))[0] == pytest.approx(65.0)
    R, z = 70.0, 95.0  # R >= local_radius(z): inside the chamfer solid
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == 0.0


# ---------------------------------------------------------------------------
# Regression (fix round 1): a non-monotonic `from_table` wall whose
# interior bulges wider than its mouths. The face-ray term must use the
# plain (unclamped) signed vertical offset to the face plane -- clamping
# it to max(|z|-h, 0) let R >= both mouth radii but R < the interior
# local_radius spuriously read back as distance 0 (as if inside the
# solid), even though `inside` correctly says False.
# ---------------------------------------------------------------------------
def test_non_monotonic_profile_face_ray_distance_is_not_clamped():
    # Both mouths (z=-50 and z=+50) have radius 30; the interior bulges out
    # to radius 80 at z=0. Both wall segments have |slope| = 1 (45 degrees).
    z_table = np.array([-50.0, 0.0, 50.0])
    r_table = np.array([30.0, 80.0, 30.0])
    profile = PoreProfile.from_table(z_table, r_table, half_thickness=50.0)

    # Open bore, R between the mouth radius (30) and the interior local_radius
    # (80): before the fix this spuriously returned 0.0 (face term clamped to
    # radial-gap-only). True nearest point is on one of the 45-degree wall
    # segments; verify against the exact point-to-segment distance.
    R, z = 70.0, 0.0
    assert profile.local_radius(np.array([z]))[0] == pytest.approx(80.0)
    tan_theta = 1.0  # |dr/dz| of both wall segments
    cos_theta = 1.0 / np.sqrt(1.0 + tan_theta ** 2)
    expected = (80.0 - R) * cos_theta
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == pytest.approx(expected)
    assert got == pytest.approx(np.sqrt(50.0))  # (80-70)*cos(45 deg) = 10/sqrt(2)

    # Inside the bulging solid interior (R >= local_radius(0) == 80): zero.
    R, z = 85.0, 0.0
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == 0.0

    # Directly above the top mouth rim (R == top mouth radius exactly): the
    # radial gap is 0, so the distance is the plain vertical offset to the
    # face plane, z - h.
    R, z = 30.0, 60.0
    got = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert got == pytest.approx(10.0)
    assert got == pytest.approx(z - 50.0)


# ---------------------------------------------------------------------------
# 45-degree chamfer: euclidean interior distance is legacy / sqrt(2).
# ---------------------------------------------------------------------------
def test_45_degree_chamfer_interior_distance():
    pore_radius, h, corner = 50.0, 100.0, 20.0
    # chamfer_depth defaults to corner_radius -> slope dr/dz = corner/corner = 1 -> 45 degrees.
    profile = PoreProfile.cylindrical(pore_radius=pore_radius, half_thickness=h, corner_radius=corner)

    z = 90.0  # within the chamfer zone (80, 100]
    local_r = profile.local_radius(np.array([z]))[0]
    assert local_r == pytest.approx(60.0)
    delta = 5.0
    R = local_r - delta

    d_legacy = profile.distance(np.array([R]), np.array([z]), metric="legacy")[0]
    d_euclid = profile.distance(np.array([R]), np.array([z]), metric="euclidean")[0]
    assert d_legacy == pytest.approx(delta, abs=1e-9)
    assert d_euclid == pytest.approx(d_legacy / np.sqrt(2.0), abs=1e-9)


# ---------------------------------------------------------------------------
# local_radius must equal the legacy _compute_local_pore_radius expression,
# for both chamfer_depth < half_thickness and chamfer_depth > half_thickness.
# ---------------------------------------------------------------------------
def _legacy_compute_local_pore_radius(Z, pore_radius, membrane_half_thickness, corner_radius, chamfer_depth):
    """Verbatim port of CylindricalPore._compute_local_pore_radius, used
    here only as an independent oracle (this module must not import
    sem.pore_geometry, which depends on scipy)."""
    Z = np.asarray(Z, dtype=float)
    base_radius = np.full_like(Z, pore_radius, dtype=float)
    if corner_radius is None or corner_radius <= 0:
        return base_radius
    cd = chamfer_depth if chamfer_depth is not None else corner_radius
    if cd is None or cd <= 0:
        return base_radius
    edge_radius = pore_radius + corner_radius
    z_edge_dist = np.maximum(membrane_half_thickness - np.abs(Z), 0.0)
    in_chamfer_zone = z_edge_dist < cd
    chamfer_progress = np.zeros_like(Z, dtype=float)
    chamfer_progress[in_chamfer_zone] = np.clip(z_edge_dist[in_chamfer_zone] / cd, 0.0, 1.0)
    return np.where(in_chamfer_zone, edge_radius + (pore_radius - edge_radius) * chamfer_progress, base_radius)


@pytest.mark.parametrize("chamfer_depth", [30.0, 80.0])  # cd < h and cd > h (h=50)
def test_local_radius_matches_legacy_compute_local_pore_radius(chamfer_depth):
    pore_radius, half_thickness, corner_radius = 50.0, 50.0, 10.0
    profile = PoreProfile.cylindrical(
        pore_radius=pore_radius,
        half_thickness=half_thickness,
        corner_radius=corner_radius,
        chamfer_depth=chamfer_depth,
    )
    Z = np.linspace(-half_thickness, half_thickness, 21)
    expected = _legacy_compute_local_pore_radius(
        Z, pore_radius, half_thickness, corner_radius, chamfer_depth
    )
    np.testing.assert_allclose(profile.local_radius(Z), expected, atol=1e-9)


# ---------------------------------------------------------------------------
# from_table round-trips local_radius against the equivalent classmethod.
# ---------------------------------------------------------------------------
def test_from_table_round_trips_local_radius():
    h = 50.0
    # Wider-than-the-slab table; from_table must clip/resample to +-h.
    z = np.array([-60.0, -50.0, 0.0, 50.0, 60.0])
    r = np.array([70.0, 50.0, 25.0, 50.0, 70.0])
    profile = PoreProfile.from_table(z, r, half_thickness=h)

    assert profile.vertices[0, 1] == pytest.approx(-h)
    assert profile.vertices[-1, 1] == pytest.approx(h)
    assert profile.vertices[0, 0] == pytest.approx(50.0)
    assert profile.vertices[-1, 0] == pytest.approx(50.0)

    reference = PoreProfile.double_cone(inner_radius=25.0, outer_radius=50.0, half_thickness=h)
    z_query = np.linspace(-h, h, 11)
    np.testing.assert_allclose(profile.local_radius(z_query), reference.local_radius(z_query))


def test_from_table_rejects_non_increasing_z():
    with pytest.raises(ValueError):
        PoreProfile.from_table(np.array([-50.0, -50.0, 50.0]), np.array([50.0, 40.0, 50.0]), 50.0)


def test_from_table_rejects_table_narrower_than_slab():
    with pytest.raises(ValueError):
        PoreProfile.from_table(np.array([-30.0, 0.0, 30.0]), np.array([50.0, 25.0, 50.0]), 50.0)


# ---------------------------------------------------------------------------
# Fix round 2: a table whose ends fall a rounding-level distance short of
# +-h (e.g. from writing/reloading through limited-precision text) must be
# accepted and clamped, not rejected -- but a real, larger shortfall must
# still raise.
# ---------------------------------------------------------------------------
def test_from_table_clamps_rounding_level_shortfall_to_half_thickness():
    h = 100.90523721980028
    shortfall = 5e-7
    z = np.array([-h + shortfall, 0.0, h - shortfall])
    r = np.array([50.0, 25.0, 50.0])

    profile = PoreProfile.from_table(z, r, half_thickness=h)

    assert profile.vertices[0, 1] == -h
    assert profile.vertices[-1, 1] == h
    assert profile.half_thickness == h


def test_from_table_rejects_shortfall_larger_than_rounding_slack():
    h = 100.90523721980028
    shortfall = 1e-3
    z = np.array([-h + shortfall, 0.0, h - shortfall])
    r = np.array([50.0, 25.0, 50.0])

    with pytest.raises(ValueError):
        PoreProfile.from_table(z, r, half_thickness=h)


# ---------------------------------------------------------------------------
# distance_xyz is rotation-invariant about the z-axis (axisymmetry).
# ---------------------------------------------------------------------------
def test_distance_xyz_rotation_invariance():
    profile = PoreProfile.double_cone(inner_radius=25.0, outer_radius=50.0, half_thickness=50.0)
    rng = np.random.default_rng(1)
    n = 200
    R = rng.uniform(0, 80, size=n)
    theta = rng.uniform(0, 2 * np.pi, size=n)
    z = rng.uniform(-70, 70, size=n)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    points = np.column_stack([x, y, z])

    d_before = profile.distance_xyz(points, metric="euclidean")

    rotate_by = 1.234  # radians
    c, s = np.cos(rotate_by), np.sin(rotate_by)
    rotated = points.copy()
    rotated[:, 0] = c * x - s * y
    rotated[:, 1] = s * x + c * y

    d_after = profile.distance_xyz(rotated, metric="euclidean")
    np.testing.assert_allclose(d_before, d_after, atol=1e-9)


# ---------------------------------------------------------------------------
# Chunking must not change the result.
# ---------------------------------------------------------------------------
def test_chunked_distance_matches_unchunked():
    profile = PoreProfile.double_cone(inner_radius=25.0, outer_radius=50.0, half_thickness=50.0)
    rng = np.random.default_rng(2)
    R = rng.uniform(0, 80, size=53)
    z = rng.uniform(-70, 70, size=53)

    d_unchunked = profile.distance(R, z, metric="euclidean")
    d_chunked = profile.distance(R, z, metric="euclidean", chunk=7)
    np.testing.assert_array_equal(d_unchunked, d_chunked)


def test_euclidean_distance_to_profile_matches_pore_profile_distance():
    """Sanity check that the module-level function and the PoreProfile
    method agree (the method delegates to it for non-vertical walls)."""
    profile = PoreProfile.conical(bottom_radius=30.0, top_radius=60.0, half_thickness=50.0)
    R = np.array([10.0, 40.0, 70.0])
    z = np.array([-40.0, 0.0, 45.0])
    direct = euclidean_distance_to_profile(R, z, profile.vertices, profile.half_thickness)
    via_method = profile.distance(R, z, metric="euclidean")
    np.testing.assert_array_equal(direct, via_method)
