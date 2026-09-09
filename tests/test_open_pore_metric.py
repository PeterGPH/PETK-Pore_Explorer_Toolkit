"""
Behavioural tests for the `distance_metric` / `chamfer_depth` / `profile`
wiring at the full `VerticalMovementSEM` (FEM) level. Requires dolfinx --
skipped cleanly wherever dolfinx is unavailable (e.g. the pip-only local
dev environment).

Uses a small box (+-30 A, 2 A grid) so the open-pore solve is fast:
  - a plain cylinder (radius 8 A, membrane 20 A): legacy and euclidean wall
    distance metrics must give identical open-pore current (vertical wall,
    the two metrics are mathematically equal).
  - a double cone (inner 4 A, outer 12 A, membrane 20 A): the wall is
    sloped, so legacy and euclidean must differ.
  - a `profile` pore built from the double cone's own (z, r) vertices must
    reproduce the double_cone euclidean current (same conductivity field).
"""

import pytest

pytest.importorskip("dolfinx")

from sem.vertical_movement_sem import VerticalMovementSEM  # noqa: E402

BOX = {"x": (-30.0, 30.0), "y": (-30.0, 30.0), "z": (-30.0, 30.0)}
MEMBRANE_THICKNESS = 20.0
GRID_RESOLUTION = 2.0

COMMON_KWARGS = dict(
    moving_pdb="unused.pdb",
    membrane_thickness=MEMBRANE_THICKNESS,
    box_dimensions=BOX,
    grid_resolution=GRID_RESOLUTION,
    bulk_conductivity=11.2,
    membrane_conductivity=1e-7 * 11.2,
    use_vdw_radii=False,
    use_pdb2pqr=False,
    prepare_analyte=False,
    cleanup_temp_files=True,
)


def _build(pore_type, distance_metric, **extra):
    return VerticalMovementSEM(
        pore_type=pore_type,
        distance_metric=distance_metric,
        **COMMON_KWARGS,
        **extra,
    )


def test_cylinder_legacy_vs_euclidean_identical():
    """A plain (unchamfered) cylinder has a purely vertical wall: the
    euclidean and legacy wall-distance metrics are mathematically
    equivalent, so the open-pore current must match to solver tolerance."""
    legacy = _build("cylindrical", "legacy", pore_radius=8.0, corner_radius=0.0)
    euclid = _build("cylindrical", "euclidean", pore_radius=8.0, corner_radius=0.0)

    i_legacy = legacy.calculate_open_pore_current()
    i_euclid = euclid.calculate_open_pore_current()

    assert i_legacy == pytest.approx(i_euclid, rel=1e-9)


def test_double_cone_legacy_vs_euclidean_differ():
    """The double-cone wall is sloped, so legacy (which overestimates the
    true wall distance by 1/cos(theta)) must give a measurably different
    open-pore current from the true euclidean distance."""
    legacy = _build("double_cone", "legacy", pore_radius=4.0, outer_radius=12.0)
    euclid = _build("double_cone", "euclidean", pore_radius=4.0, outer_radius=12.0)

    i_legacy = legacy.calculate_open_pore_current()
    i_euclid = euclid.calculate_open_pore_current()

    assert i_legacy != pytest.approx(i_euclid, rel=1e-6)


def test_profile_pore_matches_double_cone_euclidean(tmp_path):
    """A `profile` table built from the double cone's own vertices
    (outer, -h) / (inner, 0) / (outer, +h) must reproduce the double_cone
    euclidean open-pore current, since it is the identical conductivity
    field evaluated through a different pore_type code path."""
    inner_radius, outer_radius, half_thickness = 4.0, 12.0, MEMBRANE_THICKNESS / 2.0

    csv_path = tmp_path / "double_cone_profile.csv"
    csv_path.write_text(
        "z,r\n"
        f"{-half_thickness},{outer_radius}\n"
        f"0.0,{inner_radius}\n"
        f"{half_thickness},{outer_radius}\n"
    )

    double_cone = _build(
        "double_cone", "euclidean", pore_radius=inner_radius, outer_radius=outer_radius
    )
    profile = _build("profile", "euclidean", profile_path=str(csv_path))

    i_double_cone = double_cone.calculate_open_pore_current()
    i_profile = profile.calculate_open_pore_current()

    assert i_profile == pytest.approx(i_double_cone, rel=1e-9)
