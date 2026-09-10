"""
Behavioural tests for `sem.arbd_export._membrane_mask`'s analytic (no
`base_dist_interp`) fallback. Numpy-only: `sem.arbd_export` defers its
dolfinx import to the one function that actually needs it
(`_sample_dolfinx_function_on_grid`), so this module -- and this test file
-- import cleanly without dolfinx.

Regression coverage for the Task-1 review ledger item: the fallback used to
have separate per-pore-type branches (cylindrical/double_cone/conical) that
ignored `chamfer_depth` entirely and had no `profile` branch at all (always
returning "no wall" for `profile`). It now builds a `PoreProfile` via
`PoreProfile.from_params`, the same wall representation every parametric
pore type shares, so `profile` and `chamfer_depth` get identical semantics.
"""

import numpy as np

from sem.arbd_export import _membrane_mask


class _StubSEM:
    """Minimal duck-typed stand-in for the attributes `_membrane_mask`
    reads off a real `VerticalMovementSEM` instance."""

    base_dist_interp = None
    corner_radius = 0.0
    chamfer_depth = None
    outer_radius = None
    top_radius = None
    bottom_radius = None
    profile_path = None

    def __init__(self, pore_type, membrane_thickness, **overrides):
        self.pore_type = pore_type
        self.membrane_thickness = membrane_thickness
        for key, value in overrides.items():
            setattr(self, key, value)


def _grid_points():
    x = np.linspace(-20.0, 20.0, 9)
    y = np.linspace(-20.0, 20.0, 9)
    z = np.linspace(-15.0, 15.0, 11)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])


def test_cylinder_and_equivalent_profile_table_agree(tmp_path):
    pore_radius, membrane_thickness = 8.0, 20.0
    points = _grid_points()

    cyl_mask = _membrane_mask(
        _StubSEM("cylindrical", membrane_thickness, pore_radius=pore_radius),
        points,
    )

    csv_path = tmp_path / "cylinder_profile.csv"
    half_t = membrane_thickness / 2.0
    csv_path.write_text(
        "z,r\n"
        f"{-half_t},{pore_radius}\n"
        f"{half_t},{pore_radius}\n"
    )
    profile_mask = _membrane_mask(
        _StubSEM("profile", membrane_thickness, profile_path=str(csv_path)),
        points,
    )

    assert cyl_mask.any()
    np.testing.assert_array_equal(cyl_mask, profile_mask)


def test_double_cone_and_equivalent_profile_table_agree(tmp_path):
    inner_radius, outer_radius, membrane_thickness = 4.0, 12.0, 20.0
    points = _grid_points()

    dc_mask = _membrane_mask(
        _StubSEM(
            "double_cone", membrane_thickness,
            pore_radius=inner_radius, outer_radius=outer_radius,
        ),
        points,
    )

    csv_path = tmp_path / "double_cone_profile.csv"
    half_t = membrane_thickness / 2.0
    csv_path.write_text(
        "z,r\n"
        f"{-half_t},{outer_radius}\n"
        f"0.0,{inner_radius}\n"
        f"{half_t},{outer_radius}\n"
    )
    profile_mask = _membrane_mask(
        _StubSEM("profile", membrane_thickness, profile_path=str(csv_path)),
        points,
    )

    assert dc_mask.any()
    np.testing.assert_array_equal(dc_mask, profile_mask)


def test_cylindrical_mask_respects_chamfer_depth():
    """Regression: the old per-type branch ignored chamfer_depth entirely
    (always used the plain `R > pore_radius` cylinder test), so a chamfered
    corner point that should read as bulk (outside the tapered wall) was
    incorrectly flagged as membrane."""
    pore_radius, corner_radius, chamfer_depth, membrane_thickness = 8.0, 4.0, 4.0, 20.0
    half_t = membrane_thickness / 2.0

    # Point just inside the top face, radially beyond pore_radius but
    # within the chamfer's widened opening (edge_radius = pore_radius +
    # corner_radius = 12.0) -- this is open bore under the chamfer, solid
    # under a plain (unchamfered) cylinder mask.
    z = half_t - 1.0  # deep in the 4 A chamfer zone
    point = np.array([[10.0, 0.0, z]])

    plain_mask = _membrane_mask(
        _StubSEM("cylindrical", membrane_thickness, pore_radius=pore_radius),
        point,
    )
    chamfered_mask = _membrane_mask(
        _StubSEM(
            "cylindrical", membrane_thickness, pore_radius=pore_radius,
            corner_radius=corner_radius, chamfer_depth=chamfer_depth,
        ),
        point,
    )

    assert bool(plain_mask[0]) is True
    assert bool(chamfered_mask[0]) is False


def test_missing_required_params_returns_all_bulk():
    """conical requires top_radius/bottom_radius; without them the
    analytic fallback must degrade to "no wall" rather than raising."""
    points = _grid_points()
    mask = _membrane_mask(_StubSEM("conical", 20.0), points)
    assert not mask.any()
