"""EllipticalPore: a == b must reproduce CylindricalPore exactly."""
import numpy as np
import pytest

from sem.pore_geometry import CylindricalPore, EllipticalPore, PoreGeometry


def _grid(n=41, half=60.0):
    ax = np.linspace(-half, half, n)
    return np.meshgrid(ax, ax, ax, indexing="ij")


def _sigma(pore, X, Y, Z):
    interp = pore.get_conductivity_interpolator()
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    return interp(pts)


@pytest.mark.parametrize("radius", [15.0, 25.0, 40.0])
def test_circular_limit_matches_cylindrical(radius):
    """The whole design rests on this: a == b is the circular pore."""
    X, Y, Z = _grid()
    kw = dict(membrane_half_thickness=25.0, bulk_conductivity=11.2,
              membrane_conductivity=1e-7)
    cyl = CylindricalPore(X, Y, Z, pore_radius=radius, **kw)
    ell = EllipticalPore(X, Y, Z, semi_axis_a=radius, semi_axis_b=radius, **kw)
    np.testing.assert_allclose(_sigma(ell, X, Y, Z), _sigma(cyl, X, Y, Z),
                               rtol=0, atol=1e-12)


def test_anisotropic_pore_is_anisotropic():
    """a != b must break the xy symmetry -- otherwise the feature does nothing."""
    X, Y, Z = _grid()
    ell = EllipticalPore(X, Y, Z, semi_axis_a=40.0, semi_axis_b=15.0,
                         membrane_half_thickness=25.0,
                         bulk_conductivity=11.2, membrane_conductivity=1e-7)
    interp = ell.get_conductivity_interpolator()
    # (25, 0, 0) is inside the long axis; (0, 25, 0) is outside the short one.
    on_a = interp(np.array([[25.0, 0.0, 0.0]]))[0]
    on_b = interp(np.array([[0.0, 25.0, 0.0]]))[0]
    assert on_a > on_b, "point inside semi-major axis should be more conductive"
    assert on_a == pytest.approx(11.2, rel=1e-3), "well inside the pore -> bulk"
    assert on_b == pytest.approx(1e-7, abs=1e-3), "inside the membrane -> insulating"


def test_wall_distance_exact_on_axes():
    """On the axes the linearisation must give the exact distance."""
    a, b = 40.0, 15.0
    X, Y, Z = (np.array([[[10.0]]]), np.array([[[0.0]]]), np.array([[[0.0]]]))
    ell = EllipticalPore(X, Y, Z, semi_axis_a=a, semi_axis_b=b,
                         membrane_half_thickness=25.0)
    assert ell._wall_distance()[0, 0, 0] == pytest.approx(a - 10.0, rel=1e-12)
    X, Y = np.array([[[0.0]]]), np.array([[[6.0]]])
    ell = EllipticalPore(X, Y, Z, semi_axis_a=a, semi_axis_b=b,
                         membrane_half_thickness=25.0)
    assert ell._wall_distance()[0, 0, 0] == pytest.approx(b - 6.0, rel=1e-12)


def test_on_axis_distance_is_min_semi_axis():
    """grad G vanishes on the axis; the fallback must be min(a, b)."""
    z = np.array([[[0.0]]])
    ell = EllipticalPore(np.array([[[0.0]]]), np.array([[[0.0]]]), z,
                         semi_axis_a=40.0, semi_axis_b=15.0,
                         membrane_half_thickness=25.0)
    assert ell._wall_distance()[0, 0, 0] == pytest.approx(15.0)


def test_factory_and_validation():
    X, Y, Z = _grid(n=9)
    p = PoreGeometry.create_pore("elliptical", X, Y, Z, semi_axis_a=30.0,
                                 semi_axis_b=20.0, membrane_half_thickness=25.0)
    assert isinstance(p, EllipticalPore)
    for a, b in ((0.0, 20.0), (30.0, -1.0)):
        with pytest.raises(ValueError):
            EllipticalPore(X, Y, Z, semi_axis_a=a, semi_axis_b=b,
                           membrane_half_thickness=25.0)
