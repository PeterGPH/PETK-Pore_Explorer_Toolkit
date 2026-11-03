"""prevent_analyte_overlap: the hard-core check must see the membrane's actual conductivity."""
from types import SimpleNamespace

import numpy as np
import pytest

from sem.conductivity_models import SimpleConductivityModel
from sem.pore_geometry import ConicalPore, CylindricalPore, DoubleConePore, EllipticalPore

BULK, MEMBRANE = 1.12, 1e-8          # PETK's scaled values for 11.2 S/m and 1e-7 S/m
KW = dict(membrane_half_thickness=25.0, bulk_conductivity=BULK, membrane_conductivity=MEMBRANE)


def _grid(n=49, half=60.0):
    ax = np.linspace(-half, half, n)
    return np.meshgrid(ax, ax, ax, indexing="ij")


PORES = {
    "cylindrical": lambda X, Y, Z: CylindricalPore(X, Y, Z, pore_radius=20.0, corner_radius=5.0, **KW),
    "double_cone": lambda X, Y, Z: DoubleConePore(X, Y, Z, inner_radius=20.0, outer_radius=40.0, **KW),
    "conical": lambda X, Y, Z: ConicalPore(X, Y, Z, top_radius=40.0, bottom_radius=20.0, **KW),
    "elliptical": lambda X, Y, Z: EllipticalPore(X, Y, Z, semi_axis_a=20.0, semi_axis_b=15.0, **KW),
}

# well inside the membrane material of every pore above
INSIDE = np.array([[50.0, 0.0, 0.0], [0.0, -52.0, 10.0], [45.0, 45.0, -15.0]])


@pytest.mark.parametrize("name", list(PORES))
def test_membrane_conductivity_is_what_the_grid_holds(name):
    X, Y, Z = _grid()
    pore = PORES[name](X, Y, Z)
    sigma = pore.get_conductivity_interpolator()(INSIDE)
    level = pore.get_membrane_conductivity()
    np.testing.assert_allclose(sigma, level, rtol=1e-9)
    # the old check compared against membrane_conductivity, which the grid never reaches
    assert (sigma > MEMBRANE).all()


def _analyte_check(points, base, positions, threshold):
    from sem.vertical_movement_sem import VerticalMovementSEM
    sem = SimpleNamespace(cutoff=5.0, conductivity_model=SimpleConductivityModel(bulk_conductivity=BULK),
                          prevent_analyte_overlap=True, membrane_conductivity=MEMBRANE)
    if threshold is not None:
        sem.membrane_overlap_conductivity = threshold
    radii = np.full(len(positions), 1.7)
    return VerticalMovementSEM.calculate_analyte_conductivity_modification(sem, points, positions, radii, base)


def test_hard_core_overlap_with_membrane_is_rejected():
    pytest.importorskip("dolfinx")
    from sem.vertical_movement_sem import AnalyteOverlapError
    X, Y, Z = _grid()
    pore = PORES["cylindrical"](X, Y, Z)
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    base = pore.get_conductivity_interpolator()(points)
    threshold = max(MEMBRANE, pore.get_membrane_conductivity())
    in_lumen = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 1.0], [0.0, 4.0, -2.0]])
    in_membrane = in_lumen + np.array([30.0, 0.0, 0.0])

    _analyte_check(points, base, in_lumen, threshold)              # no overlap: no error
    with pytest.raises(AnalyteOverlapError):
        _analyte_check(points, base, in_membrane, threshold)
    # before the fix the threshold was membrane_conductivity, so the same overlap went unnoticed
    _analyte_check(points, base, in_membrane, None)
