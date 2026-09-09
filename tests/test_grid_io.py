"""
Behavioural tests for `sem.grid_io`, the numpy-only module holding the
condfrac ramp, binary grid I/O, and the ramp-weighted-area inverse helpers.
None of this module may import scipy or dolfinx.
"""

import re
from pathlib import Path

import numpy as np
import pytest

import sem.grid_io
from sem.grid_io import (
    RAMP_MIN,
    RAMP_MAX,
    area_coefficients,
    bin_axes,
    condfrac,
    equivalent_area,
    invert_condfrac,
    radius_from_area,
    readbinGrid,
    write_binary_file,
)


def test_grid_io_module_source_has_no_scipy_or_dolfinx_import():
    """Static guard: sem.grid_io must stay numpy-only so it (and anything
    that only needs it) imports without scipy or dolfinx installed."""
    source = Path(sem.grid_io.__file__).read_text()
    assert not re.search(r"^\s*(import|from)\s+(scipy|dolfinx)\b", source, re.MULTILINE)


def test_condfrac_ramp_values():
    # At the ramp floor the linear formula lands on exactly 0.0 (the
    # min-term cancels in floating point); values strictly below RAMP_MIN
    # are what actually clip to the 1e-7 floor (see next test).
    result = condfrac(np.array([RAMP_MIN, RAMP_MAX, (RAMP_MIN + RAMP_MAX) / 2.0]))
    assert result[0] == pytest.approx(0.0, abs=1e-9)
    assert result[1] == pytest.approx(1.0)
    assert result[2] == pytest.approx(0.5)


def test_condfrac_floors_below_ramp_min():
    result = condfrac(np.array([RAMP_MIN - 1.0, RAMP_MIN - 0.01]))
    assert result[0] == pytest.approx(1e-7)
    assert result[1] == pytest.approx(1e-7)


def test_condfrac_clips_above_ramp_max():
    result = condfrac(np.array([RAMP_MAX + 5.0]))
    assert result[0] == pytest.approx(1.0)


def test_invert_condfrac_round_trip_inside_ramp():
    d = np.linspace(RAMP_MIN + 0.05, RAMP_MAX - 0.05, 9)
    fraction = condfrac(d.copy())
    recovered = invert_condfrac(fraction)
    np.testing.assert_allclose(recovered, d, atol=1e-5)


def test_invert_condfrac_floors_at_zero_outside_ramp():
    # A fraction at/below the floor inverts to distance 0 (never negative).
    recovered = invert_condfrac(np.array([1e-7, 0.0]))
    assert np.all(recovered >= 0.0)


def test_write_binary_file_read_back_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    nz, ny, nx = 4, 5, 6
    field = rng.random((nz, ny, nx)).astype(np.float32)
    origin = (-3.5, 2.0, 10.25)
    resolution = 1.5
    out_file = tmp_path / "grid.bin"

    write_binary_file(field, origin, resolution, str(out_file))
    val3d, dims, shape, metadata = readbinGrid(str(out_file), return_metadata=True)

    # readbinGrid always reports shape as [nx, ny, nz].
    assert shape == [nx, ny, nz]
    assert val3d.shape == (nx, ny, nz)

    # write_binary_file/readbinGrid round-trip transposes (z,y,x) -> (x,y,z).
    np.testing.assert_allclose(val3d, np.transpose(field, (2, 1, 0)), atol=1e-6)

    np.testing.assert_allclose(metadata["origin"], origin, atol=1e-6)
    np.testing.assert_allclose(metadata["spacing"], [resolution] * 3, atol=1e-6)

    x, y, z = bin_axes(metadata, shape)
    assert len(x) == nx
    assert len(y) == ny
    assert len(z) == nz
    np.testing.assert_allclose(x, origin[0] + resolution * np.arange(nx), atol=1e-5)
    np.testing.assert_allclose(y, origin[1] + resolution * np.arange(ny), atol=1e-5)
    np.testing.assert_allclose(z, origin[2] + resolution * np.arange(nz), atol=1e-5)


def test_readbingrid_missing_file_raises(tmp_path):
    missing = tmp_path / "does_not_exist.bin"
    with pytest.raises(FileNotFoundError):
        readbinGrid(str(missing))


def test_area_coefficients_default():
    c1, c0 = area_coefficients()
    assert c1 == pytest.approx(5.4)
    assert c0 == pytest.approx(7.9433, abs=1e-4)


@pytest.mark.parametrize("cos_theta", [1.0, 0.894])
@pytest.mark.parametrize("a", [5, 20, 49, 99])
def test_radius_from_area_inverts_equivalent_area(a, cos_theta):
    area = equivalent_area(a, cos_theta=cos_theta)
    recovered = radius_from_area(area, cos_theta=cos_theta)
    assert recovered == pytest.approx(a, rel=1e-9)
