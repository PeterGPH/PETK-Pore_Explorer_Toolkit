"""
Numpy-only binary grid I/O and condfrac-ramp helpers.

This module intentionally has **no** scipy or dolfinx dependency: it backs
`sem.pore_geometry`'s bin-file path and the standalone `gen_dist`/derivation
tooling, all of which must work in pip-only environments that never install
the DOLFINx/MPI stack.

Contents:
  - `condfrac` / `invert_condfrac`: the conductivity-fraction ramp used to
    turn a distance-to-solid-surface field into a conductivity blend, and
    its inverse (recover an approximate distance from a fraction).
  - `readbinGrid` / `write_binary_file`: the `.bin` grid file format shared
    by `gen_dist`, `sem.pore_geometry.BinFilePore`, and `resample_bin`.
  - `bin_axes`: build the 1-D coordinate axes for a grid read via
    `readbinGrid`.
  - `area_coefficients` / `equivalent_area` / `radius_from_area`: the
    ramp-weighted equivalent-area bookkeeping used to convert between a
    "ramp-weighted equivalent radius" and the true pore area under the
    condfrac ramp (ramp-weighted-area derivation tooling).
"""

import os
import sys

import numpy as np

# The condfrac ramp: conductivity fraction is 0 at RAMP_MIN (fully solid)
# and 1 at RAMP_MAX (fully bulk), linear in between.
RAMP_MIN = 1.3
RAMP_MAX = 4.1

_SLOPE = 1.0 / (RAMP_MAX - RAMP_MIN)
_INTERCEPT = -RAMP_MIN * _SLOPE


def condfrac(invec):
    """
    Convert distance-to-solid values to conductivity fractions.

    Points on the line (RAMP_MIN, 0) -> (RAMP_MAX, 1); below RAMP_MIN the
    ramp clips to a small positive floor (never a hard zero) and above
    RAMP_MAX it clips to 1.
    """
    result = _SLOPE * invec + _INTERCEPT
    result[result < 0] = 0.0000001
    result[result > 1] = 1.0
    return result


def invert_condfrac(fraction):
    """
    Invert `condfrac`: recover the distance implied by a conductivity
    fraction, assuming the fraction came from the linear ramp region.
    Never returns a negative distance.
    """
    return np.maximum((fraction - _INTERCEPT) / _SLOPE, 0.0)


def readbinGrid(name, mask_radius=-1, *, return_metadata=False):
    """
    Read binary grid file (from original code).

    Args:
        name: Path to binary grid file.
        mask_radius: Optional radius for masking values.
        return_metadata: If True, return a metadata dict with origin, spacing, and grid shape.

    Returns:
        Tuple containing the 3D values, physical dimensions, grid counts,
        and optionally metadata about the grid spacing/origin.
    """
    if not os.path.isfile(name):
        raise FileNotFoundError(f"{name} doesn't exist")

    with open(name, 'rb') as f:
        val1d = np.fromfile(f, dtype=np.float32)

    if val1d.size < 7:
        raise ValueError(f"Binary grid file {name} is too small to contain a header")

    resolution = float(val1d[6])
    if resolution <= 0:
        raise ValueError(f"Invalid grid spacing ({resolution}) recorded in {name}")

    delta = np.array([resolution, resolution, resolution], dtype=np.float32)
    origin = np.array([val1d[3], val1d[4], val1d[5]], dtype=np.float32)
    shape = (
        int(np.ceil(val1d[0])),
        int(np.ceil(val1d[1])),
        int(np.ceil(val1d[2]))
    )

    expected_values = shape[0] * shape[1] * shape[2]
    data = val1d[7:]

    if data.size != expected_values:
        raise ValueError(
            f"Binary grid {name} contains {data.size} values but expected {expected_values}"
        )

    val3d = np.reshape(data, shape, order='F')

    if mask_radius > 0:
        x_ = np.arange(origin[0], origin[0]+shape[0]*delta[0], delta[0])
        y_ = np.arange(origin[1], origin[1]+shape[1]*delta[1], delta[1])
        z_ = np.arange(origin[2], origin[2]+shape[2]*delta[2], delta[2])
        assert len(x_) == val3d.shape[0], "x is wrong size"
        assert len(y_) == val3d.shape[1], "y is wrong size"
        assert len(z_) == val3d.shape[2], "z is wrong size"
        xx, yy, zz = np.meshgrid(x_, y_, z_, indexing='ij')

        msk = xx*xx+yy*yy > mask_radius*mask_radius
        val3d[msk] = 0.00001

    L = delta[0]*shape[0]
    W = delta[1]*shape[1]
    H = delta[2]*shape[2]
    nx = int(shape[0])
    ny = int(shape[1])
    nz = int(shape[2])
    Lm = L-delta[0]
    Wm = W-delta[1]
    Hm = H-delta[2]

    if return_metadata:
        metadata = {
            "origin": origin,
            "spacing": delta,
            "grid_shape": shape,
            "resolution": resolution,
        }
        return val3d, [Lm, Wm, Hm], [nx, ny, nz], metadata

    return val3d, [Lm, Wm, Hm], [nx, ny, nz]


def bin_axes(metadata, shape):
    """
    Build the 1-D coordinate axes `(x, y, z)` for a grid read via
    `readbinGrid(..., return_metadata=True)`.

    `metadata` must carry `origin` and `spacing` (each length-3, ordered
    x/y/z), and `shape` is the `(nx, ny, nz)` grid counts (e.g. the third
    value returned by `readbinGrid`).
    """
    origin = metadata["origin"]
    spacing = metadata["spacing"]
    axes = tuple(
        origin[i] + spacing[i] * np.arange(shape[i])
        for i in range(3)
    )
    return axes


def write_binary_file(distance_field, origin, resolution, filename):
    """
    Write a `.bin` grid file: header `[nx, ny, nz, ox, oy, oz, h]` (as
    float32) followed by the data, one z-slice at a time. `distance_field`
    must be shaped `(nz, ny, nx)`.
    """
    try:
        with open(filename, 'wb') as f:
            grid_shape = distance_field.shape
            x_count = float(grid_shape[2])
            y_count = float(grid_shape[1])
            z_count = float(grid_shape[0])

            f.write(np.array([x_count], dtype=np.float32).tobytes())
            f.write(np.array([y_count], dtype=np.float32).tobytes())
            f.write(np.array([z_count], dtype=np.float32).tobytes())
            f.write(np.array(origin, dtype=np.float32).tobytes())
            f.write(np.array([resolution], dtype=np.float32).tobytes())

            for z in range(grid_shape[0]):
                slice_data = distance_field[z, :, :].astype(np.float32)
                f.write(slice_data.tobytes())

    except Exception as e:
        print(f"ERROR: Failed to write {filename}: {e}")
        sys.exit(1)


def area_coefficients(cos_theta=1.0):
    """
    Coefficients `(c1, c0)` of the ramp-weighted equivalent-area quadratic
    `equivalent_area(a) = pi * (a**2 - c1*a + c0)`, valid for a straight
    wall inclined so that `cos_theta = cos(angle from vertical)`.

    Derivation: on a wall of slope `cos_theta`, the condfrac ramp linearly
    interpolates the local radius contribution between `RAMP_MIN/cos_theta`
    and `RAMP_MAX/cos_theta`; integrating `condfrac` radially across the
    ramp band gives the two constants below (with `d = RAMP_MAX - RAMP_MIN`):

        c1 = (RAMP_MIN + RAMP_MAX) / cos_theta
        c0 = (RAMP_MAX**2 - d*RAMP_MIN - 2*d**2/3) / cos_theta**2
    """
    d = RAMP_MAX - RAMP_MIN
    c1 = (RAMP_MIN + RAMP_MAX) / cos_theta
    c0 = (RAMP_MAX**2 - d*RAMP_MIN - 2*d**2/3) / cos_theta**2
    return c1, c0


def equivalent_area(a, cos_theta=1.0):
    """
    Ramp-weighted equivalent area `A(a) = pi*(a**2 - c1*a + c0)` for a
    ramp-weighted equivalent radius `a`. Valid for `a > RAMP_MAX/cos_theta`
    (i.e. `a` beyond the outer edge of the condfrac ramp).
    """
    c1, c0 = area_coefficients(cos_theta)
    return np.pi * (a**2 - c1*a + c0)


def radius_from_area(area, cos_theta=1.0):
    """
    Exact inverse of `equivalent_area`: recover `a` from an area.

    `equivalent_area` is the quadratic `A/pi = a**2 - c1*a + c0`, i.e.
    `a**2 - c1*a + (c0 - A/pi) = 0`. Solving via the quadratic formula:

        a = c1/2 +/- sqrt(c1**2/4 - c0 + A/pi)

    The physically valid branch (a > RAMP_MAX/cos_theta, the larger root)
    takes the `+` sign, giving:

        a = c1/2 + sqrt(A/pi - c0 + c1**2/4)
    """
    c1, c0 = area_coefficients(cos_theta)
    return c1/2.0 + np.sqrt(area/np.pi - c0 + c1**2/4.0)
