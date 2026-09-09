"""
Analytic wall-distance geometry for parametric SEM pores.

This module holds ``PoreProfile``, the axisymmetric-wall representation
shared by every parametric pore type (``cylindrical``, ``double_cone``,
``conical``, and the ``profile`` radius-vs-z table), plus the two distance
metrics that can be evaluated against it:

  - ``euclidean_distance_to_profile`` (the default): the true 3-D distance
    from a point to the solid ``S = {R >= r(z), |z| <= h}``, computed as the
    2-D distance in the meridional (R, z) half-plane to the wall polyline
    (vertices ``(r_i, z_i)``) unioned with the two face rays
    ``z = +-h, R >= r(+-h)``. Exact for sloped walls (double cone, conical,
    chamfered corners); only differs from ``legacy_distance`` when the wall
    is not perfectly vertical.
  - ``legacy_distance``: the original approximation, combining an in-plane
    radial gap and a vertical overhang independently. Exact for a purely
    vertical wall (plain cylinder); overestimates the true distance on any
    sloped wall or chamfer by up to ``1 / cos(theta)``.

This module intentionally has **no** scipy or dolfinx dependency (and does
not import from ``sem.utils``), so it — and anything that only needs it —
imports cleanly in pip-only environments that never install DOLFINx/MPI.

Contents:
  - `legacy_distance`: the verbatim old radial+vertical distance formula.
  - `_segment_distance_sq`: vectorised squared point-to-segment distance.
  - `euclidean_distance_to_profile`: true wall distance for a flattened
    batch of (R, z) query points against a wall polyline.
  - `PoreProfile`: the (r, z) wall polyline + half-thickness, with
    constructors for every parametric pore shape and `distance`/
    `distance_xyz` evaluators (chunked for large grids).
"""

from dataclasses import dataclass

import numpy as np


def legacy_distance(R, abs_z, local_radius, half_thickness):
    """
    Verbatim legacy wall-distance formula (the original
    ``sem.pore_geometry._distance_to_membrane``): combines an in-plane
    radial gap and a vertical overhang independently, added in quadrature.

    Exact for a vertical wall and flat faces; on a wall inclined from
    vertical by angle ``theta`` it overestimates the true (Euclidean)
    distance by a factor of ``1 / cos(theta)``.

    ``abs_z`` must already be ``|z|`` (the caller is responsible for taking
    the absolute value); ``local_radius`` is the wall radius at each query
    point's z (e.g. from ``PoreProfile.local_radius``).
    """
    radial_term = np.maximum(local_radius - R, 0.0)
    vertical_term = np.maximum(abs_z - half_thickness, 0.0)
    return np.sqrt(radial_term ** 2 + vertical_term ** 2)


def _segment_distance_sq(Rq, Zq, r0, z0, r1, z1):
    """
    Vectorised squared Euclidean distance from each query point
    ``(Rq[i], Zq[i])`` to the line segment from ``(r0, z0)`` to
    ``(r1, z1)`` in the (R, z) meridional half-plane.
    """
    dr = r1 - r0
    dz = z1 - z0
    seg_len_sq = dr * dr + dz * dz
    if seg_len_sq == 0.0:
        # Degenerate (zero-length) segment: distance to the single point.
        pr = np.full_like(Rq, r0, dtype=float)
        pz = np.full_like(Zq, z0, dtype=float)
    else:
        t = ((Rq - r0) * dr + (Zq - z0) * dz) / seg_len_sq
        t = np.clip(t, 0.0, 1.0)
        pr = r0 + t * dr
        pz = z0 + t * dz
    dR = Rq - pr
    dZ = Zq - pz
    return dR * dR + dZ * dZ


def euclidean_distance_to_profile(R, z, vertices, half_thickness):
    """
    True Euclidean distance from flattened query points ``(R[i], z[i])`` to
    the axisymmetric solid ``S = {R' >= r(z'), |z'| <= half_thickness}``,
    whose wall is the piecewise-linear polyline ``vertices`` (rows
    ``(r_i, z_i)``, ``z`` strictly increasing, ``z[0] == -half_thickness``,
    ``z[-1] == +half_thickness``) and whose faces are the rays
    ``z = +-half_thickness, R' >= r(+-half_thickness)``.

    Points inside ``S`` get distance 0. ``R`` and ``z`` must already be
    flattened 1-D arrays of the same length (chunking is the caller's
    responsibility — see ``PoreProfile.distance``).
    """
    R = np.asarray(R, dtype=float)
    z = np.asarray(z, dtype=float)

    r_lo = vertices[0, 0]
    r_hi = vertices[-1, 0]

    local_r = np.interp(z, vertices[:, 1], vertices[:, 0])
    inside = (np.abs(z) <= half_thickness) & (R >= local_r)

    best_sq = np.full(R.shape, np.inf, dtype=float)
    for i in range(len(vertices) - 1):
        r0, z0 = vertices[i]
        r1, z1 = vertices[i + 1]
        seg_sq = _segment_distance_sq(R, z, r0, z0, r1, z1)
        best_sq = np.minimum(best_sq, seg_sq)

    # Face rays z = +-h, R' >= r(+-h). The vertical gap max(|z|-h, 0) is the
    # same for both faces; computing both and taking the min picks whichever
    # face radius gives the smaller (correct) radial gap.
    vertical_gap = np.maximum(np.abs(z) - half_thickness, 0.0)
    for r_face in (r_lo, r_hi):
        radial_gap = np.maximum(r_face - R, 0.0)
        face_sq = radial_gap ** 2 + vertical_gap ** 2
        best_sq = np.minimum(best_sq, face_sq)

    dist = np.sqrt(best_sq)
    dist = np.where(inside, 0.0, dist)
    return np.maximum(dist, 0.0)


def _load_profile_table(path):
    """
    Read a two-column ``(z, r)`` profile table from a CSV file. A header
    row is optional (any row that fails to parse as two floats is dropped);
    ``#`` starts a comment.
    """
    data = np.genfromtxt(path, delimiter=",", comments="#", dtype=float, invalid_raise=False)
    data = np.atleast_2d(data)
    if data.size == 0:
        raise ValueError(f"Profile table {path!r} contains no data rows")
    valid = ~np.isnan(data).any(axis=1)
    data = data[valid]
    if data.shape[0] == 0 or data.shape[1] < 2:
        raise ValueError(f"Profile table {path!r} must have at least 2 numeric columns (z, r)")
    return data[:, 0], data[:, 1]


@dataclass(frozen=True)
class PoreProfile:
    """
    The (r, z) wall polyline defining an axisymmetric parametric pore, plus
    its half-thickness. ``vertices`` is a ``(k, 2)`` array of ``(r_i, z_i)``
    with ``z`` strictly increasing, ``z[0] == -half_thickness`` and
    ``z[-1] == +half_thickness``. The solid region is
    ``S = {R >= r(z), |z| <= half_thickness}`` where ``r(z)`` is the
    piecewise-linear interpolation of ``vertices``.
    """

    vertices: np.ndarray
    half_thickness: float

    @classmethod
    def cylindrical(cls, pore_radius, half_thickness, corner_radius=0.0, chamfer_depth=None):
        """
        Plain cylinder (2 vertices), or a cylinder with a linearly tapered
        corner chamfer of radial extent ``corner_radius`` over an axial
        depth ``chamfer_depth`` (default: ``corner_radius``) at each face.

        Reproduces ``sem.pore_geometry.CylindricalPore._compute_local_pore_radius``
        exactly via ``local_radius`` (including the ``chamfer_depth >
        half_thickness`` edge case, where the taper never reaches
        ``pore_radius`` within the slab).
        """
        pore_radius = float(pore_radius)
        half_thickness = float(half_thickness)

        no_chamfer = (
            corner_radius is None
            or corner_radius <= 0
            or (chamfer_depth is not None and chamfer_depth <= 0)
        )
        if no_chamfer:
            vertices = np.array(
                [[pore_radius, -half_thickness], [pore_radius, half_thickness]], dtype=float
            )
            return cls(vertices=vertices, half_thickness=half_thickness)

        corner_radius = float(corner_radius)
        cd = float(chamfer_depth) if chamfer_depth is not None else corner_radius
        edge_radius = pore_radius + corner_radius

        if cd < half_thickness:
            vertices = np.array(
                [
                    [edge_radius, -half_thickness],
                    [pore_radius, -half_thickness + cd],
                    [pore_radius, half_thickness - cd],
                    [edge_radius, half_thickness],
                ],
                dtype=float,
            )
        else:
            # chamfer_depth >= half_thickness: the taper spans the whole
            # slab as a single linear segment on each side of z=0 and never
            # bottoms out at pore_radius (unless cd == half_thickness
            # exactly, where r_center reduces to pore_radius).
            r_center = edge_radius - corner_radius * (half_thickness / cd)
            vertices = np.array(
                [
                    [edge_radius, -half_thickness],
                    [r_center, 0.0],
                    [edge_radius, half_thickness],
                ],
                dtype=float,
            )
        return cls(vertices=vertices, half_thickness=half_thickness)

    @classmethod
    def double_cone(cls, inner_radius, outer_radius, half_thickness):
        """Hourglass pore: constriction ``inner_radius`` at z=0, rim
        ``outer_radius`` at ``|z| == half_thickness`` (3 vertices)."""
        half_thickness = float(half_thickness)
        vertices = np.array(
            [
                [float(outer_radius), -half_thickness],
                [float(inner_radius), 0.0],
                [float(outer_radius), half_thickness],
            ],
            dtype=float,
        )
        return cls(vertices=vertices, half_thickness=half_thickness)

    @classmethod
    def conical(cls, bottom_radius, top_radius, half_thickness):
        """Single frustum: ``bottom_radius`` at z=-half_thickness,
        ``top_radius`` at z=+half_thickness (2 vertices, signed z)."""
        half_thickness = float(half_thickness)
        vertices = np.array(
            [[float(bottom_radius), -half_thickness], [float(top_radius), half_thickness]],
            dtype=float,
        )
        return cls(vertices=vertices, half_thickness=half_thickness)

    @classmethod
    def from_table(cls, z, r, half_thickness):
        """
        Build a profile from an arbitrary radius-vs-z table (the
        ``profile`` pore type). ``z`` must be strictly increasing and span
        at least ``[-half_thickness, +half_thickness]``; the table is
        clipped/resampled to exactly ``+-half_thickness`` by linear
        interpolation. ``r`` must be positive and finite throughout.
        """
        z = np.asarray(z, dtype=float)
        r = np.asarray(r, dtype=float)
        half_thickness = float(half_thickness)

        if z.ndim != 1 or r.ndim != 1 or z.shape != r.shape:
            raise ValueError("z and r must be 1-D arrays of the same length")
        if z.size < 2:
            raise ValueError("Profile table needs at least 2 points")
        if not (np.all(np.isfinite(z)) and np.all(np.isfinite(r))):
            raise ValueError("Profile table contains non-finite values")
        if not np.all(np.diff(z) > 0):
            raise ValueError("Profile table z values must be strictly increasing")
        if not np.all(r > 0):
            raise ValueError("Profile table radii must be positive")
        if z[0] > -half_thickness or z[-1] < half_thickness:
            raise ValueError(
                f"Profile table must span at least [-{half_thickness}, {half_thickness}]; "
                f"got [{z[0]}, {z[-1]}]"
            )

        r_lo = float(np.interp(-half_thickness, z, r))
        r_hi = float(np.interp(half_thickness, z, r))
        interior = (z > -half_thickness) & (z < half_thickness)
        new_z = np.concatenate(([-half_thickness], z[interior], [half_thickness]))
        new_r = np.concatenate(([r_lo], r[interior], [r_hi]))

        vertices = np.column_stack([new_r, new_z]).astype(float)
        return cls(vertices=vertices, half_thickness=half_thickness)

    @classmethod
    def from_params(
        cls,
        pore_type,
        *,
        membrane_thickness,
        pore_radius=None,
        corner_radius=0.0,
        chamfer_depth=None,
        outer_radius=None,
        top_radius=None,
        bottom_radius=None,
        profile_path=None,
        profile_table=None,
    ):
        """
        Build a `PoreProfile` for a named parametric ``pore_type``
        ("cylindrical", "double_cone", "conical", or "profile") from the
        same keyword arguments used throughout `sem.pore_geometry` and
        `sem.vertical_movement_sem`. ``membrane_thickness`` is the full
        thickness (not half).
        """
        half_thickness = float(membrane_thickness) / 2.0
        pore_type = pore_type.lower()

        if pore_type == "cylindrical":
            if pore_radius is None:
                raise ValueError("cylindrical pore requires pore_radius")
            return cls.cylindrical(
                pore_radius, half_thickness, corner_radius=corner_radius, chamfer_depth=chamfer_depth
            )
        if pore_type == "double_cone":
            if pore_radius is None or outer_radius is None:
                raise ValueError("double_cone pore requires pore_radius (inner) and outer_radius")
            return cls.double_cone(pore_radius, outer_radius, half_thickness)
        if pore_type == "conical":
            if top_radius is None or bottom_radius is None:
                raise ValueError("conical pore requires top_radius and bottom_radius")
            return cls.conical(bottom_radius, top_radius, half_thickness)
        if pore_type == "profile":
            if profile_table is not None:
                z, r = profile_table
            elif profile_path is not None:
                z, r = _load_profile_table(profile_path)
            else:
                raise ValueError("profile pore requires profile_path or profile_table")
            return cls.from_table(z, r, half_thickness)
        raise ValueError(f"Unknown pore_type for PoreProfile.from_params: {pore_type!r}")

    def is_vertical_cylinder(self):
        """True for a plain (unchamfered) cylinder: 2 vertices at equal r.
        Used to fast-path the Euclidean metric to the (bit-identical)
        legacy formula, which is exact for a vertical wall."""
        return self.vertices.shape[0] == 2 and self.vertices[0, 0] == self.vertices[1, 0]

    def local_radius(self, z):
        """Wall radius r(z): piecewise-linear interpolation of `vertices`,
        clamped to the end vertices' radius beyond ``+-half_thickness``."""
        z = np.asarray(z, dtype=float)
        return np.interp(z, self.vertices[:, 1], self.vertices[:, 0])

    def distance(self, R, z, *, metric="euclidean", chunk=1 << 22):
        """
        Distance from each point ``(R[i], z[i])`` to the pore solid.
        ``R`` and ``z`` may be any (matching) shape; the result has the
        same shape. Evaluated in chunks of up to ``chunk`` flattened points
        so large 3-D grids (e.g. 301**3) run in bounded memory.

        ``metric="euclidean"`` (default) is the true 3-D wall distance;
        ``metric="legacy"`` is the original radial+vertical approximation.
        A plain vertical cylinder always uses the (bit-identical) legacy
        formula regardless of ``metric``, since the two are mathematically
        equal for a vertical wall.
        """
        if metric not in ("euclidean", "legacy"):
            raise ValueError(f"Unknown distance metric: {metric!r}")

        R = np.asarray(R, dtype=float)
        z = np.asarray(z, dtype=float)
        if R.shape != z.shape:
            raise ValueError("R and z must have the same shape")

        orig_shape = R.shape
        R_flat = np.ravel(R)
        z_flat = np.ravel(z)
        n = R_flat.size
        out = np.empty(n, dtype=float)

        use_legacy = metric == "legacy" or self.is_vertical_cylinder()

        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            r_chunk = R_flat[start:end]
            z_chunk = z_flat[start:end]
            if use_legacy:
                local_r = self.local_radius(z_chunk)
                out[start:end] = legacy_distance(r_chunk, np.abs(z_chunk), local_r, self.half_thickness)
            else:
                out[start:end] = euclidean_distance_to_profile(
                    r_chunk, z_chunk, self.vertices, self.half_thickness
                )
        return out.reshape(orig_shape)

    def distance_xyz(self, points, *, metric="euclidean"):
        """Distance from each Cartesian point ``points[i] = (x, y, z)`` to
        the pore solid (axisymmetric: only ``R = sqrt(x**2+y**2)`` and
        ``z`` matter, so this is rotation-invariant about the z-axis)."""
        points = np.asarray(points, dtype=float)
        R = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
        z = points[:, 2]
        return self.distance(R, z, metric=metric)
