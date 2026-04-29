"""
Utilities for rotating PDB structures around the centre of a conductivity grid.

This reproduces the orientation used by scipy.ndimage.rotate when rotating the
distance/conductivity map, enabling consistent orientation scans.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple, List, Optional

import numpy as np

logger = logging.getLogger(__name__)


def _require_biopython():
    try:
        from Bio.PDB import PDBParser, PDBIO  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for rotation support. Install with `pip install biopython`."
        ) from exc
    return PDBParser, PDBIO


def _require_griddata():
    try:
        from gridData import Grid  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "gridData is required to read DX grids for rotation. Install with `pip install gridDataFormats`."
        ) from exc
    return Grid


def _rotation_matrix(rx: float, ry: float, rz: float) -> np.ndarray:
    """Return combined rotation matrix (Rx -> Ry -> Rz)."""
    rx_rad, ry_rad, rz_rad = np.deg2rad([rx, ry, rz])

    cx, sx = np.cos(rx_rad), np.sin(rx_rad)
    cy, sy = np.cos(ry_rad), np.sin(ry_rad)
    cz, sz = np.cos(rz_rad), np.sin(rz_rad)

    rx_mat = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry_mat = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz_mat = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])

    return rz_mat @ ry_mat @ rx_mat


def _rotate_about_pivot(coords: np.ndarray, rotation: np.ndarray, pivot: np.ndarray) -> np.ndarray:
    """Rotate coordinates (N,3) around pivot using active rotation matrix."""
    return (coords - pivot) @ rotation.T + pivot


@dataclass
class RotationSpec:
    rx: float
    ry: float
    rz: float

    def label(self) -> str:
        return f"rx{self.rx:.2f}_ry{self.ry:.2f}_rz{self.rz:.2f}"


def rotate_pdb_to_grid_center(
    input_pdb: Path,
    dx_path: Optional[Path],
    rotation: RotationSpec,
    output_path: Path,
) -> Path:
    """
    Rotate `input_pdb` so that it matches the scipy.ndimage.rotate convention
    used on the conductivity/distance grid stored in `dx_path`. If `dx_path`
    is None, the structure is rotated about its centre of mass.
    """
    PDBParser, PDBIO = _require_biopython()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(input_pdb))
    atoms = list(structure.get_atoms())
    if not atoms:
        raise ValueError(f"No atoms found in {input_pdb}")

    coords = np.array([atom.coord for atom in atoms], dtype=float)

    if dx_path is not None:
        Grid = _require_griddata()
        grid = Grid(str(dx_path))
        origin = np.array(grid.origin, dtype=float)
        delta = np.array(grid.delta, dtype=float)
        nx, ny, nz = grid.grid.shape
        pivot = origin + delta * np.array([(nx - 1) / 2.0, (ny - 1) / 2.0, (nz - 1) / 2.0])
    else:
        try:
            masses = np.array([getattr(atom, "mass", 0.0) for atom in atoms], dtype=float)
        except Exception:
            masses = np.zeros(len(atoms))
        if not np.any(masses > 0):
            masses = np.ones(len(atoms))
        pivot = np.average(coords, axis=0, weights=masses)

    rotation_matrix = _rotation_matrix(rotation.rx, rotation.ry, rotation.rz)
    rotated_coords = _rotate_about_pivot(coords, rotation_matrix, pivot)

    structure_copy = structure.copy()
    for atom, new_coord in zip(structure_copy.get_atoms(), rotated_coords):
        atom.coord = new_coord

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure_copy)
    io.save(str(output_path))

    logger.debug("Saved rotated structure to %s", output_path)
    return output_path


def rotation_matrix_from_spec(rotation: RotationSpec) -> np.ndarray:
    """Return rotation matrix corresponding to the given RotationSpec."""
    return _rotation_matrix(rotation.rx, rotation.ry, rotation.rz)


def parse_angle_file(path: Path, limit: int | None = None) -> Iterable[RotationSpec]:
    """
    Parse a whitespace-separated file containing rx ry rz per line.
    Lines starting with '#' or blank lines are ignored.
    """
    rotations = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rx, ry, rz = map(float, parts[:3])
            except ValueError:
                continue
            rotations.append(RotationSpec(rx, ry, rz))
            if limit is not None and len(rotations) >= limit:
                break
    return rotations


def _quaternion_to_matrix(q: np.ndarray) -> np.ndarray:
    """Convert unit quaternion [x, y, z, w] to rotation matrix."""
    x, y, z, w = q
    return np.array([
        [1 - 2 * (y**2 + z**2), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x**2 + z**2), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x**2 + y**2)],
    ])


def _matrix_to_euler_xyz(rotation: np.ndarray) -> RotationSpec:
    """
    Convert rotation matrix (Rz @ Ry @ Rx convention) to Euler angles in degrees.
    """
    sy = -rotation[2, 0]
    cy = np.clip(np.sqrt(1 - sy**2), 1e-12, None)

    if np.isclose(cy, 0.0):
        # Gimbal lock: ry = ±90°, set rx = 0, solve for rz
        ry = np.degrees(np.arcsin(sy))
        rx = 0.0
        rz = np.degrees(np.arctan2(-rotation[0, 1], rotation[1, 1]))
    else:
        ry = np.degrees(np.arcsin(sy))
        rx = np.degrees(np.arctan2(rotation[2, 1], rotation[2, 2]))
        rz = np.degrees(np.arctan2(rotation[1, 0], rotation[0, 0]))

    return RotationSpec(rx=rx, ry=ry, rz=rz)


def random_uniform_rotations(count: int, seed: int | None = None) -> List[RotationSpec]:
    """
    Generate `count` random orientations uniformly on SO(3), returning Euler angles.
    Uses Shoemake's algorithm (uniform quaternions).
    """
    rng = np.random.default_rng(seed)
    rotations: List[RotationSpec] = []
    for _ in range(count):
        u1, u2, u3 = rng.uniform(0.0, 1.0, 3)
        q = np.array([
            np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
            np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
            np.sqrt(u1) * np.sin(2 * np.pi * u3),
            np.sqrt(u1) * np.cos(2 * np.pi * u3),
        ])
        R = _quaternion_to_matrix(q)
        rotations.append(_matrix_to_euler_xyz(R))
    return rotations


__all__ = [
    "RotationSpec",
    "rotate_pdb_to_grid_center",
    "parse_angle_file",
    "rotation_matrix_from_spec",
    "random_uniform_rotations",
]
