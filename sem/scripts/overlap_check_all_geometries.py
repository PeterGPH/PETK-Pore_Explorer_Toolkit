#!/usr/bin/env python3
"""
Standalone overlap check for all three solid-state nanopore geometries.

For a given analyte PDB and a configurable (z, rotation) pose grid, this
script reports — per geometry — which poses are physically realisable
(analyte does not collide with the membrane) and which must be skipped.

It re-implements the same analytical predicates that
``VerticalMovementSEM._assert_radius_overlap`` uses (cylindrical,
double_cone, and the new conical branch), so it does NOT need to spin up
the FEM stack, conductivity grids, or MPI. Runs in seconds even for
thousands of poses.

Outputs:
    - Per-geometry summary printed to stdout (skip count, skip fraction,
      worst clash).
    - One CSV per geometry: pose index, z, rx, ry, rz, status, min_distance,
      worst_atom_index, local_radius_at_worst, skip_reason.
    - A combined summary CSV with per-pose pass/fail across all three
      geometries.

Usage (typical):
    python -m sem.scripts.overlap_check_all_geometries \
        --analyte petk/Demo/centered_1AOI.pdb \
        --membrane-thickness 200 \
        --cyl-pore-radius 100 \
        --dc-inner-radius 60  --dc-outer-radius 120 \
        --con-bottom-radius 60 --con-top-radius 120 \
        --z-start 110 --z-end -110 --z-step 10 \
        --rotation-grid 0,90,180,270 \
        --buffer 0.0 \
        --output-dir overlap_check_results
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Path setup so the script runs both as a module and as a plain file.
# ---------------------------------------------------------------------------
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

try:
    from sem.van_der_waals import VanDerWaalsRadii
except Exception:  # pragma: no cover — fallback if VdW table is unavailable
    VanDerWaalsRadii = None  # type: ignore[assignment]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pose specification
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RotationSpec:
    rx: float
    ry: float
    rz: float

    def label(self) -> str:
        return f"rx{self.rx:.1f}_ry{self.ry:.1f}_rz{self.rz:.1f}"


@dataclass(frozen=True)
class Pose:
    z: float
    rotation: RotationSpec

    def label(self) -> str:
        return f"z{self.z:+.2f}_{self.rotation.label()}"


# ---------------------------------------------------------------------------
# Geometry specifications. One dataclass per pore type. Each implements
# local_radius(z) and reuses _distance_to_membrane.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CylindricalSpec:
    pore_radius: float          # Å
    membrane_thickness: float   # Å (full thickness, NOT half)
    name: str = "cylindrical"

    def local_radius(self, z: np.ndarray) -> np.ndarray:
        return np.full_like(z, self.pore_radius, dtype=float)


@dataclass(frozen=True)
class DoubleConeSpec:
    inner_radius: float         # Å (constriction at z=0)
    outer_radius: float         # Å (rim at |z|=half_thickness)
    membrane_thickness: float   # Å
    name: str = "double_cone"

    def local_radius(self, z: np.ndarray) -> np.ndarray:
        half = self.membrane_thickness / 2.0
        z_frac = np.clip(np.abs(z) / half, 0.0, 1.0)
        return self.inner_radius + (self.outer_radius - self.inner_radius) * z_frac


@dataclass(frozen=True)
class ConicalSpec:
    bottom_radius: float        # Å (at z = -half_thickness)
    top_radius: float           # Å (at z = +half_thickness)
    membrane_thickness: float   # Å
    name: str = "conical"

    def local_radius(self, z: np.ndarray) -> np.ndarray:
        half = self.membrane_thickness / 2.0
        thickness = 2.0 * half
        # Asymmetric, monotone in *signed* z.
        t = np.clip((z + half) / thickness, 0.0, 1.0)
        return self.bottom_radius + (self.top_radius - self.bottom_radius) * t


GeometrySpec = CylindricalSpec | DoubleConeSpec | ConicalSpec


def _distance_to_membrane(
    R: np.ndarray, abs_z: np.ndarray, local_radius: np.ndarray, half_thickness: float
) -> np.ndarray:
    """Same formula as VerticalMovementSEM._distance_to_membrane and
    pore_geometry._distance_to_membrane. Distance from each atom to the
    membrane solid; zero means the atom is inside the wall."""
    radial_term = np.maximum(local_radius - R, 0.0)
    vertical_term = np.maximum(abs_z - half_thickness, 0.0)
    return np.sqrt(radial_term ** 2 + vertical_term ** 2)


# ---------------------------------------------------------------------------
# Per-pose overlap check
# ---------------------------------------------------------------------------
@dataclass
class PoseResult:
    pose: Pose
    geometry: str
    overlap: bool
    min_distance: float
    threshold_at_min: float
    worst_atom_index: int
    worst_atom_position: Tuple[float, float, float]
    local_radius_at_worst: float
    n_overlapping: int


def check_pose(
    geometry: GeometrySpec,
    pose: Pose,
    atom_positions: np.ndarray,
    atom_radii: np.ndarray,
    buffer: float,
    fixed_threshold: Optional[float],
) -> PoseResult:
    """Run the analytical overlap check for one pose against one geometry."""
    R = np.sqrt(atom_positions[:, 0] ** 2 + atom_positions[:, 1] ** 2)
    z = atom_positions[:, 2]
    abs_z = np.abs(z)
    half = geometry.membrane_thickness / 2.0
    local_radius = geometry.local_radius(z)
    distances = _distance_to_membrane(R, abs_z, local_radius, half)
    if fixed_threshold is not None:
        thresholds = np.full_like(distances, fixed_threshold + buffer)
    else:
        thresholds = atom_radii + buffer
    overlap_mask = distances <= thresholds
    # Worst atom = smallest (distance - threshold), i.e. deepest penetration.
    margin = distances - thresholds
    worst_idx = int(np.argmin(margin))
    return PoseResult(
        pose=pose,
        geometry=geometry.name,
        overlap=bool(np.any(overlap_mask)),
        min_distance=float(distances[worst_idx]),
        threshold_at_min=float(thresholds[worst_idx]),
        worst_atom_index=worst_idx,
        worst_atom_position=tuple(map(float, atom_positions[worst_idx])),
        local_radius_at_worst=float(local_radius[worst_idx]),
        n_overlapping=int(np.sum(overlap_mask)),
    )


# ---------------------------------------------------------------------------
# Analyte loading (PDB → positions + vdW radii)
# ---------------------------------------------------------------------------
def load_analyte(pdb_path: Path, default_radius: float = 1.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (positions [N,3], radii [N], elements [N])."""
    positions: List[List[float]] = []
    elements: List[str] = []
    with pdb_path.open() as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip()
            if not element:
                # Fall back to atom name's first letter
                element = line[12:16].strip().lstrip("0123456789")[:1].upper()
            positions.append([x, y, z])
            elements.append(element.upper())
    if not positions:
        raise ValueError(f"No ATOM/HETATM records in {pdb_path}")
    pos_arr = np.asarray(positions, dtype=float)
    elem_arr = np.asarray(elements)
    radii = np.full(len(pos_arr), default_radius, dtype=float)
    if VanDerWaalsRadii is not None:
        # VdW table lookup per element, fallback to default for unknowns.
        try:
            table = getattr(VanDerWaalsRadii, "RADII", None)
            if table:
                for i, e in enumerate(elem_arr):
                    radii[i] = float(table.get(e, default_radius))
        except Exception:
            pass
    return pos_arr, radii, elem_arr


# ---------------------------------------------------------------------------
# Pose generation: rigid translation along z + Euler-angle rotation about COM
# ---------------------------------------------------------------------------
def _rotation_matrix(rx_deg: float, ry_deg: float, rz_deg: float) -> np.ndarray:
    """Same convention as sem/rotation.py: Rz @ Ry @ Rx, active rotation."""
    rx, ry, rz = np.deg2rad([rx_deg, ry_deg, rz_deg])
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def apply_pose(
    base_positions: np.ndarray, com: np.ndarray, pose: Pose
) -> np.ndarray:
    """Rotate about COM, then translate so that the COM lands at z=pose.z
    in the lab frame (x and y of COM are preserved at 0, 0 because the
    PDB is assumed to be already centred — see petk/Demo/centered_1AOI.pdb)."""
    R_mat = _rotation_matrix(pose.rotation.rx, pose.rotation.ry, pose.rotation.rz)
    rotated = (base_positions - com) @ R_mat.T + com
    rotated[:, 2] += pose.z - com[2]
    return rotated


# ---------------------------------------------------------------------------
# Sweep orchestration
# ---------------------------------------------------------------------------
@dataclass
class Sweep:
    z_values: np.ndarray
    rotations: List[RotationSpec]

    def poses(self) -> Iterable[Pose]:
        for z in self.z_values:
            for r in self.rotations:
                yield Pose(z=float(z), rotation=r)

    @property
    def n_poses(self) -> int:
        return int(len(self.z_values) * len(self.rotations))


def run_sweep(
    geometry: GeometrySpec,
    base_positions: np.ndarray,
    radii: np.ndarray,
    sweep: Sweep,
    buffer: float,
    fixed_threshold: Optional[float],
) -> List[PoseResult]:
    com = base_positions.mean(axis=0)
    results: List[PoseResult] = []
    for pose in sweep.poses():
        positions = apply_pose(base_positions, com, pose)
        result = check_pose(
            geometry=geometry,
            pose=pose,
            atom_positions=positions,
            atom_radii=radii,
            buffer=buffer,
            fixed_threshold=fixed_threshold,
        )
        results.append(result)
    return results


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def summarise(results: Sequence[PoseResult]) -> dict:
    n_total = len(results)
    n_overlap = sum(1 for r in results if r.overlap)
    n_skip = n_overlap
    n_pass = n_total - n_overlap
    if n_overlap > 0:
        worst = min((r for r in results if r.overlap), key=lambda r: r.min_distance - r.threshold_at_min)
    else:
        worst = None
    summary = {
        "geometry": results[0].geometry if results else None,
        "n_total": n_total,
        "n_pass": n_pass,
        "n_skip": n_skip,
        "skip_fraction": n_skip / n_total if n_total else 0.0,
    }
    if worst is not None:
        summary["worst_clash"] = {
            "pose": worst.pose.label(),
            "z": worst.pose.z,
            "min_distance": worst.min_distance,
            "threshold": worst.threshold_at_min,
            "penetration": worst.threshold_at_min - worst.min_distance,
            "atom_index": worst.worst_atom_index,
            "atom_position": worst.worst_atom_position,
            "local_pore_radius": worst.local_radius_at_worst,
            "n_overlapping_atoms": worst.n_overlapping,
        }
    return summary


def write_csv(results: Sequence[PoseResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "pose_index", "z", "rx", "ry", "rz", "geometry", "status",
            "min_distance_A", "threshold_A", "penetration_A",
            "worst_atom_index", "worst_atom_x", "worst_atom_y", "worst_atom_z",
            "local_pore_radius_A", "n_overlapping_atoms",
        ])
        for i, r in enumerate(results):
            writer.writerow([
                i, r.pose.z, r.pose.rotation.rx, r.pose.rotation.ry, r.pose.rotation.rz,
                r.geometry, "SKIP" if r.overlap else "PASS",
                f"{r.min_distance:.4f}",
                f"{r.threshold_at_min:.4f}",
                f"{r.threshold_at_min - r.min_distance:.4f}",
                r.worst_atom_index,
                f"{r.worst_atom_position[0]:.3f}",
                f"{r.worst_atom_position[1]:.3f}",
                f"{r.worst_atom_position[2]:.3f}",
                f"{r.local_radius_at_worst:.4f}",
                r.n_overlapping,
            ])


def write_combined_csv(
    per_geometry: dict, sweep: Sweep, output_path: Path
) -> None:
    """One row per pose with pass/fail across all geometries."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    geom_names = list(per_geometry.keys())
    poses = list(sweep.poses())
    with output_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        header = ["pose_index", "z", "rx", "ry", "rz"] + [f"{g}_status" for g in geom_names]
        writer.writerow(header)
        for i, pose in enumerate(poses):
            row = [i, pose.z, pose.rotation.rx, pose.rotation.ry, pose.rotation.rz]
            for g in geom_names:
                r = per_geometry[g][i]
                row.append("SKIP" if r.overlap else "PASS")
            writer.writerow(row)


def print_summary_table(per_geometry_summary: dict) -> None:
    print()
    print("=" * 78)
    print(" Overlap Check Summary ".center(78, "="))
    print("=" * 78)
    header = f"{'Geometry':<14} {'Total':>7} {'Pass':>7} {'Skip':>7} {'Skip %':>9} {'Worst penetration (Å)':>25}"
    print(header)
    print("-" * 78)
    for name, summary in per_geometry_summary.items():
        worst = summary.get("worst_clash")
        worst_pen = f"{worst['penetration']:.3f}" if worst else "  —  "
        print(
            f"{name:<14} {summary['n_total']:>7} {summary['n_pass']:>7} {summary['n_skip']:>7} "
            f"{100 * summary['skip_fraction']:>8.2f}% {worst_pen:>25}"
        )
    print("=" * 78)
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_rotation_grid(spec: str) -> List[RotationSpec]:
    """Parse a rotation-grid string. Two supported shapes:
        "0,90,180,270"            → identity rx, ry; rz sweeps these values
        "rx=0,90;ry=0,90;rz=0,180" → outer product
    """
    spec = spec.strip()
    if "=" in spec:
        axes = {"rx": [0.0], "ry": [0.0], "rz": [0.0]}
        for part in spec.split(";"):
            key, vals = part.split("=")
            axes[key.strip()] = [float(v) for v in vals.split(",") if v.strip()]
        rotations = [
            RotationSpec(rx, ry, rz)
            for rx in axes["rx"]
            for ry in axes["ry"]
            for rz in axes["rz"]
        ]
    else:
        rz_values = [float(v) for v in spec.split(",") if v.strip()]
        rotations = [RotationSpec(0.0, 0.0, rz) for rz in rz_values]
    if not rotations:
        rotations = [RotationSpec(0.0, 0.0, 0.0)]
    return rotations


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--analyte", type=Path, required=True, help="Path to analyte PDB (centred at origin).")
    p.add_argument("--default-radius", type=float, default=1.5, help="Fallback vdW radius for unknown elements (Å).")

    p.add_argument("--membrane-thickness", type=float, required=True, help="Full membrane thickness L (Å).")

    # Cylindrical
    p.add_argument("--cyl-pore-radius", type=float, default=None, help="Cylindrical pore radius (Å). Omit to skip cylindrical.")

    # Double cone
    p.add_argument("--dc-inner-radius", type=float, default=None, help="Double-cone constriction radius at z=0 (Å).")
    p.add_argument("--dc-outer-radius", type=float, default=None, help="Double-cone rim radius at |z|=L/2 (Å).")

    # Conical
    p.add_argument("--con-bottom-radius", type=float, default=None, help="Conical bottom-face radius at z=-L/2 (Å).")
    p.add_argument("--con-top-radius", type=float, default=None, help="Conical top-face radius at z=+L/2 (Å).")

    # Sweep
    p.add_argument("--z-start", type=float, required=True)
    p.add_argument("--z-end", type=float, required=True)
    p.add_argument("--z-step", type=float, required=True)
    p.add_argument("--rotation-grid", type=str, default="0",
                   help="Either a comma list of rz angles, or rx=...;ry=...;rz=... outer product.")

    # Predicate
    p.add_argument("--buffer", type=float, default=0.0, help="Extra clearance added to thresholds (Å).")
    p.add_argument("--fixed-threshold", type=float, default=None,
                   help="Fixed distance threshold per atom (Å). If set, overrides per-atom vdW radii.")

    p.add_argument("--output-dir", type=Path, default=Path("overlap_check_results"))
    p.add_argument("--quiet", action="store_true")

    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Build geometry specs
    geometries: List[GeometrySpec] = []
    if args.cyl_pore_radius is not None:
        geometries.append(CylindricalSpec(
            pore_radius=args.cyl_pore_radius,
            membrane_thickness=args.membrane_thickness,
        ))
    if args.dc_inner_radius is not None and args.dc_outer_radius is not None:
        geometries.append(DoubleConeSpec(
            inner_radius=args.dc_inner_radius,
            outer_radius=args.dc_outer_radius,
            membrane_thickness=args.membrane_thickness,
        ))
    if args.con_bottom_radius is not None and args.con_top_radius is not None:
        geometries.append(ConicalSpec(
            bottom_radius=args.con_bottom_radius,
            top_radius=args.con_top_radius,
            membrane_thickness=args.membrane_thickness,
        ))
    if not geometries:
        logger.error(
            "No geometries specified. Provide at least one of "
            "--cyl-pore-radius, (--dc-inner-radius and --dc-outer-radius), "
            "or (--con-bottom-radius and --con-top-radius)."
        )
        return 2

    # Load analyte
    logger.info("Loading analyte: %s", args.analyte)
    base_positions, radii, elements = load_analyte(args.analyte, default_radius=args.default_radius)
    logger.info("Loaded %d atoms (mean vdW radius %.3f Å)", len(base_positions), radii.mean())

    # Build sweep
    if args.z_step == 0:
        logger.error("--z-step must be non-zero.")
        return 2
    if (args.z_end - args.z_start) * args.z_step < 0:
        # Step has the wrong sign — auto-fix
        step = -abs(args.z_step) if args.z_start > args.z_end else abs(args.z_step)
    else:
        step = args.z_step
    n_steps = int(np.floor((args.z_end - args.z_start) / step + 1e-9)) + 1
    z_values = args.z_start + step * np.arange(n_steps)
    rotations = _parse_rotation_grid(args.rotation_grid)
    sweep = Sweep(z_values=z_values, rotations=rotations)
    logger.info(
        "Sweep: %d z-values × %d rotations = %d poses per geometry",
        len(z_values), len(rotations), sweep.n_poses,
    )

    # Run
    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_geometry_results: dict = {}
    per_geometry_summary: dict = {}
    for geom in geometries:
        logger.info("Running %s …", geom.name)
        results = run_sweep(
            geometry=geom,
            base_positions=base_positions,
            radii=radii,
            sweep=sweep,
            buffer=args.buffer,
            fixed_threshold=args.fixed_threshold,
        )
        per_geometry_results[geom.name] = results
        summary = summarise(results)
        per_geometry_summary[geom.name] = summary
        write_csv(results, args.output_dir / f"poses_{geom.name}.csv")

    # Combined CSV (one row per pose, columns per geometry)
    write_combined_csv(per_geometry_results, sweep, args.output_dir / "poses_combined.csv")

    # JSON summary
    with (args.output_dir / "summary.json").open("w") as fh:
        json.dump(per_geometry_summary, fh, indent=2, default=str)

    # Pretty print
    if not args.quiet:
        print_summary_table(per_geometry_summary)
        print(f"Per-geometry CSVs and combined CSV written to: {args.output_dir}")
        print(f"Summary JSON: {args.output_dir / 'summary.json'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
