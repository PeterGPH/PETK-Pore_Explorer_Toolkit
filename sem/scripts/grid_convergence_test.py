#!/usr/bin/env python3
"""
Utility script to run a simple grid-size convergence study by either generating
distance fields at multiple resolutions or loading existing binary grids (e.g.,
conductivity maps) and comparing sampled values.
"""

if __name__ == "__main__" and __package__ is None:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))

import argparse
from pathlib import Path
import numpy as np

try:
    from .gen_dist import generate_binary_distance_field
except ImportError:
    try:
        from sem.scripts.gen_dist import generate_binary_distance_field
    except ImportError:
        from gen_dist import generate_binary_distance_field

try:
    from ..utils import readbinGrid
except ImportError:
    try:
        from sem.utils import readbinGrid
    except ImportError:
        from utils import readbinGrid


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Compare distance/conductivity grids across resolutions by generating them from XYZ or loading existing BIN files."
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--xyz-file",
        type=Path,
        help="Input XYZ file describing the pore structure (generates fresh grids).",
    )
    source_group.add_argument(
        "--bin-files",
        nargs="+",
        type=Path,
        help="Existing BIN files to compare (skips XYZ generation).",
    )

    parser.add_argument(
        "--bounds-min",
        nargs=3,
        type=float,
        metavar=("XMIN", "YMIN", "ZMIN"),
        default=None,
        help="Minimum XYZ bounds (Å) for the distance grid (required with --xyz-file).",
    )
    parser.add_argument(
        "--bounds-max",
        nargs=3,
        type=float,
        metavar=("XMAX", "YMAX", "ZMAX"),
        default=None,
        help="Maximum XYZ bounds (Å) for the distance grid (required with --xyz-file).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Distance cutoff (Å) when generating distance fields (XYZ mode).",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=float,
        default=[1.5, 1.2, 1.0, 0.8, 0.6],
        help="List of grid resolutions (Å) to test when generating from XYZ, coarse to fine.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=200,
        help="Number of random sample points used to compare grids.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for sampling points.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("grid_convergence"),
        help="Directory to place generated binary grids (XYZ mode).",
    )
    parser.add_argument(
        "--keep-binaries",
        action="store_true",
        help="Keep generated binary files (default removes them on success, XYZ mode).",
    )

    args = parser.parse_args()

    if args.xyz_file is not None:
        if args.bounds_min is None or args.bounds_max is None:
            parser.error("--bounds-min and --bounds-max are required when using --xyz-file")
    else:
        # bin-file mode does not use these parameters
        args.bounds_min = None
        args.bounds_max = None

    return args


def _sample_points(bounds_min, bounds_max, samples, seed):
    rng = np.random.default_rng(seed)
    bounds_min = np.asarray(bounds_min, dtype=float)
    bounds_max = np.asarray(bounds_max, dtype=float)
    return rng.uniform(bounds_min, bounds_max, size=(samples, 3))


def _trilinear_sample(grid, origin, spacing, points):
    shape = grid.shape
    spacing = np.asarray(spacing, dtype=float)
    origin = np.asarray(origin, dtype=float)
    max_index = np.array(shape) - 1
    results = np.full(len(points), np.nan, dtype=np.float64)

    for idx, point in enumerate(points):
        rel = (point - origin) / spacing
        if np.any(rel < 0) or np.any(rel > max_index):
            continue

        x0 = int(np.floor(rel[0]))
        y0 = int(np.floor(rel[1]))
        z0 = int(np.floor(rel[2]))

        x1 = min(x0 + 1, max_index[0])
        y1 = min(y0 + 1, max_index[1])
        z1 = min(z0 + 1, max_index[2])

        tx = rel[0] - x0
        ty = rel[1] - y0
        tz = rel[2] - z0

        c000 = grid[x0, y0, z0]
        c100 = grid[x1, y0, z0]
        c010 = grid[x0, y1, z0]
        c110 = grid[x1, y1, z0]
        c001 = grid[x0, y0, z1]
        c101 = grid[x1, y0, z1]
        c011 = grid[x0, y1, z1]
        c111 = grid[x1, y1, z1]

        c00 = c000 * (1 - tx) + c100 * tx
        c01 = c001 * (1 - tx) + c101 * tx
        c10 = c010 * (1 - tx) + c110 * tx
        c11 = c011 * (1 - tx) + c111 * tx

        c0 = c00 * (1 - ty) + c10 * ty
        c1 = c01 * (1 - ty) + c11 * ty

        results[idx] = c0 * (1 - tz) + c1 * tz

    return results


def _generate_single_grid(args, resolution, output_dir):
    output_path = output_dir / f"{args.xyz_file.stem}_{resolution:.2f}A.bin"
    metadata = generate_binary_distance_field(
        args.xyz_file,
        args.bounds_min[0],
        args.bounds_min[1],
        args.bounds_min[2],
        args.bounds_max[0],
        args.bounds_max[1],
        args.bounds_max[2],
        resolution,
        args.cutoff,
        output_path,
    )
    grid, dims, shape, read_meta = readbinGrid(output_path, return_metadata=True)
    return {
        "resolution": read_meta["resolution"],
        "grid": grid,
        "shape": shape,
        "dimensions": dims,
        "origin": read_meta["origin"],
        "spacing": read_meta["spacing"],
        "bin_file": output_path,
        "metadata": metadata,
    }


def _load_existing_bin(bin_path):
    grid, dims, shape, read_meta = readbinGrid(bin_path, return_metadata=True)
    spacing = read_meta["spacing"]
    resolution = read_meta.get("resolution", spacing[0])
    return {
        "resolution": resolution,
        "grid": grid,
        "shape": shape,
        "dimensions": dims,
        "origin": read_meta["origin"],
        "spacing": spacing,
        "bin_file": bin_path,
        "metadata": read_meta,
    }


def run_convergence_test(args):
    datasets = []

    if args.bin_files:
        for bin_path in args.bin_files:
            datasets.append(_load_existing_bin(bin_path))

        # Sample domain from first dataset
        base = datasets[0]
        bounds_min = np.asarray(base["origin"], dtype=float)
        bounds_max = bounds_min + (np.asarray(base["shape"]) - 1) * np.asarray(base["spacing"])
        points = _sample_points(bounds_min, bounds_max, args.samples, args.seed)

    else:
        output_dir = args.output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        points = _sample_points(args.bounds_min, args.bounds_max, args.samples, args.seed)

        for res in args.resolutions:
            grid_info = _generate_single_grid(args, res, output_dir)
            datasets.append(grid_info)

    # Sample all datasets on the chosen points
    for grid_info in datasets:
        grid_info["samples"] = _trilinear_sample(
            grid_info["grid"],
            grid_info["origin"],
            grid_info["spacing"],
            points,
        )

    reference = min(datasets, key=lambda entry: entry["resolution"])
    ref_values = reference["samples"]

    print("Grid convergence summary (relative to finest grid):")
    print(f"Reference resolution: {reference['resolution']:.3f} Å")
    print(f"Sample points: {len(points)}")
    print("")

    for data in datasets:
        diff = np.abs(data["samples"] - ref_values)
        valid = np.isfinite(diff)
        mean_abs = float(np.mean(diff[valid])) if np.any(valid) else float("nan")
        max_abs = float(np.max(diff[valid])) if np.any(valid) else float("nan")
        grid_points = data["grid"].size
        print(
            f"{data['resolution']:.3f} Å grid -> "
            f"{grid_points:,} voxels | "
            f"mean|Δ|={mean_abs:.4f} Å | "
            f"max|Δ|={max_abs:.4f} Å | "
            f"bin: {data['bin_file']}"
        )

    if args.bin_files is None and not args.keep_binaries:
        for data in datasets:
            try:
                data["bin_file"].unlink()
            except OSError:
                pass


def main():
    args = _parse_args()
    run_convergence_test(args)


if __name__ == "__main__":
    main()
