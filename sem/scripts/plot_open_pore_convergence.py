#!/usr/bin/env python3
"""
Plot open-pore current versus grid size from SEM open_pore output files.

Usage:
    python plot_open_pore_convergence.py results/*.txt --output plot.png
"""

import argparse
import re
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt


class DataPoint(NamedTuple):
    source: Path
    grid_size: float
    open_pore_current: float


GRID_RE = re.compile(r"#\s*Grid\s+resolution:\s*([0-9.+-Ee]+)")
CURRENT_RE = re.compile(r"^\s*([0-9.+-Ee]+)\s*$")


def parse_result_file(path: Path) -> DataPoint:
    grid_size = None
    open_pore_current = None

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            match = GRID_RE.match(line)
            if match:
                grid_size = float(match.group(1))
                continue

            if line.startswith("#"):
                continue

            match = CURRENT_RE.match(line)
            if match:
                open_pore_current = float(match.group(1))

    if grid_size is None or open_pore_current is None:
        raise ValueError(f"Could not parse grid/current from {path}")

    return DataPoint(path, grid_size, open_pore_current)


def plot_data(points: list[DataPoint], output: Path | None):
    points.sort(key=lambda p: p.grid_size)
    grids = [p.grid_size for p in points]
    currents = [p.open_pore_current for p in points]

    plt.figure(figsize=(6, 4))
    plt.plot(grids, currents, marker="o")
    plt.xlabel("Grid size (Å)")
    plt.ylabel("Open-pore current (nA)")
    plt.title("Open-pore current vs grid size")
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    for grid, current in zip(grids, currents):
        plt.annotate(f"{current:.2f}", (grid, current), textcoords="offset points", xytext=(0, 5), ha="center")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(output, dpi=300)
        print(f"Saved plot to {output}")
    else:
        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot open-pore current vs grid size.")
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="Files or directories to search (directories are scanned recursively)."
    )
    parser.add_argument(
        "--pattern",
        default="*_open_pore_current.txt",
        help="Glob pattern to match result files inside directories (default: *_open_pore_current.txt).",
    )
    parser.add_argument("--output", type=Path, help="Optional output image file (PNG, PDF, etc.)")
    args = parser.parse_args()

    result_files: list[Path] = []
    for entry in args.paths:
        if entry.is_file():
            result_files.append(entry)
        elif entry.is_dir():
            matches = sorted(entry.rglob(args.pattern))
            if not matches:
                print(f"Warning: no files matching {args.pattern} found under {entry}")
            result_files.extend(matches)
        else:
            print(f"Warning: {entry} does not exist")

    if not result_files:
        raise SystemExit("No result files found.")

    points = [parse_result_file(path) for path in result_files]
    plot_data(points, args.output)


if __name__ == "__main__":
    main()
