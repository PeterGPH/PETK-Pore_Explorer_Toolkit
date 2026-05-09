#!/usr/bin/env python3
"""
Compare two ARBD .dx files — typically an open-pore (no-analyte) export
against a with-analyte snapshot from the same SEM run.

Outputs three things per invocation:
  1. <prefix>_slices.png — three side-by-side XZ slices through y=0:
       (a) open-pore field
       (b) with-analyte field (same colour scale as a)
       (c) difference (with - open), symmetric colour scale
  2. <prefix>_zprofile.png — two-row figure:
       (top)    1-D profile along the pore axis (x=y=0) for both fields
       (bottom) the difference along the same axis
  3. <prefix>_axis_profile.csv — z, open_value, with_value, diff per row

Usage:
  python sem/scripts/compare_arbd_dx.py \\
      --open-dx pore_cylindrical_openpore_POT.dx \\
      --with-dx pore_cylindrical_z+0000.0A_POT.dx \\
      --label "K+ potential (kcal/mol)" \\
      --output-prefix POT_TRP_compare
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional, Sequence


def _load_grid(path: Path):
    try:
        from gridData import Grid
    except ImportError:
        print("gridData not installed. pip install gridDataFormats", file=sys.stderr)
        raise
    return Grid(str(path))


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--open-dx", required=True, type=Path,
                   help="Path to the open-pore (no-analyte) .dx file.")
    p.add_argument("--with-dx", required=True, type=Path,
                   help="Path to the with-analyte .dx file.")
    p.add_argument("--label", default="value",
                   help="Quantity label for the colour bars / y-axes.")
    p.add_argument("--output-prefix", default="arbd_compare",
                   help="Prefix for output files (default: arbd_compare).")
    p.add_argument("--slice-axis", default="y", choices=["x", "y", "z"],
                   help="Axis perpendicular to the slice plane (default: y → XZ view).")
    p.add_argument("--show", action="store_true",
                   help="Display the figures interactively after saving.")
    args = p.parse_args(argv)

    if not args.open_dx.exists():
        print(f"open-dx not found: {args.open_dx}", file=sys.stderr); return 2
    if not args.with_dx.exists():
        print(f"with-dx not found: {args.with_dx}", file=sys.stderr); return 2

    import numpy as np
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g_open = _load_grid(args.open_dx)
    g_with = _load_grid(args.with_dx)

    if g_open.grid.shape != g_with.grid.shape:
        print(
            f"ERROR: shape mismatch — open {g_open.grid.shape} vs "
            f"with {g_with.grid.shape}. The two files must come from the "
            f"same SEM run with the same export grid.",
            file=sys.stderr,
        )
        return 2

    diff = g_with.grid - g_open.grid
    nx, ny, nz = g_open.grid.shape
    x0, y0, z0 = g_open.origin
    dx, dy, dz = g_open.delta

    print(f"Loaded grids: shape {(nx, ny, nz)}  "
          f"origin {(x0, y0, z0)}  delta {(dx, dy, dz)}")

    # --- Slice setup --------------------------------------------------------
    axis_idx = {"x": 0, "y": 1, "z": 2}[args.slice_axis]
    mid = [nx // 2, ny // 2, nz // 2][axis_idx]

    if args.slice_axis == "y":
        # XZ slice; rows = z, cols = x
        s_open = g_open.grid[:, mid, :].T
        s_with = g_with.grid[:, mid, :].T
        s_diff = diff[:, mid, :].T
        extent = [x0, x0 + nx * dx, z0, z0 + nz * dz]
        xlabel, ylabel = "X (Å)", "Z (Å)"
    elif args.slice_axis == "x":
        s_open = g_open.grid[mid, :, :].T
        s_with = g_with.grid[mid, :, :].T
        s_diff = diff[mid, :, :].T
        extent = [y0, y0 + ny * dy, z0, z0 + nz * dz]
        xlabel, ylabel = "Y (Å)", "Z (Å)"
    else:  # z
        s_open = g_open.grid[:, :, mid].T
        s_with = g_with.grid[:, :, mid].T
        s_diff = diff[:, :, mid].T
        extent = [x0, x0 + nx * dx, y0, y0 + ny * dy]
        xlabel, ylabel = "X (Å)", "Y (Å)"

    # Shared colour scale for the two field panels; symmetric scale for diff.
    vmin = float(np.nanmin(np.array([s_open.min(), s_with.min()])))
    vmax = float(np.nanmax(np.array([s_open.max(), s_with.max()])))
    abs_diff = float(np.nanmax(np.abs(s_diff))) if np.any(np.isfinite(s_diff)) else 1.0
    if abs_diff == 0:
        abs_diff = 1.0

    # --- Figure 1: three side-by-side slices --------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    im0 = axes[0].imshow(s_open, origin="lower", extent=extent,
                         vmin=vmin, vmax=vmax, cmap="seismic", aspect="equal")
    axes[0].set_title("Open pore (no analyte)")
    axes[0].set_xlabel(xlabel); axes[0].set_ylabel(ylabel)
    fig.colorbar(im0, ax=axes[0], label=args.label, shrink=0.8)

    im1 = axes[1].imshow(s_with, origin="lower", extent=extent,
                         vmin=vmin, vmax=vmax, cmap="seismic", aspect="equal")
    axes[1].set_title("With analyte")
    axes[1].set_xlabel(xlabel)
    fig.colorbar(im1, ax=axes[1], label=args.label, shrink=0.8)

    im2 = axes[2].imshow(s_diff, origin="lower", extent=extent,
                         vmin=-abs_diff, vmax=abs_diff, cmap="seismic",
                         aspect="equal")
    axes[2].set_title("Difference  (with − open)")
    axes[2].set_xlabel(xlabel)
    fig.colorbar(im2, ax=axes[2], label=f"Δ {args.label}", shrink=0.8)

    fig.suptitle(f"{args.open_dx.name}    vs    {args.with_dx.name}",
                 fontsize=10)
    plt.tight_layout()
    out_slices = f"{args.output_prefix}_slices.png"
    plt.savefig(out_slices, dpi=150)
    print(f"Saved: {out_slices}")

    # --- Figure 2: 1-D profile along the pore axis (x=y=0) ------------------
    x_idx = int(round((0.0 - x0) / dx))
    y_idx = int(round((0.0 - y0) / dy))
    x_idx = max(0, min(x_idx, nx - 1))
    y_idx = max(0, min(y_idx, ny - 1))
    z_axis = z0 + np.arange(nz) * dz
    prof_open = g_open.grid[x_idx, y_idx, :]
    prof_with = g_with.grid[x_idx, y_idx, :]
    prof_diff = prof_with - prof_open

    fig2, (axA, axB) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1]})
    axA.plot(z_axis, prof_open, color="C0", lw=2.0, label="open pore")
    axA.plot(z_axis, prof_with, color="C3", lw=2.0, label="with analyte")
    axA.set_ylabel(args.label, fontsize=12, weight="bold")
    axA.legend(fontsize=10, loc="best")
    axA.grid(True, alpha=0.3)
    axA.set_title(f"Profile along pore axis (X = {x0 + x_idx*dx:.1f} Å, "
                  f"Y = {y0 + y_idx*dy:.1f} Å)",
                  fontsize=11, weight="bold")

    axB.plot(z_axis, prof_diff, color="C2", lw=2.0, label="with − open")
    axB.axhline(0.0, color="gray", lw=0.6, alpha=0.7)
    axB.set_xlabel("Z (Å)", fontsize=12, weight="bold")
    axB.set_ylabel(f"Δ {args.label}", fontsize=11, weight="bold")
    axB.grid(True, alpha=0.3)

    plt.tight_layout()
    out_profile = f"{args.output_prefix}_zprofile.png"
    plt.savefig(out_profile, dpi=150)
    print(f"Saved: {out_profile}")

    # --- CSV: per-z values along the pore axis ------------------------------
    out_csv = f"{args.output_prefix}_axis_profile.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z_A", "open_pore", "with_analyte", "diff"])
        for k in range(nz):
            w.writerow([f"{z_axis[k]:.4f}",
                        f"{prof_open[k]:.6e}",
                        f"{prof_with[k]:.6e}",
                        f"{prof_diff[k]:.6e}"])
    print(f"Saved: {out_csv}")

    # Summary on stdout
    print()
    print(f"Δ on axis : min {prof_diff.min():.4e},  max {prof_diff.max():.4e},  "
          f"|max| {np.abs(prof_diff).max():.4e}")
    print(f"Δ on slice: min {s_diff.min():.4e},  max {s_diff.max():.4e},  "
          f"|max| {abs_diff:.4e}")

    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
