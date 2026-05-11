#!/usr/bin/env python3
"""
Compare two OpenDX (.dx) grids produced by the SEM ARBD-export pipeline.

Used by results/compare.sh in arbd_test runs:

    python sem/scripts/compare_arbd_dx.py \
        --open-dx     pore_*_openpore_open_pore_phi.dx \
        --with-dx     pore_*_z+0001.0A_open_pore_phi.dx \
        --label       "phi (V)" \
        --output-prefix phi_TRP_compare

Produces:
  {prefix}_slices.png         3-panel slice comparison (open / with / diff)
  {prefix}_zprofile.png       on-axis (x=y=0) profile vs Z
  {prefix}_axis_profile.csv   raw numeric profile for both grids
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def _load_dx(path):
    from gridData import Grid  # type: ignore
    g = Grid(str(path))
    return g.grid, np.array(g.delta), np.array(g.origin)


def _on_axis_profile(grid, delta, origin):
    """Profile along z at x=y=0; falls back to grid centre if x=y=0 isn't on the grid."""
    nx, ny, nz = grid.shape
    ix = int(round((0.0 - origin[0]) / delta[0]))
    iy = int(round((0.0 - origin[1]) / delta[1]))
    ix = max(0, min(nx - 1, ix))
    iy = max(0, min(ny - 1, iy))
    z = origin[2] + np.arange(nz) * delta[2]
    return z, grid[ix, iy, :]


def _slice_xz(grid, delta, origin):
    """y=0 slice for imshow."""
    nx, ny, nz = grid.shape
    iy = int(round((0.0 - origin[1]) / delta[1]))
    iy = max(0, min(ny - 1, iy))
    xz = grid[:, iy, :].T  # (nz, nx) so origin='lower' makes Z point up
    extent = [
        origin[0],
        origin[0] + nx * delta[0],
        origin[2],
        origin[2] + nz * delta[2],
    ]
    return xz, extent


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--open-dx",   required=True, help="Open-pore reference DX file")
    p.add_argument("--with-dx",   required=True, help="DX file from an analyte-included step")
    p.add_argument("--label",     default="value",
                   help='Quantity label for plot axes / colorbars (e.g. "phi (V)")')
    p.add_argument("--output-prefix", default="dx_compare",
                   help="Prefix for output files")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    open_path = Path(args.open_dx)
    with_path = Path(args.with_dx)
    if not open_path.exists():
        sys.exit(f"Open-pore DX not found: {open_path}")
    if not with_path.exists():
        sys.exit(f"With-analyte DX not found: {with_path}")

    g_open, d_open, o_open = _load_dx(open_path)
    g_with, d_with, o_with = _load_dx(with_path)

    if g_open.shape != g_with.shape:
        sys.exit(
            f"Shape mismatch: {open_path.name}={g_open.shape}  vs  "
            f"{with_path.name}={g_with.shape}"
        )
    if not (np.allclose(d_open, d_with) and np.allclose(o_open, o_with)):
        logger.warning(
            "Grid spacing/origin differ between the two DX files; comparing element-wise "
            "but x/y/z coordinates may be misaligned (delta_open=%s delta_with=%s).",
            d_open.tolist(), d_with.tolist(),
        )

    diff = g_with - g_open

    # Lazy import: matplotlib only required if we actually plot.
    import matplotlib.pyplot as plt

    # --- Slice comparison ---------------------------------------------------
    xz_open,  extent = _slice_xz(g_open,  d_open, o_open)
    xz_with,  _      = _slice_xz(g_with,  d_with, o_with)
    xz_diff,  _      = _slice_xz(diff,    d_open, o_open)

    vmax = float(np.nanmax(np.abs([xz_open, xz_with])))
    dmax = float(np.nanmax(np.abs(xz_diff))) or 1.0

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, data, title, vmag in (
        (axes[0], xz_open, f"Open pore: {open_path.name}", vmax),
        (axes[1], xz_with, f"With analyte: {with_path.name}", vmax),
        (axes[2], xz_diff, "Difference (with - open)", dmax),
    ):
        im = ax.imshow(data, origin="lower", extent=extent, cmap="seismic",
                       vmin=-vmag, vmax=vmag, aspect="equal")
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Z (A)")
        ax.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax, label=args.label)
    slices_path = f"{args.output_prefix}_slices.png"
    fig.savefig(slices_path, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", slices_path)

    # --- On-axis profile ----------------------------------------------------
    z_open, p_open = _on_axis_profile(g_open, d_open, o_open)
    z_with, p_with = _on_axis_profile(g_with, d_with, o_with)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.plot(z_open, p_open, label="Open pore", lw=2)
    ax.plot(z_with, p_with, label="With analyte", lw=2)
    ax.set_xlabel("Z (A)")
    ax.set_ylabel(args.label)
    ax.set_title("On-axis profile (x = y = 0)")
    ax.grid(alpha=0.3)
    ax.legend()
    zprof_path = f"{args.output_prefix}_zprofile.png"
    fig.savefig(zprof_path, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", zprof_path)

    csv_path = f"{args.output_prefix}_axis_profile.csv"
    np.savetxt(
        csv_path,
        np.column_stack([z_open, p_open, p_with, p_with - p_open]),
        header="z_A,open_pore,with_analyte,difference",
        delimiter=",",
        comments="",
    )
    logger.info("Wrote %s", csv_path)

    diff_stats = {
        "max_abs_diff": float(np.nanmax(np.abs(diff))),
        "rms_diff":     float(np.sqrt(np.nanmean(diff ** 2))),
        "shape":        tuple(int(s) for s in g_open.shape),
    }
    logger.info("Diff stats: %s", diff_stats)


if __name__ == "__main__":
    main()
