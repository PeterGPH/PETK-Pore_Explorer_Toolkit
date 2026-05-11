#!/usr/bin/env python3
"""
Box-size convergence plot for a padding-sweep benchmark CSV.

Input:  benchmark_padding_sweep.csv from
        `benchmark_kowalczyk_2011.py --padding-sweep ...`

What it does (different from the default plot the benchmark already makes):

  • Default plot is G(d) per padding — answers "does SEM follow the
    analytical Hall curve at a given box size".
  • This plot is G(box) per diameter — answers "for a fixed pore, how does
    the SEM result converge as the box grows, and where is the asymptote".

Two panels:

  Left  : G_SEM vs effective box edge L_box, one curve per diameter,
          with a dashed horizontal line at the analytical Hall G(d, l=l_real).
          Markers move down toward the dashed line as L_box grows.

  Right : Relative deviation eps(L) = (G_SEM(L) - G_Hall) / G_Hall vs L_box.
          Same curves; the dashed line at eps=0 is the analytical asymptote.
          eps(small L) is the box artefact magnitude in fractional terms.

Optional convergence fit (default on):

  For each diameter with >= 3 padding points, fit
      G(L) = G_inf + a / L
  by ordinary least squares on (1/L, G). G_inf is the extrapolated infinite-box
  value; the residual G_inf - G_Hall is the implementation-level error.
  The fit and its extrapolation are reported in stdout and on the plot.

Usage:
  python plot_padding_sweep_convergence.py PATH/TO/benchmark_padding_sweep.csv
  python plot_padding_sweep_convergence.py CSV --output-png convergence.png
  python plot_padding_sweep_convergence.py CSV --no-fit
  python plot_padding_sweep_convergence.py CSV --geometry cylindrical
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional, Sequence


# Kowalczyk-style empirical taper (rim_diameter - constriction_diameter)
KOWALCZYK_TAPER_NM = 20.0


def kowalczyk_widest_nm(d_nm: float, geometry: str, taper_nm: float) -> float:
    """Return the widest pore radius in nm — the radial extent that has to
    fit inside the box. Cylindrical = d; double_cone = d + taper."""
    if geometry == "double_cone":
        return d_nm + taper_nm
    return d_nm


def effective_box_edge_nm(
    pad_nm: float, d_nm: float, l_nm: float, geometry: str, taper_nm: float
) -> float:
    """Effective box edge in nm. The benchmark sizes XY = widest + 2*pad,
    Z = l + 2*pad — they are *not* the same when the pore is wide. Use the
    XY edge for convergence plotting since the access region is dominated
    by the lateral extent.
    """
    widest = kowalczyk_widest_nm(d_nm, geometry, taper_nm)
    return widest + 2.0 * pad_nm


def fit_inv_L(L_arr, G_arr):
    """Least-squares fit G = G_inf + a / L. Returns (G_inf, a, r2)."""
    import numpy as np
    L = np.asarray(L_arr, dtype=float)
    G = np.asarray(G_arr, dtype=float)
    inv_L = 1.0 / L
    A = np.column_stack([np.ones_like(inv_L), inv_L])
    sol, *_ = np.linalg.lstsq(A, G, rcond=None)
    G_inf, a = float(sol[0]), float(sol[1])
    G_pred = G_inf + a * inv_L
    ss_res = float(np.sum((G - G_pred) ** 2))
    ss_tot = float(np.sum((G - G.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return G_inf, a, r2


def load_rows(path: Path, geometry_filter: Optional[str]) -> list:
    """Read a combined benchmark_padding_sweep.csv (which has padding_nm)."""
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if geometry_filter and r["geometry"] != geometry_filter:
                continue
            try:
                if not r["G_sem_nS"]:
                    continue  # Skip rows where SEM didn't run / failed.
                rows.append({
                    "padding_nm": float(r["padding_nm"]),
                    "diameter_nm": float(r["diameter_nm"]),
                    "geometry": r["geometry"],
                    "G_sem_nS": float(r["G_sem_nS"]),
                    "G_hall_l_real_nS": float(r["G_hall_l_real_eq3_nS"]),
                    "G_hourglass_lower_nS": float(r["G_hourglass_lower_eq7_nS"]),
                    "G_hourglass_upper_nS": float(r["G_hourglass_upper_eq7_nS"]),
                })
            except (KeyError, ValueError):
                continue
    return rows


_PAD_DIR_RE = __import__("re").compile(r"pad_(\d+(?:\.\d+)?)nm$")


def load_rows_from_dir(root: Path, geometry_filter: Optional[str]) -> list:
    """Walk ``root/pad_<N>nm/benchmark_pad<N>nm.csv`` files and assemble a
    combined row list in memory. Padding is parsed from the subdir name.

    This is the recovery path for an interrupted padding sweep where the
    top-level combined CSV was never written but the per-padding CSVs were.
    """
    rows = []
    if not root.is_dir():
        raise FileNotFoundError(f"Not a directory: {root}")
    pad_dirs = sorted(p for p in root.iterdir() if p.is_dir())
    found_any = False
    for pad_dir in pad_dirs:
        m = _PAD_DIR_RE.match(pad_dir.name)
        if not m:
            continue
        pad_nm = float(m.group(1))
        # The per-padding CSV is named to match the dir.
        candidates = list(pad_dir.glob("benchmark_pad*.csv"))
        if not candidates:
            continue
        per_csv = candidates[0]
        found_any = True
        with per_csv.open() as fh:
            reader = csv.DictReader(fh)
            for r in reader:
                if geometry_filter and r["geometry"] != geometry_filter:
                    continue
                try:
                    if not r["G_sem_nS"]:
                        continue
                    rows.append({
                        "padding_nm": pad_nm,
                        "diameter_nm": float(r["diameter_nm"]),
                        "geometry": r["geometry"],
                        "G_sem_nS": float(r["G_sem_nS"]),
                        "G_hall_l_real_nS": float(r["G_hall_l_real_eq3_nS"]),
                        "G_hourglass_lower_nS": float(r["G_hourglass_lower_eq7_nS"]),
                        "G_hourglass_upper_nS": float(r["G_hourglass_upper_eq7_nS"]),
                    })
                except (KeyError, ValueError):
                    continue
    if not found_any:
        raise FileNotFoundError(
            f"No pad_<N>nm/benchmark_pad*.csv files under {root}"
        )
    return rows


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("csv_path", type=Path, nargs="?", default=None,
                   help="Path to benchmark_padding_sweep.csv (combined). "
                        "Optional if --from-dir is given.")
    p.add_argument("--from-dir", type=Path, default=None,
                   help="Directory containing pad_<N>nm/benchmark_pad<N>nm.csv "
                        "subdirs. Use this when an interrupted padding sweep "
                        "left the per-padding CSVs but never wrote the "
                        "top-level combined CSV.")
    p.add_argument("--geometry", choices=["cylindrical", "double_cone"],
                   default="cylindrical",
                   help="Which geometry to plot (default: cylindrical).")
    p.add_argument("--membrane-thickness-nm", type=float, default=20.0,
                   help="Membrane thickness used for the runs (nm).")
    p.add_argument("--taper-nm", type=float, default=KOWALCZYK_TAPER_NM,
                   help="Kowalczyk D − d taper (nm). Only matters for double_cone.")
    p.add_argument("--output-png", type=Path, default=None,
                   help="Output PNG path. Defaults to convergence.png "
                        "next to the input CSV.")
    p.add_argument("--no-fit", action="store_true",
                   help="Disable the G(L) = G_inf + a/L extrapolation fit.")
    p.add_argument("--show", action="store_true",
                   help="Display the plot interactively.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    if args.from_dir is not None:
        try:
            rows = load_rows_from_dir(args.from_dir, args.geometry)
        except FileNotFoundError as e:
            print(str(e), file=sys.stderr)
            return 2
        source_path = args.from_dir
    else:
        if args.csv_path is None:
            print("Provide either CSV_PATH or --from-dir DIR.", file=sys.stderr)
            return 2
        if not args.csv_path.exists():
            print(f"CSV not found: {args.csv_path}", file=sys.stderr)
            return 2
        rows = load_rows(args.csv_path, args.geometry)
        source_path = args.csv_path

    if not rows:
        print(f"No usable {args.geometry} rows in {source_path}", file=sys.stderr)
        return 2

    # Bucket rows by diameter
    by_d = defaultdict(list)
    for r in rows:
        by_d[r["diameter_nm"]].append(r)

    # Compute effective box edge for each row
    for r in rows:
        r["L_box_nm"] = effective_box_edge_nm(
            r["padding_nm"], r["diameter_nm"], args.membrane_thickness_nm,
            r["geometry"], args.taper_nm,
        )

    # Sort each diameter's rows by L_box
    for d_nm in by_d:
        by_d[d_nm].sort(key=lambda r: r["L_box_nm"])

    # Plot
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_abs, ax_rel) = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = plt.get_cmap("plasma")
    sorted_d = sorted(by_d.keys())
    n_d = len(sorted_d)

    print()
    print("=" * 70)
    print(" Box-size convergence per diameter ".center(70, "="))
    print("=" * 70)
    print(f"{'d (nm)':>7} {'#pts':>5} {'G_inf (nS)':>11} {'a (nS·nm)':>12} "
          f"{'r2':>6} {'G_Hall (nS)':>12} {'G_inf-Hall':>11}")
    print("-" * 70)

    for k, d_nm in enumerate(sorted_d):
        sub = by_d[d_nm]
        Ls = [r["L_box_nm"] for r in sub]
        Gs = [r["G_sem_nS"] for r in sub]
        G_hall = sub[0]["G_hall_l_real_nS"]
        color = cmap(k / max(1, n_d - 1))

        # Left panel: absolute G
        ax_abs.plot(Ls, Gs, marker="o", linestyle="-", color=color,
                    linewidth=1.2, markersize=7, label=f"d = {d_nm:.0f} nm")
        ax_abs.axhline(G_hall, linestyle="--", color=color,
                       linewidth=0.8, alpha=0.6)

        # Right panel: relative deviation
        eps = [(g - G_hall) / G_hall for g in Gs]
        ax_rel.plot(Ls, eps, marker="o", linestyle="-", color=color,
                    linewidth=1.2, markersize=7, label=f"d = {d_nm:.0f} nm")

        # Fit and extrapolation
        if not args.no_fit and len(sub) >= 3:
            G_inf, a, r2 = fit_inv_L(Ls, Gs)
            # Plot the fitted curve from the smallest L to a long extension
            import numpy as np
            L_dense = np.linspace(min(Ls) * 0.95, max(Ls) * 1.5, 80)
            G_fit = G_inf + a / L_dense
            ax_abs.plot(L_dense, G_fit, ":", color=color, linewidth=0.8, alpha=0.7)
            # Mark extrapolated G_inf at the right edge of the abs panel
            ax_abs.scatter([L_dense[-1]], [G_inf], marker="*", color=color,
                           s=80, edgecolor="black", linewidth=0.5, zorder=10)
            print(f"{d_nm:>7.1f} {len(sub):>5d} {G_inf:>11.3f} {a:>12.3f} "
                  f"{r2:>6.3f} {G_hall:>12.3f} {G_inf - G_hall:>+11.3f}")
        else:
            print(f"{d_nm:>7.1f} {len(sub):>5d} {'(no fit)':>11} "
                  f"{'':>12} {'':>6} {G_hall:>12.3f} {'':>11}")

    print("-" * 70)
    print("G_inf is the extrapolated infinite-box SEM value (G(L)=G_inf+a/L fit).")
    print("G_Hall is the analytical eq 3 with the real membrane thickness.")
    print("Residual G_inf − G_Hall ideally < a few % of G_Hall.")

    ax_abs.set_xlabel("Box edge L_box (nm)", fontsize=12, weight="bold")
    ax_abs.set_ylabel("Open-pore conductance G (nS)", fontsize=12, weight="bold")
    ax_abs.set_title(f"Box-size convergence — {args.geometry} geometry",
                     fontsize=12, weight="bold")
    ax_abs.grid(True, alpha=0.3)
    ax_abs.legend(fontsize=9, loc="best")
    ax_abs.text(
        0.99, 0.02,
        "dashed: analytical Hall asymptote (l = real)\n"
        "dotted: G(L) = G_inf + a/L fit\n"
        "stars : extrapolated G_inf",
        transform=ax_abs.transAxes,
        ha="right", va="bottom", fontsize=8,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"),
    )

    ax_rel.set_xlabel("Box edge L_box (nm)", fontsize=12, weight="bold")
    ax_rel.set_ylabel(r"(G$_{\rm SEM}$ − G$_{\rm Hall}$) / G$_{\rm Hall}$",
                      fontsize=12, weight="bold")
    ax_rel.set_title("Relative deviation from Hall (l = real) reference",
                     fontsize=12, weight="bold")
    ax_rel.axhline(0.0, linestyle="--", color="gray", linewidth=1.0, alpha=0.6)
    ax_rel.grid(True, alpha=0.3)
    ax_rel.legend(fontsize=9, loc="best")

    plt.tight_layout()

    if args.output_png is not None:
        out_png = args.output_png
    elif args.from_dir is not None:
        out_png = args.from_dir / "padding_sweep_convergence.png"
    else:
        out_png = args.csv_path.parent / "padding_sweep_convergence.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=150)
    print()
    print(f"Saved: {out_png}")
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
