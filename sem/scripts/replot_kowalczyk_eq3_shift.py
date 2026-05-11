"""
Replot a Kowalczyk-2011 benchmark CSV with optional shifts applied to the
Eq 3 (Hall cylinder) curve only.

Why:
    Eq 3 (Kowalczyk 2011) models a *perfect* D = d cylinder of length l
    embedded in an infinitely thin membrane between two large reservoirs.
    SEM-cylindrical solves the same geometry but on a discrete conductivity
    map with a finite transition shell at the wall. The transition shell
    erodes the pore inward and dilates the membrane outward, so the
    *effective* d and l the FD/FEM solver sees are not exactly the
    geometric ones. This script lets you overlay an Eq 3 curve evaluated
    at (d + dd, l + dl) so you can read off the apparent shift.

    Per the validation criterion in the parent benchmark's docstring,
    SEM-cylindrical points should track the (unshifted) Hall curve at the
    real l = 20 nm. If they instead track the curve at (d - 2, l + 2),
    that's a clean signature of the conductivity-map boundary smoothing.

What this script does NOT touch:
    - Eq 1 (naive cylinder)             ← always plotted at the geometric d, l
    - Eqs 5-7 (hourglass bracket)       ← always plotted at d, l, D = d + taper
    - The 8.6 nm "paper-fit" curve      ← always plotted at l_eff = 8.6 nm

Inputs:
    A CSV produced by `benchmark_kowalczyk_2011.py`, with these columns:
        diameter_nm, geometry,
        G_naive_eq1_nS,
        G_hall_l_real_eq3_nS, G_hall_l_eff_8.6nm_nS,
        G_hourglass_lower_eq7_nS, G_hourglass_upper_eq7_nS,
        G_sem_nS

Usage (single CSV — fixed colours per geometry):
    python -m sem.scripts.replot_kowalczyk_eq3_shift \
        --csv benchmark_kowalczyk_2011/benchmark_results.csv \
        --out benchmark_kowalczyk_2011/eq3_shift.png \
        --eq3-d-offset-nm -2.0 --eq3-l-offset-nm 2.0

Usage (padding sweep — one CSV per padding, dots coloured by padding):
    python -m sem.scripts.replot_kowalczyk_eq3_shift \
        --csv pad20/benchmark_results.csv pad50/benchmark_results.csv pad100/benchmark_results.csv \
        --padding-nm 20 50 100 \
        --out benchmark_kowalczyk_2011/eq3_shift_padsweep.png \
        --eq3-d-offset-nm -2.0 --eq3-l-offset-nm 2.0
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import List, Tuple

# ---------------------------------------------------------------------------
# Analytical formulas (copied verbatim from benchmark_kowalczyk_2011.py
# so this script has no dependency on the compiled bytecode).
# ---------------------------------------------------------------------------

EFFECTIVE_THICKNESS_NM = 8.6  # the paper's empirical cylinder fit (l_eff)


def naive_cylindrical_nS(d_nm: float, l_nm: float, sigma: float) -> float:
    """G = sigma * pi d^2 / (4 l). Eq 1, Kowalczyk 2011. Returns nS."""
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    return sigma * math.pi * d**2 / (4.0 * l) * 1e9


def hall_cylindrical_nS(d_nm: float, l_nm: float, sigma: float) -> float:
    """G = sigma * [4 l / (pi d^2) + 1/d]^-1. Eq 3, Kowalczyk 2011. Returns nS.
    Hall (planar disc) access-resistance term R_access = rho / (2 d).
    """
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    inv_G = 4.0 * l / (math.pi * d**2) + 1.0 / d
    return sigma / inv_G * 1e9


def hourglass_bracket_nS(
    d_nm: float, l_nm: float, D_nm: float, sigma: float
) -> Tuple[float, float]:
    """Hourglass / hyperboloid bracket from eqs 5–7 of Kowalczyk 2011.

    Lower G uses Hille hemispherical access (rho / (pi D)) per side;
    upper G uses Hall planar-disc access (rho / (2 D)) per side.
    """
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    D = D_nm * 1e-9
    rho = 1.0 / sigma
    if D <= d:
        G = hall_cylindrical_nS(d_nm, l_nm, sigma)
        return G, G
    sin2_a = (D**2 - d**2) / (l**2 + D**2 - d**2)
    sin_a = math.sqrt(sin2_a)
    cos_a = math.sqrt(1.0 - sin2_a)
    R_hyp = (
        (2.0 * rho / (math.pi * d))
        * (sin_a / (1.0 - cos_a))
        * math.atan(math.sqrt(D**2 - d**2) / d)
    )
    R_lower_G = R_hyp + rho / D                 # Hall access  -> upper R, lower G
    R_upper_G = R_hyp + 2.0 * rho / (math.pi * D)  # Hille access -> lower R, upper G
    return (1.0 / R_lower_G) * 1e9, (1.0 / R_upper_G) * 1e9


def kowalczyk_taper_D_nm(d_nm: float, taper_nm: float) -> float:
    return d_nm + taper_nm


# ---------------------------------------------------------------------------
# CSV reader + plotting
# ---------------------------------------------------------------------------


def _read_csv(path: Path) -> List[dict]:
    rows: List[dict] = []
    with path.open("r") as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    return rows


def _f(s: str):
    s = (s or "").strip()
    return None if s == "" else float(s)


def replot(args: argparse.Namespace) -> None:
    import matplotlib

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # ---- read every CSV, tag with its padding label (or None) -------------
    csv_paths: List[Path] = list(args.csv)
    if args.padding_nm:
        if len(args.padding_nm) != len(csv_paths):
            raise SystemExit(
                f"--padding-nm has {len(args.padding_nm)} value(s) but --csv has "
                f"{len(csv_paths)} path(s); they must match 1:1."
            )
        labels: List[float | None] = list(args.padding_nm)
    else:
        labels = [None] * len(csv_paths)

    # rows_by_padding: {padding_label_or_None: [row_dict, ...]}
    rows_by_padding: dict = {}
    all_diameters: set = set()
    for path, lab in zip(csv_paths, labels):
        rs = _read_csv(path)
        if not rs:
            raise SystemExit(f"No rows in {path}")
        rows_by_padding[lab] = rs
        all_diameters.update(float(r["diameter_nm"]) for r in rs)

    diameters_nm = sorted(all_diameters)
    d_dense = np.linspace(min(diameters_nm), max(diameters_nm), 200)
    multi = len(rows_by_padding) > 1

    # Analytical references --------------------------------------------------
    naive = [naive_cylindrical_nS(d, args.membrane_thickness_nm, args.sigma) for d in d_dense]
    hall  = [hall_cylindrical_nS(d, args.membrane_thickness_nm, args.sigma) for d in d_dense]

    # Eq 3 curve evaluated at (d + dd, l + dl) — the shifted reference.
    # Only diameters that stay positive after the offset are plotted.
    dd = args.eq3_d_offset_nm
    dl = args.eq3_l_offset_nm
    shift_active = (dd != 0.0) or (dl != 0.0)
    if shift_active:
        l_shift = args.membrane_thickness_nm + dl
        d_shift = d_dense + dd
        valid = d_shift > 0
        hall_shift = [
            hall_cylindrical_nS(float(ds), l_shift, args.sigma) if v else float("nan")
            for ds, v in zip(d_shift, valid)
        ]

    # SEM points, grouped by padding label and geometry ---------------------
    # Each entry: (label_or_None, geometry, list_of_(d, G))
    sem_series: List[Tuple[object, str, List[Tuple[float, float]]]] = []
    for lab, rs in rows_by_padding.items():
        cyl = [
            (float(r["diameter_nm"]), float(r["G_sem_nS"]))
            for r in rs
            if r["geometry"] == "cylindrical" and _f(r.get("G_sem_nS")) is not None
        ]
        dc = [
            (float(r["diameter_nm"]), float(r["G_sem_nS"]))
            for r in rs
            if r["geometry"] == "double_cone" and _f(r.get("G_sem_nS")) is not None
        ]
        if cyl:
            sem_series.append((lab, "cylindrical", cyl))
        if dc:
            sem_series.append((lab, "double_cone", dc))

    # Build the colour map keyed by padding label so each padding gets a
    # distinct shade. With one padding (or none), fall back to fixed C3/C1.
    pad_keys = [k for k in rows_by_padding.keys()]
    n_pads = max(1, len(pad_keys))
    cmap = plt.get_cmap(args.cmap)
    pad_color: dict = {}
    if multi:
        # Sort padding labels for deterministic colour order; unlabeled entries
        # (None) sort to the end and get the last colour slot.
        sorted_keys = sorted(
            pad_keys,
            key=lambda k: (k is None, k if k is not None else 0.0),
        )
        for i, k in enumerate(sorted_keys):
            # Spread across [0.1, 0.9] of the colormap so end-stops aren't
            # too pale / too dark to read against the white background.
            t = 0.1 + 0.8 * (i / max(1, n_pads - 1))
            pad_color[k] = cmap(t)

    MARKER = {"cylindrical": "o", "double_cone": "s"}
    FALLBACK_COLOR = {"cylindrical": "C3", "double_cone": "C1"}

    # Plot -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    ax.plot(d_dense, naive, "--", color="gray", linewidth=1.0,
            label=f"Eq 1 naive cyl, l = {args.membrane_thickness_nm:.0f} nm")
    ax.plot(d_dense, hall, "-", color="C0", linewidth=1.6,
            label=f"Eq 3 cyl + Hall access, l = {args.membrane_thickness_nm:.0f} nm")

    if shift_active:
        sign_d = "+" if dd >= 0 else ""
        sign_l = "+" if dl >= 0 else ""
        ax.plot(d_dense, hall_shift, "-.", color="C4", linewidth=1.8,
                label=(f"Eq 3 shifted, d' = d {sign_d}{dd:g} nm, "
                       f"l' = {args.membrane_thickness_nm:.0f} {sign_l}{dl:g} nm"))

    for lab, geom, pts in sem_series:
        xs, ys = zip(*pts)
        marker = MARKER[geom]
        if multi:
            color = pad_color[lab]
            pad_str = "unlabeled" if lab is None else f"padding = {lab:g} nm"
            geom_str = "cyl" if geom == "cylindrical" else "double_cone"
            label = f"SEM {geom_str}, {pad_str}"
        else:
            color = FALLBACK_COLOR[geom]
            geom_str = ("cylindrical" if geom == "cylindrical"
                        else f"double_cone, D = d + {args.taper_nm:.0f} nm")
            label = f"SEM {geom_str} (this work)"
        ax.plot(xs, ys, marker=marker, linestyle="", color=color,
                markersize=8, label=label, zorder=5)

    ax.set_xlabel("Pore diameter d (nm)", fontsize=12, weight="bold")
    ax.set_ylabel("Open-pore conductance G (nS)", fontsize=12, weight="bold")
    ax.set_title(
        f"PETK / SEM vs Kowalczyk 2011 — open-pore G(d), σ = {args.sigma} S/m, "
        f"V = {args.voltage_mv} mV"
        + ("  [Eq 3 shifted]" if shift_active else ""),
        fontsize=12, weight="bold",
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    ax.set_xlim(left=min(diameters_nm), right=max(diameters_nm) * 1.02)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved comparison plot: {args.out}")
    if args.show:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--csv", type=Path, nargs="+", required=True,
                   help=("One or more benchmark_results.csv files produced by "
                         "benchmark_kowalczyk_2011.py. Pass multiple files (one "
                         "per padding value) to colour SEM markers by padding."))
    p.add_argument("--padding-nm", type=float, nargs="*", default=None,
                   help=("Padding label (nm) paired 1:1 with each --csv path, "
                         "in the same order. If given, SEM markers are coloured "
                         "by padding using --cmap; if omitted, single-CSV mode "
                         "uses fixed colours (C3 cyl, C1 double_cone)."))
    p.add_argument("--cmap", default="plasma",
                   help=("Matplotlib colormap used to colour SEM markers when "
                         "multiple CSVs / paddings are given. Default 'plasma' "
                         "(matches the original benchmark's plot_padding_sweep)."))
    p.add_argument("--out", type=Path, default=Path("kowalczyk_eq3_shift.png"),
                   help="Output PNG path.")

    # Reference parameters — must match the values used to generate the CSV
    # so the analytical curves are drawn in the same regime as the SEM points.
    p.add_argument("--membrane-thickness-nm", type=float, default=20.0,
                   help="Real SiN membrane thickness (nm). Default 20 (Kowalczyk 2011).")
    p.add_argument("--sigma", type=float, default=10.5,
                   help="Bulk conductivity (S/m). 1 M KCl @ 23 °C ≈ 10.5.")
    p.add_argument("--voltage-mv", type=float, default=200.0,
                   help="Voltage label for the title (does not affect curves).")
    p.add_argument("--taper-nm", type=float, default=20.0,
                   help="Empirical D − d (nm) for the hourglass bracket. Default 20.")

    # The new shift flags — default 0, applied to Eq 3 only.
    p.add_argument("--eq3-d-offset-nm", type=float, default=0.0,
                   help=("Diameter offset (nm) applied ONLY to the Eq 3 (Hall) "
                         "reference curve. The shifted curve uses d' = d + offset. "
                         "Negative values model the conductivity-map boundary "
                         "eroding the pore inward (e.g. -2 nm)."))
    p.add_argument("--eq3-l-offset-nm", type=float, default=0.0,
                   help=("Membrane-thickness offset (nm) applied ONLY to the "
                         "Eq 3 (Hall) reference curve. The shifted curve uses "
                         "l' = l + offset. Positive values model the conductivity "
                         "transition shell extending past each membrane face "
                         "(e.g. +2 nm)."))

    p.add_argument("--show", action="store_true",
                   help="Display the plot interactively (default: save only).")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    replot(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
