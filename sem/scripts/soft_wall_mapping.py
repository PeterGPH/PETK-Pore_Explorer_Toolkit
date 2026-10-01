#!/usr/bin/env python3
"""
Reproduce the soft-wall -> sharp-wall mapping of SEM cylindrical pores.

Runs the axisymmetric reference solver (sem.axisym_conductance) for

  * kappa : the exact end-correction coefficient kappa(L/a) of a sharp-walled
            cylinder (Richardson-extrapolated over three grids);
  * ramp  : the converged SEM-ramp open-pore resistance on a (d, L) grid,
            compared with Hall's formula at the area-conserving offset
            (d - 2t, L + 2t) and with sem.open_pore_theory's full model;
  * box   : the same ramp map in the finite golden-aspect box of the paper
            (equal-area radius (d + 2p)/sqrt(pi), height
            max(L + 2p, 1.2 (d + 2p) + L), insulating side wall);

and writes CSV tables plus a summary figure.

Usage
-----
    python -m sem.scripts.soft_wall_mapping all --out validation/soft_wall_mapping --jobs 4
    python -m sem.scripts.soft_wall_mapping figure --out validation/soft_wall_mapping

The full run takes ~10 min on 4 cores. Lengths are in Angstrom inside the
solver; tables report d and L in nm.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from sem import open_pore_theory as opt
from sem.axisym_conductance import AxisymmetricPore, pore_resistance

T_DEF = opt.deficit_length()

KAPPA_LAMBDAS = [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
                 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 20.0]
RAMP_D_NM = [1.5, 2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20, 25, 30, 40, 50]
RAMP_L_NM = [1, 2, 3, 5, 10, 15, 20, 30, 40]
BOX_D_NM = [5, 10, 15, 20, 25]
BOX_PD = [0.35, 0.5, 1.0, 2.0]


# ---------------------------------------------------------------- workers
def _kappa_job(args):
    lam, level = args
    a = 100.0
    L = lam * a
    h_fine = min(0.4, 0.1 * L) / 2 ** level
    res = pore_resistance(AxisymmetricPore.cylinder(a, L), "sharp",
                          h_fine=h_fine, band=min(8.0, 0.4 * a))
    kappa = (res.resistance * math.pi * a ** 2 - L) / (2 * a)
    return dict(lam=lam, level=level, h_fine=h_fine, kappa=kappa)


def _ramp_job(args):
    d_nm, L_nm = args
    t0 = time.time()
    res = pore_resistance(AxisymmetricPore.cylinder(5 * d_nm, 10 * L_nm), "ramp", h_fine=0.1)
    return dict(d_nm=d_nm, L_nm=L_nm, R_ramp=res.resistance, seconds=time.time() - t0)


def _box_job(args):
    d_nm, pd = args
    p = pd * d_nm
    side = d_nm + 2 * p
    res = pore_resistance(AxisymmetricPore.cylinder(5 * d_nm, 200.0), "ramp", h_fine=0.1,
                          far_field="box", box_radius=10 * side / math.sqrt(math.pi),
                          box_height=10 * max(20 + 2 * p, 1.2 * side + 20))
    return dict(d_nm=d_nm, L_nm=20, p_over_d=pd, R_box=res.resistance)


def _run(job, tasks, jobs):
    if jobs <= 1:
        return [job(t) for t in tasks]
    with Pool(jobs) as pool:
        return list(pool.imap(job, tasks))


def _write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _read(path):
    with open(path) as fh:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]


# ------------------------------------------------------------------ steps
def run_kappa(out: Path, jobs: int):
    raw = _run(_kappa_job, [(lam, k) for lam in KAPPA_LAMBDAS for k in range(3)], jobs)
    rows = []
    for lam in KAPPA_LAMBDAS:
        k = [r["kappa"] for r in sorted((r for r in raw if r["lam"] == lam), key=lambda r: r["level"])]
        d1, d2 = k[1] - k[0], k[2] - k[1]
        order = math.log2(d1 / d2) if d1 * d2 > 0 else float("nan")
        p = order if 0.8 < order < 2.5 else 4.0 / 3.0   # 270-degree corner: O(h^(4/3))
        k_ext = k[2] + d2 / (2 ** p - 1)
        rows.append(dict(L_over_a=lam, kappa=k_ext, kappa_err=abs(k_ext - k[2]), order=order,
                         kappa_closed_form=float(opt.end_correction_coefficient(lam))))
    _write(out / "kappa.csv", rows)
    return rows


def run_ramp(out: Path, jobs: int):
    rows = []
    for r in _run(_ramp_job, list(itertools.product(RAMP_D_NM, RAMP_L_NM)), jobs):
        d, L = 10 * r["d_nm"], 10 * r["L_nm"]
        R = r["R_ramp"] * 1e10               # Ohm at sigma = 1 S/m
        R_paper = opt.sem_cylinder_resistance(d, L, model="paper")
        R_full = opt.sem_cylinder_resistance(d, L, model="full")
        a_e = opt.equal_area_radius(d / 2)
        R_norim = opt.sharp_cylinder_resistance(2 * a_e, L + 2 * T_DEF)
        rows.append(dict(d_nm=r["d_nm"], L_nm=r["L_nm"], R_ramp_rho_per_A=r["R_ramp"],
                         dev_paper_pct=100 * (R / R_paper - 1),
                         dev_no_rim_pct=100 * (R / R_norim - 1),
                         dev_full_pct=100 * (R / R_full - 1)))
    _write(out / "ramp_grid.csv", rows)
    return rows


def run_box(out: Path, jobs: int):
    ramp = {(r["d_nm"], r["L_nm"]): r["R_ramp_rho_per_A"] for r in _read(out / "ramp_grid.csv")}
    rows = []
    for r in _run(_box_job, list(itertools.product(BOX_D_NM, BOX_PD)), jobs):
        R_inf = ramp[(float(r["d_nm"]), 20.0)]
        d = 10 * r["d_nm"]
        R_paper = opt.sem_cylinder_resistance(d, 200.0, model="paper") / 1e10
        rows.append(dict(d_nm=r["d_nm"], p_over_d=r["p_over_d"],
                         box_vs_infinite_pct=100 * (R_inf / r["R_box"] - 1),
                         box_vs_paper_reference_pct=100 * (R_paper / r["R_box"] - 1)))
    _write(out / "golden_box.csv", rows)
    return rows


# ----------------------------------------------------------------- figure
LIGHT = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a",
             ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", mid="#f0efec", red="#e34948")


def make_figure(out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    from scipy.optimize import brentq

    c = LIGHT
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": c["ink2"], "axes.labelcolor": c["ink"],
                         "xtick.color": c["ink2"], "ytick.color": c["ink2"], "axes.linewidth": 0.6,
                         "grid.color": c["grid"], "grid.linewidth": 0.6, "legend.frameon": False})
    kap = _read(out / "kappa.csv")
    ramp = _read(out / "ramp_grid.csv")
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)

    # (A) kappa(L/a)
    ax = axs[0, 0]
    lam = np.logspace(-2.2, 1.4, 300)
    ax.axhline(math.pi / 4, color=c["orange"], lw=1.2)
    ax.axhline(opt.KAPPA_INF, color=c["ink2"], lw=0.8)
    ax.plot(lam, opt.end_correction_coefficient(lam), color=c["blue"], lw=1.5, label="closed form")
    ax.plot([r["L_over_a"] for r in kap], [r["kappa"] for r in kap], "o", ms=4.5, mfc=c["blue"],
            mec="white", mew=0.8, label="axisymmetric solve")
    ax.text(0.0075, math.pi / 4 + 0.0012, r"Hall: $\pi/4$", color=c["orange"], fontsize=8)
    ax.text(0.0075, opt.KAPPA_INF + 0.0012, r"flanged tube: 0.8216", color=c["ink2"], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(6e-3, 25)
    ax.set_ylim(0.78, 0.83)
    ax.set_xlabel("membrane thickness / pore radius,  $L/a$")
    ax.set_ylabel(r"end correction per end / $a$,  $\kappa$")
    ax.legend(loc="center right", fontsize=8)
    ax.grid(True, which="major")
    ax.set_title("A  Exact access term of a sharp cylinder", loc="left", fontsize=9)

    # (B) deviations at L = 20 nm
    ax = axs[0, 1]
    sel = sorted((r for r in ramp if r["L_nm"] == 20), key=lambda r: r["d_nm"])
    d = [r["d_nm"] for r in sel]
    ax.axhline(0, color=c["ink2"], lw=0.6)
    for key, col, lab in [("dev_paper_pct", c["orange"], "Hall at $(d-0.54, L+0.54)$"),
                          ("dev_no_rim_pct", c["aqua"], r"exact $\kappa$, no rim term"),
                          ("dev_full_pct", c["blue"], "full mapping")]:
        ax.plot(d, [r[key] for r in sel], "-o", color=col, lw=1.5, ms=4, mec="white", mew=0.6, label=lab)
    ax.set_xlabel("nominal pore diameter $d$ (nm)")
    ax.set_ylabel(r"$G_\mathrm{theory}/G_\mathrm{SEM} - 1$ (%)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True)
    ax.set_title("B  Theory vs converged SEM, $L$ = 20 nm", loc="left", fontsize=9)

    # (C) map of the current mapping's error over (d, L)
    ax = axs[1, 0]
    D = sorted({r["d_nm"] for r in ramp})
    Ls = sorted({r["L_nm"] for r in ramp})
    Z = np.array([[next(r["dev_paper_pct"] for r in ramp if r["d_nm"] == dd and r["L_nm"] == ll)
                   for dd in D] for ll in Ls])
    cmap = LinearSegmentedColormap.from_list("div", [c["blue"], c["mid"], c["red"]])
    lim = np.abs(Z).max()
    im = ax.pcolormesh(np.arange(len(D) + 1), np.arange(len(Ls) + 1), Z, cmap=cmap,
                       norm=TwoSlopeNorm(0, -lim, lim), edgecolors="white", linewidth=0.8)
    ax.set_xticks(np.arange(len(D)) + 0.5, [f"{x:g}" for x in D], fontsize=7, rotation=90)
    ax.set_yticks(np.arange(len(Ls)) + 0.5, [f"{x:g}" for x in Ls], fontsize=7)
    ax.set_xlabel("nominal pore diameter $d$ (nm)")
    ax.set_ylabel("membrane thickness $L$ (nm)")
    cb = fig.colorbar(im, ax=ax, shrink=0.9)
    cb.set_label("Hall-offset error (%)", fontsize=8)
    full_max = max(abs(r["dev_full_pct"]) for r in ramp if r["d_nm"] >= 2)
    ax.set_title(f"C  Hall-offset error; full mapping |err| ≤ {full_max:.2f}%", loc="left", fontsize=9)

    # (D) rim coefficient collapse
    ax = axs[1, 1]
    lam_e, coef = [], []
    for r in ramp:
        if r["d_nm"] < 2:
            continue
        a_e = float(opt.equal_area_radius(5 * r["d_nm"]))
        L_e0 = 10 * r["L_nm"] + 2 * T_DEF
        # thickness of the sharp-walled pore (radius a_e) that conducts like the map
        L_e = brentq(lambda x: opt.sharp_cylinder_resistance(2 * a_e, x) / 1e10
                     - r["R_ramp_rho_per_A"], 0.5 * L_e0, L_e0)
        lam_e.append(L_e0 / a_e)
        coef.append(0.5 * (L_e0 - L_e) * a_e ** (1 / 3) / T_DEF ** (4 / 3))
    ax.plot(lam_e, coef, "o", ms=4.5, mfc=c["blue"], mec="white", mew=0.8, label="from SEM solves")
    lg = np.logspace(np.log10(min(lam_e)), np.log10(max(lam_e)), 200)
    a_ref = 50.0
    c_fit = opt.rim_length_deficit(a_ref, lg * a_ref) * a_ref ** (1 / 3) / T_DEF ** (4 / 3)
    ax.plot(lg, c_fit, color=c["ink"], lw=1.2, label="closed form")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda = (L + 2t)/a_e$")
    ax.set_ylabel(r"rim deficit $\Delta\ell\, a_e^{1/3}/t^{4/3}$, $\Delta\ell = (L+2t-L_e)/2$")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_title("D  Rounded-rim correction collapses", loc="left", fontsize=9)

    for ext in ("pdf", "png"):
        fig.savefig(out / f"soft_wall_mapping.{ext}", dpi=200)
    return out / "soft_wall_mapping.png"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["kappa", "ramp", "box", "figure", "all"])
    ap.add_argument("--out", type=Path, default=Path("validation/soft_wall_mapping"))
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.step in ("kappa", "all"):
        run_kappa(args.out, args.jobs)
    if args.step in ("ramp", "all"):
        run_ramp(args.out, args.jobs)
    if args.step in ("box", "all"):
        run_box(args.out, args.jobs)
    if args.step in ("figure", "all"):
        print("wrote", make_figure(args.out))


if __name__ == "__main__":
    main()
