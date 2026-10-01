#!/usr/bin/env python3
"""
Reproduce the soft-wall -> sharp-wall mapping of SEM cylindrical pores.

Steps (each writes a CSV into --out; `report` and `figure` only read them):

  ritz      exact kappa(L/a) and rim strength k(L/a) of a sharp cylinder
            (sem.cylinder_mode_matching), the inputs of sem.open_pore_theory;
  rimshape  thick-membrane rim coefficient C for linear ramps of different
            shape (r_min/r_cut), at fixed deficit length t = 2.7 A;
  ramp      converged SEM-ramp resistance (sem.axisym_conductance,
            Richardson-extrapolated, ~1e-6) on a (d, L) grid;
  ramps     the same for other ramps (r_min, r_cut), a held-out test;
  box       the default ramp in the finite golden-aspect box of the paper
            (equal-area radius (d + 2p)/sqrt(pi), height
            max(L + 2p, 1.2 (d + 2p) + L), insulating side wall);
  report    deviations of Hall's offset formula and of
            sem.open_pore_theory from the solves (ramp_grid.csv, other_ramps.csv);
  figure    summary figure.

Usage
-----
    python -m sem.scripts.soft_wall_mapping all --out validation/soft_wall_mapping --jobs 4

`all` takes ~1 h on 4 cores. Lengths are in Angstrom inside the solvers;
tables report d and L in nm.
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
from sem.cylinder_mode_matching import sharp_cylinder_coefficients

T_DEF = opt.deficit_length()

RAMP_D_NM = [1.5, 2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20, 25, 30, 40, 50]
RAMP_L_NM = [1, 2, 3, 5, 10, 15, 20, 30, 40]
OTHER_RAMPS = [(1.3, 3.0), (1.3, 5.0), (1.3, 6.0), (0.0, 4.1), (2.5, 4.1)]
OTHER_DL = [(2, 2), (3, 10), (5, 1), (5, 20), (10, 3), (10, 20), (25, 20), (50, 5), (50, 40)]
SHAPE_RATIOS = [0.0, 0.1, 0.2, 0.3, opt.R_MIN / opt.R_CUT, 0.4, 0.5, 0.6, 0.7, 0.8]
BOX_D_NM = [5, 10, 15, 20, 25]
BOX_PD = [0.35, 0.5, 1.0, 2.0]


# ---------------------------------------------------------------- workers
def _ramp_job(args):
    d_nm, L_nm, r_min, r_cut = args
    t0 = time.time()
    res = pore_resistance(AxisymmetricPore.cylinder(5 * d_nm, 10 * L_nm), "ramp",
                          r_min=r_min, r_cut=r_cut, h_fine=min(0.1, (r_cut - r_min) / 20),
                          extrapolate=True)
    return dict(d_nm=d_nm, L_nm=L_nm, r_min=r_min, r_cut=r_cut,
                R_rho_per_A=res.resistance, seconds=round(time.time() - t0, 1))


def _shape_job(rho):
    # t fixed at 2.7 A; thick membrane (L/a = 4) so that k(lambda) = k_inf
    r_min, r_cut = 2 * T_DEF * rho / (1 + rho), 2 * T_DEF / (1 + rho)
    a, L = 100.0, 400.0
    res = pore_resistance(AxisymmetricPore.cylinder(a, L), "ramp", r_min=r_min, r_cut=r_cut,
                          h_fine=min(0.1, (r_cut - r_min) / 20), extrapolate=True)
    a_e = float(opt.equal_area_radius(a, r_min, r_cut))
    kappa, _ = sharp_cylinder_coefficients((L + 2 * T_DEF) / a_e)
    dl = 0.5 * (L + 2 * T_DEF + 2 * kappa[0] * a_e - res.resistance * math.pi * a_e ** 2)
    return dict(r_min=r_min, r_cut=r_cut, x=(r_cut - r_min) / T_DEF, delta_l=dl,
                C=dl * a_e ** (1 / 3) / T_DEF ** (4 / 3))


def _box_job(args):
    d_nm, pd = args
    p = pd * d_nm
    side = d_nm + 2 * p
    res = pore_resistance(AxisymmetricPore.cylinder(5 * d_nm, 200.0), "ramp", h_fine=0.1,
                          far_field="box", box_radius=10 * side / math.sqrt(math.pi),
                          box_height=10 * max(20 + 2 * p, 1.2 * side + 20), extrapolate=True)
    return dict(d_nm=d_nm, L_nm=20, p_over_d=pd, R_box_rho_per_A=res.resistance)


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
def run_ritz(out: Path, jobs: int):
    lam = np.logspace(np.log10(0.02), np.log10(20.0), 70)
    kappa, k_rim = sharp_cylinder_coefficients(lam)
    rows = [dict(L_over_a=x, kappa=k, kappa_closed_form=float(opt.end_correction_coefficient(x)),
                 k_rim=kr, k_rim_interp=float(opt.rim_strength(x)))
            for x, k, kr in zip(lam, kappa, k_rim)]
    _write(out / "sharp_cylinder.csv", rows)


def run_rimshape(out: Path, jobs: int):
    _write(out / "rim_shape.csv", _run(_shape_job, SHAPE_RATIOS, jobs))


def run_ramp(out: Path, jobs: int):
    tasks = [(d, L, opt.R_MIN, opt.R_CUT) for d, L in itertools.product(RAMP_D_NM, RAMP_L_NM)]
    _write(out / "ramp_grid.csv", _run(_ramp_job, tasks, jobs))


def run_ramps(out: Path, jobs: int):
    tasks = [(d, L, r0, r1) for (r0, r1) in OTHER_RAMPS for d, L in OTHER_DL]
    _write(out / "other_ramps.csv", _run(_ramp_job, tasks, jobs))


def run_box(out: Path, jobs: int):
    _write(out / "golden_box.csv", _run(_box_job, list(itertools.product(BOX_D_NM, BOX_PD)), jobs))


def _deviations(r):
    d, L = 10 * r["d_nm"], 10 * r["L_nm"]
    kw = dict(r_min=r["r_min"], r_cut=r["r_cut"])
    R = r["R_rho_per_A"] * 1e10                       # Ohm at sigma = 1 S/m
    t = opt.deficit_length(**kw)
    a_e = opt.equal_area_radius(d / 2, **kw)
    out = dict(dev_paper_pct=100 * (R / opt.sem_cylinder_resistance(d, L, model="paper", **kw) - 1),
               dev_no_rim_pct=100 * (R / opt.sharp_cylinder_resistance(2 * a_e, L + 2 * t) - 1),
               dev_full_pct=100 * (R / opt.sem_cylinder_resistance(d, L, **kw) - 1))
    return {k: round(float(v), 5) for k, v in out.items()}


def run_report(out: Path):
    for name in ("ramp_grid.csv", "other_ramps.csv"):
        path = out / name
        if not path.exists():
            continue
        keep = ("d_nm", "L_nm", "r_min", "r_cut", "R_rho_per_A", "seconds")
        rows = [{**{k: r[k] for k in keep}, **_deviations(r)} for r in _read(path)]
        _write(path, rows)
        for lo, label in ((2.0, "d >= 2 nm"), (0.0, "all d")):
            sel = [r for r in rows if r["d_nm"] >= lo]
            full = [abs(r["dev_full_pct"]) for r in sel]
            paper = [r["dev_paper_pct"] for r in sel]
            print(f"{name:16s} {label:9s}: full max {max(full):.4f}% rms "
                  f"{math.sqrt(np.mean(np.square(full))):.4f}%  |  Hall offset "
                  f"{min(paper):+.2f} .. {max(paper):+.2f}%")
    box = out / "golden_box.csv"
    if box.exists() and (out / "ramp_grid.csv").exists():
        inf = {(r["d_nm"], r["L_nm"]): r["R_rho_per_A"] for r in _read(out / "ramp_grid.csv")}
        rows = []
        for r in _read(box):
            R_inf = inf[(r["d_nm"], 20.0)]
            R_paper = opt.sem_cylinder_resistance(10 * r["d_nm"], 200.0, model="paper") / 1e10
            rows.append(dict(d_nm=r["d_nm"], p_over_d=r["p_over_d"],
                             R_box_rho_per_A=r["R_box_rho_per_A"],
                             box_vs_infinite_pct=round(100 * (R_inf / r["R_box_rho_per_A"] - 1), 4),
                             box_vs_paper_reference_pct=round(
                                 100 * (R_paper / r["R_box_rho_per_A"] - 1), 4)))
        _write(box, rows)
        v = [r["box_vs_infinite_pct"] for r in rows]
        print(f"golden box vs infinite reservoirs: {min(v):+.2f} .. {max(v):+.2f}%")


# ----------------------------------------------------------------- figure
LIGHT = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a",
             ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", mid="#f0efec", red="#e34948")


def make_figure(out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    c = LIGHT
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": c["ink2"], "axes.labelcolor": c["ink"],
                         "xtick.color": c["ink2"], "ytick.color": c["ink2"], "axes.linewidth": 0.6,
                         "grid.color": c["grid"], "grid.linewidth": 0.6, "legend.frameon": False})
    cyl = _read(out / "sharp_cylinder.csv")
    ramp = _read(out / "ramp_grid.csv")
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)

    # (A) kappa(L/a)
    ax = axs[0, 0]
    lam = np.logspace(-2.2, 1.4, 300)
    ax.axhline(math.pi / 4, color=c["orange"], lw=1.2)
    ax.axhline(opt.KAPPA_INF, color=c["ink2"], lw=0.8)
    ax.plot(lam, opt.end_correction_coefficient(lam), color=c["blue"], lw=1.5, label="closed form")
    ax.plot([r["L_over_a"] for r in cyl][::3], [r["kappa"] for r in cyl][::3], "o", ms=4.5,
            mfc=c["blue"], mec="white", mew=0.8, label="exact (mode matching)")
    ax.text(0.0075, math.pi / 4 + 0.0012, r"Hall: $\pi/4$", color=c["orange"], fontsize=8)
    ax.text(0.0075, opt.KAPPA_INF + 0.0012, "flanged tube: 0.82167", color=c["ink2"], fontsize=8)
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
        ax.plot(d, [r[key] for r in sel], "-o", color=col, lw=1.5, ms=4, mec="white", mew=0.6,
                label=lab)
    ax.set_xlabel("nominal pore diameter $d$ (nm)")
    ax.set_ylabel(r"$G_\mathrm{theory}/G_\mathrm{SEM} - 1$ (%)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True)
    ax.set_title("B  Theory vs converged SEM, $L$ = 20 nm", loc="left", fontsize=9)

    # (C) map of the Hall-offset error over (d, L)
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
    ax.set_title(f"C  Hall-offset error; full mapping |err| ≤ {full_max:.2f}%", loc="left",
                 fontsize=9)

    # (D) rim coefficient vs lambda, against the corner theory C (k/k_inf)^2
    ax = axs[1, 1]
    lam_e, coef = [], []
    for r in ramp:
        if r["d_nm"] < 2:
            continue
        a_e = float(opt.equal_area_radius(5 * r["d_nm"]))
        L_e0 = 10 * r["L_nm"] + 2 * T_DEF
        lam0 = L_e0 / a_e
        dl = 0.5 * (L_e0 + 2 * float(opt.end_correction_coefficient(lam0)) * a_e
                    - r["R_rho_per_A"] * math.pi * a_e ** 2)
        lam_e.append(lam0)
        coef.append(dl * a_e ** (1 / 3) / T_DEF ** (4 / 3))
    ax.plot(lam_e, coef, "o", ms=4.5, mfc=c["blue"], mec="white", mew=0.8, label="from SEM solves")
    lg = np.logspace(np.log10(min(lam_e)), np.log10(max(lam_e)), 200)
    ax.plot(lg, opt.rim_coefficient(lg), color=c["ink"], lw=1.2,
            label=r"corner theory $C\,(k/k_\infty)^2$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda = (L + 2t)/a_e$")
    ax.set_ylabel(r"rim deficit  $\Delta\ell\, a_e^{1/3}/t^{4/3}$")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_title("D  Rounded-rim correction", loc="left", fontsize=9)

    for ext in ("pdf", "png"):
        fig.savefig(out / f"soft_wall_mapping.{ext}", dpi=200)
    return out / "soft_wall_mapping.png"


STEPS = {"ritz": run_ritz, "rimshape": run_rimshape, "ramp": run_ramp, "ramps": run_ramps,
         "box": run_box}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=list(STEPS) + ["report", "figure", "all"])
    ap.add_argument("--out", type=Path, default=Path("validation/soft_wall_mapping"))
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, step in STEPS.items():
        if args.step in (name, "all"):
            t0 = time.time()
            step(args.out, args.jobs)
            print(f"{name}: {time.time() - t0:.0f} s", flush=True)
    if args.step in ("report", "all"):
        run_report(args.out)
    if args.step in ("figure", "all"):
        print("wrote", make_figure(args.out))


if __name__ == "__main__":
    main()
