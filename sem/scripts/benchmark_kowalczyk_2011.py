#!/usr/bin/env python3
"""
Benchmark PETK / SEM open-pore conductance against Kowalczyk et al. 2011.

Reference:
    Kowalczyk, Grosberg, Rabin, Dekker.
    "Modeling the conductance and DNA blockade of solid-state nanopores."
    Nanotechnology 22 (2011) 315101.

The paper measures G(d) at 1 M KCl for SiN nanopores 5–100 nm in diameter
through a 20 nm-thick membrane. Three analytical references are introduced
in the paper and reproduced here for comparison:

  Eq 1 (naive cylinder):
      G_naive  = sigma * pi * d^2 / (4 * l)

  Eq 3 (cylinder with Hall access resistance, R_access = rho/(2 d)):
      G_Hall   = sigma * [ 4 l / (pi d^2) + 1 / d ]^(-1)

  Eqs 5–7 (hourglass / hyperboloid pore, exact analytical solution):
      sin^2 alpha = (D^2 - d^2) / (l^2 + D^2 - d^2)
      R_hyp = (2 rho / (pi d)) * (sin a / (1 - cos a))
                                * arctan( sqrt(D^2 - d^2) / d )
      R_lower = R_hyp + 2 rho / (pi D)            # Hille hemisphere bound
      R_upper = R_hyp + rho / D                    # Hall planar disc bound
      G_lower = 1 / R_upper,  G_upper = 1 / R_lower
  with D ≈ d + 20 nm from TEM tomographs of e-beam-drilled SiN pores.

This script does three things:

  1. Computes the analytical references on a fine d grid.
  2. (Optional, --with-sem) Generates a SEM JSON config for each diameter,
     runs `python -m sem <config> open_pore`, parses the resulting
     *_open_pore_current.txt, and converts I → G = I / V.
  3. Writes a results CSV and a comparison plot.

Validation criteria (the right answer):
  • SEM cylindrical predictions should track the Hall curve (eq 3) computed
    with the *real* membrane thickness (20 nm), NOT the paper's empirical
    "effective" l_eff = 8.6 nm. The 8.6 nm is what you get when you force a
    cylinder model to fit experimental data of an actually-hourglass pore;
    your SEM is solving a true cylinder, so it should match the cylinder
    formula at l = 20 nm.
  • SEM double_cone predictions (with D = d + 20 nm) should fall *inside*
    the hourglass bracket (G_lower < G_SEM < G_upper). The bracket is
    parameter-free; staying inside it confirms SEM is solving the
    hyperboloid Laplace equation correctly.
  • At small d (≲ 15 nm) the access resistance is small and all three
    curves converge with eq 1 (naive cylinder).

Box / system size:
  By default the script uses a FIXED box padding (--box-padding-nm, default
  50 nm) for every diameter. This means the box-truncation bias on G is part
  of the result and varies systematically with d (the bias is ~0% at small
  d, growing to ~5–30% at d = 100 nm depending on padding). Use this default
  when you want to *investigate* how system size affects conductance —
  comparing your G(d) to the analytical curves at one fixed box gives a
  visible measure of the truncation, and using --padding-sweep produces
  G(d) curves at multiple box sizes overlaid for a direct comparison.

  If instead you want diameter-adaptive padding (no truncation bias visible
  in the plot), pass --box-padding-factor 3.0 (~5% bias) or 5.0 (~1% bias).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Paper constants (Kowalczyk 2011 §1, Materials and methods)
# ---------------------------------------------------------------------------
SIGMA_S_PER_M_DEFAULT = 10.5         # 1 M KCl at 23 °C (CRC Handbook 2010)
MEMBRANE_THICKNESS_NM_DEFAULT = 20.0 # SiN membrane thickness used in paper
EFFECTIVE_THICKNESS_NM = 8.6         # Empirical best-fit l for cyl+Hall (Fig 1)
DEFAULT_VOLTAGE_MV = 200.0           # ±0.2 V applied in the paper

DEFAULT_DIAMETERS_NM = [5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 65, 80, 100]


# ---------------------------------------------------------------------------
# Analytical references
# ---------------------------------------------------------------------------
def naive_cylindrical_nS(d_nm: float, l_nm: float, sigma: float) -> float:
    """G = sigma * pi d^2 / (4 l). Eq 1, Kowalczyk 2011. Returns nS."""
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    return sigma * math.pi * d ** 2 / (4.0 * l) * 1e9


def hall_cylindrical_nS(d_nm: float, l_nm: float, sigma: float) -> float:
    """G = sigma * [4 l / (pi d^2) + 1/d]^-1. Eq 3, Kowalczyk 2011. Returns nS.
    Uses the Hall (planar disc) access-resistance term R_access = rho / (2 d).
    """
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    inv_G = (4.0 * l) / (math.pi * d ** 2) + 1.0 / d
    return sigma / inv_G * 1e9


def hourglass_bracket_nS(
    d_nm: float, l_nm: float, D_nm: float, sigma: float
) -> Tuple[float, float]:
    """Hourglass / hyperboloid bracket from eqs 5–7 of Kowalczyk 2011.

    Lower G bound uses Hille hemispherical access (rho / (pi D)) per side;
    upper G bound uses Hall planar-disc access (rho / (2 D)) per side.
    """
    d = d_nm * 1e-9
    l = l_nm * 1e-9
    D = D_nm * 1e-9
    rho = 1.0 / sigma

    if D <= d:
        # Degenerate: cylinder. Return the Hall cylinder both sides.
        G = hall_cylindrical_nS(d_nm, l_nm, sigma)
        return G, G

    sin2_a = (D ** 2 - d ** 2) / (l ** 2 + D ** 2 - d ** 2)
    sin_a = math.sqrt(sin2_a)
    cos_a = math.sqrt(1.0 - sin2_a)
    R_hyp = (2.0 * rho / (math.pi * d)) \
            * (sin_a / (1.0 - cos_a)) \
            * math.atan(math.sqrt(D ** 2 - d ** 2) / d)
    R_lower_G = R_hyp + rho / D                    # Hall on each face → smaller G
    R_upper_G = R_hyp + 2.0 * rho / (math.pi * D)  # Hille on each face → larger G
    G_lower_nS = 1.0 / R_lower_G * 1e9
    G_upper_nS = 1.0 / R_upper_G * 1e9
    return G_lower_nS, G_upper_nS


def kowalczyk_taper_D_nm(d_nm: float, taper_nm: float = 20.0) -> float:
    """Empirical D ≈ d + 20 nm for e-beam-drilled SiN pores (Kowalczyk 2011).
    The actual taper depends on TEM beam conditions; 20 nm is the typical value.
    """
    return d_nm + taper_nm


# ---------------------------------------------------------------------------
# SEM config writer + runner
# ---------------------------------------------------------------------------
def make_sem_config(
    *,
    d_nm: float,
    l_nm: float,
    voltage_mv: float,
    sigma_s_per_m: float,
    geometry: str,
    grid_resolution_A: float,
    box_padding_nm: float,
    output_prefix: str,
    box_padding_factor: float = 3.0,
    D_nm: Optional[float] = None,
    fixed_box_nm: Optional[float] = None,
    fixed_z_box_nm: Optional[float] = None,
    mesh_engine: Optional[str] = None,
    gmsh_fine_size_A: Optional[float] = None,
    gmsh_coarse_size_A: Optional[float] = None,
    gmsh_fine_box_A: Optional[List[float]] = None,
    gmsh_fine_margin_nm: Optional[float] = None,
) -> dict:
    """Construct a SEM JSON config for an open-pore run.

    Two box-sizing modes:

    * ``fixed_box_nm`` is None (default) — padding-driven: the box is sized
      to ``pore_widest + 2 * effective_padding`` in XY, ``L + 2 *
      effective_padding`` in Z, where ``effective_padding = max(
      box_padding_nm, box_padding_factor * pore_widest / 2)``. With factor=0
      this reduces to a fixed padding around the pore (so total box scales
      with d); with factor>0 the padding scales with d to keep truncation
      bias below a target.

    * ``fixed_box_nm`` is a number — the total CUBIC box edge in nm,
      overriding all padding. Every dimension is exactly that size,
      regardless of pore diameter. Raises ValueError if the pore doesn't
      fit (radial: pore_widest/2 must be < edge/2; vertical: L/2 must be <
      edge/2). Use this when you want the box-truncation bias to be visible
      and to vary with d in a controlled way.
    """
    radius_A = (d_nm / 2.0) * 10.0
    thickness_A = l_nm * 10.0
    pore_widest_nm = max(d_nm, D_nm or d_nm)

    if fixed_box_nm is not None:
        half_edge_nm = fixed_box_nm / 2.0
        # Validate that the geometry fits inside the requested box.
        if pore_widest_nm / 2.0 >= half_edge_nm:
            raise ValueError(
                f"Pore radial extent {pore_widest_nm/2:.2f} nm does not fit in "
                f"fixed cubic box of edge {fixed_box_nm:.2f} nm. "
                f"Either increase --fixed-box-nm or reduce --taper-nm "
                f"(double_cone) / pore diameter."
            )
        if l_nm / 2.0 >= half_edge_nm:
            raise ValueError(
                f"Membrane half-thickness {l_nm/2:.2f} nm does not fit in "
                f"fixed cubic box of edge {fixed_box_nm:.2f} nm. "
                f"Increase --fixed-box-nm or reduce --membrane-thickness-nm."
            )
        half_box_A = half_edge_nm * 10.0
        half_box_z_A = half_edge_nm * 10.0
    else:
        effective_pad_nm = max(box_padding_nm, box_padding_factor * pore_widest_nm / 2.0)
        half_box_A = (pore_widest_nm / 2.0 + effective_pad_nm) * 10.0
        # Z follows the padding rule by default. If --fixed-z-box-nm is set,
        # override Z only and leave XY untouched — useful for isolating the
        # lateral access truncation (vary XY via padding) from the vertical
        # access truncation (held constant).
        if fixed_z_box_nm is not None:
            if l_nm >= fixed_z_box_nm:
                raise ValueError(
                    f"Membrane thickness {l_nm:.2f} nm does not fit in "
                    f"fixed Z box of {fixed_z_box_nm:.2f} nm. "
                    f"Increase --fixed-z-box-nm or reduce --membrane-thickness-nm."
                )
            half_box_z_A = (fixed_z_box_nm / 2.0) * 10.0
        else:
            half_box_z_A = (l_nm / 2.0 + effective_pad_nm) * 10.0
    geom = {"membrane_thickness": thickness_A}
    if geometry == "cylindrical":
        geom["pore_type"] = "cylindrical"
        geom["pore_radius"] = radius_A
        geom["corner_radius"] = 0.0
    elif geometry == "double_cone":
        if D_nm is None:
            raise ValueError("double_cone requires D_nm")
        geom["pore_type"] = "double_cone"
        geom["pore_radius"] = radius_A
        geom["outer_radius"] = (D_nm / 2.0) * 10.0
    else:
        raise ValueError(f"Unsupported geometry: {geometry}")
    simulation = {
        "voltage": voltage_mv,
        "bulk_conductivity": sigma_s_per_m,
        "grid_resolution": grid_resolution_A,
        "use_vdw_radii": False,
        "use_pdb2pqr": False,
        "force_field": "CHARMM",
        "default_radius": 1.5,
        "membrane_conductivity": 1e-7,
    }
    # Optional graded-mesh fields. When all three (mesh_engine, fine_size,
    # coarse_size) are present, SEM will use a fine cell size inside the
    # gmsh_fine_box around the pore and a coarse cell size outside. When
    # omitted (or only some are set), SEM falls back to its own default
    # (uniform mesh keyed off grid_resolution).
    if mesh_engine is not None:
        simulation["mesh_engine"] = mesh_engine
    if gmsh_fine_size_A is not None:
        simulation["gmsh_fine_size"] = gmsh_fine_size_A
    if gmsh_coarse_size_A is not None:
        simulation["gmsh_coarse_size"] = gmsh_coarse_size_A

    # Resolve the gmsh fine box. Three precedences:
    #   1. Explicit --gmsh-fine-box-A wins if given.
    #   2. Otherwise, --gmsh-fine-margin-nm auto-computes a per-diameter fine
    #      box of (pore_widest + 2*margin) × (pore_widest + 2*margin) ×
    #      (l + 2*margin), guaranteeing the fine region always covers the
    #      pore wall plus a buffer regardless of d.
    #   3. Otherwise no fine box is emitted (SEM falls back to its default).
    # If a caller supplies both, that's a config error — caught upstream.
    fine_box_A: Optional[List[float]] = None
    if gmsh_fine_box_A is not None:
        fine_box_A = list(gmsh_fine_box_A)
    elif gmsh_fine_margin_nm is not None:
        # widest in nm → Å, plus 2 × margin
        fine_xy_A = (pore_widest_nm + 2.0 * gmsh_fine_margin_nm) * 10.0
        fine_z_A  = (l_nm           + 2.0 * gmsh_fine_margin_nm) * 10.0
        fine_box_A = [fine_xy_A, fine_xy_A, fine_z_A]

    # Sanity: the fine box should not exceed the simulation box. If it
    # does, clamp to the simulation box edge (avoids gmsh confusion when the
    # "fine region" leaks past the domain boundary).
    if fine_box_A is not None:
        sim_box_A = [2.0 * half_box_A, 2.0 * half_box_A, 2.0 * half_box_z_A]
        fine_box_A = [min(f, b) for f, b in zip(fine_box_A, sim_box_A)]
        simulation["gmsh_fine_box"] = fine_box_A

    return {
        "metadata": {
            "generated_by": "benchmark_kowalczyk_2011.py",
            "calculation_mode": "open_pore",
            "version": "1.0",
        },
        "input": {"moving_pdb": ""},  # Not required for open_pore
        "pore_geometry": geom,
        "simulation": simulation,
        "movement": {"z_start": 0.0, "z_end": 0.0, "z_step": 1.0},
        "output": {"output_prefix": output_prefix, "preview_frames": 0},
        "box_dimensions": {
            "x": [-half_box_A, half_box_A],
            "y": [-half_box_A, half_box_A],
            "z": [-half_box_z_A, half_box_z_A],
        },
    }


_OPEN_PORE_LINE_RE = re.compile(r"^\s*([0-9.+\-eE]+)\s*$")


def parse_open_pore_result(path: Path) -> Optional[float]:
    """Read the *_open_pore_current.txt file and return I (nA) or None."""
    if not path.exists():
        return None
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = _OPEN_PORE_LINE_RE.match(line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    continue
    return None


def run_sem_open_pore(
    config_path: Path, output_dir: Path, sem_invoke: List[str], log_path: Path
) -> Optional[float]:
    """Run `python -m sem <config> open_pore` and return I (nA) or None.

    The subprocess is launched with ``cwd=output_dir`` so SEM's relative-path
    conventions for output files line up. Because of that, the config path
    handed to the subprocess MUST be absolute — otherwise SEM tries to
    resolve it relative to the new cwd and fails.
    """
    abs_config_path = config_path.resolve()
    abs_output_dir = output_dir.resolve()
    cmd = list(sem_invoke) + [str(abs_config_path), "open_pore"]
    try:
        with log_path.open("w") as logfh:
            logfh.write(f"# cmd: {' '.join(cmd)}\n")
            logfh.write(f"# cwd: {abs_output_dir}\n")
            logfh.flush()
            subprocess.run(
                cmd, cwd=str(abs_output_dir),
                stdout=logfh, stderr=subprocess.STDOUT, check=True,
            )
    except subprocess.CalledProcessError:
        return None

    # Find the open_pore_current.txt that SEM wrote (it's relative to the
    # subprocess cwd, which is output_dir).
    output_prefix = json.loads(abs_config_path.read_text())["output"]["output_prefix"]
    result_path = abs_output_dir / f"{output_prefix}_open_pore_current.txt"
    return parse_open_pore_result(result_path)


# ---------------------------------------------------------------------------
# Sweep orchestration
# ---------------------------------------------------------------------------
@dataclass
class SweepResult:
    diameter_nm: float
    geometry: str
    G_naive_nS: float
    G_hall_l_real_nS: float
    G_hall_l_eff_nS: float
    G_hourglass_lower_nS: float
    G_hourglass_upper_nS: float
    G_sem_nS: Optional[float]


def run_sweep(args: argparse.Namespace) -> List[SweepResult]:
    diameters = sorted({float(d) for d in args.diameters})
    rows: List[SweepResult] = []
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Pre-flight: --gmsh-fine-box-A and --gmsh-fine-margin-nm are mutually
    # exclusive. Both specify the fine-region size; allowing both would be
    # ambiguous.
    if args.gmsh_fine_box_A is not None and args.gmsh_fine_margin_nm is not None:
        raise ValueError(
            "Specify either --gmsh-fine-box-A (explicit edges) or "
            "--gmsh-fine-margin-nm (auto-scaled), not both."
        )

    # Pre-flight: validate --fixed-z-box-nm vs membrane thickness.
    if args.fixed_z_box_nm is not None and args.fixed_box_nm is None:
        if args.membrane_thickness_nm >= args.fixed_z_box_nm:
            print(
                f"\nERROR: --fixed-z-box-nm {args.fixed_z_box_nm:.1f} nm is not "
                f"larger than --membrane-thickness-nm "
                f"{args.membrane_thickness_nm:.1f} nm. The membrane would "
                f"fill or exceed the box vertically, leaving no bulk solvent.\n"
            )
            raise ValueError("fixed_z_box_nm must be > membrane_thickness_nm")

    # Pre-flight: when --fixed-box-nm is set, warn up front about (d, geometry)
    # combinations that won't fit. With --with-sem the ValueError will fire
    # later anyway; without --with-sem this is the only feedback the user gets.
    if args.fixed_box_nm is not None:
        half_edge_nm = args.fixed_box_nm / 2.0
        bad: List[str] = []
        if args.membrane_thickness_nm / 2.0 >= half_edge_nm:
            bad.append(
                f"  membrane half-thickness {args.membrane_thickness_nm/2:.2f} nm "
                f">= half-edge {half_edge_nm:.2f} nm"
            )
        for d_nm in diameters:
            D_nm = kowalczyk_taper_D_nm(d_nm, args.taper_nm)
            for geometry in args.geometries:
                widest = D_nm if geometry == "double_cone" else d_nm
                if widest / 2.0 >= half_edge_nm:
                    bad.append(
                        f"  d={d_nm:.1f} nm, {geometry}: widest radius "
                        f"{widest/2:.2f} nm >= half-edge {half_edge_nm:.2f} nm"
                    )
        if bad:
            print()
            print(f"WARNING: {len(bad)} (d, geometry) configuration(s) will NOT "
                  f"fit in the {args.fixed_box_nm:.1f} nm cubic box:")
            for line in bad:
                print(line)
            print(f"  Analytical references will still be computed; SEM runs "
                  f"for these would error.")
            print()

    sem_invoke: Optional[List[str]] = None
    if args.with_sem:
        if args.conda_env:
            sem_invoke = ["conda", "run", "-n", args.conda_env, "python", "-m", "sem"]
        else:
            sem_invoke = [args.python_exec, "-m", "sem"]

    for d_nm in diameters:
        D_nm = kowalczyk_taper_D_nm(d_nm, args.taper_nm)
        G_naive = naive_cylindrical_nS(d_nm, args.membrane_thickness_nm, args.sigma)
        G_hall_real = hall_cylindrical_nS(d_nm, args.membrane_thickness_nm, args.sigma)
        G_hall_eff = hall_cylindrical_nS(d_nm, EFFECTIVE_THICKNESS_NM, args.sigma)
        G_hg_lower, G_hg_upper = hourglass_bracket_nS(
            d_nm, args.membrane_thickness_nm, D_nm, args.sigma
        )

        for geometry in args.geometries:
            G_sem: Optional[float] = None
            if args.with_sem:
                run_dir = args.output_dir / f"sem_{geometry}_d{int(round(d_nm))}nm"
                run_dir.mkdir(parents=True, exist_ok=True)
                cfg = make_sem_config(
                    d_nm=d_nm,
                    l_nm=args.membrane_thickness_nm,
                    voltage_mv=args.voltage_mv,
                    sigma_s_per_m=args.sigma,
                    geometry=geometry,
                    grid_resolution_A=args.grid_resolution_A,
                    box_padding_nm=args.box_padding_nm,
                    box_padding_factor=args.box_padding_factor,
                    output_prefix=f"open_pore_{geometry}_d{int(round(d_nm))}",
                    D_nm=D_nm if geometry == "double_cone" else None,
                    fixed_box_nm=args.fixed_box_nm,
                    fixed_z_box_nm=args.fixed_z_box_nm,
                    mesh_engine=args.mesh_engine,
                    gmsh_fine_size_A=args.gmsh_fine_size_A,
                    gmsh_coarse_size_A=args.gmsh_coarse_size_A,
                    gmsh_fine_box_A=args.gmsh_fine_box_A,
                    gmsh_fine_margin_nm=args.gmsh_fine_margin_nm,
                )
                # Surface the effective box for transparency.
                bx = cfg["box_dimensions"]
                xy_box_nm = (bx["x"][1] - bx["x"][0]) / 10.0
                z_box_nm  = (bx["z"][1] - bx["z"][0]) / 10.0
                print(f"  [{geometry:>11}] d={d_nm:>5.1f} nm  "
                      f"box = {xy_box_nm:.0f} × {xy_box_nm:.0f} × {z_box_nm:.0f} nm³")
                cfg_path = run_dir / "config.json"
                cfg_path.write_text(json.dumps(cfg, indent=2))
                I_nA = run_sem_open_pore(
                    cfg_path, run_dir, sem_invoke,
                    run_dir / "sem_run.log",
                )
                if I_nA is not None and args.voltage_mv > 0:
                    # G [S] = I [A] / V [V]; with I in nA and V in mV,
                    # G [nS] = (I [nA] × 1e-9) / (V [mV] × 1e-3) × 1e9
                    #       = I [nA] / V [mV] × 1e3
                    G_sem = I_nA / args.voltage_mv * 1e3

            rows.append(SweepResult(
                diameter_nm=d_nm,
                geometry=geometry,
                G_naive_nS=G_naive,
                G_hall_l_real_nS=G_hall_real,
                G_hall_l_eff_nS=G_hall_eff,
                G_hourglass_lower_nS=G_hg_lower,
                G_hourglass_upper_nS=G_hg_upper,
                G_sem_nS=G_sem,
            ))

    return rows


def write_csv(rows: List[SweepResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "diameter_nm", "geometry",
            "G_naive_eq1_nS", "G_hall_l_real_eq3_nS",
            "G_hall_l_eff_8.6nm_nS",
            "G_hourglass_lower_eq7_nS", "G_hourglass_upper_eq7_nS",
            "G_sem_nS",
        ])
        for r in rows:
            w.writerow([
                f"{r.diameter_nm:.2f}", r.geometry,
                f"{r.G_naive_nS:.4f}", f"{r.G_hall_l_real_nS:.4f}",
                f"{r.G_hall_l_eff_nS:.4f}",
                f"{r.G_hourglass_lower_nS:.4f}", f"{r.G_hourglass_upper_nS:.4f}",
                "" if r.G_sem_nS is None else f"{r.G_sem_nS:.4f}",
            ])


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_results(rows: List[SweepResult], args: argparse.Namespace, png_path: Path) -> None:
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Dense grid for the analytical curves
    d_dense_nm = np.linspace(min(args.diameters), max(args.diameters), 200)
    naive = [naive_cylindrical_nS(d, args.membrane_thickness_nm, args.sigma) for d in d_dense_nm]
    hall_real = [hall_cylindrical_nS(d, args.membrane_thickness_nm, args.sigma) for d in d_dense_nm]
    hall_eff = [hall_cylindrical_nS(d, EFFECTIVE_THICKNESS_NM, args.sigma) for d in d_dense_nm]
    hg_lower = [hourglass_bracket_nS(
        d, args.membrane_thickness_nm, kowalczyk_taper_D_nm(d, args.taper_nm),
        args.sigma)[0] for d in d_dense_nm]
    hg_upper = [hourglass_bracket_nS(
        d, args.membrane_thickness_nm, kowalczyk_taper_D_nm(d, args.taper_nm),
        args.sigma)[1] for d in d_dense_nm]

    fig, ax = plt.subplots(figsize=(8.0, 6.0))

    # Reference curves
    ax.plot(d_dense_nm, naive, "--", color="gray", linewidth=1.0,
            label=f"Eq 1 naive cyl, l = {args.membrane_thickness_nm:.0f} nm")
    ax.plot(d_dense_nm, hall_real, "-", color="C0", linewidth=1.6,
            label=f"Eq 3 cyl + Hall access, l = {args.membrane_thickness_nm:.0f} nm")
    ax.plot(d_dense_nm, hall_eff, ":", color="C0", linewidth=1.0,
            label=f"Eq 3 cyl + Hall access, l_eff = {EFFECTIVE_THICKNESS_NM:.1f} nm  (paper's fit)")
    ax.fill_between(d_dense_nm, hg_lower, hg_upper, color="C2", alpha=0.20,
                    label=f"Eq 7 hourglass bracket, D = d + {args.taper_nm:.0f} nm")
    ax.plot(d_dense_nm, hg_lower, "-", color="C2", linewidth=1.0, alpha=0.6)
    ax.plot(d_dense_nm, hg_upper, "-", color="C2", linewidth=1.0, alpha=0.6)

    # SEM points
    cyl_pts = [(r.diameter_nm, r.G_sem_nS) for r in rows
               if r.geometry == "cylindrical" and r.G_sem_nS is not None]
    dc_pts = [(r.diameter_nm, r.G_sem_nS) for r in rows
              if r.geometry == "double_cone" and r.G_sem_nS is not None]
    if cyl_pts:
        xs, ys = zip(*cyl_pts)
        ax.plot(xs, ys, "o", color="C3", markersize=8,
                label="SEM cylindrical (this work)", zorder=5)
    if dc_pts:
        xs, ys = zip(*dc_pts)
        ax.plot(xs, ys, "s", color="C1", markersize=8,
                label=f"SEM double_cone, D = d + {args.taper_nm:.0f} nm (this work)",
                zorder=5)

    ax.set_xlabel("Pore diameter d (nm)", fontsize=12, weight="bold")
    ax.set_ylabel("Open-pore conductance G (nS)", fontsize=12, weight="bold")
    ax.set_title(
        f"PETK / SEM vs Kowalczyk 2011 — open-pore G(d), "
        f"σ = {args.sigma:.1f} S/m, V = {args.voltage_mv:.0f} mV",
        fontsize=12, weight="bold",
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    ax.set_xlim(0, max(args.diameters) * 1.02)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    if args.show:
        plt.show()
    print(f"Saved comparison plot: {png_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--diameters", type=float, nargs="+",
                   default=DEFAULT_DIAMETERS_NM,
                   help="Pore diameters in nm.")
    p.add_argument("--membrane-thickness-nm", type=float,
                   default=MEMBRANE_THICKNESS_NM_DEFAULT,
                   help="Real SiN membrane thickness (nm). Paper used 20.")
    p.add_argument("--sigma", type=float, default=SIGMA_S_PER_M_DEFAULT,
                   help="Bulk conductivity (S/m). 1 M KCl @ 23 °C ≈ 10.5.")
    p.add_argument("--voltage-mv", type=float, default=DEFAULT_VOLTAGE_MV,
                   help="Applied voltage in mV. Paper used ±200 mV.")
    p.add_argument("--taper-nm", type=float, default=20.0,
                   help="Empirical D − d (nm) for the hourglass bracket.")
    p.add_argument("--geometries", nargs="+",
                   choices=["cylindrical", "double_cone"],
                   default=["cylindrical", "double_cone"],
                   help="Which SEM geometries to run.")
    p.add_argument("--with-sem", action="store_true",
                   help="Actually run SEM open_pore for each diameter.")
    p.add_argument("--conda-env", default="",
                   help="Conda env to invoke (e.g. sem-dolfinx). "
                        "If empty, uses --python-exec directly.")
    p.add_argument("--python-exec", default=sys.executable,
                   help="Python executable to invoke as `<exe> -m sem`.")
    p.add_argument("--grid-resolution-A", type=float, default=2.0,
                   help="SEM uniform-mesh grid resolution (Å). Used when no "
                        "graded-mesh flags are set, or as the conductivity-grid "
                        "spacing when they are.")
    p.add_argument("--mesh-engine", default=None, choices=["gmsh"],
                   help="Mesh engine. Set to 'gmsh' to enable graded meshing "
                        "via the --gmsh-* flags below. Omit for SEM's default "
                        "uniform mesh keyed off --grid-resolution-A.")
    p.add_argument("--gmsh-fine-size-A", type=float, default=None,
                   help="Cell size (Å) inside --gmsh-fine-box-A. Smaller than "
                        "--gmsh-coarse-size-A. Requires --mesh-engine gmsh.")
    p.add_argument("--gmsh-coarse-size-A", type=float, default=None,
                   help="Cell size (Å) outside --gmsh-fine-box-A in the bulk. "
                        "Larger than --gmsh-fine-size-A. Requires "
                        "--mesh-engine gmsh.")
    p.add_argument("--gmsh-fine-box-A", type=float, nargs=3, default=None,
                   metavar=("LX_A", "LY_A", "LZ_A"),
                   help="Three numbers (Å): explicit edge lengths of the "
                        "centered fine box around the pore. Example: "
                        "--gmsh-fine-box-A 100 100 100 for a 10×10×10 nm "
                        "fine region. Use this when you want the same fine "
                        "box at every diameter. Mutually exclusive with "
                        "--gmsh-fine-margin-nm.")
    p.add_argument("--gmsh-fine-margin-nm", type=float, default=None,
                   help="Auto-scale the fine box per diameter: it becomes "
                        "(pore_widest + 2*margin) × (pore_widest + 2*margin) "
                        "× (l + 2*margin). Guarantees the fine mesh always "
                        "covers the pore wall plus a buffer of <margin> nm, "
                        "regardless of d. Recommended for diameter sweeps. "
                        "Example: --gmsh-fine-margin-nm 5 (5 nm buffer). "
                        "Mutually exclusive with --gmsh-fine-box-A.")
    p.add_argument("--box-padding-nm", type=float, default=50.0,
                   help="Padding from pore edge to box boundary (nm). "
                        "Default 50 nm matches PETK's GUI auto-box default. "
                        "With --box-padding-factor 0 (the default), this is the "
                        "actual padding used for every diameter; with "
                        "--box-padding-factor > 0 it's a floor and the actual "
                        "padding becomes max(this, factor * d / 2).")
    p.add_argument("--box-padding-factor", type=float, default=0.0,
                   help="OPTIONAL diameter-scaled padding multiplier. Default 0 "
                        "means use the fixed --box-padding-nm for every "
                        "diameter — recommended when you want the box-truncation "
                        "effect to be a visible part of the result. Set to 3.0 "
                        "to push truncation bias on G below ~5% at the cost of "
                        "a much larger box at large d, or 5.0 for ~1%.")
    p.add_argument("--padding-sweep", type=float, nargs="+", default=None,
                   help="Optional list of paddings (nm). If given, the full "
                        "d-sweep runs once for each padding value, producing an "
                        "extra plot G(d) per box size — the explicit "
                        "investigation of how system size affects conductance. "
                        "Overrides --box-padding-nm. Example: "
                        "--padding-sweep 20 50 100 200")
    p.add_argument("--fixed-box-nm", type=float, default=None,
                   help="Total CUBIC box edge in nm. Overrides padding logic: "
                        "every diameter uses exactly this box. The script "
                        "errors out if the pore doesn't fit. Example: "
                        "--fixed-box-nm 30 with --diameters 5 10 15 20 25 "
                        "for a 30 nm cubic box across a 5–25 nm cylindrical "
                        "diameter sweep.")
    p.add_argument("--fixed-z-box-nm", type=float, default=None,
                   help="Pin the Z dimension of the box to this value (nm), "
                        "while XY still follows the padding logic. Useful "
                        "for isolating the lateral access truncation from "
                        "the vertical one — combine with --padding-sweep so "
                        "only the XY box grows. Ignored if --fixed-box-nm "
                        "is also set (cubic mode wins). Example: "
                        "--fixed-z-box-nm 30 --padding-sweep 5 10 20 50 .")
    p.add_argument("--output-dir", type=Path,
                   default=Path("benchmark_kowalczyk_2011"),
                   help="Where to write configs, SEM runs, CSV and PNG.")
    p.add_argument("--show", action="store_true",
                   help="Display the plot interactively (default: save only).")
    return p.parse_args(argv)


def plot_padding_sweep(
    rows_by_padding: dict, args: argparse.Namespace, png_path: Path
) -> None:
    """For each fixed box size, plot G(d) overlaid with the analytical Hall
    curve so the user can read off the box-truncation bias visually.

    Same axes as plot_results, but each padding value is a different colour
    series of SEM markers. Hourglass bracket and Hall curves drawn once.
    """
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    d_dense = np.linspace(min(args.diameters), max(args.diameters), 200)
    hall_real = [hall_cylindrical_nS(d, args.membrane_thickness_nm, args.sigma)
                 for d in d_dense]
    hg_lower = [hourglass_bracket_nS(
        d, args.membrane_thickness_nm,
        kowalczyk_taper_D_nm(d, args.taper_nm), args.sigma)[0]
        for d in d_dense]
    hg_upper = [hourglass_bracket_nS(
        d, args.membrane_thickness_nm,
        kowalczyk_taper_D_nm(d, args.taper_nm), args.sigma)[1]
        for d in d_dense]

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    ax.plot(d_dense, hall_real, "-", color="C0", linewidth=1.5, alpha=0.7,
            label=f"Eq 3 cyl + Hall, l = {args.membrane_thickness_nm:.0f} nm")
    ax.fill_between(d_dense, hg_lower, hg_upper, color="C2", alpha=0.18,
                    label=f"Eq 7 hourglass bracket, D = d + {args.taper_nm:.0f} nm")

    cmap = plt.get_cmap("plasma")
    paddings = sorted(rows_by_padding.keys())
    for k, pad in enumerate(paddings):
        rows = rows_by_padding[pad]
        color = cmap(k / max(1, len(paddings) - 1))
        for geometry, marker in (("cylindrical", "o"), ("double_cone", "s")):
            pts = [(r.diameter_nm, r.G_sem_nS) for r in rows
                   if r.geometry == geometry and r.G_sem_nS is not None]
            if not pts:
                continue
            xs, ys = zip(*pts)
            label = (f"SEM {geometry}, padding = {pad:.0f} nm"
                     if geometry == "cylindrical" else None)
            ax.plot(xs, ys, marker=marker, linestyle="-", color=color,
                    markersize=7, linewidth=1.0, alpha=0.85, label=label)

    ax.set_xlabel("Pore diameter d (nm)", fontsize=12, weight="bold")
    ax.set_ylabel("Open-pore conductance G (nS)", fontsize=12, weight="bold")
    ax.set_title("System-size effect on SEM conductance — padding sweep",
                 fontsize=12, weight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    ax.set_xlim(0, max(args.diameters) * 1.02)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    if args.show:
        plt.show()
    print(f"Saved padding-sweep plot: {png_path}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    print("=" * 72)
    print(" Kowalczyk 2011 — open-pore conductance benchmark ".center(72, "="))
    print("=" * 72)
    print(f"Diameters [nm]      : {sorted(args.diameters)}")
    print(f"Membrane l [nm]     : {args.membrane_thickness_nm}")
    print(f"Conductivity [S/m]  : {args.sigma}")
    print(f"Voltage [mV]        : {args.voltage_mv}")
    print(f"Hourglass D − d [nm]: {args.taper_nm}")
    print(f"With SEM            : {args.with_sem}")
    print(f"Geometries          : {args.geometries}")
    print(f"Output directory    : {args.output_dir.resolve()}")
    if args.fixed_box_nm is not None:
        print(f"Fixed cubic box [nm]: {args.fixed_box_nm} "
              f"(overrides padding)")
    elif args.padding_sweep:
        print(f"Padding sweep [nm]  : {sorted(args.padding_sweep)}")
    else:
        print(f"Box padding [nm]    : {args.box_padding_nm} "
              f"(factor {args.box_padding_factor})")
    if args.fixed_z_box_nm is not None and args.fixed_box_nm is None:
        print(f"Fixed Z box [nm]    : {args.fixed_z_box_nm} "
              f"(XY still follows padding)")
    if args.mesh_engine == "gmsh" and (
        args.gmsh_fine_size_A or args.gmsh_coarse_size_A
        or args.gmsh_fine_box_A or args.gmsh_fine_margin_nm is not None
    ):
        print(f"Mesh mode           : graded (gmsh)")
        if args.gmsh_fine_size_A is not None:
            print(f"  fine cell size [Å]: {args.gmsh_fine_size_A}")
        if args.gmsh_coarse_size_A is not None:
            print(f"  coarse cell size [Å]: {args.gmsh_coarse_size_A}")
        if args.gmsh_fine_box_A is not None:
            print(f"  fine box edges [Å]: {args.gmsh_fine_box_A}  (explicit)")
        if args.gmsh_fine_margin_nm is not None:
            print(f"  fine box margin [nm]: {args.gmsh_fine_margin_nm}  "
                  f"(auto-scaled per diameter)")
    else:
        print(f"Mesh mode           : uniform "
              f"(grid_resolution = {args.grid_resolution_A} Å)")
    print("-" * 72)

    if args.padding_sweep:
        rows_by_padding: dict = {}
        for pad in sorted(args.padding_sweep):
            print(f"\n=== Padding = {pad:.1f} nm ===")
            sub_args = argparse.Namespace(**vars(args))
            sub_args.box_padding_nm = pad
            sub_args.box_padding_factor = 0.0
            sub_args.output_dir = args.output_dir / f"pad_{int(round(pad))}nm"
            sub_rows = run_sweep(sub_args)
            rows_by_padding[pad] = sub_rows
            sub_csv = sub_args.output_dir / f"benchmark_pad{int(round(pad))}nm.csv"
            write_csv(sub_rows, sub_csv)
            print(f"Wrote {sub_csv}")

        # Combined CSV across all paddings
        combined_csv = args.output_dir / "benchmark_padding_sweep.csv"
        combined_csv.parent.mkdir(parents=True, exist_ok=True)
        with combined_csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([
                "padding_nm", "diameter_nm", "geometry",
                "G_naive_eq1_nS", "G_hall_l_real_eq3_nS",
                "G_hall_l_eff_8.6nm_nS",
                "G_hourglass_lower_eq7_nS", "G_hourglass_upper_eq7_nS",
                "G_sem_nS",
            ])
            for pad, rows in rows_by_padding.items():
                for r in rows:
                    w.writerow([
                        f"{pad:.1f}", f"{r.diameter_nm:.2f}", r.geometry,
                        f"{r.G_naive_nS:.4f}", f"{r.G_hall_l_real_nS:.4f}",
                        f"{r.G_hall_l_eff_nS:.4f}",
                        f"{r.G_hourglass_lower_nS:.4f}",
                        f"{r.G_hourglass_upper_nS:.4f}",
                        "" if r.G_sem_nS is None else f"{r.G_sem_nS:.4f}",
                    ])
        print(f"Wrote combined {combined_csv}")
        plot_padding_sweep(
            rows_by_padding, args,
            args.output_dir / "benchmark_padding_sweep.png",
        )
        return 0

    # Single-padding mode
    rows = run_sweep(args)
    csv_path = args.output_dir / "benchmark_kowalczyk_2011.csv"
    png_path = args.output_dir / "benchmark_kowalczyk_2011.png"
    write_csv(rows, csv_path)
    print(f"Wrote {csv_path}")
    plot_results(rows, args, png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
