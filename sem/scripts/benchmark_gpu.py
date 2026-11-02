#!/usr/bin/env python3
"""
Benchmark the GPU backend against the DOLFINx (FEM) backend on one SEM configuration.

Runs the same workload end to end through the CLI with `--backend fem` and `--backend gpu`
and compares the results file by file:

  open_pore      `python -m sem open_pore CONFIG`
  run            `python -m sem run CONFIG`             (the z scan of the config)
  rotation_scan  `python -m sem rotation_scan CONFIG . ANGLES --mode run`
                 (N random orientations at z = 0)

Default workload: petk/Demo/config.json (the 1AOI nucleosome through a 10 nm cylindrical pore,
30 nm box, 5 A grid, 23 z positions). The script writes its own copy of the config into the
output folder, with VdW radii instead of pdb2pqr (no external tools needed), an absolute
analyte path, and the optional --grid override. Your config is not changed.

Both backends must import: run it from an environment with DOLFINx and torch, or give each
backend its own interpreter with --fem-python / --gpu-python.

Example:
  python -m sem.scripts.benchmark_gpu --out gpu_benchmark
  python -m sem.scripts.benchmark_gpu --grid 2.5 --orientations 32 --gpu-compile
  python -m sem.scripts.benchmark_gpu --skip-fem          # GPU only (no comparison)

Writes <out>/benchmark_summary.md and <out>/benchmark_summary.json.
"""

import argparse
import copy
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO / "petk" / "Demo" / "config.json"
MODES = ("open_pore", "run", "rotation_scan")


def resolve_pdb(config_path, pdb):
    """The analyte path as given (relative to the config's folder), or the same file name next to it."""
    p = Path(pdb)
    for candidate in (p if p.is_absolute() else config_path.parent / p, config_path.parent / p.name):
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"analyte {pdb} not found relative to {config_path.parent}")


def prepare_config(args, out):
    cfg = json.load(open(args.config))
    cfg = copy.deepcopy(cfg)
    cfg["input"]["moving_pdb"] = str(resolve_pdb(args.config, cfg["input"]["moving_pdb"]))
    sim = cfg["simulation"]
    if not args.keep_radii:
        sim["use_pdb2pqr"] = False
        sim["use_vdw_radii"] = True
    if args.grid:
        sim["grid_resolution"] = float(args.grid)
    sim["cleanup_temp_files"] = True
    return cfg


def write_config(cfg, path, prefix, **movement):
    c = copy.deepcopy(cfg)
    c["output"] = dict(c.get("output", {}), output_prefix=str(prefix))
    c["movement"] = dict(c["movement"], **movement)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(c, indent=2))
    return path


def run_sem(python, argv, workdir, log):
    t = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run([python, "-m", "sem"] + argv, cwd=workdir, stdout=fh,
                            stderr=subprocess.STDOUT).returncode
    if rc:
        tail = Path(log).read_text()[-2000:]
        raise SystemExit(f"`sem {' '.join(argv)}` failed (exit {rc}); see {log}\n{tail}")
    return time.time() - t


def read_open_pore(prefix):
    return float(Path(f"{prefix}_open_pore_current.txt").read_text().split()[-1])


def read_run(prefix):
    data = np.loadtxt(f"{prefix}_results.txt", comments="#", ndmin=2)
    return data[:, 0], data[:, 1], data[:, 3]     # z, current, blockage


def read_rotations(folder):
    with open(Path(folder) / "rotation_results.csv") as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r["index"]))
    return np.array([float(r["current_nA"]) for r in rows])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="SEM config (default: petk/Demo/config.json)")
    ap.add_argument("--out", type=Path, default=Path("gpu_benchmark"), help="output folder (default: gpu_benchmark)")
    ap.add_argument("--grid", type=float, default=None, help="override simulation.grid_resolution (A)")
    ap.add_argument("--orientations", type=int, default=16, help="orientations for rotation_scan (default: 16)")
    ap.add_argument("--seed", type=int, default=42, help="seed for the random orientations (default: 42)")
    ap.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    ap.add_argument("--skip-fem", action="store_true", help="run the GPU backend only")
    ap.add_argument("--keep-radii", action="store_true", help="keep the config's radius settings (pdb2pqr etc.)")
    ap.add_argument("--gpu-compile", action="store_true", help="pass --gpu-compile (one-off compile per grid)")
    ap.add_argument("--gpu-batch", type=int, default=64)
    ap.add_argument("--gpu-current", choices=["flux", "power"], default=None)
    ap.add_argument("--fem-python", default=sys.executable, help="interpreter for --backend fem")
    ap.add_argument("--gpu-python", default=sys.executable, help="interpreter for --backend gpu")
    args = ap.parse_args()
    args.config = args.config.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    cfg = prepare_config(args, out)
    from sem.rotation import random_uniform_rotations   # noqa: E402  (after argument parsing)
    angles = out / "angles.txt"
    angles.write_text("".join(f"{r.rx:.12f} {r.ry:.12f} {r.rz:.12f}\n"
                              for r in random_uniform_rotations(args.orientations, seed=args.seed)))

    backends = ["gpu"] if args.skip_fem else ["fem", "gpu"]
    gpu_flags = ["--gpu-batch", str(args.gpu_batch)] + (["--gpu-compile"] if args.gpu_compile else []) \
        + (["--gpu-current", args.gpu_current] if args.gpu_current else [])
    results = {}
    for backend in backends:
        python = args.fem_python if backend == "fem" else args.gpu_python
        flags = ["--backend", backend] + (gpu_flags if backend == "gpu" else [])
        res = results[backend] = {}
        for mode in args.modes:
            d = out / backend / mode
            prefix = d / "bench"
            if mode == "rotation_scan":
                c = write_config(cfg, d / "config.json", prefix, z_start=0.0, z_end=0.0)
                argv = ["rotation_scan", str(c), ".", str(angles), "--mode", "run", "--reuse-mesh",
                        "--reuse-open-pore", "--output-dir", str(d / "out")]
            else:
                c = write_config(cfg, d / "config.json", prefix)
                argv = [mode, str(c)]
            print(f"[{backend}] {mode} ...", flush=True)
            wall = run_sem(python, argv + flags, d, d / "log.txt")
            if mode == "open_pore":
                res[mode] = {"wall_s": wall, "cases": 1, "I0": read_open_pore(prefix)}
            elif mode == "run":
                z, cur, blk = read_run(prefix)
                res[mode] = {"wall_s": wall, "cases": len(z), "z": z.tolist(), "I": cur.tolist(), "blockage": blk.tolist()}
            else:
                cur = read_rotations(d / "out")
                res[mode] = {"wall_s": wall, "cases": len(cur), "I": cur.tolist()}
            print(f"[{backend}] {mode}: {wall:.1f} s for {res[mode]['cases']} case(s)", flush=True)

    sim = cfg["simulation"]
    lines = [
        "# GPU vs FEM benchmark",
        "",
        f"Config: `{args.config}` with analyte `{Path(cfg['input']['moving_pdb']).name}`, "
        f"grid {sim['grid_resolution']} A, box {cfg.get('box_dimensions', 'auto')}, "
        f"{'VdW radii' if not args.keep_radii else 'radii as configured'}.",
        f"GPU options: {' '.join(gpu_flags)}",
        "",
        "Wall-clock times include process start-up, structure preparation and the open-pore solve.",
        "",
        "| mode | cases | FEM s | GPU s | speed-up | max abs(dI)/I | max abs(d blockade) (points) |",
        "|---|---|---|---|---|---|---|",
    ]
    for mode in args.modes:
        g = results["gpu"][mode]
        f = results.get("fem", {}).get(mode)
        if f is None:
            lines.append(f"| {mode} | {g['cases']} | - | {g['wall_s']:.1f} | - | - | - |")
            continue
        if mode == "open_pore":
            dI, dB = abs(g["I0"] / f["I0"] - 1), float("nan")
        else:
            gi, fi = np.array(g["I"]), np.array(f["I"])
            ok = np.isfinite(gi) & np.isfinite(fi)
            dI = float(np.abs(gi[ok] / fi[ok] - 1).max()) if ok.any() else float("nan")
            if mode == "run":
                gb, fb = np.array(g["blockage"]), np.array(f["blockage"])
                dB = float(np.abs(gb[ok] - fb[ok]).max()) if ok.any() else float("nan")
            else:
                g0, f0 = results["gpu"].get("open_pore", {}).get("I0"), results["fem"].get("open_pore", {}).get("I0")
                dB = float(np.abs(100 * (1 - gi[ok] / g0) - 100 * (1 - fi[ok] / f0)).max()) if g0 and f0 else float("nan")
        lines.append(f"| {mode} | {g['cases']} | {f['wall_s']:.1f} | {g['wall_s']:.1f} | "
                     f"{f['wall_s'] / g['wall_s']:.1f}x | {dI:.1e} | {dB:.1e} |")
    text = "\n".join(lines) + "\n"
    (out / "benchmark_summary.md").write_text(text)
    (out / "benchmark_summary.json").write_text(json.dumps({"config": str(args.config), "simulation": sim,
                                                            "gpu_flags": gpu_flags, "results": results}, indent=1))
    print("\n" + text)
    print(f"Written: {out / 'benchmark_summary.md'}")


if __name__ == "__main__":
    main()
