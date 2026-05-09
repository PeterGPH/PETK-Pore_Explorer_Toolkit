#!/usr/bin/env python3
"""
Comprehensive integration test for SEM's ARBD potential-grid export.

What it covers:
  - open_pore mode for cylindrical and double_cone geometries
  - uniform vs multigrid (gmsh) mesh modes
  - run mode (with-analyte) with stride 0 (final z only) and stride 2
  - custom ion species and custom export resolution
  - the new "one-shot no-analyte + per-z with-analyte" behaviour

Each test case:
  1. Generates a config.json under <output-dir>/<case_name>/
  2. (Optional) invokes `python -m sem <config> <mode>` via the conda env
  3. Checks the expected list of .dx files exists
  4. Validates DX header + non-trivial size for each file

Usage:
  # Dry run — just emit configs, do not run SEM
  python sem/scripts/test_arbd_export_integration.py --dry-run \\
      --output-dir test_runs

  # Full run on a workstation that has the dolfinx env
  python sem/scripts/test_arbd_export_integration.py \\
      --conda-env dolfinx-test-fem-setup \\
      --analyte-pdb petk/analytes/centered_amino_acids/centered_GLY.pdb \\
      --output-dir test_runs

  # Run only the cheap open_pore cases
  python sem/scripts/test_arbd_export_integration.py \\
      --conda-env dolfinx-test-fem-setup \\
      --filter open_pore --output-dir test_runs

Exit code is 0 iff every selected case passed.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Test-case definition
# ---------------------------------------------------------------------------
@dataclass
class TestCase:
    name: str
    description: str
    mode: str                                     # "open_pore" or "run"
    pore_type: str                                # cylindrical / double_cone
    geometry_params: Dict[str, float]             # pore_radius, outer_radius...
    membrane_thickness_A: float = 200.0
    voltage_mv: float = 200.0
    sigma_S_per_m: float = 10.5
    grid_resolution_A: float = 5.0
    box_edge_A: float = 600.0
    z_start_A: float = 100.0
    z_end_A: float = -100.0
    z_step_A: float = 50.0                        # coarse for speed
    arbd_ions: List[List] = field(
        default_factory=lambda: [["POT", +1], ["CLA", -1]]
    )
    arbd_resolution_A: Optional[float] = None
    arbd_stride: int = 0
    arbd_wall_height: float = 100.0
    arbd_temperature_K: float = 295.0
    mesh_mode: str = "uniform"                    # uniform / multigrid
    gmsh_fine_size_A: Optional[float] = None
    gmsh_coarse_size_A: Optional[float] = None
    gmsh_fine_box_A: Optional[List[float]] = None
    output_prefix_override: Optional[str] = None  # else "test_<name>"

    @property
    def output_prefix(self) -> str:
        return self.output_prefix_override or f"test_{self.name}"


# ---------------------------------------------------------------------------
# The test cases
# ---------------------------------------------------------------------------
def build_test_cases() -> List[TestCase]:
    """The matrix of cases. Tweak values to balance coverage vs runtime."""
    cyl = {"pore_radius": 50.0, "corner_radius": 0.0}      # 5 nm pore
    dc  = {"pore_radius": 30.0, "outer_radius": 70.0}      # double-cone
    return [
        TestCase(
            name="open_pore_cyl_uniform",
            description="Open-pore cylindrical, uniform mesh, default ions.",
            mode="open_pore",
            pore_type="cylindrical",
            geometry_params=cyl,
        ),
        TestCase(
            name="open_pore_dc_uniform",
            description="Open-pore double-cone, uniform mesh, default ions.",
            mode="open_pore",
            pore_type="double_cone",
            geometry_params=dc,
        ),
        TestCase(
            name="open_pore_cyl_multigrid",
            description="Open-pore cylindrical, gmsh graded mesh.",
            mode="open_pore",
            pore_type="cylindrical",
            geometry_params=cyl,
            mesh_mode="multigrid",
            gmsh_fine_size_A=2.0,
            gmsh_coarse_size_A=10.0,
            gmsh_fine_box_A=[200.0, 200.0, 200.0],
        ),
        TestCase(
            name="open_pore_cyl_custom_ions",
            description="Open-pore cylindrical with three ion species.",
            mode="open_pore",
            pore_type="cylindrical",
            geometry_params=cyl,
            arbd_ions=[["NA", +1], ["KK", +1], ["CL", -1]],
        ),
        TestCase(
            name="open_pore_cyl_custom_resolution",
            description="Open-pore cylindrical, coarse 8 Å ARBD resolution.",
            mode="open_pore",
            pore_type="cylindrical",
            geometry_params=cyl,
            arbd_resolution_A=8.0,
        ),
        TestCase(
            name="run_cyl_stride0",
            description=(
                "Translocation cylindrical, stride 0 → no-analyte set + "
                "with-analyte snapshot at the FINAL z only."
            ),
            mode="run",
            pore_type="cylindrical",
            geometry_params=cyl,
            z_start_A=100.0, z_end_A=-100.0, z_step_A=50.0,  # 5 z steps
            arbd_stride=0,
        ),
        TestCase(
            name="run_cyl_stride2",
            description=(
                "Translocation cylindrical, stride 2 → no-analyte set + "
                "with-analyte snapshots at every-other z plus final."
            ),
            mode="run",
            pore_type="cylindrical",
            geometry_params=cyl,
            z_start_A=100.0, z_end_A=-100.0, z_step_A=50.0,  # 5 z steps
            arbd_stride=2,
        ),
    ]


# ---------------------------------------------------------------------------
# Config writer
# ---------------------------------------------------------------------------
def make_config(tc: TestCase, analyte_pdb: Optional[str]) -> dict:
    """Build a SEM config dict for one test case."""
    half = tc.box_edge_A / 2.0
    cfg = {
        "metadata": {
            "generated_by": "test_arbd_export_integration.py",
            "case_name": tc.name,
            "description": tc.description,
        },
        "input": {
            "moving_pdb": (str(analyte_pdb) if (tc.mode == "run" and analyte_pdb) else "")
        },
        "pore_geometry": {
            "pore_type": tc.pore_type,
            "membrane_thickness": tc.membrane_thickness_A,
            **tc.geometry_params,
        },
        "simulation": {
            "voltage": tc.voltage_mv,
            "bulk_conductivity": tc.sigma_S_per_m,
            "grid_resolution": tc.grid_resolution_A,
            "use_vdw_radii": False,
            "use_pdb2pqr": False,
            "force_field": "CHARMM",
            "default_radius": 1.5,
            "membrane_conductivity": 1e-7,
        },
        "movement": {
            "z_start": tc.z_start_A if tc.mode == "run" else 0.0,
            "z_end":   tc.z_end_A   if tc.mode == "run" else 0.0,
            "z_step":  tc.z_step_A,
        },
        "output": {
            "output_prefix": tc.output_prefix,
            "preview_frames": 0,
        },
        "box_dimensions": {
            "x": [-half, half],
            "y": [-half, half],
            "z": [-half, half],
        },
    }
    if tc.mesh_mode == "multigrid":
        cfg["simulation"]["mesh_engine"] = "gmsh"
        if tc.gmsh_fine_size_A is not None:
            cfg["simulation"]["gmsh_fine_size"] = tc.gmsh_fine_size_A
        if tc.gmsh_coarse_size_A is not None:
            cfg["simulation"]["gmsh_coarse_size"] = tc.gmsh_coarse_size_A
        if tc.gmsh_fine_box_A is not None:
            cfg["simulation"]["gmsh_fine_box"] = list(tc.gmsh_fine_box_A)

    arbd = {
        "ions": [list(pair) for pair in tc.arbd_ions],
        "stride": tc.arbd_stride,
        "wall_height": tc.arbd_wall_height,
        "temperature_K": tc.arbd_temperature_K,
    }
    if tc.arbd_resolution_A is not None:
        arbd["resolution"] = tc.arbd_resolution_A
    cfg["output"]["arbd_export"] = arbd
    return cfg


# ---------------------------------------------------------------------------
# Expected-files map
# ---------------------------------------------------------------------------
def _z_positions(tc: TestCase) -> List[float]:
    """Replicate the z grid that SEM's run() uses."""
    if tc.mode != "run":
        return []
    n = int(abs(tc.z_end_A - tc.z_start_A) / tc.z_step_A) + 1
    direction = 1.0 if tc.z_end_A > tc.z_start_A else -1.0
    return [tc.z_start_A + direction * tc.z_step_A * i for i in range(n)]


def _z_tag(z_A: float) -> str:
    """Mirror visualization.py / run() filename convention: '_z+020.0A'."""
    return f"_z{z_A:+07.1f}A".replace(" ", "0")


def expected_files(tc: TestCase) -> List[str]:
    """Filenames (basenames) that should exist in the run output directory."""
    files: List[str] = []
    base = f"{tc.output_prefix}_{tc.pore_type}"
    ion_names = [name for name, _ in tc.arbd_ions]
    if tc.mode == "open_pore":
        # Single set, no z-tag, no _openpore suffix.
        files.append(f"{base}_steric.dx")
        files.append(f"{base}_open_pore_phi.dx")
        for ion in ion_names:
            files.append(f"{base}_{ion}.dx")
    else:  # run
        # 1) One-shot no-analyte snapshot.
        files.append(f"{base}_openpore_steric.dx")
        files.append(f"{base}_openpore_open_pore_phi.dx")
        for ion in ion_names:
            files.append(f"{base}_openpore_{ion}.dx")
        # 2) Per-z exports.
        zs = _z_positions(tc)
        for i, z in enumerate(zs):
            is_last = (i == len(zs) - 1)
            should = is_last or (tc.arbd_stride > 0 and i % tc.arbd_stride == 0)
            if not should:
                continue
            tag = _z_tag(z)
            files.append(f"{base}{tag}_steric.dx")
            files.append(f"{base}{tag}_open_pore_phi.dx")
            for ion in ion_names:
                files.append(f"{base}{tag}_{ion}.dx")
    return files


# ---------------------------------------------------------------------------
# DX-file validator
# ---------------------------------------------------------------------------
_DX_HEADER_RE = re.compile(r"^\s*object\s+\d", re.IGNORECASE | re.MULTILINE)


def validate_dx(path: Path) -> Tuple[bool, str]:
    """Lightweight check: DX file exists, has header, has at least one numeric line."""
    if not path.exists():
        return False, "missing"
    size = path.stat().st_size
    if size < 200:
        return False, f"too small ({size} B)"
    head = path.read_text(errors="ignore")[:4000]
    if not _DX_HEADER_RE.search(head):
        return False, "no DX 'object N …' header"
    # Spot-check that there is at least one numeric data row.
    found_numeric = False
    with path.open("r", errors="ignore") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "object", "attribute", "component")):
                continue
            try:
                float(stripped.split()[0])
                found_numeric = True
                break
            except (ValueError, IndexError):
                continue
    if not found_numeric:
        return False, "no numeric data lines"
    return True, "ok"


# ---------------------------------------------------------------------------
# SEM invocation
# ---------------------------------------------------------------------------
def run_sem(config_path: Path, mode: str, conda_env: str,
            python_exec: str, output_dir: Path,
            log_path: Path) -> Tuple[bool, str]:
    """Invoke `python -m sem <config> <mode>` from output_dir."""
    abs_config = config_path.resolve()
    abs_outdir = output_dir.resolve()
    if conda_env:
        cmd = ["conda", "run", "-n", conda_env, "python", "-m", "sem",
               str(abs_config), mode]
    else:
        cmd = [python_exec, "-m", "sem", str(abs_config), mode]
    try:
        with log_path.open("w") as logfh:
            logfh.write(f"# cmd: {' '.join(cmd)}\n")
            logfh.write(f"# cwd: {abs_outdir}\n")
            logfh.flush()
            t0 = time.time()
            subprocess.run(
                cmd, cwd=str(abs_outdir),
                stdout=logfh, stderr=subprocess.STDOUT, check=True,
                timeout=3600,
            )
            elapsed = time.time() - t0
        return True, f"ok ({elapsed:.1f}s)"
    except subprocess.CalledProcessError as e:
        return False, f"non-zero exit {e.returncode} (see {log_path.name})"
    except subprocess.TimeoutExpired:
        return False, "timed out (>1h)"
    except FileNotFoundError as e:
        return False, f"command not found: {e}"


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------
@dataclass
class CaseResult:
    case: TestCase
    ran: bool
    run_msg: str
    expected: List[str]
    missing: List[str]
    invalid: List[Tuple[str, str]]   # (file, reason)
    elapsed_s: float = 0.0

    @property
    def passed(self) -> bool:
        return self.ran and not self.missing and not self.invalid


def write_case_config(tc: TestCase, root_dir: Path,
                      analyte_pdb: Optional[Path]) -> Path:
    case_dir = root_dir / tc.name
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg = make_config(tc, str(analyte_pdb) if analyte_pdb else None)
    cfg_path = case_dir / "config.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))
    return cfg_path


def execute_case(tc: TestCase, root_dir: Path, args) -> CaseResult:
    case_dir = root_dir / tc.name
    cfg_path = write_case_config(tc, root_dir, args.analyte_pdb)
    expected = expected_files(tc)

    if args.dry_run:
        return CaseResult(case=tc, ran=False, run_msg="dry-run",
                          expected=expected, missing=[], invalid=[])

    log_path = case_dir / "sem_run.log"
    t0 = time.time()
    ok, run_msg = run_sem(cfg_path, tc.mode, args.conda_env,
                          args.python_exec, case_dir, log_path)
    elapsed = time.time() - t0
    if not ok:
        return CaseResult(case=tc, ran=False, run_msg=run_msg,
                          expected=expected, missing=expected,
                          invalid=[], elapsed_s=elapsed)

    missing: List[str] = []
    invalid: List[Tuple[str, str]] = []
    for fname in expected:
        path = case_dir / fname
        if not path.exists():
            missing.append(fname)
            continue
        valid, reason = validate_dx(path)
        if not valid:
            invalid.append((fname, reason))
    return CaseResult(case=tc, ran=True, run_msg=run_msg,
                      expected=expected, missing=missing,
                      invalid=invalid, elapsed_s=elapsed)


def print_report(results: List[CaseResult], dry_run: bool) -> int:
    n_total = len(results)
    n_pass = sum(1 for r in results if r.passed)
    n_run_failed = sum(1 for r in results if not r.ran and r.run_msg != "dry-run")

    print()
    print("=" * 78)
    title = " ARBD-export integration test report (DRY RUN) " if dry_run else " ARBD-export integration test report "
    print(title.center(78, "="))
    print("=" * 78)

    for r in results:
        tc = r.case
        if dry_run:
            status = "CFG"
            extra = f"  expected files: {len(r.expected)}"
        else:
            status = "PASS" if r.passed else ("FAIL_RUN" if not r.ran else "FAIL_VERIFY")
            extra = f"  ({r.run_msg}, expected={len(r.expected)}, " \
                    f"missing={len(r.missing)}, invalid={len(r.invalid)})"
        print(f"  [{status:<11}] {tc.name:<32} {tc.description}")
        print(f"               {extra}")
        if r.missing:
            print(f"               missing files:")
            for m in r.missing[:5]:
                print(f"                 - {m}")
            if len(r.missing) > 5:
                print(f"                 ... and {len(r.missing) - 5} more")
        if r.invalid:
            print(f"               invalid files:")
            for fn, reason in r.invalid[:5]:
                print(f"                 - {fn}: {reason}")
            if len(r.invalid) > 5:
                print(f"                 ... and {len(r.invalid) - 5} more")

    print("-" * 78)
    if dry_run:
        print(f"  {n_total} configs generated.")
        return 0
    print(f"  {n_pass}/{n_total} cases passed "
          f"({n_run_failed} failed to run, "
          f"{n_total - n_pass - n_run_failed} ran but failed verification)")
    print("=" * 78)
    return 0 if n_pass == n_total else 1


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--output-dir", type=Path,
                   default=Path("arbd_test_runs"),
                   help="Per-case subdirs land under this folder (default: arbd_test_runs).")
    p.add_argument("--analyte-pdb", type=Path, default=None,
                   help="PDB file used by 'run' mode cases. Required for those "
                        "cases unless --filter excludes them. Tip: use a small "
                        "centered amino acid like centered_GLY.pdb for speed.")
    p.add_argument("--conda-env", default="",
                   help="Conda env to invoke for SEM (e.g. dolfinx-test-fem-setup). "
                        "Empty → use --python-exec directly.")
    p.add_argument("--python-exec", default=sys.executable,
                   help="Python executable for `python -m sem` (used when "
                        "--conda-env is empty).")
    p.add_argument("--filter", default=None,
                   help="Substring filter on case name (e.g. 'open_pore' to "
                        "skip the heavier run-mode cases).")
    p.add_argument("--dry-run", action="store_true",
                   help="Generate configs only; do not invoke SEM.")
    p.add_argument("--list", action="store_true",
                   help="Print the test case names and exit.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    cases = build_test_cases()
    if args.filter:
        cases = [c for c in cases if args.filter in c.name]

    if args.list:
        for c in cases:
            print(f"{c.name:<32} {c.mode:<10} {c.description}")
        return 0

    if not cases:
        print("No matching test cases.", file=sys.stderr)
        return 2

    needs_analyte = any(c.mode == "run" for c in cases)
    if needs_analyte and not args.dry_run:
        if args.analyte_pdb is None:
            print(
                "ERROR: --analyte-pdb is required when run-mode cases are selected.\n"
                "       Either pass --analyte-pdb <path>, or filter to open_pore "
                "cases via --filter open_pore.",
                file=sys.stderr,
            )
            return 2
        if not args.analyte_pdb.exists():
            print(f"ERROR: --analyte-pdb file does not exist: {args.analyte_pdb}",
                  file=sys.stderr)
            return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Test root : {args.output_dir.resolve()}")
    print(f"Mode      : {'DRY RUN' if args.dry_run else 'FULL RUN'}")
    if not args.dry_run:
        print(f"Conda env : {args.conda_env or '(none, using --python-exec)'}")
        if args.analyte_pdb:
            print(f"Analyte   : {args.analyte_pdb}")
    print(f"Cases     : {len(cases)}")
    for c in cases:
        print(f"  - {c.name}  ({c.mode}, {c.pore_type}, {c.mesh_mode})")
    print()

    results: List[CaseResult] = []
    for tc in cases:
        print(f"=== {tc.name} ===")
        r = execute_case(tc, args.output_dir, args)
        results.append(r)
        if args.dry_run:
            print(f"  config written, expects {len(r.expected)} dx files")
        else:
            tag = "PASS" if r.passed else ("FAIL_RUN" if not r.ran else "FAIL_VERIFY")
            print(f"  → {tag}  ({r.run_msg})  "
                  f"missing={len(r.missing)}  invalid={len(r.invalid)}  "
                  f"elapsed={r.elapsed_s:.1f}s")

    return print_report(results, args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())
