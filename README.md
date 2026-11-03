# PoreExplorer-PETK

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/status-beta-orange.svg)]()

**PoreExplorer Toolkit (PETK)** is a VMD-based GUI plus a Python package
(`nanopore-sem`) for setting up and running **Steric Exclusion Model (SEM)**
calculations on nanopore systems. It pairs an interactive Tcl/Tk front-end
inside [VMD](https://www.ks.uiuc.edu/Research/vmd/) with a DOLFINx-powered
finite-element solver that runs in parallel via MPI.

Use it to:

- Build cylindrical, conical, double-cone, elliptical, biological, or arbitrary-geometry pores.
- Position analyte structures (PDB) inside the pore and translate them along z.
- Solve the Laplace problem for ionic current with steric exclusion at every frame.
- Produce ionic-current vs. translocation-distance traces for benchmarking
  against experimental nanopore data.

## Repository layout

```
.
├── petk/                  # VMD GUI (Tcl/Tk) — petk_gui.tcl + tab1/2/3
│   ├── analytes/          # Example analyte PDBs
│   ├── shapes/            # Pore-geometry templates
│   └── Demo/              # Demo inputs and reference outputs
├── sem/                   # Python package: nanopore-sem
│   ├── cli.py             # `sem` entry point (works in serial and under mpirun)
│   ├── pore_geometry.py   # Geometry builders and distance fields
│   ├── vertical_movement_sem.py  # Translocation driver
│   ├── conductivity_models.py
│   ├── structure_preparation.py
│   ├── van_der_waals.py
│   ├── visualization.py
│   ├── gpu_solver.py      # `--backend gpu`: the same FEM, batched on a GPU (PyTorch)
│   └── scripts/           # pdb2xyz, gen_dist, resample_bin, benchmark_gpu
├── 1AOI.pdb / centered_1AOI.pdb   # Example nucleosome analyte
├── config.json            # Example run config produced by the GUI
├── environment.yml        # Conda env (DOLFINx + MPI + Python deps)
├── setup.py               # `pip install -e .`
└── LICENSE
```

## Prerequisites

- [VMD](https://www.ks.uiuc.edu/Research/vmd/) (1.9.3 or newer) for the GUI.
- [Anaconda or Miniconda](https://docs.conda.io/en/latest/miniconda.html) for
  the Python solver environment.
- A working MPI implementation; `environment.yml` pulls `mpich` from
  conda-forge by default.

## Installation

```bash
# 1. Clone
git clone https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit.git
cd PETK-Pore_Explorer_Toolkit

# 2. Create the conda env (DOLFINx + MPI + numerical core)
conda env create -f environment.yml
conda activate sem-env

# 3. Install the Python package in editable mode
pip install -e .

# 4. (Optional) verify the install
sem --help
```

## Quick start (GUI)

Launch VMD and open the **Tk Console** from the *Extensions* menu:

```tcl
source /path/to/PoreExplorer-PETK/petk/petk_gui.tcl
::PETK::gui::petk_gui
```

Or load it as a package:

```tcl
set petk_root "/path/to/PoreExplorer-PETK/petk"
lappend auto_path $petk_root
package require petk_gui
::PETK::gui::petk_gui
```

The GUI walks you through three tabs: **Nanopore** (build geometry),
**Analyte** (load and align a PDB), and **SEM** (configure voltage, grid
resolution, z-range, and run). It writes a `config.json` and invokes the
`sem` CLI for you.

## Quick start (CLI / headless)

Run a SEM calculation directly from a config:

```bash
sem run config.json
```

Run in parallel — `mpi4py` inside `sem.cli:main` picks up the MPI rank
automatically, so the same `sem` entry point works under `mpirun`:

```bash
mpirun -n 8 sem run config.json
```

Other available subcommands: `sem open_pore config.json` (open-pore current
only), `sem preview_only config.json`, `sem rotation_scan config.json`, and
`sem create_config <pore_type>` to write an example config file. Run
`sem --help` for the full list.

A minimal `config.json` is included at the repo root and reproduces a
1AOI nucleosome translocating through a 100 Å cylindrical pore.

## GPU backend

Every solving command takes `--backend gpu`. It solves the same P1 finite-element problem as the
DOLFINx backend on a GPU, for many z positions or orientations at once, and writes the same
output files:

```bash
sem run config.json --backend gpu            # z scan: <prefix>_results.txt (+ ARBD grids if configured)
sem open_pore config.json --backend gpu      # <prefix>_open_pore_current.txt
sem rotation_scan config.json . angles.txt --mode run --backend gpu   # rotation_results.csv, hybrid_currents.csv
sem rotation_scan config.json . angles.txt --mode open_pore --backend gpu
```

`preview_only` solves nothing, so `--backend` does not change it.

**How it matches the FEM.** DOLFINx `create_box` splits every cube into the 6 Kuhn tetrahedra
around its main diagonal. On that mesh the P1 stiffness matrix couples each vertex only to its 6
axis neighbours: an axis edge gets (area / length) / 6 times the sum of the conductivity of the 6
tetrahedra around it, and the diagonal edges get exactly zero. The GPU builds that stencil
directly, samples the conductivity at the tetrahedron centroids with the FEM's own rule (pore
interpolator, nearest atom centre within `cutoff`, conductivity model), and reports the same
top-facet flux. It solves with multigrid-preconditioned conjugate gradients from the open-pore
solution.

**Install.** Add PyTorch with the CUDA build for your GPU to the solver environment, e.g.
`pip install torch --index-url https://download.pytorch.org/whl/cu128` (the `gpu` extra in
`setup.py` lists it). If PyTorch then fails to import NumPy with a `GLIBCXX` error, put the conda
environment's `lib` folder first: `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`.

**Options.**

| option | default | meaning |
|---|---|---|
| `--gpu-method` | `fem` | `fem`: the FEM above. `fv`: cell-centred finite volumes (faster, but a different discretisation: 0.2-0.45 blockade points above the FEM at 1 Å) |
| `--gpu-current` | `flux` | `flux`: the FEM's top-facet flux (float64, 14 iterations). `power`: dissipated power / voltage (float32, 8 iterations, 1.7-2.5x faster, within ~0.02 points of the flux) |
| `--gpu-batch` | 64 | z positions or orientations per GPU batch |
| `--gpu-compile` | off | `torch.compile` the solve: 2.4-3.6x faster per case after a one-off compile of 2-6 minutes per grid size |
| `--gpu-iters`, `--gpu-dtype`, `--gpu-device` | | override the defaults above |

**Accuracy.** On 2,880 test cases (6 proteins, 4 pores, 1.0 / 1.5 / 2.0 Å grids, cutoffs 5 and
6 Å) the GPU blockade agrees with the DOLFINx backend to 7e-5 points on average and 0.009 points at
worst, and the open-pore current to 2e-7. `tests/test_gpu_solver.py` checks the stencil against an
element-by-element assembly, the potential interpolation, and the conductivity rule.

**Speed.** On one A100 (batch 64, compiled, flux current) a 12 nm box takes 15 ms per orientation
at 1 Å, 5 ms at 1.5 Å and 2.8 ms at 2 Å; the DOLFINx backend takes about 30-60 s and 3.5-7 s per
orientation at 1 Å and 2 Å on one core.

**Benchmark.** `sem/scripts/benchmark_gpu.py` runs `open_pore`, `run` and `rotation_scan` on both
backends through the CLI and compares the outputs (default workload: `petk/Demo/config.json`):

```bash
python -m sem.scripts.benchmark_gpu --out gpu_benchmark               # demo grid (5 Å)
python -m sem.scripts.benchmark_gpu --grid 2.5 --out gpu_benchmark_2p5
```

Reference output (one A100; the DOLFINx backend as one process on the same node; wall-clock
times including process start-up; 1AOI nucleosome, 30 nm box):

| grid | mode | cases | DOLFINx | GPU | speed-up | max abs(d blockade) |
|---|---|---|---|---|---|---|
| 5 Å (demo) | `run` (z scan) | 23 | 103 s | 7.5 s | 14x | 0.001 points |
| 5 Å (demo) | `rotation_scan` | 16 | 75 s | 7.5 s | 10x | 0.0006 points |
| 2.5 Å | `run` (z scan) | 23 | 930 s | 10.2 s | 91x | 0.006 points |
| 2.5 Å | `rotation_scan` | 16 | 650 s | 9.9 s | 66x | 0.005 points |

Most of the GPU time here is start-up and setup (about 7 s), so the speed-up grows with the grid
size and with the number of z positions or orientations.

**Not supported on the GPU:** gmsh meshes (`mesh_engine: gmsh`), `use_charges`, MPI (run one
process per GPU), and ARBD grids from `rotation_scan` (use `run` for one orientation). The GPU
backend applies `prevent_analyte_overlap` exactly as the FEM backend does.

## Citation

If you use this software in academic work, please cite it via the
`CITATION.cff` at the repo root, or the GitHub "Cite this repository"
button.

## License

Released under the [MIT License](LICENSE).
