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

- Build cylindrical, double-cone, biological, or arbitrary-geometry pores.
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
│   ├── cli.py             # `sem` and `sem-mpi` entry points
│   ├── pore_geometry.py   # Geometry builders and distance fields
│   ├── vertical_movement_sem.py  # Translocation driver
│   ├── conductivity_models.py
│   ├── structure_preparation.py
│   ├── van_der_waals.py
│   ├── visualization.py
│   └── scripts/           # pdb2xyz, gen_dist, resample_bin
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
sem run --config config.json
```

Run in parallel:

```bash
mpirun -n 8 sem-mpi run --config config.json
```

A minimal `config.json` is included at the repo root and reproduces a
1AOI nucleosome translocating through a 100 Å cylindrical pore.

## Citation

If you use this software in academic work, please cite it via the
`CITATION.cff` at the repo root, or the GitHub "Cite this repository"
button.

## License

Released under the [MIT License](LICENSE).
