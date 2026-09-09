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
│   ├── cli.py             # `sem` entry point (works in serial and under mpirun)
│   ├── pore_geometry.py   # Geometry builders and distance fields
│   ├── geometry_profiles.py  # Analytic wall geometry + distance metrics (no dolfinx)
│   ├── grid_io.py         # Binary grid I/O + condfrac ramp (no dolfinx)
│   ├── provenance.py      # git_commit() helper (no dolfinx)
│   ├── vertical_movement_sem.py  # Translocation driver
│   ├── conductivity_models.py
│   ├── structure_preparation.py
│   ├── van_der_waals.py
│   ├── visualization.py
│   └── scripts/           # pdb2xyz, gen_dist, derive_geometry, make_validation_table, resample_bin
├── tests/                 # pytest suite; `python3 -m pytest -q` (no dolfinx required)
├── validation/            # All-atom-vs-parametric validation table inputs
├── 1AOI.pdb / centered_1AOI.pdb   # Example nucleosome analyte
├── config.json            # Example run config produced by the GUI
├── environment.yml        # Conda env (DOLFINx + MPI + Python deps)
├── setup.py               # `pip install -e .`
├── CHANGELOG.md
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
only), `sem preview_only config.json`, `sem rotation_scan config.json`,
`sem create_config <pore_type>` to write an example config file, and
`sem derive_geometry <bin>` to fit a parametric `pore_geometry` block (or a
`profile` radius-vs-z table) from an all-atom conductivity/distance map
(see "All-atom workflow" below). Run `sem --help` for the full list.

A minimal `config.json` is included at the repo root and reproduces a
1AOI nucleosome translocating through a 100 Å cylindrical pore.

`pore_geometry.pore_type` (config key, or `sem create_config <pore_type>`)
is one of: `cylindrical`, `double_cone`, `conical`, `profile`, `biological`,
`bin_file`.

## Pore geometry

### Radius convention

For every parametric pore type (`cylindrical`, `double_cone`, `conical`,
`profile`), the configured radius (`pore_radius`, `outer_radius`,
`top_radius`/`bottom_radius`, or a `profile` table's `r` column) is the
radius of the **solid surface** -- the distance-0 point the `condfrac` ramp
starts from, not an atom's center. If you drill a real (all-atom) membrane
to a given atom-center radius, the resulting solid (van der Waals) surface
sits roughly 1 Å inside that atom-center radius, so an all-atom pore and a
parametric pore configured with the *same* radius are not quite the same
pore. Don't guess an offset by hand: use `sem derive_geometry` (below) to
fit the parametric radius that actually reproduces a given all-atom map.

`corner_radius` (cylindrical pores only) is a **linear chamfer**: it cuts a
straight bevel of radial width `corner_radius` starting at the mouth of the
pore, over an axial depth of `chamfer_depth`. If `chamfer_depth` is not
set (`null`), it defaults to `corner_radius` (a 45° chamfer, matching
pre-0.2.0 behavior); set it independently to chamfer over a different axial
depth than the radial cut width.

### Distance metric

`pore_geometry.distance_metric` selects how the wall-to-query-point
distance is computed for parametric pores: `"euclidean"` (the default) is
the true 3-D distance to the solid wall, matching what `gen_dist` measures
for an all-atom membrane; `"legacy"` reproduces the pre-0.2.0 approximation
(an in-plane radial term and a vertical term combined independently, which
overestimates the true distance on any sloped or chamfered wall).

- **Plain (unchamfered) cylinders**: the wall is purely vertical, so the two
  metrics give bit-identical results -- no need to change anything.
- **Chamfered cylinders, `double_cone`, `conical`, and `profile` pores**:
  the two metrics can differ by up to a few percent in open-pore current,
  since the wall is sloped. Use `"legacy"` only to reproduce results
  computed before this change; new work should use the `"euclidean"`
  default.

### All-atom workflow

To simulate against (or fit a parametric approximation to) a real drilled
all-atom membrane, first build its distance-to-solid-surface map:

```bash
# 1. Embed van der Waals radii into the B-factor column of a drilled
#    membrane PDB. Default table: SI 2.10 A, N 1.55 A (sem.van_der_waals);
#    pass radius_table={'SI': 1.9, 'N': 1.4} to reproduce the paper's
#    Si3N4 membranes. Membrane/pore atoms must be HETATM records -- ATOM
#    records get protein-style atom-name radii instead of element radii.
python -c "from sem.structure_preparation import pdb_add_radii_in_bfactor as f; f('membrane.pdb', 'membrane_radii.pdb')"

# 2. Convert to XYZ (x y z radius per line); pdb2xyz reads the radius
#    from the tempfactor column written above.
pdb2xyz membrane_radii.pdb membrane.xyz

# 3. Generate the distance field. Argv order is MAX then MIN:
#    gen_dist XYZ  MaxX MaxY MaxZ  MinX MinY MinZ  resolution cutoff OUT
gen_dist membrane.xyz 75 75 75 -75 -75 -75 1.0 5.0 membrane.bin
```

From `membrane.bin`, either:

1. **Simulate directly against the all-atom map** (exact; no parametric
   fit) with a `bin_file` pore:

   ```json
   "pore_geometry": {
     "pore_type": "bin_file",
     "bin_file_path": "membrane.bin",
     "bin_file_units": "distance",
     "membrane_thickness": 0.0
   }
   ```

2. **Fit a parametric (or `profile`) approximation** with
   `sem derive_geometry`, so downstream runs (and other tools that expect a
   simple radius/thickness geometry) get a lightweight equivalent:

   ```bash
   sem derive_geometry membrane.bin --pore-type cylindrical \
       --template config.json --write-config config_derived.json
   ```

   which writes a `pore_geometry` block into `config_derived.json` shaped
   like (illustrative values -- the fitted numbers depend on `membrane.bin`):

   ```json
   "pore_geometry": {
     "pore_type": "cylindrical",
     "pore_radius": 49.1,
     "corner_radius": 0.0,
     "chamfer_depth": null,
     "membrane_thickness": 100.2,
     "distance_metric": "euclidean"
   }
   ```

   or, with `--pore-type profile`, a `profile` block pointing at the
   fitted radius-vs-z table it also writes out
   (`{output-prefix}_profile_pore.csv`):

   ```json
   "pore_geometry": {
     "pore_type": "profile",
     "profile_path": "derived_profile_pore.csv",
     "membrane_thickness": 100.2,
     "distance_metric": "euclidean"
   }
   ```

   `--pore-type` also accepts `double_cone` and `conical`; add
   `--fit-corner` to also fit `corner_radius`/`chamfer_depth` for a
   cylindrical pore. Run `sem derive_geometry --help` for the full option
   list, or see `sem/scripts/derive_geometry.py`'s module docstring for the
   fitting approach.

## Citation

If you use this software in academic work, please cite it via the
`CITATION.cff` at the repo root, or the GitHub "Cite this repository"
button.

## License

Released under the [MIT License](LICENSE).
