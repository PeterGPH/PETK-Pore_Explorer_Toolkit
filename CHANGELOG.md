# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-09

### Changed

- **Euclidean wall distance is now the default** for every parametric pore
  type (`cylindrical`, `double_cone`, `conical`, `profile`). Previously the
  parametric path combined an in-plane radial term and a vertical term
  independently (exact only for a purely vertical wall); it overestimated
  the true 3-D distance on any sloped wall or chamfer by up to `1/cos(theta)`
  (about 12% at slope 0.5, up to 41% on the default 45 degree chamfer) and
  rounded the pore mouth spuriously. The new default matches the true
  Euclidean distance the all-atom (`gen_dist`) path measures, so parametric
  and all-atom pores are now built on the same geometric definition.
  - For a plain, unchamfered cylinder the two metrics are mathematically
    identical (a purely vertical wall), so open-pore current is bit-for-bit
    unchanged.
  - For a chamfered cylinder, `double_cone`, or `conical` pore, results
    change by up to a few percent; see `README.md`'s "Distance metric"
    section.
- `setup.py`'s package version is now read from `sem/__init__.py`
  (`sem.__version__`) at build time instead of being hardcoded, so
  `sem.__version__`, `setup.py --version`, and `CITATION.cff`'s `version:`
  field cannot drift out of sync. All three now read `0.2.0`.
- `sem/__init__.py` restores PEP 562 lazy attribute loading
  (`_LAZY` + `__getattr__`): `import sem` (and `sem.grid_io`,
  `sem.pore_geometry`, `sem.config`, `sem.cli`,
  `sem.scripts.derive_geometry`) no longer eagerly imports `dolfinx`.
  Dolfinx-backed symbols (`sem.VerticalMovementSEM`, `sem.PoreGeometry`,
  ...) are imported on first access instead, with a clean `ImportError` if
  dolfinx genuinely is not installed. This restores behaviour lost when a
  prior eager-import regression (`13ec995`) undid the original lazy loader
  (`b248de5`).
- The numpy-only binary-grid I/O and `condfrac` ramp helpers shared by
  `sem.pore_geometry`, `sem.utils`, and `gen_dist` now live in a new
  `sem.grid_io` module (no scipy, no dolfinx); the analytic parametric-wall
  geometry (`PoreProfile` and both distance metrics) lives in a new
  `sem.geometry_profiles` module (no scipy, no dolfinx either). This lets
  `sem.pore_geometry` itself import cleanly without dolfinx (it now gets
  `condfrac`/`readbinGrid` from `sem.grid_io` instead of `sem.utils`, which
  still imports dolfinx eagerly and re-exports them for backwards
  compatibility).

### Added

- `pore_geometry.distance_metric` config key (`"euclidean"` | `"legacy"`,
  default `"euclidean"`). Set `"legacy"` to reproduce pre-0.2.0 results
  exactly.
- `pore_geometry.chamfer_depth` config key: the axial depth of a
  cylindrical pore's mouth chamfer, independent of its radial width
  (`corner_radius`). Defaults to `corner_radius` when unset (matching the
  pre-0.2.0 45 degree chamfer), and is validated in `sem/config.py`.
- New `profile` pore type: a pore wall defined by an arbitrary
  radius-vs-z table (a two-column `z,r` CSV), for pores whose wall does not
  fit a simple analytic shape. Config key `pore_geometry.profile_path`.
  Wired through `sem/config.py`, `sem/cli.py`
  (`create_sem_from_config`/`create_example_config`), and
  `sem/scripts/benchmark_kowalczyk_2011.py`.
- `sem derive_geometry` (`sem/scripts/derive_geometry.py`, also runnable as
  `python -m sem.scripts.derive_geometry`): fits a parametric
  `pore_geometry` block (or a `profile` radius-vs-z table) to an all-atom
  `gen_dist` conductivity/distance map, so a drilled all-atom membrane and
  its parametric approximation can be matched by tool rather than by hand.
  Numpy/scipy only -- no dolfinx dependency. See `README.md`'s "All-atom
  workflow" section for the full command chain
  (`pdb_add_radii_in_bfactor` -> `pdb2xyz` -> `gen_dist` ->
  `sem derive_geometry`).
- `sem/scripts/make_validation_table.py`
  (`python -m sem.scripts.make_validation_table results.csv --output rows.tex`):
  generates LaTeX table rows directly from a validation `results.csv`
  (currents rounded to 2 decimals, deviation recomputed from full
  precision with an explicit sign), replacing hand-transcription of the
  paper's all-atom-vs-parametric comparison table. See
  `validation/all_atom_match/README.md`.
- `validation/all_atom_match/`: describes the six paper validation systems
  and the cluster paths to their all-atom `.bin` maps, with an
  empty-header `results.csv` to be filled in by the cluster validation
  rerun.
- Provenance lines in `{prefix}_open_pore_current.txt` (both the file and
  the console summary) and in `derive_geometry`'s `*_derivation.json`
  (`provenance` object): `sem_version`, `git_commit` (short SHA, or
  `"unknown"` if it cannot be determined -- never raises), and `ramp`
  (the `condfrac` ramp bounds, `RAMP_MIN`/`RAMP_MAX`). All new lines in the
  result file are `#`-prefixed comments, so existing parsers
  (`benchmark_kowalczyk_2011.parse_open_pore_result`,
  `plot_open_pore_convergence.py`) are unaffected. Implemented via a shared
  helper, `sem.provenance.git_commit`, and a pure formatting function,
  `sem.cli._open_pore_result_lines`, both unit-tested without dolfinx.
- `sem create_config conical` now works (`conical` was missing from the
  `create_config` pore-type choices, though `sem.config.create_example_config`
  already had a `conical` template), and `create_example_config` gained a
  new `profile` template.
- `sem/scripts/benchmark_kowalczyk_2011.py` (the analytical-validation
  benchmark) gained a `--distance-metric` flag and two new `--geometries`
  choices, `cylindrical_corner` (a cylinder with a rounded/chamfered
  corner) and `conical`, alongside its existing `cylindrical`/
  `double_cone` geometries.
- `tests/`: a pytest suite (`pytest.ini` sets `pythonpath = .`) covering
  config validation, the Euclidean/legacy distance metrics, `derive_geometry`
  round-trips, the provenance helpers, `make_validation_table`, and the
  dolfinx-free import surface -- runnable with `python3 -m pytest -q` in a
  pip-only environment with no dolfinx installed (dolfinx-only tests skip
  cleanly via `pytest.importorskip("dolfinx")`).
- `tests/test_import_surface.py`: extended to also import `sem.cli` and
  `sem.scripts.derive_geometry` with dolfinx blocked, and to check
  `sem.__version__` against both `setup.py --version` and `CITATION.cff`.

### Fixed

- `config.json`, `petk/config.json`, `petk/Demo/config.json` now set
  `pore_geometry.distance_metric` (and `chamfer_depth` for the cylindrical
  templates) explicitly, rather than relying on the implicit default, so
  the example configs document the geometry convention they run under.
