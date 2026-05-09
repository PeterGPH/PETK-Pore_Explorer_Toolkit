# Contributing to PoreExplorer-PETK

Thanks for your interest in contributing. This is a research codebase that
pairs a VMD/Tcl GUI (`petk/`) with a Python finite-element solver (`sem/`),
so changes can land in either layer.

## Branching model

We use a simple trunk-based flow: `main` is always green, all work happens
on short-lived feature branches that land via pull request.

```
main                           ← always green, PR-only
  ↑ pull request, CI green
feature/<short-slug>           ← one feature per branch
fix/<short-slug>               ← bug fixes
docs/<short-slug>              ← docs-only changes
chore/<short-slug>             ← build, CI, dependencies
refactor/<short-slug>          ← restructuring without behavior change
```

Branch names use the prefix that matches the kind of change. Slugs are
2–4 hyphenated words.

## Setting up a development environment

```bash
git clone https://github.com/PeterGPH/PETK-Pore_Explorer_Toolkit.git
cd PETK-Pore_Explorer_Toolkit
conda env create -f environment.yml
conda activate sem-env
pip install -e .[dev]   # editable install + dev extras (pytest, flake8)
```

Verify the install:

```bash
sem --help
pytest tests/
```

## Workflow

1. **Sync `main`.** `git checkout main && git pull origin main`.
2. **Create a feature branch.** `git checkout -b feature/<slug>`.
3. **Make focused commits.** One logical change per commit; commit message
   summarizes *what* changed in the imperative ("Add golden-aspect rule for
   FEM box sizing"), with a body explaining *why* if not obvious.
4. **Run tests locally.** `pytest tests/` — must pass.
5. **Push and open a PR.** `git push -u origin feature/<slug>`, then open
   a pull request against `main` from the GitHub UI. CI will run
   automatically.
6. **Address review comments.** Push additional commits to the same
   branch — they show up in the PR.
7. **Merge.** Once CI is green and review is resolved, the PR can be
   squash-merged. The feature branch should then be deleted.

## Commit messages

We don't enforce Conventional Commits, but follow a similar pattern:

```
<short summary in imperative mood, ~60 chars>

Optional longer paragraph explaining the why, side effects, or
follow-up work needed. Wrap at 72 chars. Reference issues with #N.
```

Avoid commits like `fix stuff`, `wip`, or `address review` — squash them
locally before pushing or rely on the squash-merge to consolidate.

## Tests

- Tests live under `tests/`. Add a test alongside any behavior change.
- Mark tests that need DOLFINx with a clear marker so CI's pip-only
  smoke job can skip them:

  ```python
  pytest.importorskip("dolfinx")
  ```

- Long-running benchmarks belong in `sem/scripts/`, not in `tests/`.

## Code style

- Python: follow PEP 8 within reason; we run `flake8` in CI for
  syntax/undefined-name checks (`E9, F63, F7, F82`).
- Tcl: match the surrounding `petk/*.tcl` style — 4-space indent,
  `::PETK::gui::` namespace.
- Type hints encouraged for new Python code, not required for older code.

## Reporting issues

Use the GitHub issue tracker. Useful issues include: a minimal reproducer,
the SEM config (or a redacted version), the relevant log tail, and the
SEM/DOLFINx versions you're running on.
