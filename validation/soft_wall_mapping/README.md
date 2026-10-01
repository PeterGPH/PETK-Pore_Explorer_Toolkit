# Soft-wall → sharp-wall mapping of SEM solid-state pores

The SEM conductivity map has no sharp wall: σ ramps linearly from 0 at
`r_min = 1.3 Å` to bulk at `r_cut = 4.1 Å` in the Euclidean distance to the
membrane (`sem.utils.condfrac`). Comparing an SEM open-pore conductance with
sharp-walled theory therefore needs two things:

1. the sharp geometry that is electrically equivalent to the ramped map, and
2. a sharp-walled formula that is accurate in its own right.

The manuscript draft did (1) with the area-conserving offset
`(d − 0.54 nm, L + 0.54 nm)`, and (2) with Hall/Kowalczyk,
`R = ρ (4L/πd² + 1/d)`. The "required offset" then drifted from 0.52 nm at
d = 5 nm to 0.95 nm at d = 25 nm (L = 20 nm). This directory shows where that
drift comes from and replaces the mapping with one that holds to 0.03 %.

## Findings

All numbers come from `sem.axisym_conductance`, an axisymmetric Q1 finite-element
solver with unbounded reservoirs (exact far-field Robin condition). Its conductance
is grid-converged to ~1e-5, and it reproduces the thin-aperture limit (R = ρ/2a)
and the flanged-tube end correction (0.8216 a).

| Piece of the pore | What the ramp does | Treatment |
|---|---|---|
| Bore | Field is exactly axial, so the layers conduct in parallel | Equal-area radius `a_e` is **exact** (`a − a_e` = 2.685–2.700 Å) |
| Faces / access region (the radial "dr" part) | Current runs tangentially along the face | Same deficit length t = (r_min + r_cut)/2 = 2.7 Å; L → L + 2t, exact by translation invariance |
| Rim | Euclidean distance rounds the rim (radius ≈ t) | End-length deficit `Δℓ = c(λ) t^{4/3} a_e^{-1/3}`; c ≈ 0.53 for L ≳ a; 0.71 Å per end at d = 5 nm |
| **Theory** | Hall's access term ρ/(4a) per side is the **thin-plate** limit | Exact end correction `κ(L/a)·a`: π/4 at L→0, 0.82155 for L ≳ a, so Hall's access term is **4.6 % small** |

The drift seen against Hall/Kowalczyk is therefore Hall's access error, not the
soft wall. The access term's share of the resistance grows with d, so the error
does too: Hall overestimates G by 0.8 % at L = 4d, 2.3 % at L = 0.8d, and at most
3.2 % at L ≈ 0.2d.

Combined (`sem.open_pore_theory.sem_cylinder_resistance`):

```
R_SEM = ρ [L_e + 2 κ(L_e/a_e) a_e] / (π a_e²),   L_e = L + 2t − 2Δℓ
```

| Reference vs converged SEM ramp (G_ref/G_SEM − 1) | d = 2–50 nm, L = 1–40 nm | L = 20 nm, d = 5 → 25 nm |
|---|---|---|
| Hall at `(d − 0.54, L + 0.54)` (draft) | −6.5 … +3.1 % | −0.04 → +2.02 % |
| Exact κ, equal-area `a_e`, `L + 2t`, no rim term | −7.4 … −0.08 % | −0.59 → −0.20 % |
| **Full mapping** | **max 0.032 %, rms 0.007 %** | **≤ 0.005 %** |

Additional checks:

- **Other ramps.** The rim coefficient depends on the ramp only through
  `(r_cut − r_min)/t`; it is calibrated on one thick geometry per ramp shape.
  On (d, L) grids that were not part of any fit, for r_cut = 3–6 Å and
  r_min = 0–2.5 Å, the full mapping holds to ≤ 0.015 % for d ≥ 3 nm.
- **d = 1.5 nm** (outside the fit): 0.18 %.
- **Golden-aspect box** (paper's domain rule, equal-area radius `(d + 2p)/√π`):
  the box changes G by −0.24 % to +0.54 % from the infinite-reservoir value for
  p/d = 0.35–2 and d = 5–25 nm. That is why the draft's SEM sits 2.2 % below
  Hall-offset at d = 25 nm, p/d = 0.35 (draft: "worst case 2.4 %").

## Files

| File | Content |
|---|---|
| `kappa.csv` | κ(L/a) of a sharp cylinder (Richardson-extrapolated, error estimate, closed form) |
| `ramp_grid.csv` | Converged SEM-ramp resistance on the (d, L) grid and the deviation of each reference |
| `golden_box.csv` | Finite golden-aspect box vs infinite reservoirs |
| `soft_wall_mapping.{pdf,png}` | Summary figure (A: κ; B: deviations at L = 20 nm; C: Hall-offset error map; D: rim collapse) |
| `manuscript_soft_wall_mapping.tex` | Drop-in replacement for Appendix C, the Results 3.1 validation paragraph, and a list of other draft sentences this changes |

## Reproduce

```bash
python -m sem.scripts.soft_wall_mapping all --out validation/soft_wall_mapping --jobs 4   # ~10 min
pytest tests/test_open_pore_theory.py                                                     # ~10 s
```

## Using the mapping

```python
from sem import open_pore_theory as opt
R = opt.sem_cylinder_resistance(d=100.0, L=200.0, sigma=10.5)   # Å, S/m -> Ohm
d_e, L_e = opt.effective_cylinder(100.0, 200.0)                  # equivalent sharp pore (Å)
d_sem, L_sem = opt.sem_cylinder_for_sharp_pore(100.0, 200.0)     # SEM inputs that conduct like a sharp 10 nm pore
```

Scope: straight cylinders with the Euclidean wall distance. The same
construction applies to tapered walls, but each needs its own exact reference.
`sem.axisym_conductance` solves any piecewise-linear axisymmetric wall, with
either the Euclidean or the horizontal (`a(z) − r`) distance used by the
tapered pores on `main`.
