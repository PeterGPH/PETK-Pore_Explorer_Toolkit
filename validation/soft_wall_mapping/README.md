# Soft-wall → sharp-wall mapping of SEM solid-state pores

The SEM conductivity map has no sharp wall: σ ramps linearly from 0 at
`r_min = 1.3 Å` to bulk at `r_cut = 4.1 Å` in the Euclidean distance to the
membrane (`sem.utils.condfrac`). Comparing an SEM open-pore conductance with
sharp-walled theory therefore needs two things:

1. the sharp geometry that is electrically equivalent to the ramped map, and
2. a sharp-walled formula that is accurate in its own right.

The manuscript draft did (1) with the area-conserving offset
`(d − 0.54 nm, L + 0.54 nm)`, and (2) with Hall/Kowalczyk,
`R = ρ (4L/πd² + 1/d)`. The offset needed for agreement then drifted from
0.52 nm at d = 5 nm to 0.94 nm at d = 25 nm (L = 20 nm). This directory
shows where that drift comes from and replaces the mapping with one that has
no parameter fitted to SEM solutions.

## Findings

Three independent references back the numbers below:

- `sem.axisym_conductance`: an axisymmetric Q1 finite-element solver with
  unbounded reservoirs. It is Richardson-extrapolated to ~1e-6 and agrees with
  an independent finite-volume solver to ≤ 8e-7.
- `sem.cylinder_mode_matching`: an exact Ritz/mode-matching solution of the
  sharp cylinder (variational, converged to 1e-9).
- A matched-asymptotics analysis of the ramp at walls and corners.

| Piece of the pore | What the ramp does | Treatment |
|---|---|---|
| Bore | Away from the ends the field is axial, so the layers conduct in parallel | Equal-area radius `a_e² = (a − t)² + Δr²/12` is **exact**: t is the ramp mean, Δr²/12 its variance |
| Faces / access region (the radial "dr" part) | Current runs tangentially along the face | The same deficit length t = (r_min + r_cut)/2 = 2.7 Å. No correction at second order for any ramp profile; the first "series" (normal) term is third order. L → L + 2t |
| Rim | The Euclidean distance rounds the rim (radius ≈ t) | End-length deficit `Δℓ = C (k/k∞)² t^{4/3} a_e^{−1/3}`. This is the 270° corner law: C = −2πβk∞² = 0.535 from a local corner constant β, and k is the rim strength of the sharp cylinder. Δℓ is 0.71 Å per end at d = 5 nm, 0.40 Å at 25 nm |
| **Theory** | Hall's access term ρ/(4a) per side is the **thin-plate** limit | Exact end correction `κ(L/a)·a`: π/4 as L → 0, rising to **0.82167** (flanged tube) by L ≈ a. Hall's access term is therefore **4.6 % small** |

The drift seen against Hall/Kowalczyk is Hall's access error, not the soft wall.
The access term's share of the resistance grows with d, so the error grows too.
For a sharp pore, Hall overestimates G by 0.8 % at L = 4d, 2.3 % at L = 0.8d,
and at most 3.2 % at L ≈ 0.2d.

At d = 5 nm the draft's ~0 % agreement is a cancellation of three effects:

| Effect | Contribution to G |
|---|---|
| Hall access | −0.67 % |
| Equal-area variance term | +0.12 % |
| Rim | +0.59 % |
| **Net** | **−0.04 %** |

How the "required offset" grows with d depends on the convention.

- **Continuum, infinite reservoirs:** it reaches 0.77 nm at d = 25 nm when
  applied to both d and L, or 0.86 nm when applied to d alone.
- **With the golden-aspect box:** the box adds ≈ 0.04 nm, giving ≈ 0.90 nm,
  which accounts for most of the reported 0.94 nm.

Combined (`sem.open_pore_theory.sem_cylinder_resistance`), with λ = (L + 2t)/a_e:

```
R_SEM = ρ [L + 2t + 2 κ(λ) a_e − 2 Δℓ] / (π a_e²)
```

Accuracy against the converged SEM solves (G_ref/G_SEM − 1, max |error|):

| Reference | d = 2–50 nm, L = 1–40 nm | L ≥ 20 nm | L = 20 nm, d = 5 → 25 nm |
|---|---|---|---|
| Hall at `(d − 0.54, L + 0.54)` (draft) | −6.5 … +3.1 % | −1.9 … +2.8 % | −0.04 → +2.02 % |
| Exact κ, equal-area `a_e`, `L + 2t`, no rim term | −7.4 … −0.08 % | −0.94 … −0.08 % | −0.59 → −0.20 % |
| **Full mapping (no parameter fitted to SEM)** | **0.12 %** (rms 0.016 %) | **0.014 %** | **≤ 0.0015 %** |

The full mapping's residual sits in the thinnest membranes (L = 1–3 nm: 0.12 %;
L = 5–15 nm: 0.05 %), where L + 2t is only a few t and the corner expansion is
strained. At d = 1.5 nm, outside the asymptotic regime, it is 0.30 %.

Additional checks:

- **Other ramps** (`other_ramps.csv`; r_cut = 3–6 Å, r_min = 0–2.5 Å). The only
  ramp-shape input is the thick-membrane rim constant C(Δr/t), which is
  calibrated on one geometry per shape. On these (d, L) points the mapping holds
  to ≤ 0.03 % for d ≥ 3 nm. Hall's offset errs by −11 … +2.9 % on the same points.
- **Corner theory vs solve.** The thick-membrane solve gives C = 0.5349 for the
  default ramp; the corner theory, −2πβk∞² with β = −0.1494 and k∞ = 0.75485,
  gives 0.5349.
- **Golden-aspect box** (`golden_box.csv`; the paper's domain rule, equal-area
  radius `(d + 2p)/√π`). The box changes G by −0.24 % to +0.54 % from the
  infinite-reservoir value for p/d = 0.35–2 and d = 5–25 nm. At d = 25 nm,
  p/d = 0.35 the draft's SEM therefore sits 2.2 % below Hall-offset (draft:
  "worst case 2.4 %").

## Files

| File | Content |
|---|---|
| `sharp_cylinder.csv` | Exact κ(L/a) and rim strength k(L/a) of a sharp cylinder, with the closed forms used by `open_pore_theory` |
| `rim_shape.csv` | Thick-membrane rim coefficient C for linear ramps of different shape |
| `ramp_grid.csv` | Converged SEM-ramp resistance on the (d, L) grid, and the deviation of each reference |
| `other_ramps.csv` | The same for r_cut = 3–6 Å, r_min = 0–2.5 Å (held out) |
| `golden_box.csv` | Finite golden-aspect box vs infinite reservoirs |
| `soft_wall_mapping.{pdf,png}` | Summary figure (A: κ; B: deviations at L = 20 nm; C: Hall-offset error map; D: rim term vs corner theory) |
| `manuscript_soft_wall_mapping.tex` | Drop-in replacement for Appendix C, the Results 3.1 validation paragraph, and a list of other draft sentences this changes |

## Reproduce

```bash
python -m sem.scripts.soft_wall_mapping all --out validation/soft_wall_mapping --jobs 4   # ~1 h
pytest tests/test_open_pore_theory.py                                                     # ~20 s
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
