# Overlap-Skipping Test Protocol for Solid-State Nanopores

**Scope.** A standardized, geometry-agnostic test suite that quantifies how well a
solid-state nanopore (cylindrical, double-cone / hourglass, conical) "skips"
overlap from analytes during translocation sensing. The protocol is built on
the PETK / SEM stack, but the metrics are pore- and tool-independent.

The protocol mirrors the analogous tests used for biological nanopores
(MspA, α-hemolysin, CsgG, aerolysin) so results from solid-state and biological
pores can be compared on the same axes.

---

## 1. Operational definition of "overlap"

Three regimes need to be tested separately. They have different physical origins
and different mitigation strategies.

| # | Name | What overlaps | Biological-pore analogue | Primary lever |
|---|------|---------------|--------------------------|---------------|
| A | Spatial (read-head) overlap | Multiple sub-features of a single analyte (k-mer effect) | MspA constriction (0.6 nm) vs α-HL barrel (5 nm) | Pore geometry, sensing-zone length L_s |
| B | Temporal (coincidence) overlap | Two or more analytes in the pore at the same time | Event pile-up at high concentration | Concentration, voltage, deconvolution algorithm |
| C | Signal-class overlap | Distinct analytes producing indistinguishable blockades | Homopolymer ambiguity in DNA seq | Discrimination feature space |

A "complete" test reports each separately and gives a single composite score.

---

## 2. Pore geometries under test

The PETK shape primitives needed:

- **Cylindrical** — radius `r`, length `L` (membrane thickness). One
  geometric scale.
- **Conical** — tip radius `r_tip`, base radius `r_base`, length `L`. Half-angle
  `θ = atan((r_base − r_tip) / L)`.
- **Double-cone (hourglass)** — same as two conical pores joined at the
  constriction `r_c`, with full length `L = 2 L_half`. Symmetric around z = 0.

All three are parameterised in `petk/config.json` via `pore_type`,
`pore_radius`, `corner_radius`, `membrane_thickness`. For double-cone you'll
need to add a `constriction_radius` field and a `cone_angle` (or two-radii)
field — note this in the changelog if you introduce it.

For each geometry, run the **same nominal minimum diameter** (e.g., 2 r_min =
2.0 nm) so the open-pore conductance is similar and the comparison is fair.

---

## 3. Pre-test characterisation (must run first)

These four numbers are reported for every pore before any analyte is used.
Without them you cannot interpret the overlap-skip metrics.

1. **Open-pore conductance G₀** at the test voltage (e.g., 200 mV). Sanity
   check: cross-check against the analytical Hall formula
   G₀ = σ · (4 L / (π d²) + 1 / d)⁻¹ for the cylindrical case. Use this to
   calibrate `bulk_conductivity` in the SEM config.
2. **Sensing-zone length L_s** — defined as the FWHM of the conductance
   sensitivity profile S(z) = |dG/dz| obtained by translating a hard, neutral
   point sphere (r = 1 Å) along the central z-axis from −L to +L. This is the
   single most important geometric number for spatial-overlap behaviour.
3. **Sensing-zone asymmetry** — for conical pores, report
   skew = (z_peak − z_center) / L. Cylindrical → 0, conical → ~±0.3, double-cone
   → 0 by construction.
4. **Edge-effect length** — the z-range over which |dG/dz| > 5% of its
   maximum, outside the membrane. Reports how far the analyte "is felt"
   before/after entering.

Implementation note: PETK already produces `sem_results_plot_conductance.png`
from a z-scan; the same scan with a 1 Å probe gives S(z). Add a
`scripts/compute_sensitivity_profile.py` if it doesn't exist.

---

## 4. Test A — Spatial (read-head) overlap

### A.1 Point-Spread Function

A single hard probe of effective radius r_p (run for r_p ∈ {1, 2, 4, 8} Å) is
translated along z. Record ΔG(z) = G₀ − G(z). Fit a Gaussian to extract:

- mean position μ
- width σ (FWHM = 2.355 σ)
- amplitude A

**Metric A.1**: PSF FWHM, normalised by minimum pore diameter:
`L_s* = FWHM / d_min`. Smaller is better. MspA reference value: ~0.6.
α-HL reference value: ~3.0.

### A.2 Two-Point Rayleigh Resolution

Two identical probes are placed on the z-axis with separation Δz, both
translated rigidly. Record combined ΔG(z; Δz). Define the resolved threshold
using the Sparrow criterion (the smallest Δz at which the second derivative
of the central minimum is zero) — Rayleigh is fine too but Sparrow is
geometry-tolerant.

**Metric A.2**: minimum resolvable separation R_min, in Å, at the midplane.
Also report R_min(z) at z = ±L/4 to capture conical asymmetry.

Expected ordering for a 2 nm minimum diameter, 30 nm membrane:
double-cone < conical < cylindrical (smaller R_min = better overlap-skip).

### A.3 Tandem-Bead k-mer Test

A linear chain of 12 distinguishable beads (alternating "heavy" / "light",
encoded by van der Waals radius differences of 0.5 Å) is dragged through the
pore at constant velocity. Record the current trace.

**Metric A.3**: effective k-mer N_eff = L_s / L_bead, where L_bead is the
bead-bead spacing. The fraction of the sequence that can be unambiguously
basecalled (with a simple Viterbi decoder, see §7) is reported as the k-mer
recovery rate.

This is the direct analogue of the biological k-mer benchmark and is the
metric that should be highlighted in any cross-pore comparison plot.

---

## 5. Test B — Temporal (coincidence) overlap

### B.1 Pile-up Probability

This is partly geometry-dependent (residence time τ scales with L and
inversely with E-field gradient) and partly concentration-dependent. From the
Smoluchowski capture rate

  R_capture = 2 π D r_eff c

and the residence time τ ≈ L_eff / v_drift, the dimensionless pile-up
parameter is

  Π = R_capture · τ

**Metric B.1**: Π at a reference concentration (say c = 1 µM) and reference
voltage (200 mV). Lower is better. Conical pores typically have shorter
effective residence (because the drag drops away from the tip) and hence
lower Π for the same throughput.

### B.2 Two-Event Deconvolution Test

Synthetic traces are generated by superposing two single-analyte events with
controlled inter-arrival time Δt_a. For each Δt_a, run a deconvolution
algorithm (matched filter, then Viterbi using the PSF from §A.1 as the kernel)
and measure detection accuracy.

**Metric B.2**: t* = smallest Δt_a at which both events are recovered with
≥95% accuracy. Report normalised by the open-pore RC time: t*/τ_RC.

### B.3 Concentration sweep

Repeat the live-translocation simulation at c ∈ {0.1, 1, 10} µM and report
the recovered event count vs the true count. The slope of the
"missed-events fraction vs Π" curve is the practical pile-up resilience.

---

## 6. Test C — Signal-class (discrimination) overlap

This is conceptually a classification problem, but it interacts with geometry:
a longer L_s averages out features that distinguish classes.

### C.1 Class library

A standard library of 10 analytes spanning the dimensions you care about. For
peptide-style work (1AOI is in your demo, so you're already on this track),
use the 20 canonical amino acids each as a free monomer. For nucleotide work,
use the 4 dNMPs.

For each analyte, run 100 random orientations through the pore at the same
z-trajectory and compute the (ΔG_max, dwell-time, ΔG-area) feature triple.

### C.2 Class-overlap metric

Use the Davies–Bouldin index, or equivalently the silhouette score, on the
3-D feature cloud. Lower DB = better class separation = lower signal-class
overlap.

**Metric C.1**: DB index across the 20-AA cloud.
**Metric C.2**: pairwise confusion matrix (20×20). Report the maximum
off-diagonal entry — the worst-confused pair.

### C.3 Orientation-averaged sensitivity

Because real translocation samples orientations stochastically, also compute
the standard deviation of ΔG_max across orientations for a single analyte.
This is an irreducible noise floor the algorithm has to fight.

---

## 7. Algorithmic toolkit (borrowed from biological pores)

These are the deconvolution and classification methods that biological-pore
people have validated, and that map cleanly onto solid-state data.

| Problem | Biological-pore tool | Solid-state port |
|---|---|---|
| k-mer / spatial overlap | nanopolish HMM, scrappie/bonito CNN-RNN | 1-D CNN on simulated tandem-bead traces; ground truth from PETK |
| coincidence overlap | matched-filter + Viterbi (Pedone et al.) | same; PSF kernel from §A.1 |
| signal-class overlap | random-forest on event features (Wendell et al.) | same |
| translocation-speed control | helicase, Hel308, T4 gp43 | not directly portable — substitute pressure/dielectrophoretic trapping (Yusko, Wanunu) |

The single most useful imported idea is the **k-mer HMM**: treat the current
trace as emissions from a hidden state that is the analyte's z-coordinate
(or sequence index), with emission means given by the PSF convolution. The
Viterbi backtrace gives the most-likely position sequence, and the per-state
posterior gives a confidence. This is exactly the framework that lets MspA
basecall with k=3–5 even though no individual base is delta-resolved.

For a solid-state pore this means: even a long-L_s cylindrical pore can be
made to "skip overlap" *algorithmically*, provided the PSF is well-known and
stationary. The point of the geometric metrics in §4 is to tell you how much
work the algorithm has to do.

---

## 8. Reporting template

Every pore tested should produce one row of the following table (use the
existing `petk/Demo/results/` layout).

```
geometry        d_min   L      L_s/d   R_min    N_eff   Π       t*/τ_RC   DB
cylindrical     2.0 nm  30 nm  3.1     1.5 nm   8.0     0.18    1.4       2.1
double-cone     2.0 nm  30 nm  0.9     0.5 nm   2.4     0.07    0.8       0.9
conical         2.0 nm  30 nm  1.7     0.9 nm   4.3     0.10    1.0       1.4
α-HL (ref)      1.4 nm  10 nm  3.6     1.4 nm   5.1     —       —         —
MspA (ref)      1.2 nm  9.6 nm 0.6     0.4 nm   1.8     —       —         —
```

A single composite score can be computed as

  S_overlap = w₁ · (1/L_s*) + w₂ · (1/R_min) + w₃ · (1/N_eff)
            + w₄ · (1/Π)   + w₅ · (1/DB)

with weights chosen for the application (sequencing → emphasise A; counting
→ emphasise B; identification → emphasise C). State the weights explicitly.

---

## 9. Validation against published data

Before claiming a geometry "skips overlap" well, verify the protocol on cases
with a known answer:

1. Re-derive Hall's open-pore formula for a cylinder (Test 3.1).
2. Reproduce the experimental L_s ≈ 0.6 nm for MspA from Manrao et al. 2012
   (Test A.1 on `bio_pore/centered_3X2R.pdb` if MspA is added — currently
   the bio_pore folder has α-HL 7AHL, 1UUN, 6MRT, 3X2R; check which of those
   is MspA before running).
3. Reproduce Plesa et al. 2013 pile-up curves on a SiN cylindrical pore at
   r = 5 nm, 200 mV, c = 1 µM (Test B.3).

Publish (in supplementary) the residual between PETK predictions and these
three benchmarks. If any residual exceeds 20% the protocol is not yet
calibrated and conclusions about geometry comparisons should be deferred.

---

## 10. Implementation checklist (PETK-specific)

The following scripts need to exist (some already do, some are new):

- [x] `sem/__main__.py` — z-scan driver (already there)
- [ ] `sem/scripts/compute_sensitivity_profile.py` — S(z), L_s
- [ ] `sem/scripts/two_point_resolution.py` — R_min sweep
- [ ] `sem/scripts/tandem_bead_translocation.py` — k-mer trace + Viterbi
- [ ] `sem/scripts/coincidence_pileup.py` — Π and t*
- [ ] `sem/scripts/class_overlap.py` — DB, confusion matrix
- [ ] `sem/scripts/run_overlap_battery.py` — orchestrator that runs A+B+C
       and emits the §8 reporting table as JSON + Markdown

The `ChargeAwareConductivityModel` already in `sem/conductivity_models.py`
should be exercised in addition to `SimpleConductivityModel` for §A.1 and
§A.2 — Debye screening expands the effective L_s and the screening
length must therefore be reported alongside every L_s value.

---

## 11. Caveats specific to solid-state pores

1. **Translocation control is the missing piece.** Biological pores have
   helicases that ratchet DNA at ~450 bp/s; solid-state pores translocate
   at 1–10 µs/base. The k-mer HMM assumes a roughly known velocity, so any
   solid-state port has to either (a) include velocity as a hidden variable,
   or (b) introduce a separate slowing mechanism (pressure, double-pore
   trapping, DNA origami capture). State which assumption the test makes.
2. **Drift.** Solid-state pores enlarge under high voltage and current; L_s
   drifts. The full protocol must be run at the *start* and *end* of an
   experimental session (in simulation, repeat with a 10% pore-radius
   inflation) and the metrics quoted as bands, not single numbers.
3. **Surface charge.** The PSF for a charged pore depends on ionic strength.
   Run §4 at three ionic strengths (e.g., 0.1, 1, 3 M KCl) and report L_s
   for each. This is the single biggest difference between a SiN/SiO₂ pore
   and a (largely uncharged) MspA constriction.
4. **Asymmetry.** Conical pores rectify; the "approach" and "exit" PSFs
   differ. Report S(z) in both directions of translocation.

---

*Last revision: 2026-05-06.*
