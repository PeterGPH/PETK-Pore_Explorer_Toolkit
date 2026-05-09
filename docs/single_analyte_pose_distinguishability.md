# Single-Analyte Pose Distinguishability Test

**Question.** For one analyte translocating and rotating through a solid-state
nanopore (cylindrical / double-cone / conical), is the steric overlap between
analyte and pore distinguishable from the ionic-current signal alone? In
formal terms: is the map

  Φ : (z, R) ∈ ℝ × SO(3)  →  G ∈ ℝ

injective above the measurement noise floor σ_G, and how does the answer
depend on pore geometry?

This complements the multi-analyte protocol in
`overlap_skip_test_methodology.md`. There the question was about resolving
multiple molecules; here the question is about resolving the *pose* of one.

---

## 1. Variables

- z ∈ [−L, L] — centre-of-mass z position of the analyte (1 Å step is fine
  for a 30 Å analyte; coarser for larger).
- R ∈ SO(3) — analyte orientation, parameterised in PETK by
  `RotationSpec(rx, ry, rz)` in degrees (Rx → Ry → Rz convention; rotations
  about the grid centre).
- V_overlap(z, R) — geometric steric overlap between the analyte's vdW
  volume and the pore's solid-volume complement (this is the *ground truth*
  that the current is supposed to encode).
- G(z, R) — predicted conductance from SEM (`SimpleConductivityModel` or
  `ChargeAwareConductivityModel`).
- σ_G — noise floor. From either the grid-convergence residual already
  produced by `scripts/plot_open_pore_convergence.py` or an experimental
  RMS noise level (typical SiN, 100 kHz: σ_G ≈ 50–200 pS).

The distinguishability question is whether two distinct poses (z₁, R₁) ≠
(z₂, R₂) ever satisfy |G(z₁, R₁) − G(z₂, R₂)| < σ_G — i.e., are *confusable*
— and how many such confusable pairs exist for each geometry.

---

## 2. Pre-test: define the steric-overlap ground truth

For every (z, R) sample, compute V_overlap by intersecting the vdW
representation of the analyte with the pore solid. This is independent of
the conductivity model and gives the *physical* overlap that the test is
trying to recover.

Implementation (drop into `sem/scripts/compute_steric_overlap.py`):

```python
def compute_steric_overlap(analyte_xyz, vdw_radii, pore_grid, dx_step):
    """
    Returns V_overlap in Å^3.
    pore_grid: 3-D boolean grid, True inside the pore solid (membrane).
    """
    overlap = 0.0
    for atom_xyz, r in zip(analyte_xyz, vdw_radii):
        # Voxelise the atom's vdW sphere on the pore grid.
        atom_mask = sphere_mask(atom_xyz, r, pore_grid.shape, dx_step)
        overlap += np.sum(atom_mask & pore_grid) * dx_step ** 3
    return overlap
```

The protocol is to compute V_overlap **for the same (z, R) grid** as the
SEM run, and store both arrays alongside each other. Every plot that uses
G should also be plottable against V_overlap as the x-axis.

---

## 3. Five core tests

### Test 1 — Pure translocation, fixed orientation

Hold R = R₀ (a canonical orientation: e.g., principal moment of inertia
aligned with z). Sweep z from −L to +L at the chosen step. Record G(z; R₀)
and V_overlap(z; R₀).

**Metrics:**
- ΔG_max = G₀ − min G(z) — depth of blockade
- z_resolution = number of z values for which G(z) is uniquely recoverable
  given σ_G. Compute as the number of distinct bins after applying a
  σ_G/|dG/dz| Voronoi smoothing.
- Monotonicity test: is G(z) a monotone-then-monotone function (single
  minimum) for this geometry? Cylindrical → flat plateau interior, sharp
  edges; conical → asymmetric single minimum at apex; double-cone →
  symmetric single minimum at constriction.

**Expected ordering of z-distinguishability** (best → worst):
double-cone > conical > cylindrical, because the sensitivity peak is
sharpest where the geometry has the steepest radius gradient.

### Test 2 — Pure rotation, fixed z

Hold z = z₀ at three reference positions: z₀ ∈ {−L/4, 0, +L/4}. For each
z₀, sample SO(3) on a uniform grid (e.g., 18°×18°×18° → 8000 orientations,
or use Sobol sequences on SO(3) for ~500 quasi-random orientations if
runtime is tight).

Record G(R; z₀) and V_overlap(R; z₀).

**Metrics:**
- σ_R(G) = standard deviation of G over orientation. This is the
  *orientation noise* that any classification has to fight.
- Rotation-resolved bin count: number of distinguishable G values at the
  given z₀, given σ_G.
- Symmetry residuals:
  - Cylindrical: G should be invariant under rotation about z. Compute
    σ_φ(G) where φ = rz (the axial angle). Should be ≤ σ_G if the
    geometry is correctly cylindrically symmetric. Any larger residual
    indicates either (a) grid anisotropy, or (b) a genuinely
    non-axisymmetric analyte response that the cylinder still can't
    distinguish.
  - Double-cone: G(R) = G(P R) where P is z-mirror at z = 0. Test by
    comparing the rotation cube at +L/4 with that at −L/4 with rx → −rx.
    Residuals quantify any unintended mirror-asymmetry.
  - Conical: no symmetry constraint, but the (+L/4, R) cube and (−L/4, R)
    cube should differ significantly. Compute the Frobenius distance and
    report it.

### Test 3 — Joint (z, R) injectivity

This is the actual distinguishability test. Sample (z, R) jointly, e.g.,
21 z-values × 500 Sobol-on-SO(3) orientations = 10500 poses per pore.

Build the empirical distribution of pairwise distances |ΔG_ij| over all
i < j pairs. Define the **confusion fraction**:

  C(σ_G) = #{(i, j) : i < j, |G_i − G_j| < σ_G} / N_pairs

This is a single scalar between 0 and 1. C = 0 means the pose is perfectly
recoverable; C = 0.5 means half of all pose pairs are indistinguishable.

Plot C(σ_G) vs σ_G on a log-log axis for each geometry. The slope and the
intercept (at the experimental σ_G) constitute the geometry's
distinguishability fingerprint.

**Metric:** C(σ_G = 100 pS) — single number, lower is better. Expected
ordering: double-cone < conical < cylindrical.

### Test 4 — Information capacity (mutual information)

Treat z, R, V_overlap each as a hidden variable and G as the observed
variable. Estimate

  I(G ; z),  I(G ; R),  I(G ; V_overlap)

using the Kraskov k-NN estimator (`sklearn.feature_selection.mutual_info_regression`
works for 1-D z; for R use the Kozachenko–Leonenko on the 3-sphere).

**Metrics (all in bits):**
- I(G; z) — translocation-position capacity
- I(G; R) — orientation capacity (this is the part that biological
  pores generally lose because of rotational averaging during fast
  translocation; solid-state pores can in principle recover it if the
  analyte is trapped)
- I(G; V_overlap) — steric-overlap recoverability. If this is close to
  H(V_overlap) then the conductance is essentially a function of the
  overlap volume alone, and pose-distinguishability is bounded by how
  many distinct (z, R) map to the same V_overlap.

A pore "distinguishes overlap well" if I(G; V_overlap) is high *and*
I(G; z) + I(G; R) > I(G; V_overlap) — i.e., the conductance encodes more
than just the scalar overlap volume; it encodes geometric details of *how*
the analyte overlaps.

### Test 5 — Linearised diagnostic: G vs V_overlap

Plot G against V_overlap across all (z, R) samples. Compute:
- Pearson r — how scalarly the conductance tracks the overlap.
- Local Lipschitz: for every pair within ΔV_overlap < δ, what is the
  spread of ΔG? If wide, the conductance carries shape information beyond
  the volume.
- Residual after a 1-D fit G = a + b · V_overlap + ε. The residual ε(z, R)
  is the part of the signal that *is* pose information beyond overlap.

A geometry that distinguishes overlap *well* shows tight scatter around
the G(V_overlap) trend (high Pearson r) with informative residuals — the
trend tells you the volume, the residuals tell you the pose.

A geometry that distinguishes overlap *poorly* either shows large scatter
(volume itself isn't recovered) or has residuals indistinguishable from
σ_G (no pose information beyond volume).

---

## 4. Reporting template

One row per geometry × analyte:

```
geometry        analyte   ΔG_max  z_bins  σ_R(G)/G₀  C(100pS)  I(G;z)   I(G;R)   I(G;V)   r(G,V)
cylindrical     1AOI      18%     12      0.3%       0.41      2.1 b    0.4 b    3.0 b    0.91
double-cone     1AOI      32%     28      4.1%       0.08      4.6 b    2.2 b    3.5 b    0.78
conical         1AOI      24%     19      2.7%       0.17      3.4 b    1.6 b    3.3 b    0.82
```

The interpretation:

- A high I(G; R) on the double-cone with low Pearson r means the pore is
  responding to *shape* of overlap, not just *volume* — it can
  distinguish poses that share V_overlap but differ in atomic distribution.
- A high Pearson r with high I(G; V) but low I(G; R) (cylindrical) means
  the pore is essentially a Coulter counter — overlap volume in, signal
  out, no pose information.

---

## 5. Sampling guidance

**Don't sample SO(3) by Euler grids.** Uniform Euler grids over-sample the
poles and under-sample the equator. Two acceptable alternatives:

1. **Quasi-uniform on SO(3) via SU(2) → unit quaternions, then Sobol on
   the 3-sphere.** ~500 poses gives ~12° mean nearest-neighbour spacing.
2. **Hopf fibration / HEALPix-on-SO(3) (Yershova 2010).** Produces an
   exactly uniform grid. Use ~648 (level-1) or ~5184 (level-2) poses.

For a peptide like 1AOI with an obvious long axis, also include a "dense
strip" of orientations near the canonical axial alignment ±20°, where the
signal varies fastest with orientation.

For z, sample more densely in the entrance/exit edge region (where
|dG/dz| is largest) and coarsely in the open-pore plateau.

---

## 6. Symmetry sanity checks

Before reporting any C(σ_G), verify the following:

| Pore | Invariance | Expected residual |
|---|---|---|
| Cylindrical | G(rz) = G(0) ∀ rz | ≤ σ_G |
| Cylindrical | G(z, rx, ry) = G(z, −rx, −ry) | ≤ σ_G |
| Double-cone | G(z) = G(−z) (with R fixed) | ≤ σ_G |
| Conical | (no useful symmetry) | — |

Failures here usually mean grid anisotropy from the
`scipy.ndimage.rotate`-based pipeline. Fix by aligning the grid axis to
the pore axis exactly, using a rotation-invariant analyte representation
(spherical-harmonic expansion of the vdW indicator) for the symmetry
checks, or by averaging G over a small azimuthal range.

---

## 7. Distinguishing "overlap distinguishable" from "overlap recoverable"

These are not the same and the test should report both.

- **Distinguishable**: two different poses give different G. Reported by
  the confusion fraction C(σ_G) (Test 3).
- **Recoverable**: from G alone, you can *infer* the pose. This is
  stronger and is reported by the inverse map's accuracy. To test, train
  a small regressor (random forest or MLP) on (G_features → pose) on 80%
  of the (z, R) cube and report the test-set RMSE on the held-out 20%.
  Use only G-derived features (G itself, dG/dt if you have a translocation
  trajectory, blockade duration, area). Report:
  - RMSE_z (Å)
  - RMSE_R (degrees, geodesic distance on SO(3))

Distinguishability without recoverability is common (the signal contains
the information but it is not extractable without prior knowledge of the
analyte structure). Both numbers belong in the report.

---

## 8. PETK implementation checklist

- [x] z-scan driver (existing `__main__.py`)
- [x] rotation utilities (`sem/rotation.py`, `RotationSpec`)
- [ ] `sem/scripts/compute_steric_overlap.py` — V_overlap ground truth
- [ ] `sem/scripts/sample_so3.py` — Sobol/HEALPix on SO(3); emits a
      list of `RotationSpec`
- [ ] `sem/scripts/joint_pose_scan.py` — runs the (z, R) cube, stores
      a single HDF5 with z, R, V_overlap, G arrays
- [ ] `sem/scripts/pose_distinguishability.py` — computes ΔG_max,
      z_bins, σ_R(G), C(σ_G), I(G; ·), Pearson r, RMSE_z, RMSE_R from
      the HDF5
- [ ] `sem/scripts/symmetry_audit.py` — runs §6 invariance checks and
      reports residuals as a pass/fail table
- [ ] `sem/scripts/run_pose_battery.py` — orchestrator that runs the
      full battery for one (geometry, analyte) pair and emits the §4
      table row

Each script should accept a single config JSON (same schema as the
existing PETK config, plus `rotation_spec_list` and `z_list`) so the
whole battery can be re-run with one command.

---

## 9. Caveats

1. **Charge.** With `ChargeAwareConductivityModel`, V_overlap is no
   longer the only physical predictor; the analyte's charge distribution
   adds a second channel. Run the battery once with the simple model
   (geometry-only) and once with the charge-aware model, and report the
   *gain* in I(G; R) from charge. That gain is the signature of charge-
   sensitive sensing.
2. **Conformational flexibility.** A real analyte is not rigid. The test
   above holds the structure rigid and varies pose. To extend, re-run
   Test 2 with N conformers from an MD ensemble at fixed z₀; the spread
   of G then has two components, σ_R (rotation) and σ_C (conformation).
   Report both.
3. **Translocation kinetics.** Fast translocation averages over R; slow
   (or trapped) translocation does not. Whether I(G; R) is *useful* in
   practice depends on the trapping time vs the rotational correlation
   time. Quote the I(G; R) numbers with the rotational diffusion time
   τ_rot of the analyte estimated from its hydrodynamic radius; if τ_rot
   << dwell time, the orientation channel is effectively averaged out
   and I(G; R) overstates the practically recoverable information.
4. **Grid resolution.** The σ_R(G) residual on a cylindrical pore is the
   tightest sanity check on grid anisotropy. If it doesn't drop with
   `grid_resolution` halving, the rotation pipeline (`scipy.ndimage.rotate`
   in `rotation.py`) has a discretisation artefact that needs fixing
   before any C(σ_G) is trustworthy.

---

*Last revision: 2026-05-06.*
