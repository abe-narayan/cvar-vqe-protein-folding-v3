# Sprint 14 — ENER findings

The controlled comparison of the two physical energy models (Legacy 11-term knowledge-based;
AMBER ff14SB/GBn2 all-atom via OpenMM) and the structure of their landscapes.

Instrument confirmed at start: `shipped 3.4540004952559396, pool_best 1.7108244199364904,
top75_best 2.3061526409453816, synthesis_fit 3.2040761603809194, n_zero_recall 18`.

Every arm that reads a native quantity is labelled ORACLE. `rmsd` in the enumerated cache is
a stored ORACLE label used only for post-hoc scoring.

Modules: `s14/ener_lib.py` (shared), `s14/ener_subset.py`, `s14/ener_matrix.py`,
`s14/ener_decoy.py`, `s14/ener_complement.py`, `s14/ener_invisible.py`, `s14/ener_refine.py`.

---

## E0. THE CACHED AMBER SUBSET IS NOT A UNIFORM SAMPLE — every previous pooled AMBER
## number is conditioned on an ORACLE selection

**Status: DEMONSTRATED.** `s14/ener_subset.py`, `s14/results/ener_subset.json`.

`s13/results/qarch_enum_<PDB>.npz` carries genuine AMBER on 2,955 of 262,144 configurations
(1.13% of the space). Verified against the data, that subset is three strata:

| `amber_kind` | stratum | n/target | selection | mean RMSD vs full space |
|---|---|---|---|---|
| 0 | uniform | ~1,194 | uniform over the space | **+0.008 A** (KS 0.026) |
| 1 | prior | ~1,150 | sampled from the leakage-safe 1-local prior | — |
| 2 | ORACLE band | 600 | drawn from the **1.00% lowest TRUE RMSD** | selected on the label |

Pooled over the nine targets the **whole subset is 0.401 A better in RMSD than the space it
came from** (KS statistic 0.206 against the full population). The band's threshold is
1.757–2.427 A per target — i.e. exactly the best 1.00% of each enumeration.

**Operating rule for this sprint, pinned:** any AMBER statistic meant to describe the
objective on a native-free population must be computed on `amber_kind == 0` only
(`E.enum(pdb).uniform_mask`). The uniform stratum is clean: mean RMSD bias +0.008 A,
KS 0.026. Statistics on the full 2,955 are conditioned on an oracle and must be labelled
**ORACLE DIAGNOSTIC**. The band stratum remains useful — it is the only near-native AMBER
population that exists — but only as an in-band diagnostic, never as evidence about the
objective's behaviour on the space.

This does not overturn the Sprint 13 AMBER table, which already stratified correctly
(`s13/qarch_FINDINGS.md` section 3c reports per-population rows). It pins the rule so that
the pooled figure is never quoted.

---

## E1. THE CONTROLLED MATRIX — arms A-E, H, I on one fixed ensemble

**Status: DEMONSTRATED** (ranking metrics; the RMSD label is post-hoc scoring only).
`s14/ener_matrix.py`, `s14/results/ener_matrix.json`.

Ensemble: the uniform AMBER stratum, ~1,194 configurations per target on the nine
fully-enumerated n=9 k=4 targets. No configuration was selected by any energy or by the
label, so all arms see exactly the same conformations. Pooled references:
**random draw 3.789 A, pool best 1.295 A, ORACLE snap 1.921 A.**

| arm | objective | rho | **rho in low-E decile** | pair acc | selected RMSD | vs random draw | W/L |
|---|---|---|---|---|---|---|---|
| A | Legacy total | +0.093 | **+0.079** | 0.529 | 3.493 | -0.296 [-0.845,+0.283] | 6/3 |
| B | **AMBER total** | +0.085 | **-0.075** | 0.529 | 3.834 | **+0.045** [-0.540,+0.549] | 3/6 |
| E | 1-local prior | +0.013 | -0.031 | 0.504 | 4.334 | **+0.545** [+0.113,+1.018] | 1/8 |
| E | radius of gyration | -0.186 | -0.373 | 0.422 | 5.531 | +1.742 [+0.861,+2.555] | 1/8 |
| H | rank(Leg)+rank(AMB) | +0.106 | +0.035 | 0.535 | 3.508 | -0.281 [-0.601,+0.044] | 5/4 |
| H | robust-z variant | +0.110 | +0.038 | 0.536 | 3.493 | -0.296 [-0.845,+0.283] | 6/3 |
| — | Legacy minus steric | -0.047 | +0.052 | 0.483 | 3.941 | +0.152 | 6/3 |
| — | `leg_steric` alone | +0.154 | **undefined** | 0.591 | 3.702 | -0.087 [-0.189,+0.016] | 6/3 |
| — | `leg_torsion` alone | -0.001 | +0.010 | 0.500 | **2.954** | **-0.835** [-1.607,-0.105] | 6/3 |
| — | `amb_nonbonded` | +0.091 | +0.021 | 0.532 | 4.117 | +0.328 [+0.040,+0.667] | 3/6 |

**Not one CI on a physical energy excludes zero.** Every lambda in the sweep
`rank(Legacy) + lambda * rank(AMBER)` for lambda in {0.25, 0.5, 1, 2, 4} lands between
-0.091 and -0.281 A with a CI straddling zero; the combination never beats Legacy alone
(3.493). H buys nothing.

**Cascades (arms C and D), each against a RANDOM keep-set of the same size:**

| arm | selected RMSD | matched random-prefilter NULL | value of the prefilter |
|---|---|---|---|
| C Legacy(top 10%) then AMBER | 3.630 | 3.682 | **-0.052 A** |
| D AMBER(top 10%) then Legacy | 3.719 | 3.641 | **+0.078 A — worse than random** |

Both cascades are inside their own null. **Arm D is the sharper statement: pre-filtering
with AMBER and then ranking with Legacy is worse than pre-filtering at random.**

**Arm I (post-hoc validation only — the guard vetoes, the prior selects):** every veto arm
is at or worse than the unvetoed prior. `I_prior_veto_leg_steric_10` is +0.757 A worse than
random. A 50% veto by `leg_steric` gives -0.133 A [-0.608,+0.357], the best veto arm and
still inside noise.

### E1a. `leg_steric` is a BINARY GATE, not a ranking function — this is the mechanism

**Status: DEMONSTRATED.** Degeneracy diagnostic in `ener_matrix.py`.

Per target, over the uniform ensemble:

| objective | fraction tied at the minimum | distinct values | **fraction tied inside the lowest-energy decile** |
|---|---|---|---|
| `leg_steric` | **0.735** | 174 | **1.0000** |
| `leg_compactness` | 0.270 | 842 | **1.0000** |
| `amb_bond` | **1.0000** | **1** | 1.0000 |
| `amb_angle` | 0.538 | 5 | 0.874 |
| `leg_coop_sheet` | 0.446 | 2 | 0.457 |
| `legacy` / `amber` | 0.0008 | ~1,180 | 0.008 |

`leg_steric` is exactly zero on **73.5%** of the space and is **constant across the whole
lowest-energy decile**, which is why its decile rho is undefined rather than small. Its
global rho of +0.154 and pair accuracy of 0.591 are generated entirely by the 26.5% of the
space where it fires. That is the mechanism behind "Legacy's entire skill is clash
rejection": the term is a pass/fail gate with a 73.5% pass rate, and **inside the pass set
it carries no ordering information at all**. `amb_bond` is literally constant (one distinct
value across 262,144 ideal-geometry configurations) and `amb_angle` takes five values —
under ideal backbone geometry AMBER's bonded terms cannot vary, so AMBER's five-term
decomposition is effectively three terms, and one of those dominates (E2).

---

## E2. THE DECOY DISCRIMINATION THRESHOLD — the single most important number

**Status: DEMONSTRATED** (the pairwise-accuracy formulation; the anchor test is
ORACLE DIAGNOSTIC because it selects anchors on the label).
`s14/ener_decoy.py`, `s14/results/ener_decoy.json`, `s14/results/ener_decoy_amber.log`.

### E2a. The obvious version of the experiment is confounded, and the confound is the answer

Pairwise decision accuracy `P(sign(dE) == sign(dRMSD))` conditioned on the pair's mutual
CA-RMSD looked like a clean resolution curve — Legacy 0.521 -> 0.578 and AMBER 0.505 ->
0.605 as separation goes from 1 A to 5 A. **It is not.** The mean quality gap `|dRMSD|`
rises with separation in lockstep (0.103, 0.248, 0.348, 0.486, 0.721, 1.123, 1.807, 2.531
across the eight separation bins). Conditioning on BOTH (`per_objective_2d`) shows the
separation dependence essentially vanishes at fixed gap — Legacy at gap < 0.25 A reads
0.505 / 0.525 / 0.509 / 0.504 / 0.506 / 0.482 / 0.470 across separation, i.e. flat and at
chance:

> **Neither energy has structural resolution in the geometric sense. What both have is a
> QUALITY-GAP threshold: they can order two structures only when the structures differ
> enough in RMSD, and it does not matter how far apart the two structures are.**

### E2b. The threshold. |dRMSD| needed to reach a given pairwise decision accuracy

Pooled over the nine targets, 3.6 M sampled pairs, uniform ensemble:

| objective | acc at gap < 0.25 A | acc at 0.5 A | **gap for 55%** | **gap for 60%** |
|---|---|---|---|---|
| `leg_steric` | 0.511 | 0.532 | **0.58 A** | 1.24 A |
| `leg_electrostatic` | 0.506 | 0.519 | 1.00 A | 2.30 A |
| **Legacy total** | 0.504 | 0.511 | **1.31 A** | 2.03 A |
| `amb_solvation` | 0.504 | 0.514 | 1.38 A | 2.93 A |
| `amb_nonbonded` | 0.502 | 0.508 | 1.50 A | 3.24 A |
| **AMBER total** | 0.501 | 0.507 | **1.58 A** | **3.87 A** |
| rank(Leg)+rank(AMB) | 0.503 | 0.511 | 1.21 A | 2.01 A |
| 1-local prior | 0.500 | 0.501 | **never** | never |
| `leg_torsion` | 0.500 | 0.499 | **never** | never |
| Legacy minus steric | 0.501 | 0.499 | **never** | never |
| `amb_bond`, `amb_angle`, `leg_aromatic`, `leg_hbond_local` | ~0.50 | ~0.50 | never | never |
| `leg_solvation`, `leg_hbond_longrange`, `leg_compactness`, `leg_coop_sheet` | <0.50 | <0.50 | never (ANTI-ranking) | never |

**THE NUMBER. At a quality gap below 0.25 A no objective in this project — physical or
otherwise — exceeds 0.511 pairwise accuracy. Legacy needs a 1.31 A gap to reach a coin-flip
plus five points; AMBER needs 1.58 A, and 3.87 A to reach sixty percent.** The pool best on
this ensemble is 1.295 A against a mean of 3.789 A, so the entire useful range of the
search is about 2.5 A wide: an objective with 1.3-1.6 A resolution can resolve roughly ONE
bit of that range. That is the architecture's behaviour in one number, and it explains
everything downstream — why the certified global optimum is worse than random sampling, why
optimising harder makes structures worse, and why no reranker has ever transferred.

The per-target CIs on the crossing point are wide ([0.55, 1.70] for Legacy, [0.41, 0.93]
for AMBER) because a per-target crossing is a noisy statistic; the pooled curve is the
reliable object and the ORDERING (`leg_steric` < Legacy < AMBER, with the 1-local
objectives never crossing) is stable across all nine targets.

### E2c. Legacy actively PREFERS decoys to good structures

**Status: ORACLE DIAGNOSTIC** (anchors selected on true RMSD). Decoys built by mutating m
of the anchor's 9 torsion states, m = 1..9, 60 per cell, 8 anchors per target. Reported as
`P(E_anchor < E_decoy)` minus the SAME construction around RANDOM anchors — the null that
removes anchor-typicality:

| objective | 1.0 A | 1.5 A | 2.0 A | 3.0 A | 4.0 A | 5.0 A |
|---|---|---|---|---|---|---|
| **Legacy total** | **-0.121** | -0.043 | **-0.189** | **-0.167** | -0.121 | -0.031 |
| Legacy minus steric | -0.110 | -0.060 | -0.206 | -0.177 | -0.160 | -0.115 |
| `leg_solvation` | **-0.247** | -0.178 | -0.209 | -0.201 | -0.222 | **-0.272** |
| `leg_compactness` | -0.028 | -0.166 | -0.160 | -0.127 | -0.185 | -0.175 |
| **AMBER total** | n/a | n/a | **+0.056** | +0.037 | **+0.113** | +0.089 |
| `leg_steric` | +0.049 | +0.030 | -0.003 | -0.008 | +0.049 | **+0.114** |
| `leg_torsion` | +0.011 | **+0.161** | **+0.221** | **+0.214** | +0.213 | +0.162 |
| 1-local prior | +0.046 | **+0.142** | **+0.168** | +0.148 | +0.125 | +0.118 |

**Legacy's excess AUC is negative at every separation.** Presented with a good structure and
a decoy built from it, Legacy prefers the decoy more often than a null anchor would predict,
and the effect is largest (-0.19) at 2 A — exactly the band where the answer lives.
`leg_solvation` is the worst single offender at -0.25. This is a stronger statement than
"Legacy does not rank": it is a **measured anti-preference for near-native geometry**.

**AMBER is the only physical objective with a positive excess AUC** (+0.04 to +0.11, best at
4-5 A). That is genuine, small, and consistent with E2b: AMBER discriminates, but only once
the decoy is bad enough. **This is the strongest positive result for AMBER in this sprint
and it is a validation signal, not a ranking signal.**

The only objectives with a large positive excess are `leg_torsion` (+0.16 to +0.22) and the
1-local prior (+0.12 to +0.17) — both Ramachandran-shaped and neither a physical energy;
see E7 for their null control.

---

## E3. COMPONENT DECOMPOSITION — both totals ARE a steric potential, confirmed exactly

**Status: DEMONSTRATED, exact.** `s14/ener_complement.py`, `s14/results/ener_complement.json`.
Legacy on the FULL 262,144-configuration enumeration (no sampling error at all); AMBER on
the uniform stratum.

| Legacy term | weight | **covariance share** | abs(w)*sd | rho with total | rho with RMSD | frac exactly 0 |
|---|---|---|---|---|---|---|
| **steric** | 4.0 | **0.988** | **12.766** | +0.326 | +0.156 | **0.741** |
| contact | 1.0 | 0.013 | 0.466 | -0.029 | +0.040 | 0.019 |
| electrostatic | 1.0 | 0.007 | 0.252 | +0.264 | +0.126 | 0.444 |
| coop_helix | 2.0 | 0.005 | 0.507 | +0.305 | +0.000 | 0.951 |
| hbond_local | 1.0 | 0.003 | 0.558 | +0.406 | -0.031 | 0.501 |
| torsion | 0.15 | 0.001 | 0.249 | +0.017 | +0.009 | 0.000 |
| aromatic | 0.8 | -0.000 | 0.126 | +0.365 | -0.018 | 0.556 |
| **coop_sheet** | 2.0 | -0.000 | 0.034 | -0.037 | -0.025 | **1.000** |
| compactness | 0.4 | -0.002 | 0.534 | +0.143 | -0.177 | 0.273 |
| solvation | 0.5 | -0.006 | 0.745 | +0.448 | -0.128 | 0.000 |
| hbond_longrange | 3.0 | -0.007 | 0.315 | -0.095 | -0.071 | 0.985 |

| AMBER term | covariance share (raw) | share after monotone conditioning | rho with total | distinct values |
|---|---|---|---|---|
| **nonbonded** | **1.000** | **0.995** | +0.995 | 1,174 |
| solvation | 0.000 | 0.275 | +0.275 | 1,194 |
| torsion | -0.000 | -0.030 | -0.030 | 1,192 |
| angle | -0.000 | 0.022 | +0.041 | **5** |
| bond | 0.000 | 0.000 | undefined | **1** |

**CONFIRMED, and more strongly than the brief states it.** Legacy `steric` share **0.988**
(the brief recorded 0.955) and AMBER `nonbonded` **1.000** raw, **0.995** after monotone
rank conditioning — so the AMBER figure is not the delta-spike artefact it could have been.
abs(w)*sd makes the mechanism unmissable: Legacy's steric term contributes 12.766 units of
weighted spread against 0.745 for the next largest. The other ten terms are not
contributing signal that the total drowns out; they are contributing **1.2% of the
variance**, and three of them (`coop_sheet` identically zero everywhere,
`hbond_longrange` zero on 98.5%, `coop_helix` zero on 95.1%) are essentially inert on this
space. Under ideal backbone geometry `amb_bond` is CONSTANT and `amb_angle` takes five
values, so AMBER's five-term decomposition is really three terms with one of them at 99.5%.

### E3a. Reweighting the terms does not recover anything — REFUTED

**Status: REFUTED (my own hypothesis).** Leave-one-target-out linear reweighting of all
eleven Legacy terms on rank-normalised features:

| | rho | decile rho | selected RMSD |
|---|---|---|---|
| in-sample (the CEILING of any linear reweighting) | +0.209 | +0.039 | 3.274 |
| **leave-one-target-out (honest)** | **-0.031** | **-0.087** | **4.036** |
| Legacy total, same ensemble | +0.093 | +0.079 | 3.493 |

LOTO vs random draw **+0.247 A [-0.280, +0.776], 3W/6L**; LOTO vs the Legacy total
**+0.543 A [-0.316, +1.443], 3W/6L**. The in-sample rho of +0.209 is the entire size of
the effect and **none of it survives to a held-out target**. There is no reweighting of
Legacy's physical terms that generalises. The default weights are not destroying signal:
there is no signal in the components to destroy.

---

## E4. COMPLEMENTARITY (H7) — real but priced out, exactly as the project's arithmetic predicts

**Status: HYPOTHESIS H7 REFUTED as an exploitable route; the decorrelation itself is
DEMONSTRATED.** `s14/ener_complement.py`.

Truth-partialled error correlation (each channel's rank residualised on rank(RMSD), then
correlated) versus the raw rank correlation, uniform ensemble, nine targets:

| pair | raw | **truth-partialled error corr** |
|---|---|---|
| Legacy / AMBER | +0.077 | **+0.096** |
| Legacy / `leg_steric` | +0.333 | +0.371 |
| AMBER / `leg_steric` | +0.425 | +0.373 |
| AMBER / `amb_nonbonded` | +0.995 | +0.995 |
| AMBER / rg | -0.467 | -0.399 |
| prior / `leg_torsion` | +0.680 | +0.684 |
| Legacy / prior | -0.011 | -0.016 |

**Legacy and AMBER genuinely are near-decorrelated** — partialled error correlation +0.096.
H7's premise is correct. Partialling truth barely changes the numbers here (+0.077 ->
+0.096), so this is not the "both tracking compactness" artefact the brief warns about.

**But the decorrelation cannot be spent.** Channel skill on the same ensemble is Legacy
rho +0.093 and AMBER rho +0.085. Pricing the optimal linear fusion by the multiple-
correlation identity `R^2 = (r_a^2 + r_b^2 - 2 r_a r_b r_ab)/(1 - r_ab^2)`, fitted PER
TARGET (an oracle ceiling), the top of the whole table is:

| pair | rho_fused | **gain over the stronger channel** |
|---|---|---|
| Legacy / Legacy-minus-steric | +0.389 | +0.0832 |
| Legacy / `leg_steric` | +0.400 | +0.0490 |
| **Legacy / AMBER** | +0.325 | **+0.0324** |
| Legacy / `amb_nonbonded` | +0.316 | +0.0317 |
| Legacy / prior | +0.194 | +0.0046 |
| AMBER / prior | +0.239 | +0.0020 |

The two largest entries are Legacy fused with pieces of ITSELF — a reweighting, which E3a
shows does not survive to a held-out target. **The genuine Legacy-AMBER fusion ceiling is
+0.032 rho, with a per-target-optimal weight, and it is unattainable:** the identical
machinery (E3a) loses 0.240 rho going from in-sample to leave-one-target-out. The gain is
second order in the weak channel's unique skill (`r_b - r_ab*r_a` = +0.055 for
Legacy/AMBER), which is the project's established arithmetic reproduced independently here.
Measured directly in the matrix, `rank(Legacy)+rank(AMBER)` selects 3.508 A against Legacy
alone at 3.493 A — the fusion is 0.015 A **worse**, and every lambda in the sweep is inside
the noise.

**Verdict on H7: the information is there and is worth at most +0.032 rho under an oracle
weight; there is no route to spending it.**

### E4a. Pareto structure — there is not enough conflict to justify it, REFUTED

**Status: REFUTED.** Six axes (1-local prior, Legacy, AMBER, radius of gyration, hydrogen
bonding, helix fraction), uniform ensemble:

* mean absolute pairwise rank correlation **0.222**; only **4.9 of 15** axis pairs are
  negatively correlated at all, and the strongest conflicts are AMBER vs rg (-0.467) and
  prior vs hbond (-0.228) — both conflicts between an energy and a shape descriptor, not
  between two claims about nativeness;
* the non-dominated frontier is **14.6% of the pool (174 of 1,194)** and its best member is
  **1.665 A against 1.713 A for a RANDOM subset of the same size** — a 0.048 A difference,
  i.e. the frontier is worth almost exactly what its cardinality is worth;
* the frontier **retains only 28.3% of the truly-best 1%** of the pool. A multi-objective
  formulation over these axes DISCARDS 72% of the good structures.

A scalar energy is not discarding structurally useful states that a Pareto formulation
would keep. **Do not build the multi-objective machinery.**

---

## E5. NORMALISATION — verified, and it is provably the wrong lever

**Status: DEMONSTRATED.** `s14/ener_norm.py`, constants exported to
`s14/cache/ener_norm.json` for the other workstreams.

Derived constants, pooled (Legacy and its terms exactly, from the full enumeration):

| objective | sd | MAD | **grad_sd** (exact mean single-flip sd) | sd/grad_sd | log10 range |
|---|---|---|---|---|---|
| legacy | 12.78 | 0.913 | **7.847** | 1.63 | 6.81 |
| `leg_steric` | 3.191 | **0.000** | 1.957 | 1.63 | 8.44 |
| legacy_nosteric | 1.539 | 0.825 | 0.718 | 2.14 | 5.02 |
| prior | 2.257 | 1.562 | 0.752 | 3.00 | 0.59 |
| `leg_torsion` | 1.662 | 1.148 | 0.554 | 3.00 | 4.39 |
| **amber (raw)** | **2.27e+17** | 8.68e+04 | n/a | n/a | **16.32** |
| amber (softcore) | 8.183 | 5.903 | n/a | n/a | 1.47 |

`grad_sd` is the exact mean single-flip standard deviation — the scale a local search
actually experiences, and the correct normaliser, computed with no sampling error. Raw
AMBER's sd is 2.27e+17 and is **not an estimate of anything** (Sprint 13: 40-68% of
ideal-geometry configurations exceed 1e4 kcal/mol); use rank normalisation for AMBER.
`leg_steric`'s MAD is **exactly 0**, which is E1a's gate structure showing up in the
normalisation constant itself.

**Monotone invariance, verified rather than assumed, on every metric this sprint uses:**

| AMBER transform | rho | decile rho | pair acc | selected RMSD | log10 range |
|---|---|---|---|---|---|
| raw | +0.0847 | -0.0754 | 0.5293 | 3.834 | **16.32** |
| softcore (sign-preserving log) | +0.0847 | -0.0754 | 0.5293 | 3.834 | **1.47** |
| rank | +0.0847 | -0.0754 | 0.5293 | 3.834 | 3.08 |
| robust z (median/MAD) | +0.0847 | -0.0754 | 0.5293 | 3.834 | 15.98 |
| affine 3.7x - 11 | +0.0847 | -0.0754 | 0.5293 | 3.834 | 16.67 |
| clip at p99 (NOT monotone) | +0.0847 | -0.0754 | 0.5293 | 3.834 | 12.17 |
| clip at p90 (NOT monotone) | +0.0842 | -0.0754 | 0.5295 | 3.834 | 8.19 |

**Bit-identical across every monotone transform**, while the dynamic range collapses from
16.32 decades to 1.47. Even the non-monotone caps move rho by 0.0005. This confirms the
brief on a second, independent instrument and closes the question:
**every scale-based intervention is exactly a no-op on ranking. The damage is in the
ordering. Any proposal to "fix AMBER by normalising it" is treating the wrong problem and
can be rejected without running it.**
---

## E6. ARMS F, G, J — SEARCH DESTROYS STRUCTURE; REFINEMENT DOES NOTHING

**Status: DEMONSTRATED.** `s14/ener_refine.py`, `s14/results/ener_refine.json`,
`s14/results/ener_refine.log`. Nine targets, five starts each = 45 (target, start) cells,
all four stages measured on the SAME structure so the refinement effect is attributable.

| start | raw | geometric cleanup | **F: Legacy descent** | **G: AMBER descent** | **J: AMBER refine** | J - raw | CI | W/L |
|---|---|---|---|---|---|---|---|---|
| random | 3.795 | 3.801 | 3.658 | 3.625 | 3.807 | +0.012 | [-0.008,+0.030] | 3/6 |
| prior argmin | 4.334 | 4.334 | 3.909 | 3.980 | 4.337 | +0.003 | [-0.051,+0.057] | 3/6 |
| Legacy argmin | 3.493 | 3.495 | 3.789 | 3.704 | 3.493 | +0.000 | [-0.017,+0.017] | 4/5 |
| AMBER argmin | 3.834 | 3.834 | 3.779 | 4.251 | 3.862 | **+0.028** | **[+0.010,+0.047]** | 1/8 |
| **ORACLE pool best** | **1.295** | 1.293 | **3.684** | **3.685** | **1.321** | +0.026 | [-0.006,+0.058] | 2/7 |

### E6a. Both energies drive a 1.3 A structure to 3.7 A. This is the study's worst number.

Started from the best structure in the pool at **1.295 A**, coordinate descent on Legacy
returns **3.684 A** and coordinate descent on genuine AMBER returns **3.685 A**. Neither is
noise, neither is a local-minimum artefact of one target, and the two energies agree almost
exactly on where to go. That is the anti-ranking measured directly rather than inferred
from a correlation: **the descent direction of both physical energies points away from the
native, from a near-native start, on nine of nine targets.**

Pooled over all 45 cells:

| arm | mean RMSD | vs raw | CI | W/L |
|---|---|---|---|---|
| geometric cleanup only | 3.352 | +0.001 | [-0.001,+0.004] | 26/19 |
| **J: AMBER refinement** (`refine_coords` k=10 steps=0) | 3.364 | **+0.014** | [-0.001,+0.028] | 13/32 |
| **F: Legacy descent** | 3.764 | **+0.414** | **[+0.038,+0.796]** | 18/26 |
| **G: AMBER descent** | 3.849 | **+0.499** | **[+0.069,+0.929]** | 17/27 |

Arm F at its limit — Legacy's **CERTIFIED global optimum** over all 262,144 configurations —
gives **3.920 A against a whole-space mean of 3.781 A** and a space best of 0.969 A. That
reproduces the Sprint 13 figure (+0.139 A worse than random) exactly, on an independently
written implementation, and it now has a mechanism: E2b's 1.31 A resolution limit means the
last 2.5 A of the search is below the objective's noise floor, so the descent is a random
walk with a systematic drift set by whatever the steric gate happens to prefer.

### E6b. AMBER refinement is a no-op on CA-RMSD — the hypothesis is REFUTED

**Status: REFUTED (my own hypothesis, and the sprint brief's most plausible rescue for
AMBER).** Genuine restrained ff14SB/GBn2 `refine_coords(k_restraint=10, steps=0)` moves the
CA trace by **+0.014 A [-0.001, +0.028], 13W/32L** over 45 cells — statistically
indistinguishable from zero and pointing very slightly the wrong way. On the AMBER-argmin
start it is significantly harmful (+0.028 [+0.010, +0.047], 1W/8L).

The two important sub-results:

* **Refinement does not rescue a bad structure.** From the prior argmin at 4.334 A it
  returns 4.337 A. There is no cell in the matrix where refinement turns a poor structure
  into a good one.
* **Refinement does not destroy a good one either.** From the ORACLE pool best at 1.295 A
  it returns 1.321 A. This is a genuine, if modest, positive: **AMBER refinement is safe.**
  A pipeline may run it for all-atom physical validity (which is what it is for) without
  paying a CA-RMSD penalty. It just must not be credited with accuracy.

The measured cost is 1.6-12.1 s per structure (mean about 4 s at `steps=0` with the memo
warm), against 0.3 ms for a Legacy evaluation.

**Attribution, stated explicitly as the brief requires:** every RMSD in the J column is
within 0.03 A of its own raw start. **All of the variation across the J column comes from
the start, none of it from the refinement.**

### E6c. Geometric cleanup is exactly a no-op here, by construction

`I.project` moves the structures by +0.001 A [-0.001, +0.004]. This is expected and is
reported so it is not mistaken for a result: the discrete torsion configurations are BUILT
on the ideal-geometry manifold, so the projection has nothing to do. It is not evidence
that projection is useless in the pipeline, where the input is an arbitrary CA trace.

---

## E7. THE ONE POSITIVE RESULT AND ITS NULLS

> **SUPERSEDED — see E7-REFUTED immediately below. The measurement below stands; its
> interpretation does not. `leg_torsion` was tested at n=126 and is the worst arm tried.** — `leg_torsion` survives the helix control
## but fails the concentration check

**Status: HYPOTHESIS, not DEMONSTRATED.** `s14/ener_nulls.py`, `s14/results/ener_nulls.json`.

`leg_torsion` (the exactly-1-local Ramachandran penalty) selected 2.954 A against a random
draw of 3.789 — **-0.835 A [-1.607, -0.105], 6W/3L**, the only CI in the whole E1 matrix
that excluded zero on the good side. The project's standing warning is that this shape of
result is a helix artefact. The controls:

| arm | selected RMSD | vs random draw | CI | W/L |
|---|---|---|---|---|
| **`leg_torsion` (the claim)** | **2.954** | **-0.835** | **[-1.607,-0.105]** | 6/3 |
| NULL: pure helix fraction, zero physics | 3.819 | +0.030 | [-0.300,+0.389] | 5/4 |
| NULL: constant alpha-helix point (-63, -42) | 4.011 | +0.222 | [-0.100,+0.547] | 3/6 |
| `leg_torsion` with helix fraction partialled out | 3.280 | -0.509 | [-1.740,+0.628] | 4/5 |
| 1-local prior | 4.334 | +0.545 | [+0.113,+1.018] | 1/8 |
| prior with helix fraction partialled out | 3.503 | -0.286 | [-1.098,+0.467] | 5/4 |

Against the **constant-helix baseline** (the baseline the project says must be cleared):
`leg_torsion` **-1.057 [-2.033, -0.194], 6W/3L**; Legacy total **-0.518 [-0.903, -0.067],
8W/1L**; AMBER total -0.177 [-0.752, +0.322], 5W/4L; pure helix fraction -0.192, 4W/5L.

Split by the TARGET's own native helicity, `leg_torsion`'s advantage is **larger on the
non-helical targets (-0.995) than on the helical ones (-0.274)** — the opposite sign of a
helix artefact.

**But two caveats that keep this at HYPOTHESIS:**

1. **The null is uninformative on this ensemble.** rho(helix fraction, RMSD) is **+0.003**
   here — the helix confound channel is *dead* on these nine targets (seven of nine have
   native helix fraction below 0.5). So the control cannot detect the artefact even if it
   exists; it only establishes that the artefact is not driving *this* number.
2. **The result fails the concentration check.** Per-target deltas vs the pool mean are
   [+0.82, -2.36, -0.13, -1.12, -0.92, +0.37, -2.74, +0.06, -1.49]. **Two targets (7N2I,
   2MK7) carry 68% of the effect;** drop-top-1 is -0.596, drop-top-2 **-0.344**, drop-top-3
   **-0.153**. By this project's own rule a result carried by two of nine targets is not a
   result. n = 9 is too small for the mandated drop-top-10 check at all, and that is a hard
   limit of the enumerated instrument.

**The same check applied to Legacy is worse:** Legacy's -0.296 A advantage has a top-2 share
of **1.02** — drop-top-2 leaves **+0.009**, i.e. **Legacy's entire apparent advantage over a
random draw is two targets.** AMBER's +0.045 is likewise two targets, in the other
direction.

**Conclusion: `leg_torsion` is the most promising residue in the matrix and is worth a
confirmation on a larger target set, but nothing here supports deploying it. Note also that
it is not a physical energy — it is a Ramachandran prior, and E2b shows it never reaches
0.55 pairwise accuracy at any quality gap, so whatever it does is an argmin effect on a
1-local objective, not ranking.**
---

## E7-REFUTED. `leg_torsion` does not replicate at n=126. The concentration warning was right.

**Status: REFUTED.** Coordinator's `s14/hamil3.py`, mirrored here per the rule that a
correction is written next to what it replaces. **E7 above is kept in place unaltered** —
its measurement stands, its interpretation does not.

E7 filed `leg_torsion` as the one surviving positive: 2.954 A selected, -0.835 A
[-1.607, -0.105], 6W/3L on the nine enumerated targets, clearing the constant-alpha-helix
baseline (-1.057, 6W/3L) and surviving the helix controls. It was filed as **HYPOTHESIS**,
not DEMONSTRATED, on two explicit grounds: two of nine targets carried 68% of the effect
(drop-top-2 = -0.344, drop-top-3 = -0.153), and the helix null could not discriminate
because rho(helix fraction, RMSD) was +0.003 on that ensemble.

Tested at n=126, 3,000 native-free configurations per target, both a mixed and a strictly
uniform proposal:

| arm | rho all | **rho in-decile** | argmin RMSD | vs incumbent |
|---|---|---|---|---|
| prior + distogram | 0.388 | +0.367 | 3.571 | +0.367 [+0.251,+0.487] |
| distogram only | 0.539 | +0.095 | 3.633 | +0.429 [+0.283,+0.580] |
| prior only | 0.304 | +0.356 | 3.836 | +0.632 [+0.388,+0.885] |
| `leg_torsion` + distogram | 0.331 | +0.068 | 3.910 | +0.706 [+0.516,+0.913] |
| `leg_torsion` + prior | 0.269 | +0.140 | 4.306 | +1.102 [+0.800,+1.415] |
| **`leg_torsion` ONLY** | 0.230 | **-0.073** | **4.933** | **+1.729 [+1.419,+2.052]** |

Under a strictly uniform proposal: global rho **-0.073**, argmin 5.195, +1.991 against the
incumbent. **It is the worst arm tested, it anti-ranks in the low-energy decile under both
proposals, and it poisons every combination it enters** — added to prior+distogram it moves
the argmin from 3.571 to 3.670; added to the distogram alone, 3.633 to 3.910. There is no
weight at which it contributes, which kills the "weak but complementary" reading as well as
the "strong alone" one.

**I have no methodological objection to the n=126 test and do not want it rerun.** It is
consistent with everything my own instruments already said about this term, and in
retrospect three of my own measurements predicted it:

* **E2b**: `leg_torsion` **never reaches 0.55 pairwise decision accuracy at any quality
  gap** — it was in the "never" row of the resolution table alongside the 1-local prior and
  Legacy-minus-steric.
* **E1**: its global rho is -0.001 and its in-decile rho +0.010 on the uniform ensemble.
  The 2.954 A was **always** a pure argmin effect on a 1-local objective with no ranking
  content, and I labelled it as such.
* **E3**: it carries **0.001 of Legacy's covariance share**, and E10b shows 97.7% of an
  interior residue's variance sits on the two residues CA-RMSD cannot see.

**What this establishes methodologically, and why it belongs in this file rather than only
in the coordinator's:** the lead looked like the best native-free selection result the
project had produced, on a nine-target instrument with complete enumeration and no sampling
error whatsoever. Exhaustive enumeration does not protect against target concentration. The
only thing that caught it was the mandatory per-target concentration check, on n=9, where
the brief's own drop-top-10 is not even computable. **Report drop-top-k at whatever k the
sample allows, and treat a top-2 share above ~0.5 as disqualifying regardless of the CI.**

By the same rule, applied to my own headline arms and recorded here so the verdict is not
read more strongly than it should be: **Legacy's -0.296 A selection advantage over a random
draw has a top-2 share of 1.02** (drop-top-2 leaves +0.009), and AMBER's +0.045 A is
likewise two targets. Neither was ever claimed as a result — the E11 verdict rests on the
filter curve (E9), the clash audit (E8) and the resolution threshold (E2b), all of which
hold 6-9 targets of 9 or are exact over the full enumeration — but the concentration figures
belong next to them.

**The E11 verdict is unchanged.** `leg_torsion` was listed under "still open", not under
either model's legitimate role, and it is not a physical energy. It moves to REFUTED.

---

## E8. GEOMETRIC VALIDITY — the builder is rigid; the optimiser has no artefact to exploit

**Status: DEMONSTRATED.** `s14/ener_geom.py`, `s14/results/ener_geom.json`. Audited on the
uniform ensemble AND adversarially on each energy's own lowest decile, which is where an
exploit would appear and a uniform sample would miss it.

| quantity | value | sd | verdict |
|---|---|---|---|
| N-CA / CA-C / C-N / C=O | 1.4580 / 1.5250 / 1.3290 / 1.2310 A | **0.00000** | rigid |
| N-CA-C / CA-C-N / C-N-CA | 111.000 / 116.200 / 121.700 deg | **0.00000** | rigid |
| omega (folded; the raw dihedral is circular and its arithmetic mean is meaningless) | **180.0000 deg** | **0.00000** | trans, rigid |
| CA-CA virtual bond | 3.8040 A | **0.00000** | rigid, no chain breaks |
| N-CA-C-CB improper | -122.5863 deg, **100.00% L** | 0.00000 | correct chirality, rigid |

**Every geometric degree of freedom except phi and psi is bit-constant.** There is no
geometric artefact for an optimiser to exploit: it cannot stretch a bond, open an angle,
flip a peptide to cis, invert a centre, or break the chain. This closes the adversarial
question for the discrete torsion space.

**Non-bonded clash (heavy atoms at least two residues apart — including bonded and 1-3
pairs reports the C=O bond length of 1.231 A as "the minimum distance" and calls every
structure clashed, which was a real error in the first version of this audit):**

| population | min heavy-atom separation | mean | **fraction with a clash < 2.0 A** |
|---|---|---|---|
| uniform ensemble | **0.2277 A** | 2.670 | **13.07%** |
| Legacy's lowest decile | 1.6205 A | 2.772 | **0.56%** |
| **AMBER's lowest decile** | 1.9830 A | **3.023** | **0.09%** |

**This is the clean quantification of what both models actually do, and it is the strongest
evidence in the sprint for AMBER's legitimate role.** The space contains 13.07% clashed
configurations with separations down to 0.23 A. Legacy's low decile removes 96% of them;
**AMBER's low decile removes 99.3% of them and has the larger minimum separation.** AMBER is
a strictly better clash detector than Legacy — which is exactly what an all-atom force field
should be, and exactly consistent with E2c (AMBER is the only physical objective with a
positive excess decoy AUC) and E6b (AMBER refinement is safe).

**Ramachandran validity is NOT enriched by either energy**: 91.29% allowed in the uniform
ensemble, 91.47% in Legacy's low decile, 91.14% in AMBER's low decile. The torsion library
already supplies Ramachandran plausibility; the energies add none.

### E8a. AMBER single points are not rigid-motion invariant to better than 2.4e-4 relative

| check | Legacy | AMBER |
|---|---|---|
| rigid rotation + 50 A translation, max relative change | **4.20e-13** | **2.44e-04** |
| the same, max absolute change | 1.42e-13 | **1.61e+14 kcal/mol** |
| 1e-9 relative coordinate perturbation, max relative change | 1.02e-06 | — |
| cached enumeration vs recomputation, max absolute | 3.66e-05 (float32 storage) | — |

Legacy is invariant to machine precision. AMBER is invariant only to **2.4e-4 relative**,
which on the 40-68% of configurations whose energy exceeds 1e4 kcal/mol is an absolute
change of up to 1.6e14 kcal/mol. **On clashing ideal-geometry configurations an AMBER single
point is not a reproducible number in absolute terms.** It is still reproducible enough to
rank (a 2.4e-4 relative jitter cannot reorder configurations separated by decades), so this
does not invalidate the ranking results — but it does mean any AMBER-in-the-loop optimiser
must never difference two raw AMBER energies, and it independently reinforces E5: use rank
normalisation.

Legacy's 1e-9 -> 1.02e-6 amplification is a factor of ~1000, which is the steric term's
stiffness showing up as condition number. It is well behaved but it is not a smooth
objective near a contact.
---

## E9. THE VALIDATION ROLE, PRICED — AMBER filters without discarding good structures;
## Legacy filters by discarding them

**Status: DEMONSTRATED.** `s14/ener_filter.py`, `s14/results/ener_filter.json`.

The project's established terminal-operator law is `d_out = 1.16*d_set_mean +
0.04*d_set_best` (R^2 0.893), so a filter's job is to move the surviving set's MEAN. Each
objective removes the worst x% of the pool; the null removes x% AT RANDOM, which is the only
control for the fact that shrinking a set changes its statistics.

**Remove the worst 10%** (random null: mean 3.789, best 1.321, operator 4.448):

| guard | surviving mean | vs random removal | W/L | surviving **best** | vs null | operator delta |
|---|---|---|---|---|---|---|
| **`leg_steric`** | 3.731 | **-0.058 [-0.104,-0.010]** | **7/2** | 1.295 | **-0.026** | **-0.069 [-0.122,-0.012]** |
| **Legacy total** | 3.742 | **-0.048 [-0.085,-0.009]** | **7/2** | 1.392 | **+0.071** | -0.052 [-0.094,-0.008] |
| **AMBER total** | 3.768 | **-0.021 [-0.041,-0.002]** | **6/3** | 1.295 | **-0.026** | **-0.026 [-0.049,-0.003]** |
| `amb_nonbonded` | 3.768 | **-0.021 [-0.041,-0.002]** | 6/3 | 1.295 | -0.026 | -0.025 [-0.049,-0.003] |
| 1-local prior | 3.781 | -0.008 [-0.019,+0.002] | 6/3 | 1.295 | -0.026 | -0.011 [-0.024,+0.002] |
| Legacy minus steric | 3.806 | +0.017 [-0.040,+0.079] | 5/4 | 1.552 | +0.231 | +0.029 [-0.046,+0.111] |

**These are the only CIs in the entire ENER matrix that exclude zero on the good side for a
physical energy.** They are small — 0.021 to 0.058 A on the pool mean, 0.026 to 0.069 A
through the operator — but they are real, they hold 6-7 targets of 9, and they are exactly
the size a 13.07% clash rate should be worth.

**The discriminating fact between the two models is in the `best` column, not the `mean`
column.** Carried to larger removal fractions:

| guard | surviving best at 25% removed | at 50% | at 75% | at 90% |
|---|---|---|---|---|
| random null | 1.345 | 1.433 | 1.592 | 1.801 |
| **Legacy total** | 1.564 (**+0.219**) | 1.729 (**+0.296**) | 1.977 (**+0.385**) | 2.181 (**+0.380**) |
| Legacy minus steric | 1.638 (+0.293) | 1.760 (+0.327) | 2.001 (+0.409) | 2.276 (+0.475) |
| **AMBER total** | 1.295 (**-0.050**) | 1.417 (-0.016) | 1.672 (+0.079) | 1.769 (**-0.033**) |
| **`amb_nonbonded`** | 1.295 (**-0.050**) | 1.433 (+0.000) | 1.645 (+0.053) | 1.755 (**-0.046**) |
| `leg_steric` | 1.432 (+0.088) | 1.432 (-0.001) | 1.634 (+0.042) | 1.800 (-0.001) |

**AMBER throws away 90% of the pool and still holds the best structure in it.** Legacy
throws away 25% and has already lost 0.219 A of the best; at 75% it has lost 0.385 A —
worse than removing 75% at random. **`legacy_nosteric` is the worst filter at every
fraction, which is the third independent confirmation that Legacy's only usable content is
its steric gate.**

That is the decision-grade separation of the two models' roles: **as a filter both improve
the mean, but only AMBER (and Legacy's steric term in isolation) does so without destroying
the top of the pool. The Legacy TOTAL is an actively unsafe filter above about 10%
removal.**

### E9a. Remaining matrix metrics, for completeness

| objective | top-10 recovery of the truly-best 1% (chance = 0.010) | mean RMSD of the energy's own best 1% | calibration slope on rank-normalised E | Pearson r |
|---|---|---|---|---|
| Legacy total | 0.011 | 3.602 | +0.234 | +0.115 |
| AMBER total | 0.022 | 3.624 | +0.344 | +0.095 |
| 1-local prior | 0.022 | 3.985 | +0.036 | +0.017 |
| `leg_steric` | 0.000 | 3.683 | +0.740 | +0.179 |
| `leg_torsion` | 0.067 | 3.669 | +0.024 | -0.002 |
| `amb_nonbonded` | 0.022 | 3.645 | +0.362 | +0.102 |
| rank(Leg)+rank(AMB) | 0.000 | 3.685 | +0.334 | +0.125 |
| radius of gyration | 0.000 | 5.082 | -0.725 | -0.190 |

Top-10 recovery is at chance for every objective. The calibration slope is the total RMSD
swing across the objective's ENTIRE energy range: Legacy +0.234 A and AMBER +0.344 A against
a within-pool RMSD spread of about 1.2 A, i.e. **going from the best-ranked to the
worst-ranked configuration changes the expected RMSD by less than a third of a standard
deviation.**

**ORACLE-snap percentile** (where the native's in-space representative sits on each energy's
axis; 0 = the energy puts the native at the very bottom, 0.5 = no information):

| | pooled | per target |
|---|---|---|
| Legacy | **0.623** | 0.443, 0.911, 0.745, 0.546, 0.900, 0.008, 0.962, 0.823, 0.271 |
| AMBER | **0.457** | 0.369, 0.001, 0.003, 0.599, 0.836, 0.568, 0.291, 0.638, 0.813 |

**On this uniform ensemble Legacy places the native at the 62nd percentile — WORSE than a
coin flip and worse than the 32nd-40th percentile the brief records; AMBER places it at the
46th, i.e. exactly at no information.** The discrepancy with the brief's figure is a
population difference and is important: the Sprint 13 number was measured on populations
that included garbage configurations, where separating the native from garbage pulls the
percentile down. **On a population that is already Ramachandran-plausible — which is what a
sequence-conditioned torsion library always produces, and therefore what a real search
always sees — Legacy's native percentile is on the wrong side of 0.5.** Legacy's own value
of 0.008 on 6S0N against 0.962 on 7N2I shows how target-dependent this is; it is not a
stable property of the energy.
---

## E10. THE INVISIBLE DEGREES OF FREEDOM — a CORRECTION to the stated defect,
## and the variance each objective spends on them

**Status: DEMONSTRATED, exact.** `s14/ener_invisible.py`,
`s14/results/ener_invisible.json`, `s14/results/ener_invisible_amber.log`.

### E10a. CORRECTION: FOUR torsions are inert for the CA trace, not three

The sprint brief records that `core/project.py` leaves `phi[0]`, `psi[n-1]` and `phi[n-1]`
inert for the CA trace, so that residue n-1 is entirely invisible to CA-RMSD. Measured
directly, by perturbing one torsion at a time in the production ideal-geometry builder
(`I.build_ca` -> `core.project.build_ca_exact`) on 1CS9:

    phi[0] 0.000e+00   psi[0] 1.124e-07
    phi[1] 8.192e-02   psi[1] 2.480e-01
    ...
    phi[7] 3.553e-01   psi[7] 1.244e-01
    phi[8] 0.000e+00   psi[8] 0.000e+00

**`psi[0]` is inert too.** Both terminal residues are entirely invisible, not one. Confirmed
independently and exactly by the total Sobol index over the full enumeration: CA-RMSD's
total Sobol index is **0.00000** (machine zero, max 2.2e-16 across all nine targets) on BOTH
residue 0 and residue n-1, and 0.389 on an average interior residue.

The reason is geometric and general, not a bug in one builder: an n-residue CA trace has
n-2 virtual bond ANGLES (at residues 1..n-2) and n-3 virtual DIHEDRALS, and the standard
ideal-geometry map writes each of those in terms of `phi[i], psi[i]` for interior i only.
The four terminal torsions never enter.

**Consequence, which corrects every qubit count in the sprint.** At k=4 and n=9 the
brief's `n_res*log2(k)` = 18 nominal qubits, and the brief's correction gives 16 live. The
measurement gives **14 live qubits: 4 are dead, not 2, and 2 of 9 residues (22%) are
decision-free for CA-RMSD.** At n=9 that is a large fraction; it shrinks as 2/n and is
therefore worst on exactly the short peptides this project works on.

### E10b. How much of each objective rides on the states CA-RMSD cannot see

Exact total Sobol index over all 262,144 configurations, pooled over the nine targets
(a total Sobol index is the share of the objective's variance that vanishes when that
residue is averaged out: main effect plus every interaction involving it):

| objective | T[residue 0] | T[residue n-1] | T[mean interior] | **per-invisible-residue share of an interior residue** |
|---|---|---|---|---|
| **CA-RMSD (the metric)** | **0.00000** | **-0.00000** | 0.38929 | **0.000** |
| Legacy total | 0.01308 | 0.02484 | 0.49208 | **0.039** |
| `leg_steric` | 0.01283 | 0.02488 | 0.49155 | 0.038 |
| Legacy minus steric | 0.01592 | 0.01750 | 0.27805 | 0.060 |
| `leg_hbond_longrange` | 0.00000 | **0.19535** | 0.59463 | 0.164 |
| `leg_solvation` | 0.04347 | 0.02016 | 0.25160 | 0.126 |
| **1-local prior** | **0.08433** | **0.11932** | 0.11377 | **0.895** |
| **`leg_torsion`** | **0.09511** | **0.12345** | 0.11163 | **0.977** |
| **AMBER total** (softcore-conditioned, by k-state sweep) | **0.0349** | **0.1273** | 0.4610 | **0.176** |
| AMBER total (raw, uninterpretable) | 0.3817 | 0.2812 | 1.1337 | — |

Three results:

1. **Legacy is NOT heavily loaded on the invisible states — my own hypothesis is REFUTED.**
   Legacy spends 3.8% as much variance per invisible residue as per interior residue;
   together the two invisible residues carry **3.8% of Legacy's total decomposed variance**.
   That is not the mechanism of Legacy's energy-structure decoupling. The decoupling is the
   1.31 A resolution limit (E2b) and the steric gate (E1a), not an invisible-DOF leak.

2. **AMBER is 4.5x more exposed than Legacy.** After the mandatory monotone conditioning
   (the raw variance share exceeds 1 for an interior residue, which is exactly the
   delta-spike artefact the brief warns about and proves the raw number is unusable), AMBER
   spends **17.6% as much per invisible residue as per interior residue** — 0.176 against
   Legacy's 0.039. That is expected and correct physics: `phi[n-1]` genuinely places C, CB
   and O, and AMBER scores those atoms while the CA metric cannot. It means roughly
   **10% of AMBER's decision variance is spent on chemistry that CA-RMSD is blind to.**
   Some of that is real all-atom quality, but from the point of view of the CA-RMSD
   objective it is pure noise, and it is a genuine, quantified contribution to the
   energy-structure decoupling — reported here for the first time.

3. **The 1-local objectives are catastrophically exposed, and this is the finding with the
   most immediate consequence.** `leg_torsion` (0.977) and the empirical prior (0.895) put
   essentially the SAME variance on the invisible residues as on any other — necessarily,
   since a 1-local objective is a sum of independent per-residue terms. **At n=9,
   2/9 = 22% of everything a 1-local torsion prior decides is invisible to CA-RMSD.**
   Any predictor-quality knob, any coverage figure and any per-residue torsion accuracy
   measured on all n residues is therefore **overstated by about 22% at this chain length**,
   and a prior that simply declined to model the two terminal residues would lose nothing.
   This is directly relevant to the torsion-restraint coverage results and to the
   OBJ and SHIFT workstreams; the corresponding correction is to score per-residue channels
   on residues 1..n-2 only.

`leg_hbond_longrange` is the single most exposed Legacy term (0.195 on residue n-1, 16.4% of
an interior residue's share) — and it is also one of the two terms whose selected-RMSD CI
excludes zero on the BAD side. Those two facts are consistent.---

## E13. TAIL-RESTRICTED DISCRIMINATION, and the FILTERING/ORDERING decomposition that
## explains why good bulk discrimination coexists with a bad argmin

**Status: DEMONSTRATED, exact** (full 262,144-configuration enumeration for Legacy and
every `leg_*`; AMBER on the uniform stratum per E0). `s14/ener_tail.py`,
`s14/results/ener_tail.json`. Machinery exported as `ener_lib.tail_accuracy` and
`ener_lib.selection_decomposition` so other objectives compose with it.

Requested by the coordinator after observing that the pure distogram has the best bulk
pairwise accuracy of anything measured (0.654 below a 0.25 A gap, where every physical
energy sits at 0.51) and the **worst argmin of its family** (3.237 against the 1-local
prior's 2.878, whose bulk accuracy is chance-level). Neither `pair_accuracy` over uniform
pairs nor `decile_rho` at a fixed 10% predicts what an argmin returns, and this sprint
leaned on both.

### E13a. The prediction was stated first and holds: raw tail accuracy is confounded

E2a established that accuracy is a function of the QUALITY GAP. Inside a tail the available
gaps are small by construction, so **every objective's tail accuracy must fall toward chance
as q shrinks whether or not it has lost skill.** Measured, the mean `|dRMSD|` available
inside Legacy's own tail falls from 0.920 A (whole space) to **0.416 A at q = 0.2%**. A raw
tail-accuracy table would therefore rank objectives by how heterogeneous their tails happen
to be. Two controls are carried throughout: a **random-tail null** (the same objective
scored on a random subset of the same size) and a **gap-matched** accuracy restricted to
pairs with `|dRMSD|` in a fixed [0.5, 1.5) A window.

**The random-tail null is 0.524-0.527 at every q, not 0.500.** Every number below must be
read against that, not against a coin flip.

### E13b. An energy is BLINDEST inside its own tail

Gap-matched pairwise accuracy, pairs with `|dRMSD|` in [0.5, 1.5) A, nine targets:

| objective | q=0.2% | 0.5% | **1%** | 2% | 5% | 10% | 25% | whole space |
|---|---|---|---|---|---|---|---|---|
| **random-tail NULL** | 0.527 | 0.523 | **0.523** | 0.526 | 0.525 | 0.524 | 0.524 | 0.524 |
| **Legacy total** | 0.475 | 0.448 | **0.462** | 0.526 | 0.530 | 0.534 | 0.505 | 0.524 |
| Legacy minus steric | 0.460 | 0.472 | 0.486 | 0.515 | 0.509 | 0.531 | 0.503 | 0.490 |
| 1-local prior | 0.488 | 0.497 | 0.489 | 0.491 | 0.497 | 0.495 | 0.494 | 0.502 |
| `leg_torsion` | 0.522 | 0.514 | 0.510 | 0.508 | 0.501 | 0.499 | 0.497 | 0.500 |
| `leg_solvation` | 0.528 | 0.472 | 0.452 | 0.454 | 0.475 | 0.473 | 0.469 | 0.464 |
| AMBER total | — | — | — | — | 0.488 | 0.462 | 0.494 | 0.522 |
| `amb_nonbonded` | — | — | — | — | 0.506 | 0.502 | 0.505 | 0.525 |
| `leg_steric` | — | — | — | — | — | — | — | 0.581 |

(Cells with fewer than 30 distinct tail members are suppressed. A tail of a dozen
configurations still supplies 200,000 sampled pairs, and the first version of this table
reported 0.000 and 1.000 from exactly that pseudo-replication before the guard was added.)

**Legacy scores 0.462 inside its own lowest 1%, against 0.523 for Legacy on a random subset
of the same size** — and gap-matching does not rescue it, so this is not the shrinking-gap
artefact. `leg_steric` is entirely tied inside every tail (E1a) and cannot be scored at all
below q = 100%. **An energy is at its blindest precisely on the configurations it itself
selects.** That is a sharper statement than "it does not rank", and it is the correct
generalisation of the E2b resolution floor to the regime a search occupies.

### E13c. THE DECOMPOSITION. Filtering helps a little; ordering hurts more.

    sel_rmsd = pool_mean + (tail_mean - pool_mean) + (sel_rmsd - tail_mean)
                          \____ FILTERING ____/   \____ ORDERING ____/

Exact and additive. FILTERING is the value of being in the tail at all; ORDERING is the
value of the objective's ranking WITHIN its own tail. At q = 1%, nine targets, pool mean
3.781 A:

| objective | tail mean | selected | **FILTERING** | **ORDERING** | CI (ordering) | W/L |
|---|---|---|---|---|---|---|
| **Legacy total** | 3.657 | 3.920 | **-0.124** | **+0.263** | **[+0.046,+0.484]** | 4/5 |
| **`leg_steric`** | 3.712 | 3.691 | **-0.068** | **-0.021** | [-0.125,+0.082] | 5/4 |
| **AMBER total** | 3.648 | 3.834 | **-0.141** | +0.186 | [-0.206,+0.567] | 3/6 |
| `amb_nonbonded` | 3.626 | 4.117 | **-0.163** | **+0.491** | **[+0.095,+0.918]** | 2/7 |
| Legacy minus steric | 3.831 | 4.907 | +0.050 | **+1.077** | **[+0.383,+1.842]** | 1/8 |
| `leg_solvation` | 4.241 | 4.455 | **+0.460** | +0.213 | [-0.114,+0.542] | 5/4 |
| `leg_compactness` | 4.070 | 4.158 | +0.289 | **+0.088** | **[+0.037,+0.138]** | 1/8 |
| `leg_contact` | 3.806 | 3.950 | +0.025 | +0.143 | [-0.165,+0.497] | 6/3 |
| 1-local prior | 3.895 | 3.636 | +0.114 | -0.259 | [-0.642,+0.126] | 6/3 |

**Legacy's certified global optimum being +0.139 A worse than a random draw — the sprint's
central negative, first measured in Sprint 13 and independently reproduced in E6a —
decomposes exactly: FILTERING -0.124, ORDERING +0.263, net +0.139.** The additive identity
closes on the number. The energy's tail IS better than the pool; its ranking inside that
tail then throws the gain away and more, and the ordering CI excludes zero.

Three consequences:

1. **This resolves the coordinator's paradox as a general mechanism, not a curiosity.**
   Bulk discrimination and argmin quality measure DIFFERENT terms of an additive identity.
   An objective can have excellent bulk accuracy (a good FILTERING term) and a bad argmin
   (a positive ORDERING term) with no contradiction. The two should never be reported as
   though one predicted the other, and the decomposition should be reported instead of
   either. **This applies to the distogram and to any learned objective, and is offered to
   OBJ and to the coordinator for exactly that use.**
2. **`leg_steric` is confirmed a pure filter for the third independent time.** Its ORDERING
   term is **-0.021 [-0.125, +0.082]** — indistinguishable from zero — while its FILTERING
   term is -0.068. After the degeneracy diagnostic (E1a) and the filter curve (E9), this is
   the cleanest statement of it: the term has a filtering component and *no ordering
   component at all*.
3. **AMBER has the best FILTERING term of any objective measured** (-0.141 total, -0.163 for
   `amb_nonbonded`), better than Legacy's -0.124, and a harmful ordering term. That is the
   same conclusion E8 and E9 reached from clash rates and the filter curve, now derived a
   third way from an exactly additive decomposition. **The E11 verdict — AMBER is a
   validator, not a ranker — is now supported by three independent instruments.**

### E13d. Recommended reporting standard, and a caution about this sprint's own statistics

Neither global rho, nor in-decile rho, nor bulk pairwise accuracy should be used alone to
predict selection quality. **Report the FILTERING/ORDERING decomposition at the q an argmin
actually samples, with the random-tail null beside it.** Where a single number is wanted,
the ORDERING term is the one that decides whether optimising the objective harder helps or
hurts — it is the sign of "optimise harder, get worse" measured directly.

I note against my own work that E1's headline metrics (global rho, decile rho, selected
RMSD) are the statistics this section shows to be insufficient in isolation. The E11 verdict
does not change — it rests on E2b, E8 and E9, and E13 strengthens it — but E1's table should
be read with E13c beside it.

### E13e. On the coordinator's companion result

The finding that a learned structural objective (distogram) reaches 0.654 pairwise accuracy
below a 0.25 A quality gap, where every physical energy sits at 0.51, is **consistent with
E2b and does not contradict it.** E2b was scoped, deliberately, to the physical objectives
and to the components of this project's two energy models; the floor was always a claim
about Legacy and AMBER, not about the problem. That scoping is now load-bearing and I am
glad it was written that way. Combined with E13c the picture is coherent: the physical
energies supply a FILTERING term and no usable ORDERING term, so **something other than a
physical energy has to do the ranking** — which is the E11 verdict stated from the other
side.

Recorded for the corrections list, at the coordinator's request: they published a two-target
version of that discrimination table (1CS9 and 2MK7) showing the opposite conclusion, which
reversed completely on nine targets and is retracted in `s14/coord_FINDINGS.md` C15 with the
original left in place.

**One of those two overlaps with my own concentration failure, and I am stating the overlap
precisely rather than rounding it into a tidier claim.** The two targets carrying 68% of the
E7 `leg_torsion` effect were **7N2I (-2.74) and 2MK7 (-2.36)**, not 1CS9 and 2MK7 — the full
per-target deltas are [1CS9 +0.82, 2MK7 -2.36, 2P5H -0.13, 6EY3 -1.12, 6F3V -0.92,
6S0N +0.37, 7N2I -2.74, 8IS3 +0.06, 9UV5 -1.49]. So **2MK7 is the single target implicated
in both failures**; 1CS9 is implicated only in theirs and is in fact the target on which
`leg_torsion` did WORST. The defensible statement is therefore the narrow one: **2MK7 is a
recurring outlier on this instrument, and n=2 samples from these nine targets have now
produced a reversed conclusion twice.** The general rule stands on its own without needing
a two-target villain: **on a nine-target instrument, report drop-top-k at whatever k the
sample allows and treat a top-2 share above ~0.5 as disqualifying regardless of the CI.**
---

## E14. CAN PHYSICS REPLACE THE IDEAL-GEOMETRY PROJECTION? No — it costs exactly the same,
## and the 0.157 A is the price of un-shrinking an average, not the price of ideal geometry

**Status: DEMONSTRATED, n=126, the full development instrument.** `s14/ener_avgrefine.py`,
`s14/results/ener_avgrefine.json`, `s14/results/ener_avgrefine.log`. Proposed by the
coordinator off `s14/avgspace.py`, which measured that projecting the coordinate average
onto the ideal-geometry manifold costs **+0.157 A [+0.124, +0.191], 18W/108L**. The
proposal: a force field does not need ideal geometry, so let restrained AMBER recover
validity instead and see whether it can be had for less.

### E14a. The test as posed was not runnable, and the fix is recorded

`I.coordinate_average` averages `u["W"]`, which is a **CA trace only**;
`core.amber.refine_coords` requires N, CA, C, O and CB. The only route in this repository
from a bare CA trace to an all-atom structure is `core.amber.refine_ca`, which **fits the
discrete torsion states to the CA trace first** — and that fit IS a projection onto the
ideal-geometry manifold, the exact step under test. Run that way the experiment measures
projection-then-AMBER and reports it as AMBER-instead-of-projection.

The fix: every window carries its own PHI/PSI, so each has a full ideal-geometry backbone.
Superpose each window's full backbone onto that window's own CA trace **in the medoid
frame** — the frame `coordinate_average` works in — then average every atom. No projection
anywhere. **The reconstruction is asserted at run time, not assumed**, and it validates:

| | coordinator's `avgspace.py` | this module | agreement |
|---|---|---|---|
| arm A, coordinate average -> project | 3.205 | **3.205** | exact |
| arm B, coordinate average raw | 3.048 | **3.050** | 0.002 A |

**A correction to my own first implementation, kept here because the assertion is what
caught it.** `geo.build_backbone_batch(PHI, PSI)` returns each window in the *builder's*
reference frame, not `W`'s. Applying the W-derived superposition to it superposes nothing
and averages structures in unrelated frames. The identity check fired at **16.1 A** and the
"-86% bond contraction" that version produced is **RETRACTED as an artefact of my bug**. The
corrected figure is below and is a quarter of that.

### E14b. The result: AMBER ties the projection, to within 0.001 A

126 targets, same top-75 windows in every arm:

| arm | RMSD | FAIL18 | other-108 | vs arm A | CI | W/L | geom deviation | Ramachandran | clashed |
|---|---|---|---|---|---|---|---|---|---|
| **A projected** (incumbent) | **3.205** | 6.034 | 2.734 | — | — | — | 0 by construction | 0.734 | 0.000 |
| B raw average | 3.050 | 5.846 | 2.584 | **-0.155** | [-0.191,-0.122] | 99/27 | **0.2480** | 0.836 | **0.452** |
| **B + AMBER k=10** | **3.207** | 6.029 | 2.736 | **+0.001** | **[-0.014,+0.016]** | 53/73 | **0.0134** | **0.860** | **0.000** |
| B + AMBER k=2 | 3.248 | 6.068 | 2.778 | +0.043 | [+0.020,+0.066] | 51/75 | 0.0152 | 0.863 | 0.000 |
| B + AMBER k=0 (free) | 3.477 | 6.026 | 3.053 | +0.272 | [+0.177,+0.370] | 31/95 | 0.0170 | 0.868 | 0.000 |

> **SUPERSEDED IN PART by E14e below: k=10 is ONE POINT on a frontier that does not
> turn. At k=30 -- geometry at 1.13x the force field's own equilibrium -- AMBER beats
> the projection by -0.022 A [-0.036, -0.009] on 126 targets. The k=10 measurement
> below stands; the conclusion drawn from it alone does not.**

**This is the coordinator's outcome 2, measured to a CI 0.03 A wide.** Restrained AMBER at
the production restraint returns **3.207 A against the projection's 3.205 A**: a difference
of +0.001 A with a confidence interval that excludes anything larger than 0.016 A in either
direction. **Ideal geometry was not the binding constraint.** Replacing a geometric
projection with a genuine all-atom force field, which is a strictly weaker constraint,
recovers exactly the same amount of accuracy — none.

Arm B's -0.155 A reproduces the coordinator's +0.157 A to within 0.002 A and its 99W/27L
matches the 18W/108L, so the instrument is sound and the two measurements compose.

**The validity is real, not nominal.** The raw average is **45.2% clashed** with a geometry
deviation of 0.248; after AMBER it is **0% clashed** at 0.0134, which is experimental-
structure quality. And AMBER's output is **more Ramachandran-valid than the incumbent
projection's** — 0.860 against 0.734 — at identical CA-RMSD. So the correct statement is
not that AMBER is useless here: **at equal accuracy AMBER delivers a better structure than
the projection does. It just does not deliver a better number.**

### E14c. The mechanism, predicted before the run and confirmed at n=126

The prediction was written into the module docstring before any OpenMM time was spent:
averaging superposed structures CONTRACTS them, because the mean of a set of unit vectors is
shorter than a unit vector. If that is the mechanism, the +0.157 A is the price of
RE-EXPANDING a shrunken structure, *any* operation restoring valid bond lengths must pay it,
and AMBER will land near 3.2 rather than near 3.05.

Measured over 126 targets:

| bond | in the averaged backbone | ideal | deviation |
|---|---|---|---|
| N-CA | 1.1201 | 1.4580 | **-23.18%** |
| CA-C | 1.1801 | 1.5250 | **-22.62%** |
| C-N | 0.8642 | 1.3290 | **-34.97%** |
| CA-CA virtual | 2.9488 | 3.8040 | **-22.48%** |
| C=O | 0.7863 | 1.2310 | **-36.12%** |

**Mean contraction -25.81%**, per-target range [-56.85%, -0.49%], and
**rho(contraction, per-target projection cost) = -0.552** — the more contracted the average,
the more the projection costs it. The prediction is confirmed and the mechanism is
identified:

> **The 0.157 A is not the price of IDEAL geometry. It is the price of restoring ANY valid
> bond length to a structure that coordinate averaging shrank by a quarter. The coordinate
> average is 3.048 A because it is not a molecule; the moment it is required to be one, by
> either route, it costs the same 0.157 A.**

This also explains why the free relaxation (k=0) is worst at +0.272: with no restraint the
force field re-expands the structure to its own preferred geometry and drifts a further
1.526 A of CA movement away from the average, discarding the ensemble information that made
arm B good in the first place.

I record one caution against my own analysis: this correlation was -0.829 at n=6 and -0.224
at n=10 before settling at -0.552 at n=126. Only the n=126 figure should be quoted, and the
intermediate values are noted here as a reminder that a correlation on nine or fewer targets
is not evidence — the same lesson as E7-REFUTED.

### E14d. A correction to an inference drawn from my own work

The coordinator wrote that E6b ("AMBER refinement is safe but inert, +0.014 A") implied
AMBER could not recover the 0.157 A either. **E6b does not license that inference**: every
structure in E6b was already ON the ideal-geometry manifold, where AMBER has nothing to fix.
E14 is the case where the input is off-manifold, and it needed its own measurement. The
inference happens to be correct — but it was not established until now, and the two results
are independent.

Together E6b and E14 do now make the complete statement: **AMBER refinement changes CA-RMSD
by +0.014 A [-0.001, +0.028] on structures that are already valid, and by +0.001 A
[-0.014, +0.016] relative to the projection on structures that are not. In both regimes it
is accuracy-neutral and validity-positive.** That is a clean and complete characterisation
of what all-atom physical realism costs and what it buys in this system, and it is the same
verdict E8, E9 and E13c reached from three other directions.

### E14e. THE RESTRAINT FRONTIER — and a correction to my own concentration verdict

**Status: DEMONSTRATED (the frontier); my intermediate "FAILS concentration" verdict is
RETRACTED as a statistical error of mine.** `s14/ener_avgrefine.py --check`,
`s14/results/ener_avgrefine_k30.log`, `ener_avgrefine_tight50.json`.

E14b reported k=10 and concluded AMBER ties the projection. **That was one point on a
frontier, and the frontier does not turn.** The restraint constant trades RMSD against
geometric strain continuously: as k grows the relaxation is held ever closer to the raw
average, so RMSD approaches arm B's 3.050 and the geometry approaches arm B's broken 0.248.
The question is therefore not "does AMBER tie" but **"how much accuracy is available at
geometry that is still valid"** — and validity needs a measured yardstick, not a judgement.

**The yardstick, measured rather than asserted: free relaxation (k=0) IS ff14SB's own
equilibrium geometry, and it sits at 0.0170 RMS relative deviation from the ideal-builder
constants.** Any arm at or below 0.0170 is as well-formed as the force field itself wants;
anything above is strained by the force field's own standard.

| arm | RMSD | vs A | CI | W/L | geom dev | **x ff14SB equilibrium** | clashed | Rama |
|---|---|---|---|---|---|---|---|---|
| A projected | 3.205 | — | — | — | 0 by construction | — | 0.000 | 0.734 |
| B raw average | 3.050 | -0.155 | [-0.191,-0.122] | 99/27 | 0.2480 | 14.6x | **0.452** | 0.836 |
| B+AMBER k=10 | 3.207 | +0.001 | [-0.014,+0.016] | 53/73 | 0.0134 | **0.79x** | 0.000 | 0.860 |
| **B+AMBER k=30** (n=126) | **3.183** | **-0.022** | **[-0.036,-0.009]** | 67/59 | 0.0192 | **1.13x** | 0.000 | 0.874 |
| B+AMBER k=100 (n=50) | 3.101 | -0.070 | [-0.098,-0.043] | 35/15 | 0.0374 | 2.20x | 0.000 | 0.894 |
| B+AMBER k=300 (n=50) | 3.066 | -0.105 | [-0.140,-0.073] | 40/10 | 0.0665 | **3.91x** | 0.000 | 0.857 |

**So a better number IS available, and the size of it is set by how much strain is accepted.**
k=100 and k=300 recover 45% and 68% of the projection cost — but at 2.2x and 3.9x the force
field's own equilibrium strain, which is buying RMSD with distortion. **k=30 is the arm that
matters: 1.13x equilibrium, zero clashes, and -0.022 A [-0.036, -0.009] on all 126 targets.**

### E14f. RETRACTED: my "FAILS the concentration check" verdict was a statistical error

I reported to the coordinator, and they published, that k=30 fails the mandated
concentration check. **That verdict was wrong and the error was mine.** It rested on a
threshold I invented on the spot — "drop-top-10 and drop-top-20 must keep the sign and the
top-10 share must be below ~0.5" — which I never calibrated.

**A raw drop-top threshold is not a valid concentration test.** When an effect's mean is
small relative to the per-target spread, discarding the ten most favourable targets
mechanically removes a large share of the total **even if every target carries an identical
effect**. The observed statistic must be compared with its own null: a uniform effect of the
same mean and the same per-target sd (`concentration_null`, 4,000 simulations).

k=30, n=126, mean -0.0221, per-target sd 0.0789 (mean/sd = **0.28**):

| statistic | observed | **NULL: a UNIFORM effect of the same mean and sd** | percentile |
|---|---|---|---|
| drop-top-10 | -0.0080 | -0.0094 [-0.0234, +0.0048] | 0.62 |
| drop-top-20 | +0.0028 | **+0.0006 [-0.0136, +0.0150]** | **0.62** |
| top-10 share | 0.668 | **0.698 [0.389, 1.543]** | **0.62** |

**Every statistic sits at the 62nd percentile of its own null.** A perfectly uniform effect
of this size and noise would produce a drop-top-20 that reverses sign and a top-10 share of
0.70. There is **no evidence of concentration in any arm**, including the n=50 case on which
I based the retraction (observed share 1.163 against a null mean of **1.410**, i.e. *less*
concentrated than a uniform effect would typically look):

| arm | mean/sd | share observed / null | drop-top-20 observed / null | verdict |
|---|---|---|---|---|
| k=30, n=126 | 0.28 | 0.668 / 0.698 | +0.0028 / +0.0006 | not concentrated |
| k=30, n=50 | 0.33 | 1.163 / 1.410 | +0.0243 / +0.0236 | not concentrated |
| k=100, n=50 | 0.70 | 0.678 / 0.615 | -0.0007 / -0.0061 | not concentrated |
| k=300, n=50 | 0.85 | 0.553 / 0.535 | -0.0163 / -0.0263 | not concentrated |

**The symmetric caution, which matters as much as the retraction.** "Not concentrated" here
does **not** mean "demonstrated uniform". At mean/sd = 0.28 the drop-top test has almost no
power and **cannot settle the question either way**. What actually carries the result is the
CI on the mean and the per-fold consistency — and the per-fold table is the genuinely
informative one:

| fold | n | A | B+AMBER k=30 | diff | W/L |
|---|---|---|---|---|---|
| 0 | 25 | 2.991 | 2.969 | -0.0213 | 15/10 |
| 1 | 23 | 3.169 | 3.155 | -0.0138 | 12/11 |
| 2 | 25 | 3.116 | 3.100 | -0.0165 | 12/13 |
| 3 | 23 | 3.472 | 3.452 | -0.0195 | 11/12 |
| 4 | 30 | 3.282 | 3.246 | -0.0357 | 17/13 |

**All five pinned folds carry the same sign with only a 2.6x range** (-0.014 to -0.036), and
FAIL18 moves with the rest (-0.011 on 18, -0.024 on the other 108). That is the real evidence
for uniformity, and it is independent of the drop-top machinery. At n=50 the same table
showed a 25x fold range; the spread narrowed as n grew, which is what noise does.

**Two corrections carried forward from this episode, both mine:**

1. **The W/L is not the discriminating statistic — but neither is a raw drop-top number.**
   I first read 26W/24L as reassurance (wrong: the median-vs-mean gap already said the mean
   was tail-driven), then read an uncalibrated drop-top curve as disqualifying (also wrong:
   that is what a uniform effect looks like at this signal-to-noise). **The correct test is
   the drop-top statistic against its own uniform-effect null, and it should be reported
   with mean/sd beside it so a reader can see when the test has no power.**
2. This project's mandated concentration check, applied literally, would have discarded a
   real effect. It caught `leg_torsion` (mean/sd was large and the effect genuinely sat on
   two of nine targets) and it misfired here. **The check is necessary and not sufficient;
   it needs the null.**

### E14g. What the k=30 result is, and what it is not

**Baseline check (this comparison is the CONSERVATIVE one).** `I.project` returns two arms:
`fit_ca` (lam=0) at 3.2052 and `ca` (lam=0.3, Ramachandran-penalised) at 3.2126. Arm A above
uses `fit_ca`, the BETTER of the two, matching `avgspace.py`. Against the lam=0.3 arm the
k=30 gain is larger: **-0.0295 A [-0.0465, -0.0122], 70W/56L**. So -0.022 A is the
conservative figure, not a favourable choice of baseline.

**IS:** a native-free, post-processing change to the incumbent pipeline that emits
**3.183 A against 3.205 A on the 126-target development instrument, -0.022 A
[-0.036, -0.009]**, consistent in sign across all five pinned folds, with no evidence of
target concentration, and with **strictly better stereochemistry than the incumbent** — zero
clashes against the raw average's 45.2%, geometry at 1.13x the force field's own equilibrium,
and Ramachandran validity **0.874 against the projection's 0.734**.

**IS NOT:**

* **not benchmark-validated.** This is `tuning126`, the development instrument. Nothing here
  has been near the held-out benchmark, and the project's history is that no dev improvement
  has ever transferred (the last attempt: +0.0103 A, CI [-0.1596, +0.1803]).
* **not free of a selection concern, and the concern is WORSE than I first stated it.**
  I told the coordinator that the rule used to pick k=30 -- match the force field's own
  equilibrium strain -- is native-free and would have selected k=30 blind. **Checked, it
  would not: the rule as stated is DEGENERATE.** k=0 has strain 0.0170 = exactly 1.00x
  equilibrium *by definition*, because k=0 IS the equilibrium, so "closest to equilibrium"
  trivially selects free relaxation -- which emits 3.477 A, the worst arm on the frontier.

  | k | strain | x equilibrium | \|dist from 1.0x\| | RMSD |
  |---|---|---|---|---|
  | 0 | 0.0170 | **1.00x (circular)** | **0.00** | 3.477 |
  | 2 | 0.0152 | 0.89x | 0.11 | 3.248 |
  | 10 | 0.0134 | 0.79x | 0.21 | 3.207 |
  | **30** | 0.0192 | 1.13x | 0.13 | **3.183** |
  | 100 | 0.0374 | 2.20x | 1.20 | 3.101 |
  | 300 | 0.0665 | 3.91x | 2.91 | 3.066 |

  The variant "the tightest restraint whose strain does not EXCEED equilibrium" is
  well-formed and non-circular, but it selects **k=10**, which emits +0.001 A -- no gain at
  all. Selecting k=30 requires a threshold of ~1.15x equilibrium, and 1.15 was chosen with
  the frontier visible. **So: the strain yardstick legitimately RULES OUT k=100 and k=300
  as strained, but it does not select k=30 out of {k=2, k=10, k=30}. k=30 was chosen on
  development-set RMSD from within the geometrically acceptable set.** That is dev-set
  selection over a small grid. It is mild, but it is real, and it means a pre-registered
  k could land anywhere from +0.001 (k=10) to -0.022 (k=30). **A confirmation must
  pre-register the restraint constant; there is no native-free rule in hand that picks it.**
* **not large.** 0.022 A is 0.7% of a 3.205 A baseline, and the brief's own standard is that
  0.02 A is not an improvement without statistical support. It has that support here; it is
  still 0.022 A, and it is a post-processing change, not an architectural one.
* **not a contradiction of E6b.** E6b measured AMBER refinement on structures already ON the
  ideal-geometry manifold, where it is inert (+0.014 A). E14 measures it on structures OFF
  the manifold, where the projection is the alternative. Both are true and they are
  different experiments.

**The mechanism remains E14c's.** All of the accuracy on this frontier comes from *not
fully re-expanding* a coordinate average that shrank by 25.8%. AMBER's advantage over the
projection is that it can hold a structure at 1.13x its own equilibrium strain, whereas the
projection must go all the way to exact ideal geometry. **The gain is a strain allowance,
not physics finding a better structure** — which is why it grows monotonically with the
restraint and why the whole frontier is bounded by arm B's 3.050 A.

---

# E11. VERDICT — the legitimate role of each model

Every number below is on the nine fully-enumerated n=9, k=4 targets, on a uniform ensemble
of ~1,194 configurations per target unless stated. **Scope limit, stated once and meant:
n=9 at k=4. Nothing here establishes how either model behaves at larger n or larger k, and
nine targets cannot support the drop-top-10 concentration check this project mandates.**

## LEGACY — a clash gate with a 1.31 A resolution limit. Not a ranking function, and an
## actively unsafe filter above 10% removal.

* Its variance is **98.8% one term** (`steric`), which is **exactly zero on 73.5% of the
  space and constant across its entire lowest-energy decile** — a binary gate with a 73.5%
  pass rate, carrying no ordering information inside the pass set (E1a, E3).
* It needs a **1.31 A quality gap** to reach 55% pairwise decision accuracy and reaches
  0.504 at gaps below 0.25 A (E2b).
* It places the ORACLE snap at the **62nd percentile** of its own energy on a
  Ramachandran-plausible population — the wrong side of a coin flip (E9a).
* Presented with a good structure and a decoy made from it, it **prefers the decoy** at
  every separation, worst (-0.19 excess AUC) at 2 A (E2c).
* Coordinate descent on it takes a **1.295 A structure to 3.684 A** (E6a); its certified
  global optimum is 3.920 A against a whole-space mean of 3.781 A -- which decomposes
  exactly as FILTERING -0.124 plus ORDERING **+0.263 [+0.046, +0.484]** (E13c). Its tail is
  better than the pool; its ranking inside that tail throws the gain away and more.
* It is **blindest inside its own tail**: 0.462 gap-matched pairwise accuracy in its lowest
  1%, against 0.523 for Legacy on a random subset of the same size (E13b).
* Reweighting its eleven terms: in-sample rho +0.209, **leave-one-target-out -0.031** (E3a).
* As a 10% filter it is worth **-0.048 A [-0.085, -0.009]** on the pool mean, but it costs
  **+0.219 A on the pool best by 25% removal** and +0.385 A by 75% (E9).
* Its apparent -0.296 A selection advantage over a random draw has a **top-2 share of
  1.02** — drop the two best targets and it is +0.009 (E7).

**Legitimate role: a 10%-removal geometric veto, and nothing else.** It must not be a
search objective, a ranking objective, a reranker, or a filter above ~10%. Its steric term
used alone is a slightly better version of the same gate (-0.058 A [-0.104, -0.010], 7W/2L,
and it does not damage the pool best), so if Legacy is kept for the veto, **keep only the
steric term**: `legacy_nosteric` is the worst filter at every removal fraction tested.

## AMBER — a genuine physical validator. The best clash detector measured, safe to refine
## with, and useless as a ranking or search objective.

The positives, all measured and all real:

* **The best clash detector in the study.** The uniform ensemble is 13.07% clashed with
  heavy-atom separations down to 0.23 A; AMBER's lowest decile is **0.09% clashed with a
  minimum separation of 1.98 A**, against Legacy's 0.56% / 1.62 A (E8).
* **The only physical objective with a positive excess decoy AUC** against a random-anchor
  null: +0.056 at 2 A rising to +0.113 at 4 A (E2c).
* **The only filter that never damages the top of the pool.** It can discard 90% of the
  ensemble and still hold the best structure in it (-0.033 A vs the random-removal null at
  90%), while Legacy has lost 0.380 A of the best by then (E9).
* **The best FILTERING term of any objective measured** (-0.141, and -0.163 for
  `amb_nonbonded`, against Legacy's -0.124) in the exactly-additive selection decomposition
  (E13c) -- a third independent instrument reaching the validator conclusion, after the
  clash audit (E8) and the filter curve (E9).
* **Refinement is safe.** `refine_coords(k=10, steps=0)` moves CA-RMSD by
  **+0.014 A [-0.001, +0.028]** over 45 (target, start) cells, and from the ORACLE pool best
  at 1.295 A returns 1.321 A (E6b).
* **It can replace the ideal-geometry projection and beat it, at valid geometry**
  (E14e, n=126, the full development instrument): restrained AMBER at k=30 on the raw
  coordinate average gives **3.183 A against the projection's 3.205 A, -0.022 A
  [-0.036, -0.009]**, same sign in all five pinned folds, no evidence of target
  concentration against the proper null, **0% clashes against the raw average's 45.2%**,
  geometry at **1.13x the force field's own equilibrium strain**, and Ramachandran validity
  **0.874 against the projection's 0.734**. **The one candidate deployable result in this
  file: if a physically valid structure is wanted, AMBER-relax the coordinate average
  rather than projecting it** -- better accuracy AND strictly better stereochemistry,
  ~12 s per target. Caveats in E14g: dev-instrument only, never benchmark-validated, 0.7%
  of the baseline, and the restraint constant was chosen with the frontier visible even
  though the rule for choosing it is native-free.

The negatives, equally real:

* Its variance is **100.0% one term** (`nonbonded`; 99.5% after monotone conditioning).
  Under ideal backbone geometry `amb_bond` is CONSTANT and `amb_angle` takes five values, so
  the five-term decomposition is three terms with one at 99.5% (E3).
* It needs a **1.58 A quality gap** for 55% accuracy and **3.87 A for 60%** — the coarsest
  resolution of any objective that resolves at all (E2b).
* It places the ORACLE snap at the **46th percentile** — exactly no information (E9a).
* Coordinate descent on it takes a **1.295 A structure to 3.685 A** (E6a).
* **Refinement never rescues a bad structure**: from the prior argmin at 4.334 A it returns
  4.337 A. All variation in the refinement column comes from the start (E6b).
* AMBER pre-filtering followed by Legacy ranking is **worse than random pre-filtering**
  (3.719 vs 3.641), and rank(Legacy)+rank(AMBER) is 0.015 A worse than Legacy alone (E1).
* Single points are reproducible only to **2.4e-4 relative** under a rigid motion — up to
  1.6e14 kcal/mol absolute on clashing configurations (E8a).
* It spends **~4.5x more of its variance than Legacy on the two residues CA-RMSD cannot
  see** (0.176 vs 0.039 per invisible residue relative to an interior one) (E10b).

**Legitimate role: post-hoc all-atom physical validation, and refinement for physical
validity — never for accuracy, never in a search loop, never as a reranker.** As a validation
filter it is worth -0.021 A [-0.041, -0.002] on the pool mean at 10% removal and is safe at
any removal fraction. At 28 ms per distinct single point and ~4 s per refinement (against
0.3 ms for Legacy) it should run **once, at the end, on the candidates that are already
selected**, and its output should be a pass/fail flag and a refined structure, not a score.

## THE TWO TOGETHER

* **They are genuinely near-decorrelated** (truth-partialled error correlation +0.096) —
  hypothesis H7's premise is CORRECT.
* **The decorrelation is worth at most +0.032 rho** under a per-target-optimal weight, which
  is an oracle ceiling, and the identical machinery loses 0.240 rho going from in-sample to
  leave-one-target-out. In the actual matrix every combination is inside its own noise and
  the best is 0.015 A worse than Legacy alone. **H7 is REFUTED as an exploitable route.**
* **There is no Pareto conflict worth exploiting**: mean |axis correlation| 0.222, the
  frontier is 14.6% of the pool and beats a same-size random subset by 0.048 A while
  discarding 72% of the truly-best 1%. **Do not build multi-objective machinery.**
* **All scale-based interventions are provably no-ops.** Rank, softcore, robust-z and affine
  transforms give bit-identical rho, decile rho, pair accuracy and selected RMSD while
  collapsing AMBER's range from 16.32 decades to 1.47. The damage is in the ordering.

## WHAT PHYSICAL REALISM COSTS

The projection onto ideal geometry costs 0.157 A (coordinator, `avgspace.py`), and E14 shows
that is **not the price of ideal geometry** but the price of restoring ANY valid bond length
to a coordinate average that coordinate averaging shrank by **25.8%** --
rho(contraction, per-target projection cost) = -0.552, n=126. A genuine force field, a
strictly weaker constraint than the ideal-geometry manifold, pays exactly the same 0.157 A.
Every torsion-space method pays it before it starts, because a torsion representation IS
that manifold. **But the price is not fixed: it is a function of how much geometric strain
is allowed** (E14e). The projection must go all the way to exact ideal geometry; a force
field can stop at its own equilibrium, and that strain allowance is worth -0.022 A at 1.13x
equilibrium, -0.070 A at 2.2x and -0.105 A at 3.9x. The frontier is bounded by the raw
average's 3.050 A, and **all of it is the un-expansion of a shrunken structure, not physics
locating a better one.**

## THE ONE NUMBER

**No objective in this project — physical or otherwise — exceeds 0.511 pairwise decision
accuracy when two structures differ by less than 0.25 A in RMSD. Legacy needs a 1.31 A gap
to reach 0.55; AMBER needs 1.58 A and 3.87 A to reach 0.60.** The useful range of the search
on this instrument is about 2.5 A wide (pool mean 3.789, pool best 1.295), so both physical
energies can resolve roughly one bit of it. That single threshold explains the certified
global optimum being worse than random, the descent from 1.295 A to 3.68 A, and every
failed reranker in the project's history.

## WHAT I REFUTED, INCLUDING MY OWN HYPOTHESES

1. **REFUTED (mine): "the totals destroy signal the components carry."** LOTO reweighting of
   all eleven Legacy terms goes from in-sample rho +0.209 to **-0.031** held out, and selects
   4.036 A against Legacy's own 3.493 (E3a).
2. **REFUTED (mine): "AMBER refinement turns poor structures into good ones."**
   +0.014 A [-0.001, +0.028], 13W/32L, and all variation comes from the start (E6b).
3. **REFUTED (mine): "Legacy's decoupling is driven by the invisible terminal DOFs."**
   Legacy spends 3.8% of an interior residue's variance per invisible residue (E10b).
4. **REFUTED: H7 as an exploitable route.** Decorrelation is real (+0.096 partialled error
   correlation) and worth +0.032 rho under an oracle weight; unattainable (E4).
5. **REFUTED: the Pareto hypothesis.** Not enough conflict; the frontier discards 72% of the
   good structures (E4a).
6. **REFUTED: the naive decoy experiment.** Discrimination is NOT a function of structural
   separation; it is a function of the quality gap, and the apparent separation curve is the
   gap confound (E2a).
7. **REFUTED (mine, twice, in my own instruments):** the first fusion pricing used the
   truth-partialled correlation in the multiple-correlation identity, which needs the raw
   correlation, and produced fictitious rho_fused values up to +0.867; the first geometry
   audit measured omega as a raw circular mean (reading 38 deg for a rigid 180 deg trans
   bond) and counted bonded pairs as clashes (reporting 100% of structures clashed). Both
   are corrected in place and both corrections are recorded here.
8. **CORRECTED, in the sprint brief:** the CA trace has **four** inert torsions, not three —
   `psi[0]` is inert as well, so **residue 0 is entirely invisible too** and live qubits at
   n=9, k=4 are **14, not 16** (E10a).
9. **CORRECTED, in the sprint brief:** Legacy's `steric` covariance share is **0.988**, not
   0.955, and AMBER's `nonbonded` share survives monotone conditioning at 0.995 (E3).
10. **PINNED, not previously stated:** the cached AMBER subset is 40% ORACLE-conditioned;
    only `amber_kind == 0` may be used for native-free statistics (E0).

## WHAT IS STILL OPEN

* ~~`leg_torsion`~~ **CLOSED — REFUTED at n=126, see E7-REFUTED.** It was the only open
  lead; it is the worst arm the coordinator tested, anti-ranks in the low-energy decile
  under both proposals, and poisons every combination it enters. My own concentration
  warning (two of nine targets, 68% of the effect) called it correctly, and E2b had already
  placed it in the "never reaches 0.55 pairwise accuracy at any gap" row. **Nothing in the
  ENER matrix is now open on the positive side: no physical energy, and no component or
  combination of one, ranks structures in this space.**
* Everything here is n=9, k=4. Whether the 1.31 / 1.58 A resolution thresholds scale with n
  is untested and is the most valuable follow-up.
* AMBER *with minimisation* as a ranking objective (not `single_point`) remains untested as
  a ranker; E6b tests it only as a refiner. At ~4 s per structure it cannot be a search
  objective regardless of how it ranks.

---

# E12. REPRODUCTION

Seed 20260905 everywhere; deterministic given the seed. Run from the repository root with
`PYTHONIOENCODING=utf-8`. Thread caps are set at import by `s13/qarch_lib.py`.

    python -m s12.instrument          # confirm the pinned constants first
    python -m s14.ener_subset         # E0   AMBER subset provenance          ~10 s
    python -m s14.ener_matrix         # E1   arms A-E, H, I + degeneracy      ~3 min
    python -m s14.ener_decoy          # E2   decoy resolution, no AMBER       ~6 min
    python -m s14.ener_decoy --amber  # E2   with genuine AMBER               ~25 min
    python -m s14.ener_complement     # E3/E4 decomposition, H7, Pareto       ~5 min
    python -m s14.ener_norm           # E5   normalisation constants          ~2 min
    python -m s14.ener_refine         # E6   arms F, G, J                     ~20 min
    python -m s14.ener_nulls          # E7   null controls                    ~1 min
    python -m s14.ener_geom           # E8   geometric audit                  ~4 min
    python -m s14.ener_filter         # E9   filter curve                     ~1 min
    python -m s14.ener_invisible --amber  # E10 invisible-DOF Sobol           ~12 min
    python -m s14.ener_tail           # E13  tail-restricted discrimination   ~4 min
    python -m s14.ener_avgrefine      # E14  AMBER vs the projection, n=126   ~100 min

Only `ener_decoy --amber`, `ener_refine`, `ener_geom` and `ener_invisible --amber` touch
OpenMM; each gates on `s13.qarch_lib.wait_for_memory(1.5)` and none of them catches a
`MemoryError`. Everything else runs on the cached enumeration at zero AMBER cost.

## Artefacts for the other workstreams

* **`s14/ener_lib.py`** — `enum(pdb)` (cached enumerated target with the AMBER strata),
  `decile_rho`, `pair_accuracy`, `argmin_rmsd` (**tie-averaged** — use this, not
  `np.argmin`), `topk_recovery`, `calibration`, `rank_norm`, `robust_z`, `softcore`,
  `pairwise_rmsd`.
* **`s14/ener_decoy.py`** — `decoys(pdb, oracle_anchor=True|False)` returns the cached
  controlled decoy sets: 8 anchors x 9 mutation counts x 60 decoys with the true mutual
  CA-RMSD and the true RMSD of every member, plus the matched random-anchor NULL set.
  Genuine AMBER on a deterministic 1-in-3 stride is cached alongside. **OBJ, VQE and SHIFT
  should train and evaluate resolution on these rather than building their own.**
* **`s14/cache/ener_norm.json`** — derived normalisation constants. Use `grad_sd` (the exact
  mean single-flip standard deviation) to put terms on the scale a local search experiences.
  **Never** use the sd of raw AMBER. For AMBER use rank normalisation.
* **`s14/cache/decoy_*.npz`, `s14/cache/decoy_amber_*.npz`** — the decoy sets and their
  AMBER energies.

## Warnings for anyone building on this

1. **Use `amber_kind == 0`** for any native-free AMBER statistic. The full 2,955-config
   subset is 0.401 A better in RMSD than the space it came from.
2. **Tie-average every argmin.** `leg_steric` is tied at its minimum on 73.5% of the
   ensemble and constant across its whole lowest-energy decile; `amb_bond` is constant;
   `amb_angle` takes five values. `np.argmin` on any of these reads the array order.
3. **Score per-residue channels on residues 1..n-2.** Both terminal residues are entirely
   invisible to CA-RMSD; at n=9 that is 22% of a 1-local objective's decisions.
4. **Do not difference two raw AMBER energies.** They are reproducible only to 2.4e-4
   relative.
5. **Do not propose a normalisation fix for AMBER.** Every monotone transform is bit-identical
   on every ranking metric; this is measured, not argued.
