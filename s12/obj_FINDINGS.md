# S12 — OBJECTIVE AND ERROR STRUCTURE

**Question.** Is the CA–CA distance matrix the right objective, and is scalar accuracy (MAE)
or ERROR STRUCTURE what determines the structure the pipeline emits?

**Instrument.** tuning126 only. benchmark60 never read; dev24 never run. `rr`/`nat_ca` used for
evaluation and (where stated) as leave-fold-out labels only. All corrupted-objective arms are
labelled ORACLE/DIAGNOSTIC — they start from the native distance matrix and are not deployable.

Code: `s12/obj_common.py` (pool pack, scorers, corruptors, the path), `s12/obj_run.py`
(parallel driver), `s12/obj_errstruct.py` (§1), `s12/obj_iso.py` (§2 analysis),
`s12/obj_alt.py` (§4 scorers + LFO fits), `s12/obj_arm4.py` (§4 driver),
`s12/obj_profile.py` (§6), `s12/obj_suff.py` (§5, §5b, §7, §7b, §8), `s12/obj_report.py`.
Results: `obj_repro`, `obj_errstruct`, `obj_iso_cheap`(+`_summary`), `obj_iso2`, `obj_arm4`,
`obj_arm4b`, `obj_profile`, `obj_suff`, `obj_confirm` — all under `s12/results/`.

---

## 0. Instrument reproduction and the OBJECTIVE CEILING (`obj_repro.json`)

All 126 targets through the production terminal path (score K=500 → top-75 → superpose on
medoid → coordinate average → L-BFGS projection onto ideal geometry).

| objective | MAE | r | emitted (λ=0 synthesis) | FAIL18 | other-108 | frac<2Å | top-75 best | argmin |
|---|---|---|---|---|---|---|---|---|
| shipped Bayes risk (baseline) | 2.339 | 0.696 | **3.205** (record 3.204) | 6.034 | 2.734 | 0.28 | 2.306 (record 2.306) | **3.454** (record 3.454) |
| L1 vs E[d] (matched scorer) | 2.339 | 0.696 | 3.243 | 5.966 | 2.789 | 0.26 | 2.324 | 3.509 |
| **ORACLE: L1 vs the TRUE matrix** | 0 | 1 | **2.395** | 3.652 | 2.185 | 0.38 | 1.711 | 1.994 |

Reproduction is exact. Two things follow immediately:

* **The scoring FORM is worth nothing.** Bayes-risk vs plain L1 against the same expectation:
  Δ = +0.038 Å [−0.005, +0.083], 53W/73L. The whole 17-bin distribution buys ≈ 0. All
  corruption arms below therefore use plain L1 without prejudice.
* **A PERFECT CA–CA distance objective emits 2.395 Å, not 2.0 Å.** Knowing the native
  distance matrix exactly and scoring the shipped K=500 pool with it is worth −0.811 Å
  [−0.994, −0.640], 116W/10L, drop-top-10 −0.597 — large, uniform, real. But it stops
  0.395 Å short of the 2.0 Å goal, and 0.68 Å short of the top-75 best member it itself
  selects (1.711). **The distance-matrix objective is a binding channel worth ≈0.81 Å, and
  it is NOT sufficient**: even at zero error the pool + terminal operator lose the rest.
* On FAIL18 a perfect objective is worth −2.38 Å (6.03 → 3.65) but still leaves them at
  3.65 Å against a 2.28 Å pool best. FAIL18 is only ~2/3 an objective problem.

### 0b. The coordinate average is a near-exact proxy for the emitted structure

Validated so the corruption sweeps below can be run without paying 2.5 s of L-BFGS per cell:

| condition | avg-RMSD | projected fit-RMSD | Pearson r (126 targets) |
|---|---|---|---|
| bayes | 3.048 | 3.205 | 0.9946 |
| L1 vs E[d] | 3.078 | 3.243 | 0.9953 |
| oracle | 2.218 | 2.395 | 0.9887 |
| Δ(oracle − bayes) | −0.830 | −0.811 | 0.9847 |

Projection adds a near-constant +0.16 Å. Sweeps are therefore reported on the coordinate
average (`avg_rmsd`); every headline row is confirmed through the real projection.

---

## 1. Error structure of the shipped distogram (`obj_errstruct.py`, no projection)

126 targets, pairs with |i−j| ≥ 2, e = expected − true.

| arm (all matched to the distogram's per-target MAE) | MAE | r | slope t~p | slope p~t | **slope_sep** | **r_sep** | eff.rank | top1 | scale_frac | sign_agree | MDS RMSD |
|---|---|---|---|---|---|---|---|---|---|---|---|
| shipped distogram | 2.339 | 0.696 | 0.724 | 0.727 | **0.198** | **0.196** | 8.16 | 0.405 | 0.157 | 0.607 | 4.621 |
| null: native + iid Gaussian | 2.339 | 0.768 | 0.626 | 0.988 | 0.277 | 0.484 | 10.34 | 0.256 | 0.018 | 0.500 | 4.602 |
| null: native + random sign, fixed magnitude | 2.339 | 0.824 | 0.708 | 0.996 | 0.363 | 0.559 | 10.54 | 0.241 | 0.024 | 0.503 | 3.978 |
| null: native + displacement field L=3.0 | 2.339 | 0.802 | 0.723 | 0.927 | 0.377 | 0.517 | 5.79 | 0.523 | 0.090 | 0.690 | 3.307 |
| null: native + displacement field L=0.5 | 2.339 | 0.753 | 0.678 | 0.858 | 0.294 | 0.461 | 7.93 | 0.431 | 0.052 | 0.618 | 3.875 |
| null: native + rank-1 pair error | 2.339 | 0.721 | 0.577 | 0.956 | 0.224 | 0.418 | 5.36 | 0.714 | 0.020 | 0.499 | 4.607 |
| null: real error marginal, pairs shuffled | 2.339 | 0.784 | 0.644 | 0.999 | 0.274 | 0.484 | 9.94 | 0.322 | 0.018 | 0.582 | 4.563 |

Definitions. `slope_sep` / `r_sep`: calibration slope and correlation of prediction against truth
**after removing everything predictable from sequence separation |i−j| alone** (per-shell means
subtracted from both). `eff.rank`: entropy participation rank of the eigenvalues of the symmetric
n×n error matrix (iid ≈ full rank, displacement field ≈ low). `top1`: Frobenius share of the
leading eigenvector. `scale_frac`: fraction of the squared error explained by a single global
"expand/contract all distances" mode. `sign_agree`: fraction of pairs SHARING A RESIDUE whose
errors agree in sign (chance = 0.50). `MDS RMSD`: CA-RMSD to native of the classical-MDS 3-D
embedding of the predicted matrix — i.e. **the objective's own geometric optimum**.

Readings:

1. **The distogram carries almost no pair-specific information.** Once |i−j| is partialled out,
   r_sep = 0.196 and slope_sep = 0.198. Every matched-MAE corruption of the NATIVE matrix keeps
   2–3× more within-shell signal (0.42–0.56) at the same MAE. MAE is therefore not even measuring
   the same thing across these arms — most of the distogram's apparent accuracy is the sequence-
   separation prior, which every candidate window in the pool already satisfies.
2. **The error is structured, and it looks like a displacement field.** sign_agree 0.607 vs 0.500
   for iid noise, eff.rank 8.16 vs 10.34 for iid, top1 0.405 vs 0.256. A coordinate displacement
   field with correlation length ≈ 0.5–3 residues reproduces those numbers (0.618–0.690, 5.8–7.9,
   0.43–0.52). 15.7% of the squared error is a single global scaling mode (iid: 1.8%).
3. **The objective's own optimum is 4.62 Å.** The best 3-D structure consistent with the
   predicted distance matrix is worse than what the pipeline actually emits (3.20 Å). The pipeline
   beats its own objective because retrieval restricts it to real peptide geometry.

The record's "calibration slope +0.376" is reproduced exactly as the **within-shell, pooled-across-
targets** slope (`obj_alt.fit_lfo`, per fold: 0.29–0.54, mean ≈ 0.40, `obj_lfo.npz`). Within a
single target it is 0.198. Both say the same thing: **across targets the model does track overall
compactness; within a target it barely orders the pairs.**

---

## 2. The ISO-MAE surface (`obj_iso_cheap`, `obj_iso2`; 126 targets, coordinate-average arm)

Native distance matrix corrupted by parameterised noise, amplitude bisected to a requested MAE,
pushed through the real path. 3 seeds per cell. ORACLE/DIAGNOSTIC arms throughout.

Emitted RMSD (coordinate average) vs requested MAE:

| family | 0.25 | 0.5 | 1.0 | 1.5 | 2.339 | 3.5 |
|---|---|---|---|---|---|---|
| iid Gaussian | 2.217 | 2.220 | 2.222 | 2.235 | **2.258** | 2.315 |
| random sign, fixed magnitude | 2.221 | 2.217 | 2.230 | 2.218 | **2.268** | 2.426 |
| heteroscedastic (σ ∝ d) | 2.221 | 2.221 | 2.226 | 2.237 | **2.271** | 2.331 |
| rank-1 pair error | 2.219 | 2.220 | 2.221 | 2.234 | **2.240** | 2.291 |
| rank-3 pair error | 2.217 | 2.223 | 2.221 | 2.242 | **2.266** | 2.310 |
| noise on \|i−j\| ≤ 4 only | 2.225 | 2.223 | 2.238 | 2.254 | **2.261** | 2.271 |
| noise on \|i−j\| > 4 only | 2.222 | 2.220 | 2.234 | 2.228 | **2.271** | 2.321 |
| real error marginal, pairs shuffled | 2.213 | 2.208 | 2.206 | 2.219 | **2.287** | 2.429 |
| displacement field L=0.25 | 2.222 | 2.215 | 2.244 | 2.240 | **2.358** | 2.597 |
| displacement field L=1 | 2.222 | 2.224 | 2.225 | 2.258 | **2.430** | 2.811 |
| displacement field L=3 | 2.226 | 2.235 | 2.263 | 2.332 | **2.589** | 3.210 |
| displacement field L=8 | 2.230 | 2.246 | 2.289 | 2.379 | **2.715** | 3.462 |
| global expansion | 2.203 | 2.205 | 2.219 | 2.269 | **2.508** | 3.045 |
| global compression | 2.250 | 2.272 | 2.345 | 2.469 | **2.743** | 3.184 |
| shrink to the target's mean distance (slope .376) | — | — | — | — | **2.846** | — |
| REFERENCE: oracle (MAE 0) | | | | | **2.218** | |
| REFERENCE: **the real distogram** | | | | | **3.048** | |
| REFERENCE: random score | | | | | **3.445** | |

**The headline.** At the shipped MAE of 2.339 Å the emitted RMSD ranges from 2.24 to 2.85 Å
depending only on the SHAPE of the error — a 0.61 Å iso-MAE spread — and **the real distogram sits
at 3.05 Å, outside and worse than the entire synthetic family.** MAE does not price this objective:
2.34 Å of i.i.d. distance error costs 0.04 Å of structure; 2.34 Å of the distogram's own error
costs 0.83 Å, a factor of 20.

Corollary: the "native + random-sign residuals beats the real distogram at the same MAE" clue is
confirmed and generalised — *every* symmetric, pair-local corruption at the same MAE beats it. The
distinguishing property is not sign and not magnitude; it is coherence along sequence separation.

---

## 3. WHICH property of the error is decisive (`obj_iso2`)

The real distogram's error e = E[d] − d_true is split into
`e = e_sep + e_res`, where e_sep is its per-\|i−j\| mean and e_res the within-shell residual, and
each component is added back to the NATIVE matrix on its own (ORACLE/DIAGNOSTIC):

| arm | MAE | emitted (avg) | FAIL18 | other-108 |
|---|---|---|---|---|
| native (oracle objective) | 0.000 | 2.218 | 3.381 | 2.024 |
| native + e_res only (pair-specific error) | 1.532 | **2.402** | 3.715 | 2.183 |
| native + e_res shuffled inside its shell | 1.521 | 2.223 | 3.389 | 2.028 |
| native + e_sep only (separation-profile error) | 1.679 | **2.846** | 5.267 | 2.443 |
| native + e_sep with shells shuffled | 2.364 | 2.563 | 4.265 | 2.279 |
| native + whole error, pairs shuffled | 2.279 | 2.356 | 3.735 | 2.126 |
| native + whole error (= the real prediction) | 2.339 | 3.078 | 5.756 | 2.632 |
| shrink to per-shell mean, slope 0.376, no noise | 0.450 | 2.285 | 3.462 | 2.089 |
| the prediction's shell profile only, pair detail deleted | 2.178 | 3.203 | 5.675 | 2.791 |

Readings:

1. **73 % of the damage is the separation-profile component.** e_sep alone (MAE 1.68) costs
   0.63 Å; e_res alone (MAE 1.53, a LARGER perturbation of the same order) costs 0.18 Å.
2. **Shuffling destroys the damage.** The identical error marginal reassigned to random pairs costs
   0.14 Å instead of 0.86 Å — **84 % of the damage is in WHICH pair gets WHICH error**, not in how
   big the errors are. Even e_res, shuffled inside its own shell, becomes free (2.223 vs 2.402).
3. **Within-shell shrinkage is not the defect.** Regressing every distance to its own shell mean
   with slope 0.376 — the exact calibration pathology named in the record — costs 0.07 Å.
   The pathology that matters is getting the shell MEANS wrong, not the within-shell slope.
4. The prediction's entire pair-specific detail is worth 0.155 Å (3.203 → 3.048): the shipped
   distogram is, operationally, a predictor of a 12-to-14-number curve.

---

## 4. ARM 4 (coordinator request) and the alternative-objective screen (`obj_arm4`, `obj_arm4b`)

All 126 targets, coordinate-average arm. `rho` = Spearman(score, true RMSD) over the K=500 pool;
`rhoBnd` = the same inside the near-native band; `natPct` = percentile of the NATIVE structure under
the score (0 = the score's argmin, 0.5 = no better than chance); `recall` = fraction of the band
kept by the top-75.

| scorer | emitted | FAIL18 | other | top-75 best | argmin | rho | rhoBnd | natPct | recall |
|---|---|---|---|---|---|---|---|---|---|
| ORACLE full L1 | **2.218** | 3.381 | 2.024 | 1.717 | 1.994 | 0.840 | 0.579 | 0.000 | 0.560 |
| ORACLE shell-profile only (n−2 numbers) | **2.299** | 3.494 | 2.100 | 1.739 | 2.505 | 0.787 | 0.350 | 0.000 | 0.494 |
| ORACLE separation-partialled deviation | **2.307** | 3.541 | 2.102 | 1.733 | 2.498 | 0.758 | 0.330 | 0.000 | 0.500 |
| ORACLE within-shell residual only (arm-4 rule, true matrix) | 3.030 | 4.989 | 2.704 | 1.874 | 2.520 | 0.472 | 0.509 | 0.000 | 0.361 |
| shipped Bayes risk | 3.048 | 5.832 | 2.584 | 2.306 | 3.454 | 0.568 | 0.126 | 0.368 | 0.302 |
| 1/sd-weighted L1 | 3.047 | 5.824 | 2.585 | 2.289 | 3.466 | 0.564 | 0.128 | 0.369 | 0.302 |
| LFO per-shell weighted L1 | 3.058 | 5.726 | 2.614 | 2.261 | 3.490 | 0.541 | 0.100 | 0.364 | 0.293 |
| plain L1 vs E[d] | 3.078 | 5.756 | 2.632 | 2.324 | 3.509 | 0.559 | 0.080 | 0.384 | 0.286 |
| scale-invariant L1 | 3.104 | 5.615 | 2.686 | 2.323 | 3.594 | 0.505 | 0.112 | 0.404 | 0.288 |
| global-offset-removed L1 | 3.121 | 5.698 | 2.692 | 2.373 | 3.484 | 0.564 | 0.076 | 0.412 | 0.282 |
| **separation-component only** | 3.148 | 5.531 | 2.750 | 2.459 | 3.901 | 0.461 | −0.059 | 0.414 | 0.235 |
| shell-profile only (deployable) | 3.163 | 5.608 | 2.755 | 2.425 | 3.893 | 0.480 | −0.056 | 0.459 | 0.231 |
| z(L1)+z(sep-residual) mix | 3.165 | 5.827 | 2.721 | 2.275 | 3.531 | 0.497 | 0.130 | 0.378 | 0.283 |
| z(L1)+z(contact CE) | 3.182 | 5.832 | 2.740 | 2.332 | 3.513 | 0.485 | 0.093 | 0.398 | 0.269 |
| LFO per-shell recalibrated L1 | 3.262 | 5.618 | 2.870 | 2.517 | 3.617 | 0.482 | −0.020 | 0.483 | 0.218 |
| **LFO Mahalanobis under the error covariance** | 3.302 | 5.831 | 2.880 | 2.225 | 3.984 | 0.349 | 0.111 | 0.385 | 0.261 |
| per-shell standardised L1 | 3.382 | 5.660 | 3.002 | 2.541 | 3.648 | 0.444 | −0.033 | 0.494 | 0.196 |
| Gram top-3 eigenvector alignment | 3.441 | 5.357 | 3.122 | 2.362 | 4.522 | 0.110 | 0.016 | 0.478 | 0.184 |
| **RANDOM score** | **3.445** | 5.286 | 3.138 | 2.076 | 4.536 | −0.005 | 0.016 | 0.566 | 0.149 |
| **ARM 4: separation-partialled deviation (deployable)** | **3.497** | 5.968 | 3.085 | 2.410 | 3.950 | 0.346 | **0.170** | 0.381 | 0.243 |
| contact-map cross-entropy alone | 3.572 | 5.985 | 3.170 | 2.334 | 3.893 | 0.309 | 0.120 | 0.384 | 0.259 |

**Arm 4 is a clean negative, and it falsifies the mechanism behind it.** Removing the separation
component from the score costs +0.42 Å and lands the pipeline *below a random score* (3.497 vs
3.445). The premise — "every pool window already satisfies the separation prior, so that component
is dead weight" — is wrong: the K=500 windows are real fragments spanning helix to extended, so
their shell profiles differ enormously and the profile is exactly what discriminates them. The
diagnostic confirms it from the other side: with the TRUE matrix, the shell profile alone reaches
2.299 (90 % of the full oracle's gain) while the within-shell residual alone reaches only 3.030
(barely better than the shipped score, 0.42 Å better than random).

Arm 4 does earn one thing: it has the **best in-band rank correlation of any deployable scorer**
(+0.170 vs +0.126 shipped) while being much worse overall. Consistent with the record's C2/in-band
work, in-band skill and pool-level skill are different quantities; the residual channel has a
little of the former and none of the latter.

Every other alternative representation loses too: contact map (alone or fused), Gram/eigenvector
subspace, scale-invariant scoring, per-shell standardisation, LFO recalibration, LFO per-shell
weights, and the LFO error-covariance Mahalanobis score. **None beats the shipped Bayes risk.**
The Mahalanobis arm is the cheapest test of "structure beats magnitude" as an inference-time fix
and it fails at −0.254 Å: knowing the error covariance does not let you undo the error.

---

## 5. What the objective's information actually is

Per-shell profile q_s = mean CA–CA distance at \|i−j\| = s, across the 126 targets:

| s | true mean | pred mean | bias | sd(true) | sd(pred) | r(pred, true) | r(**pool mean**, true) |
|---|---|---|---|---|---|---|---|
| 2 | 5.946 | 5.994 | +0.048 | 0.413 | 0.385 | 0.574 | 0.566 |
| 3 | 7.153 | 7.347 | +0.194 | 1.487 | 1.278 | 0.602 | 0.608 |
| 4 | 8.557 | 8.838 | +0.281 | 2.031 | 1.674 | 0.588 | 0.556 |
| 5 | 10.210 | 10.534 | +0.325 | 2.243 | 1.721 | 0.455 | 0.445 |
| 6 | 11.300 | 11.776 | +0.476 | 2.844 | 1.970 | 0.479 | 0.435 |
| 7 | 12.032 | 12.468 | +0.436 | 3.671 | 2.270 | 0.406 | 0.473 |
| 8 | 12.976 | 13.253 | +0.277 | 4.473 | 2.618 | 0.367 | 0.514 |

The model over-predicts every shell (peptides predicted too extended), its spread is compressed by
1.4–1.7×, and its correlation with the truth decays with separation. **The K=500 pool's own mean
profile — a native-free quantity involving no model at all — is as good a predictor of the true
profile as the distogram is, and better at s ≥ 7.**

How low-dimensional is the usable signal? (ORACLE arms, emitted avg-RMSD; random = 3.445)

| objective given to the scorer | emitted |
|---|---|
| true full distance matrix | 2.218 |
| true full shell profile (≈12–14 numbers) | 2.299 |
| true profile, shells 2–3 only | 2.871 |
| true \|i−j\|=2 mean only (1 number) | 3.055 |
| true radius of gyration only (1 number) | 2.948 |
| true maximum CA–CA distance only (1 number) | 2.946 |
| random | 3.41–3.45 (two seeds) |

One oracle number (Rg) is worth 0.50 Å; the full profile is worth 1.15 Å; all remaining pair
detail is worth 0.08 Å. The channel is a **low-dimensional shape descriptor, and it is genuinely
multi-dimensional** — no single scalar substitutes for the curve.

Shell ablation (swap one shell of the shipped profile for the truth, or vice versa) is flat:
adding one true shell buys 0.07–0.08 Å for every s; removing one shell from the true profile costs
0.00–0.05 Å. **The requirement is distributed over the whole curve, not carried by a few shells.**

Pair-level sufficiency (replace a fraction of pairs by their true distances, ORACLE selection):

| pairs replaced | 5 % | 10 % | 25 % | 50 % | 75 % | 100 % |
|---|---|---|---|---|---|---|
| chosen by largest \|error\| | 2.814 | 2.608 | **2.294** | 2.212 | 2.213 | 2.218 |
| chosen at random | 3.018 | 2.956 | 2.783 | 2.455 | 2.275 | 2.218 |
| **chosen by longest separation** | 2.878 | 2.720 | **2.443** | 2.260 | 2.216 | 2.218 |
| chosen by shortest separation | 3.075 | 3.070 | 3.031 | 2.861 | 2.668 | 2.218 |

Fixing the **longest-separation quarter of the pairs recovers 74 % of the whole oracle gain**;
fixing the shortest-separation quarter recovers 5 %. The missing information is long-range.

### 5b. Is the CA–CA distance matrix the RIGHT representation? (ORACLE comparison)

Local geometry expressed at the CA level — virtual bond angles and virtual dihedrals of the CA
trace, i.e. the CA-only analogue of a (φ,ψ)/orientation objective — scored against the NATIVE
values (ORACLE), against the oracle distance matrix on the same pools:

| oracle objective | emitted | FAIL18 | other-108 | frac < 2 Å |
|---|---|---|---|---|
| full CA–CA distance matrix | **2.218** | 3.381 | 2.024 | 0.42 |
| distance + local geometry (z-sum) | 2.416 | 3.737 | 2.196 | 0.35 |
| virtual dihedrals | 2.783 | 4.557 | 2.488 | 0.34 |
| virtual dihedrals + virtual angles | 2.806 | 4.486 | 2.526 | 0.32 |
| virtual bond angles | 2.990 | 4.681 | 2.708 | 0.29 |

**The distance matrix dominates.** Perfect local geometry is worth 0.44–0.66 Å less than perfect
distances, and *adding* perfect local geometry to a perfect distance matrix makes the emitted
structure worse (2.218 → 2.416) — the terminal operator averages 75 structures, and a local-
geometry criterion selects sets that agree locally while disagreeing globally. Caveat: this is a
CA-only proxy for the trRosetta ω/θ/φ family; no CB exists in the window universe, so a true
residue-frame orientation objective was not testable here. The CA-level result gives no reason to
expect one to beat distances.

---

## 6. Can anything deployable predict the profile better? (`obj_profile.json`) — NO

Only the profile is swapped; the distogram's own pair-specific residual is kept. All learned arms
are leave-fold-out (fitted on the other four folds' targets; asserted in `fit_all`).

| profile arm | profile MAE | pair MAE | emitted | FAIL18 | other-108 |
|---|---|---|---|---|---|
| ORACLE true profile | 0.000 | 1.533 | **2.402** | 3.715 | 2.183 |
| LFO per-shell debias | 2.458 | 2.342 | 3.071 | 5.776 | 2.621 |
| variance restoration ×1.3 about the LFO shell mean | 2.620 | 2.449 | 3.074 | 5.812 | 2.617 |
| variance restoration ×1.7 | 3.051 | 2.727 | 3.078 | 5.832 | 2.619 |
| **shipped distogram profile (baseline)** | 2.458 | 2.339 | **3.078** | 5.756 | 2.632 |
| LFO affine recal. with the pool profile as 2nd regressor | **2.394** | 2.335 | 3.089 | 5.614 | 2.668 |
| variance restoration ×2.2 | 3.767 | 3.193 | 3.127 | 5.767 | 2.687 |
| LFO per-shell affine recalibration | 2.508 | 2.353 | 3.131 | 5.692 | 2.704 |
| LFO ridge, ESM+composition+disto+pool profile | 2.960 | 2.797 | 3.309 | 5.793 | 2.895 |
| LFO ridge, ESM-2 PCA-32 pooled + composition | 2.923 | 2.759 | 3.339 | 5.733 | 2.940 |
| pool mean profile (native-free typicality) | 2.717 | 2.490 | 3.383 | 5.662 | 3.004 |

Nothing wins. The arm with the **best profile MAE** (cal2, 2.394) emits **worse** than the baseline
— the project's MAE lesson again, now at the level of the profile. Variance restoration cannot
work here for a principled reason: with r ≈ 0.4–0.6 between predicted and true profile, expanding
the prediction amplifies error faster than signal.

---

## 7. Direct verdict on the coordinator's arms 2/3/5 (retraining on the separation residual)

Simulated ceilings, no training required (ORACLE/DIAGNOSTIC; "const profile" = the training-fold
per-shell mean, i.e. exactly what arm 2's offset would supply at inference):

| arm | pair MAE | emitted | FAIL18 | other-108 |
|---|---|---|---|---|
| full oracle matrix | 0.000 | 2.218 | 3.381 | 2.024 |
| true profile + shipped residual | 1.533 | 2.402 | 3.715 | 2.183 |
| **shipped profile + PERFECT residual** (ceiling of arms 2/3/5) | 1.679 | **2.846** | 5.267 | 2.443 |
| shipped (baseline) | 2.339 | 3.078 | 5.756 | 2.632 |
| **const profile + PERFECT residual** (best case of arm 2 as specified) | 2.021 | **3.162** | 4.899 | 2.873 |
| const profile + shipped residual | 2.636 | 3.526 | 5.508 | 3.196 |
| const profile, no residual at all | 2.480 | 3.817 | 5.684 | 3.506 |

* A **perfect** within-shell residual model on top of the shipped profile is worth **−0.232 Å**.
  That is the entire ceiling of arms 2, 3 and 5 in their residual sense.
* Arm 2 **as specified is worse than doing nothing even with a perfect residual head**: subtracting
  a training-fold conditional mean and adding it back at inference replaces the per-target profile
  with a constant, which costs more (+0.084 Å at the ceiling, +0.448 Å with the residual head the
  model would actually learn) than the residual head can ever return. Arm 2 must not be run in
  that form. If it is run, the offset must be a *per-target predicted* profile, at which point the
  experiment is "predict the profile better" — section 6, which fails.
* Arm 3 (separation-balanced loss) is aimed at the right place — long-range shells are where the
  damage is (section 5) — and its ceiling is the profile channel's −0.676 Å, not the residual's
  −0.232 Å. But it only pays if reweighting actually raises r(pred, true) for the long-range
  profile, and section 6 shows every reweighting/recalibration/feature set tried leaves it at
  r ≈ 0.4–0.5, where the pool's own mean profile already sits.
* **Training was NOT run.** Two reasons, both hard: `core.predict.train_fold` writes checkpoints
  into the root `distogram_models*/` directory, which rule 5 forbids me to modify; and
  `free_gb()` read 1.46 GB during this sprint (17 concurrent python processes), below the
  brief's 1.5 GB launch floor. A retrain needs its own MODEL_DIR and a quiet box; I recommend the
  coordinator schedule arm 3 only, alone, and judge it on the per-target long-range profile
  correlation first — if r(s≥6) does not move above ≈0.6, the emitted RMSD will not move.

---

## 7b. Sufficiency: how good must the objective be to reach 2.0 Å?

Per-target, using the profile-error amplitude `a` as the knob (a = 1 is exactly the shipped
model's profile error; a = 0 is a perfect profile). Thresholds are on the coordinate average;
projection adds ≈ +0.16 Å, so "projected < 2.0 Å" ≈ "avg < 1.84 Å".

| | frac of targets reaching it with a PERFECT objective | FAIL18 | other-108 |
|---|---|---|---|
| avg < 2.50 | 0.643 | 0.333 | 0.694 |
| avg < 2.00 | 0.421 | 0.111 | 0.472 |
| **avg < 1.84 (≈ projected 2.0 Å)** | **0.357** | 0.056 | 0.407 |

Among the targets that *can* reach projected 2.0 Å with a perfect objective, **73 % already reach
it with the shipped profile error** (they tolerate a ≥ 1). The objective almost never decides a
single target's 2.0 Å; it decides the MEAN, by inflating the targets that are already bad
(distribution of emitted RMSD, perfect objective: mean 2.218, median 2.303, max 5.50; shipped
profile error: mean 2.846, median 2.643, max 7.41). This is the same shape as the record's
"the mean is set by which targets the prior FAILS on".

**Therefore: no objective upgrade of any kind reaches a 2.0 Å mean.** A perfect CA–CA distance
matrix emits 2.395 Å through the production path (2.218 before projection); 64 % of targets
cannot reach 2.0 Å even then, because their K=500 pool + top-75 average cannot express it.

---

## 8. Null controls and leakage audit

* Every corruption comparison is **paired within target** at matched MAE, so target difficulty
  cannot generate the iso-MAE spread; the spread is between arms on the same targets.
* Across targets, profile MAE correlates r = 0.785 with emitted RMSD and r = 0.708 with the gain
  a perfect objective would give. **After partialling out pool best and pool mean** those become
  0.663 and 0.586 — the relation is not a restatement of target difficulty. (Pair MAE: 0.837
  raw, 0.672 partialled — so MAE *is* a valid within-predictor difficulty index across targets;
  it is only worthless for comparing different objectives at fixed MAE.)
* A **random score** reference (3.445) is reported in every ranking table, so "better than nothing"
  claims are checkable; two arms (arm 4 at 3.497, contact-only at 3.572) fall below it.
* Leakage: deployable arms read only `exp/sd/prob/risk` (the shipped LFO distogram for the target's
  own sequence) and the pool's own coordinates. Learned arms (`fit_lfo`, `obj_profile.fit_all`)
  are fitted strictly on targets with `fold != fold(T)`. Every arm that touches `dtrue`/`nat_ca`
  is labelled ORACLE or DIAGNOSTIC in the tables and in the code. benchmark60 and dev24 untouched.
* Concentration: the oracle-objective gain is not concentrated — 116W/10L, drop-top-10 −0.597 of
  a −0.811 mean. The FAIL18 carry a much larger share (−2.041) than the other 108 (−0.449), but
  the effect survives dropping them.

---

## 9. Confirmation through the REAL projection (`obj_confirm.json`)

Every headline arm re-run through the full production terminal path (L-BFGS projection, λ=0.3
path, `fit_ca` reported). Paired against the shipped Bayes-risk baseline (3.205), 126 targets.

| arm | MAE | r | emitted | FAIL18 | other | <2Å | Δ vs shipped | CI95 | W/L | drop-10 |
|---|---|---|---|---|---|---|---|---|---|---|
| ORACLE full matrix | 0.000 | 1.000 | **2.395** | 3.652 | 2.185 | .38 | −0.811 | [−0.994,−0.640] | 116/10 | −0.597 |
| native + e_res shuffled in shell | 1.524 | 0.836 | 2.425 | 3.697 | 2.213 | .37 | −0.780 | [−0.964,−0.615] | 115/11 | −0.567 |
| **native + i.i.d. noise at the shipped MAE** | 2.339 | 0.751 | **2.443** | 3.611 | 2.248 | .37 | **−0.762** | [−0.957,−0.584] | 100/26 | −0.539 |
| shrink to per-shell mean, slope .376 | 0.694 | 0.952 | 2.471 | 3.754 | 2.257 | .37 | −0.734 | [−0.918,−0.566] | 111/15 | −0.521 |
| ORACLE shell-profile-only scorer | — | — | 2.490 | 3.784 | 2.275 | .33 | −0.715 | [−0.901,−0.545] | 106/20 | −0.500 |
| **native + real error, pairs shuffled** | 2.278 | 0.786 | **2.499** | 3.926 | 2.261 | .38 | **−0.706** | [−0.874,−0.551] | 115/11 | −0.501 |
| ORACLE separation-component scorer | — | — | 2.503 | 3.828 | 2.283 | .32 | −0.702 | [−0.885,−0.533] | 103/23 | −0.487 |
| native + e_res only (pair-specific) | 1.532 | 0.765 | 2.600 | 4.011 | 2.365 | .33 | −0.605 | [−0.786,−0.444] | 108/18 | −0.389 |
| native + e_sep only (profile error) | 1.679 | 0.894 | 3.002 | 5.502 | 2.585 | .29 | −0.203 | [−0.291,−0.119] | 96/30 | −0.100 |
| native + e_res shuffled + e_sep | 2.333 | 0.775 | 3.067 | 5.656 | 2.635 | .29 | −0.138 | [−0.219,−0.064] | 89/37 | −0.047 |
| ORACLE within-shell-residual scorer | — | — | 3.174 | 5.136 | 2.847 | .30 | −0.031 | [−0.225,+0.169] | 78/48 | +0.186 |
| **shipped Bayes risk** | 2.339 | 0.696 | **3.205** | 6.034 | 2.734 | .28 | 0 | — | — | — |
| separation-component-only scorer | — | — | 3.365 | 5.788 | 2.961 | .21 | +0.160 | [+0.043,+0.281] | 48/78 | +0.268 |
| prediction's shell profile, pair detail deleted | 2.178 | 0.713 | 3.367 | 5.859 | 2.952 | .23 | +0.162 | [+0.057,+0.270] | 46/80 | +0.243 |
| shell-profile-only scorer (deployable) | — | — | 3.367 | 5.850 | 2.953 | .21 | +0.162 | [+0.063,+0.265] | 55/71 | +0.243 |
| LFO Mahalanobis | — | — | 3.479 | 6.059 | 3.049 | .26 | +0.274 | [+0.124,+0.431] | 55/71 | +0.374 |
| **RANDOM score** | — | — | **3.601** | 5.511 | 3.282 | .24 | +0.395 | [+0.221,+0.573] | 47/79 | +0.539 |
| **ARM 4 (deployable sep-residual scorer)** | — | — | **3.642** | 6.127 | 3.227 | .27 | **+0.436** | [+0.276,+0.627] | 48/78 | +0.560 |
| contact-map cross-entropy alone | — | — | 3.738 | 6.153 | 3.335 | .24 | +0.533 | [+0.346,+0.730] | 43/83 | +0.678 |

The projected numbers reproduce the coordinate-average sweep with a near-constant +0.16–0.20 Å
offset and identical ordering. Restated on the projected scale, of the 0.811 Å that a perfect
objective is worth:

* **75 % (0.607 Å) is the separation-profile component** of the error, 25 % (0.206 Å) the
  pair-specific component;
* **87 % (0.706 Å) survives shuffling the error across pairs** — i.e. is destroyed by shuffling —
  so it is carried by WHICH pair gets WHICH error;
* **94 % (0.762 Å) is recovered by replacing the real error with i.i.d. Gaussian noise of the very
  same MAE** (100W/26L, CI excluding zero). This is the cleanest statement of the whole report:
  *at matched MAE the shipped distogram is 0.76 Å worse than white noise on the native matrix.*
* Arm 4 (deployable separation-residual scoring) is +0.436 Å [+0.276,+0.627] worse than shipped
  and +0.041 Å worse than a random score. Confirmed negative.

---

## 10. Answers, and what the coordinator should do next

**Is the objective the binding channel?** Partly, and with a hard ceiling. Replacing the shipped
distogram by the exact native distance matrix is worth **−0.811 Å** (3.205 → 2.395 projected,
116W/10L, CI [−0.994,−0.640]) — the largest single deployable-in-principle lever measured in this
sprint. But it stops at 2.395 Å. **The objective cannot deliver 2.0 Å.** Above ≈2.4 Å the binding
constraint moves to the K=500 pool and the top-75 coordinate average (pool best 1.711, oracle
re-weighting of the top-75 2.044).

**Magnitude or structure?** Structure, overwhelmingly. Through the real projection: i.i.d. Gaussian
noise on the native matrix at the *identical* MAE of 2.339 Å emits 2.443 Å — **0.762 Å better than
the shipped distogram** (100W/26L, CI [−0.957,−0.584]). The shipped error marginal reassigned to
random pairs emits 2.499 Å, so **87 % of the damage is the assignment of errors to pairs.** The
decisive structural property is coherence along sequence separation — 75 % of the damage
(0.607 of 0.811 Å) is carried by the error in the per-\|i−j\| mean profile
(a ~13-number curve), and it is concentrated in the long-range shells (fixing the longest-
separation quarter of pairs recovers 74 % of the whole oracle gain; the shortest quarter, 5 %).
MAE remains a decent *within-predictor* difficulty index across targets (r = 0.84 with emitted
RMSD, 0.67 after partialling out pool best/mean) but is worthless for comparing objectives.

**Which alternative representation is most promising?** None of the ones tested. Contact map,
Gram/eigenvector, per-shell standardisation, scale-invariance, LFO recalibration, LFO per-shell
weights, LFO error-covariance Mahalanobis, and the coordinator's separation-residual scoring all
lose to the shipped Bayes risk, two of them by falling below a random score. An ORACLE local-
geometry objective (CA virtual angles/dihedrals) also loses to the distance matrix by 0.57–0.77 Å
and degrades it when fused. **The CA–CA distance matrix is the right representation; the shipped
scoring form (Bayes risk vs plain L1) is worth 0.038 Å, i.e. nothing.**

**What would have to be true for a predictor upgrade to reach 2.0 Å?** Nothing that a predictor
can supply. Concretely, the ladder is:

| upgrade | emitted (projected) |
|---|---|
| shipped | 3.205 |
| perfect within-shell residual, shipped profile | 3.002 (−0.203) |
| perfect per-\|i−j\| profile, shipped residual | 2.600 (−0.605) |
| perfect full distance matrix | 2.395 (−0.811) |
| oracle re-weighting of the top-75 (record) | 2.044 |
| pool best member (record) | 1.711 |

To move the objective at all, a predictor would have to raise r(predicted, true) for the LONG-RANGE
per-separation mean distance from ≈0.37–0.48 to ≈0.8+, per target. Today the K=500 pool's own mean
profile — no model, no labels — matches or beats the distogram on exactly those shells, which says
the distogram has extracted the sequence-conditional compactness prior and nothing beyond it.

### Ranked recommendations

1. **Stop optimising the objective's MAE, and stop optimising the pair-specific residual.** The
   residual channel's total ceiling is −0.203 Å [−0.291,−0.119] and its deployable form is worse
   than random (+0.436 Å vs shipped). Scored *on its own* with ORACLE knowledge it is
   statistically indistinguishable from the shipped score (3.174 vs 3.205, CI [−0.225,+0.169]).
   Do NOT run arm 2 in the form specified (a training-fold shell offset destroys more than a
   perfect residual head returns; simulated at +0.084 Å even in the best case).
2. If any retrain is done, do **arm 3 (separation-balanced loss)** alone on a quiet box with its
   own MODEL_DIR, and gate it on an intermediate metric that costs no structure runs:
   **r(predicted, true) of the per-shell mean distance for s ≥ 6, per target, LFO.** If that does
   not exceed ≈0.6 the emitted RMSD will not move, and it can be checked in minutes.
3. Treat the long-range compactness curve as its own prediction problem with its own features
   (ESM contact-map statistics, predicted SS run lengths, a learned per-target Rg) rather than as
   a by-product of a per-pair distogram — but note that the pooled ESM/composition ridge tried
   here already fails, so a genuinely new information source is required, not a new head.
4. **The larger prize is not the objective.** Even perfect distances give 2.395 Å; 64 % of targets
   cannot reach 2.0 Å from their current pool + operator. Sprint effort is better spent on pool
   coverage (the record's open BREADTH question) and on the terminal operator's 0.68 Å loss
   between the top-75 best member (1.711) and the emitted average (2.395 at oracle).
5. A confirmation pass on dev24 is NOT warranted: there is no positive result to confirm. Every
   deployable arm in this report is a negative. The positive results are all ORACLE ceilings, which
   are properties of the instrument, not of a method.
