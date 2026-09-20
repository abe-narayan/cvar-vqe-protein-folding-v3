# PREREG — S30 lane P: **CAN ANY NATIVE-FREE FEATURE SET PRICE THE ORACLE ERROR, OUT OF FOLD?**

Committed **before the first arm of this experiment produces a number** (contract rule 12).
Lane P, 2026-09-20.

---

## 0. WHAT THIS REPLACES, AND WHAT ALREADY EXISTS

My original brief was L3: *can the distogram's shape error be predicted and corrected from
information demonstrably not already in it?* — i.e. hunt for a decorrelated source, then build a
corrector on it. The coordinator's mid-task message reformulates it, correctly, into **one
estimable number** using lane T's restatement of assumption B2:

> **ρ_max = √(R²(e ~ S))** — e the ORACLE error (displacement from production's answer to the
> native), S any native-free signal.

```
ρ = 0.14   (B2's ceiling, best field ever built)   <->  R²  =  1.96 %
ρ = 0.358  (3.00 Å, charter primary)               <->  R²  = 12.82 %
ρ = 0.628  (2.50 Å, charter ambitious)             <->  R²  = 39.44 %
```

**What is already measured and is NOT re-run here.** S29's 21-field survey and lane D's S30-L5/L6
correction (11 of 21 fields significant, ORACLE *global* combination ρ = 0.1693, leave-fold-out
0.0124, Gram stable rank 2.057). Lane D fitted **one global weight per field**. This experiment
differs in three ways that are the entire reason to run it:

1. it works in a **per-target native-free basis** (the pool's own leading deviation directions),
   which lane T measured to concentrate the ORACLE error ~8× above isotropy — so the bar inside it
   is far lower;
2. it lets the coefficients be **functions of native-free target features**, i.e. it supplies the
   per-target sign/magnitude at inference, which project memory
   (`in-band-ordering-is-per-target`) names as the only remaining leverage;
3. it throws **everything** at `e` at once — distogram outputs, pool statistics, geometry,
   sequence, provenance — rather than enumerating sources one at a time.

**Instrument probe already run (declared, not an arm).** Before this registration I ran a 3-target
probe that reproduced `rmsd_avg` to 1e-4 from `avg_ca` and measured the ORACLE subspace capture
fraction against a random-subspace control. It reproduces lane T's known concentration result and
contains **no arm of this experiment** (no out-of-fold fit, no feature, no prediction). It is
recorded here for chronology, not claimed as a finding.

---

## 1. THE OBJECT

Per target t (all 126 tuning targets, pinned folds, pinned benchmark — contract rule 10):

- `C0` = `shipped_record(pdb)["avg_ca"]`, production's emitted **point cloud** (mean 3.0483 Å).
  **This is an INTERMEDIATE basis, labelled as such in every sentence (contract rule 1).** The
  endpoint is built-chain Cα RMSD (3.2041 Å for `fit_ca`); any headline is converted and labelled.
- `e` = `remove_rigid(oracle_direction(C0, nat), C0)`, flattened to R^{3n}. **ORACLE.**
  `‖e‖²/n = RMSD(C0)²`, verified.
- The **native-free basis** `U_t`, K = 12 orthonormal directions, rigid-body components removed,
  built in this fixed order and QR-orthonormalised in it:
  - `Z`   — the pure **scale** direction `(C0 − centroid)` (expansion/contraction)
  - `U_0` — the **consensus** direction, mean deviation of production's top-75 from C0
  - `U_1..U_10` — PCs of the centred top-75 deviations
- **Shape/scale split** (coordinator point 5, lane F's 83 % result): `c_scale = <e, Z>`;
  `e_shape = e − c_scale·Z`. Reported separately throughout.

## 2. THE MODEL, AND WHY IT IS SIGN-EQUIVARIANT

A PC's sign is arbitrary, so a naive regression of `c_k = <e, U_k>` on target features is
ill-posed. Every predictor here is therefore either **signed** (flips with `U_k`) or **unsigned**,
and the model is linear in the signed ones:

    ĉ_k  =  σ_t · [ Σ_a θ_a s_ka  +  Σ_{a,b} θ_ab · s_ka · z_kb ]

`s` = signed direction features, `z` = unsigned direction/target features, `σ_t` = a **native-free**
magnitude scale (rms pool deviation × √3n). Sign-flipping `U_k` flips `c_k` and every `s_ka`, so
the fit is equivariant by construction. Ridge, **leave-fold-out over the 5 pinned folds**, no
exceptions (contract rule 9 and the S29 memorisation audit: the distogram memorises its training
peptides by 8×, so no in-sample number is admissible).

**Signed features** (all native-free): projections onto `U_k` of — the top-75 mean deviation; the
best-scoring pool member; the pool medoid; the **distogram least-squares descent direction** at C0;
the pure-scale direction; the EXPAND field (C0 rescaled to the pool's mean Rg); the
**peptide-origin minus fragment-origin mean deviation** (lane L's E3 provenance contrast, lifted to
the field level); the member whose Rg best matches the distogram-predicted Rg. Plus the **skewness**
of the 500 pool projections onto `U_k`, and the signed correlation across the pool between
projection on `U_k` and (shipped score / member Rg / BLOSUM sim / member distogram residual).

**Unsigned features**: eigenvalue share, projection sd, direction index, and target-level — n,
pool score mean/sd/gap, mean pairwise RMSD of the top-75, Rg(C0), pool mean Rg, distogram-predicted
Rg, **Rg disagreement** (predicted − realised; project memory's one demonstrated-but-unpriced
native-free signal), mean |d(C0) − expected|, posterior sd, posterior entropy, n_distinct, mean
BLOSUM sim, pool-deviation stable rank and PC1 share.

## 3. THE CONTROLS (contract rules 6, 7, and the coordinator's point 1)

Reported **beside every number**, never after it:

- **C-DIM**: matched-dimension random basis — a random orthonormal K-frame drawn in the rigid-free
  space, features recomputed for it, model refit. At n = 9..16 the rigid-free space is 21–42
  dimensional and the 75 pool deviations span **all of it**, so an isotropic frame is *inside the
  operator's own space* and this is not the `control-must-match-the-operators-space` error.
  Lane D's cautionary case is the reason: an ORACLE per-target fit reached ρ = 0.949 against a
  matched random-subspace control of 0.809.
- **C-PERM**: target-label permutation within length strata — destroys the feature↔target link,
  keeps dimension, basis and response marginals.
- **C-ZERO**: the plausible zero-information field (constant contraction/expansion to the pool's
  mean Rg), *not* an isotropic random field (`zero-information-control-must-be-plausible`).

**The number that counts is EXCESS over C-DIM**, not raw R².

## 4. THE BARS, REGISTERED WITH THEIR CONCLUSIONS

Primary statistic: **out-of-fold R²(e) = 1 − Σ_t‖e_t − ê_t‖² / Σ_t‖e_t‖²**, excess over C-DIM,
fold-clustered bootstrap CI, per-comparison MDE = 2.8016·SE reported beside it.

| band | ρ | verdict I commit to now |
|---|---|---|
| **R² < 1.96 %** | < 0.14 | **B2 confirmed as a measurement, not an enumeration.** Every native-free feature the repository can produce, thrown at `e` at once in the most favourable basis, is at or below the best single field ever built. L3 is closed, E3's last form is closed with it, and the deliverable is the ceiling argument. |
| **1.96–12.82 %** | 0.14–0.358 | Real signal above B2, short of 3.00 Å. Report the implied endpoint `3.0483·√(1−ρ²)` (point cloud, intermediate) and the built-chain conversion, and say plainly it does not reach the charter's primary target. |
| **≥ 12.82 %** | ≥ 0.358 | **Primary target live.** Escalate to the coordinator immediately, before any further analysis, and build the deployable arm. |
| **≥ 39.44 %** | ≥ 0.628 | Ambitious target live. |

**Ablation gate on any positive (coordinator point 3, and the charter's named trap).** A result at
or above 1.96 % must be re-run with **every distogram-derived feature removed** (the four
distogram-descent/residual/posterior/predicted-Rg families and the score-correlation feature). If
the excess collapses, what was measured is the distogram predicting its own error —
`error-coherence-decides-correctors`, recorded twice in project memory — and it is reported as that
and not as a correction channel.

**Stratum arms, registered as secondary:** FAIL18 vs the other 108, plus lane F's two
filter-independent tails (worst-18 by pool mean, worst-18 by ORACLE best-in-pool), because FAIL18
is defined by production and is circular for anything downstream of the filter.

**Multiplicity (contract rule 26).** Arms: {full, no-distogram, scale-only, shape-only} ×
{weighted, unweighted} × {126, FAIL18, 108, tail-A, tail-B} plus 3 controls. ≈ 43 comparisons. A
single 1× MDE positive is expected. **Nothing below 2× MDE is reported as a positive**, above the
contract's own 0.7×/1.0× floor.

## 5. MY OWN REGISTERED PREDICTIONS (so the null is as pre-committed as the positive)

- **P1 (headline):** out-of-fold excess R² over C-DIM is **< 1.96 %** for the full feature set.
  Falsified if the fold-clustered CI excludes 1.96 % from above.
- **P2:** the **scale** coefficient of `e` is more predictable out of fold than the shape part —
  averaging is known to contract the backbone 25.8 % and that is a single global mode. Note this is
  *not* in tension with lane F: lane F refuted scale as the mechanism of **in-pool ranking skill**,
  a different quantity.
- **P3:** C-DIM accounts for **> 70 %** of any raw per-target R².
- **P4:** FAIL18 shows **no higher** out-of-fold R² than the other 108. If it does, the tail has
  exploitable structure and that is the sprint's finding.
- **P5:** removing the distogram-derived features changes the excess by **< 1 percentage point**,
  because the excess will already be ≈ 0.

If P1 survives, the honest deliverable is a well-evidenced ceiling, which the charter values above
a fragile improvement. If P1 fails I say so in the same sentence I report the number.
