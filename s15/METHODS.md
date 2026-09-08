# METHODS

Sprint 15, quantum peptide folding research programme. This document is written to be sufficient
for an independent reimplementation, and to be defensible under hostile review. Every constant that
is load-bearing is named, its sensitivity is stated, and where it was chosen on the same data it is
reported on, that is said plainly rather than omitted.

Status of this document: **living**. Sections marked *(pending)* await a run that is in flight and
will be completed before the protocol freeze. Nothing in this file may be softened after the freeze;
corrections are appended, never applied in place.

---

## 1. PROBLEM STATEMENT

Predict the backbone Cα trace of a peptide of 9–16 residues from its amino-acid sequence alone, and
report full-chain Cα-RMSD to the experimentally determined structure.

The regime is deliberately narrow and it is not the regime most structure-prediction literature
addresses. Peptides of this length are conformationally heterogeneous, frequently lack a single
dominant fold, and are short enough that the statistical potentials calibrated on globular proteins
sit outside their fitted domain — a point Roget et al. (arXiv:2606.21241) make quantitatively for
Miyazawa–Jernigan contact potentials, and which we confirm is *not* inherited by the
distance-restraint objective used here (see RESULTS, K6).

---

## 2. DATA

### 2.1 Targets

**126 tuning targets.** Every experiment in this document is run on these unless explicitly stated.
They are enumerated in pinned, PDB-sorted order by `s12.instrument.targets()`, which reads
`s8/generate_univ/*.npz`. Each carries `pdb`, `n` (chain length), `fold` (0–4), and `seq`.

**60 protected benchmark targets.** Held out completely. **No experiment in this sprint has read
them**, and the protocol is frozen before they are touched. This is not a convention but the central
methodological commitment of the programme: a result that required looking at the benchmark to
obtain is not a prediction.

**Five pinned cross-validation folds.** Every learned component — the distogram, the bias
correction, any tuned hyperparameter — is trained on four folds and applied to the held-out fifth.
The folds were pinned after an earlier correction to identity clustering silently moved 13 targets
and invalidated every fold model then in use; they must not be recomputed.

### 2.2 The retrieval library and window universes

For each target, `s8/generate_univ/<pdb>.npz` caches **every length-n window** of a leakage-safe
library: out-of-fold peptides plus this fold's protein fragments. Per window it stores Cα
coordinates `W (nw, n, 3)`, torsions `PHI`/`PSI (nw, n)`, residue codes `S`, a provenance flag `org`
(True = peptide database, False = protein fragment), the BLOSUM62 similarity `sim` to the target,
and a stable ordering `order`.

Two arrays in these files are **labels, never inputs**: `rr` (each window's Cα-RMSD to the native)
and `nat_ca` (the native trace). They are used for evaluation and for leave-fold-out training labels
only. Any quantity derived from them is marked ORACLE throughout this programme, and the predictive
inference path is structurally incapable of reading them.

### 2.3 Provenance and reproducibility

Verified in the Phase 0 audit (`s15/audit_FINDINGS.md`), rebuilding from primary inputs:

- `peptide_db.npz` rebuilt from 1,524 PDB files: **787/787 records bit-identical**, `ca`/`phi`/`psi`
  maximum difference exactly 0.0.
- All 126 window universes: **bit-identical**, 9 arrays × 126 targets, 2.35 M windows.
- `pool_best = 1.7108244199364904` reproduced **from the PDB files**, not from cache.
- `shipped = 3.4540004952559396` reproduced at torch 1, 2 and 8 threads with **zero argmin flips**.
- The AMBER golden value −489.9138948277905 is **bit-exact**, with ff14SB/GBn2 parameters verified
  against an independently constructed ForceField at maximum difference 0.0.

**A third qualification, and the one that took longest to find.** Every multi-start module in this
sprint originally seeded its RNG with `hash(pdb)`, and **Python's string hash is salted per
process** unless `PYTHONHASHSEED` is set — which the environment lock records as `null`. Two
consecutive interpreters return `hash('1A13') = 217586290314588545` and `2408026022170661001`, so
every run drew a different set of multi-start initialisations. Within a single run all arms shared
one start set, so every paired comparison inside a result file remains start-matched and valid; but
cross-module figures compared different draws, and no multi-start constant was bit-reproducible on
re-running. `s15/seed.py` now supplies a `blake2b`-based stable seed, verified identical across two
interpreters, and the seven affected modules use it. **Exploratory results predating the fix are
correct but not bit-reproducible; the frozen protocol is re-run under stable seeding before the
benchmark is touched, and only those numbers are quoted as reproducible constants.**

This defect was invisible to every determinism check already in place, because all of them ran
within a single process or against a cached artefact. A reproducibility audit that never starts a
second interpreter cannot detect the most common source of irreproducibility in Python.

Two honest qualifications. `python -m s12.instrument` is a **cache read, not a reproduction** — in
its 2.2 s it recomputes only the Kabsch step and asserts only three of its five constants. And the
**distogram is not bit-stable under `torch.set_num_threads`**: 1, 2 and 8 threads give three byte
patterns (`prob` up to 2.0e-6, `risk` up to 1.3e-4). This moves no headline constant, but any claim
of bit-identical reproduction must state the thread count. All results here are produced at 2
threads.

### 2.4 Leakage: measured, disclosed, not removed

16 of 126 window universes contain a window that is verbatim an own-fold database member's full
sequence, and on **4 targets the top BLOSUM hit is the target's own sequence**. Measured impact on
`pool_best`: **exactly zero on all four** — they are different structural determinations of the same
sequence, lying 0.59–4.13 Å from the native.

These targets are **disclosed and retained**. Removing them would improve no number and would
constitute the selective target removal the programme forbids.

A related caution: an identity threshold of 0.6 against the fragment bank is **at the null for
peptides of this length** — random sequences score 0.56–0.63 — so containment-style filters are not
a meaningful leakage control here. Identity or verbatim-substring tests are.

---

## 3. EVALUATION

### 3.1 The RMSD definition, frozen

Independently reimplemented in Phase 0 using **Horn's quaternion method** (no SVD; a proper rotation
by construction). Agreement with the production `kabsch_rmsd_batch`: **2.04e-13 Å over 63,000 real
pool structures**. All 15 edge cases pass, including the mirror test, where both implementations
return 2.9602 Å while the same code without the determinant fix returns 0.0.

> **Frozen definition.** Cα atoms only. Every residue, including both termini. No trimming. Atom
> correspondence by residue index. Uniform weights. **Proper rotations only — reflections
> forbidden.** Ångström. Model 1 of the native. Chain breaks are scored, not rejected. Missing
> atoms are unsupported.

This is not a formality. Full-chain versus `[1:-1]` differs by up to **2.4 Å** on one structure in
this set. Any paper reporting peptide RMSD without stating which it used is not comparable, and this
particular ambiguity was left unstated in the closest recent preprint to our work.

### 3.2 Statistics

Every comparison is **paired across targets** and reported with a bootstrap 95% confidence interval
via `s12.instrument.paired` (4,000 resamples, fold-aware where folds apply). Bare means are never
reported alone. Every table carries **n**, and the **median beside the mean** — a
median-versus-mean divergence is the free diagnostic for a concentrated effect.

Additional discipline adopted after specific failures in earlier sprints:

- **A near-even win/loss record with a CI excluding zero** indicates a concentrated effect and
  triggers a null-calibrated check before the result is quoted.
- **A raw drop-top-k threshold is not a valid test of concentration.** It requires a uniform-effect
  null simulation to interpret, because dropping the largest effects from any distribution shrinks
  the mean.
- **All pre-specified statistics are reported as a single PASS/FAIL verdict**, never as fields a
  reader may select from after seeing them.
- `n ≤ 4` samples from the nine enumerated targets have reversed a conclusion four times in this
  project's history. Small-n reads are labelled as such and are never the basis of a claim.

### 3.3 The FAIL18 column, and what it is not

Tables report a `FAIL18` column for continuity with every prior sprint. It must be read with the
following caveat, established in Phase 0:

`I.FAIL18` is a **threshold artefact**. The band constant `BAND = 1.5 Å` has no derivation anywhere
in the repository, and the zero-recall count runs **45 → 2** across BAND 0.5 → 3.0. On a
99-combination sweep, **exactly one of the 18 targets (9KAR) appears in every version of the set**;
mean Jaccard similarity between versions is 0.498. Three sprints have treated this set as an object.

The column stays, because it is comparable across sprints. **No argument may rest on set
membership**, and the set is never described as "the hard targets" without this caveat attached.

### 3.4 The ORACLE discipline

Any quantity computed from native coordinates inherits **ORACLE = TRUE** and propagates it. ORACLE
arms appear throughout this work as **ceilings and diagnostics** — they answer "what would be
achievable if this component were perfect", which is the only way to attribute a deficit to a
component. They may **never** appear in a predictive headline table, and they are labelled ORACLE in
every row and in prose.

The programme's numerical goals are stated against native-free results only: **< 2.0 Å ideal,
< 2.5 Å a major success, 3.0–3.2 Å not a structural improvement** absent a major mechanistic
finding.

---

## 4. THE INCUMBENT PIPELINE

The comparator every arm is measured against. It emits **3.204 Å** mean full-chain Cα-RMSD on the
126 tuning targets (`synthesis_fit`).

1. **Retrieval.** Score every library window by BLOSUM62 similarity to the target sequence; take the
   top **K = 500**.
2. **Filtering.** Score the pool with a learned leave-fold-out distogram consistency term; keep the
   top **m = 75**.
3. **Aggregation.** Coordinate-average the 75 survivors after superposition.
4. **Projection.** Project the average onto the ideal-geometry manifold (`core.project.lam_path`,
   `ramah` penalty at λ = 0.3, multi-start, exact gradient).
5. **Relaxation.** AMBER ff14SB with GBn2 implicit solvent.

Sensitivity of the two fitted constants, both **fitted on the same 126 targets they are reported
on** — this is a real limitation and is stated rather than buried:

| constant | what it moves | range |
|---|---|---|
| `K = 500` | `pool_best` | 1.970 → 1.504 across K = 100 → 2000 |
| `K = 500` | `shipped` | 3.425 → 3.520 — **flat; not load-bearing** |
| `m = 75` | `top75_best` | 2.609 → 2.106 across m = 25 → 150 |

Undocumented and load-bearing constants inherited by this pipeline, listed so a reviewer does not
have to find them: `min_sep = 2` (no sensitivity sweep exists), `torsion_window = 8`, the AMBER
convergence tolerance of 1.0, and the `bond + angle > 1000` strain gate. There is also a
pre-registration mismatch: `s8/project_devarm.json` pins `phip@0.03` while production ships
`ramah@0.3`. It is worth 0.004 Å on RMSD, but every backbone-plausibility claim rests on the
post-hoc arm.

---

## 5. THE FOUR MANDATED SCIENTIFIC COMPONENTS

These are required variables of the programme, not decorations. Each has a measured role and each
role is defended by a control.

### 5.1 Legacy statistical potential

Eleven components combined by `DEFAULT_WEIGHTS`. The total is verified to be the weighted sum of its
components to within 2.7e-5 over 1.28e7 structures.

**The weights were never fitted.** They are documented as "a variance-balanced starting point, not a
fit." This is a strong claim to have to defend, and it is disclosed rather than presented as a
design choice. Legacy's only measured-honest role in this programme is as a **clash gate**.

### 5.2 AMBER ff14SB with GBn2 implicit solvent

Parameters verified against an independently constructed ForceField at maximum difference 0.0;
golden energy bit-exact.

**Cost, corrected in Phase 0.** With memoisation defeated and hits counted, a single-point evaluation
costs **8.3 ms at 127 atoms rising to 23.3 ms at 237 atoms** — so **8–14 ms** on the n = 9–10
enumerated targets, not the 28 ms previously assumed. A repeated structure costs 0.33 ms, a 33×
memoisation speedup. Sprint 14's budget matching therefore **over-charged AMBER by 2–3×**, meaning
AMBER received *less* search than a fair match and its negative result is, if anything, understated.

**A numerical caution that must accompany any Walsh/Fourier analysis of this energy:** 37.6% of
cached single points exceed 1e6 kcal/mol, 3.2% exceed 1e12, and the maximum is 2.6e20. The 99th
percentile is itself 1e12–1e15, so winsorising at the 99th percentile is **not** sufficient
conditioning; a Pauli spectrum of the unconditioned energy measures its worst steric clash and
nothing else.

**Oracle contamination of the cached AMBER subset, and the binding rule.** Of the labelled subset,
21.7% is oracle-conditioned (`amber_kind == 2`), 40% is unbiased (`kind == 0`), and 39% is
prior-conditioned. `kind == 0` is unbiased **on the mean** (+0.0003 Å) but **not in the tails**,
because the index array force-includes the ORACLE snap index, which lands in the "unbiased" stratum
on 16 of 19 files — where it is the minimum-RMSD member on 6/19 and inside the < 1.5 Å in-band set
on 14/19.

> **Binding rule: every tail, in-band, argmin or top-k statistic uses
> `amber_kind == 0 AND amber_idx != snap_index`, with the per-target n printed beside it.**

### 5.3 VQE

Statevector and exact-MPS simulation (`core/quantum.py`: `StatevectorCircuit`, `MPSAnsatz`, `Adam`).
Configurations are discrete torsion-library states; with a per-residue library of size k and chain
length n the nominal register is `n · log₂ k` qubits.

**The live count is `(n − 2) · log₂ k`.** Four torsions per chain — `phi[0]`, `psi[0]`, `phi[n−1]`,
`psi[n−1]` — are **inert for the Cα trace**, confirmed by two independent routes: direct
perturbation of ±0.7 rad moves the Kabsch RMSD by 0 to 1.4e-07 Å, and the analytic gradient at those
entries is 0, 2.98e-14, 0 and 0 against a maximum of 3.31e+02 and an interior median of 6.76e+01 — a
ratio of 9e-17. Qubit accounting that ignores this overstates the register by four.

### 5.4 CVaR

Verified against the Rockafellar–Uryasev definition: **maximum error 4.2e-14 over 400 cases**, and
`dCVaR/dp` has a median finite-difference error of 0.0.

Three defects in the *estimator as deployed* are carried as findings rather than fixed silently; the
most consequential is that the gradient baseline is centred on the tail only, giving a cosine of
**−0.023** with the true gradient. Any negative result obtained with a defective gradient estimator
is weaker evidence than one obtained with a sound one, so arms with a corrected baseline are run
alongside.

---

## 6. THE RESTRAINT-GENERATION ARCHITECTURE

The methodological departure of this sprint: treat structure prediction as **solving for a
restraint-consistent conformation** rather than selecting among retrieved candidates.

### 6.1 Conformational representation

Continuous backbone torsions `(φ, ψ) ∈ R^{2n}`, with ideal bond geometry. The Cα trace is built by
`core.project.build_ca_exact`, the **bit-exact builder the production projection uses**, so fitted
structures lie on the same manifold the pipeline emits.

### 6.2 The analytic gradient

`core.project.frames` with `_torsion_grad` gives an exact **O(n) reverse-mode** gradient of any
Cα-function with respect to `(φ, ψ)`. Verified against central differences on a real target:
**cosine 1.000000000000000**, with disagreement confined to entries that are numerically zero.

Every objective below is differentiated through this path, so all fits are gradient-based rather
than sampled.

### 6.3 Objectives

For pairs with `|i − j| ≥ 2`:

**Weighted least squares.** `f = Σ w_ij (d_ij(φ,ψ) − d̂_ij)²`, with `w_ij = 1/sd_ij²` from the
distogram's own per-pair uncertainty. The `sd` had never been used for anything in fourteen sprints
and is worth real accuracy.

**Maximum likelihood.** `f = −Σ log P_ij(d_ij(φ,ψ))` under the distogram's full 17-bin per-pair
distribution — strictly more informative than least squares, of which it is the Gaussian special
case.

> **Correction, 2026-09-05.** The likelihood table's bin lookup was initially written as
> `floor((d − c₀)/(c₁ − c₀))`, a **uniform-grid** interpolation. The distogram's bin centres are
> **not uniform**: they run 0.5 Å apart at short range and up to 4.0 Å apart at long range. Every
> query above ~5 Å therefore read the wrong bin (a query at 8.50 Å read the 7.25 Å bin; 15.00 Å read
> the 17.50 Å bin; everything above 15.25 Å collapsed into one). The lookup is now a `searchsorted`
> on the real centres with per-interval widths, and the derivative is re-verified against central
> differences at **cosine 1.000000000000000, maximum absolute error 5e-09**. All maximum-likelihood
> results predating this correction are void and are recorded as such.

**Torsion prior regularisation.** A von Mises log-density pulling toward retrieval-conditioned
circular-mean torsions — a genuine probabilistic prior with an exact gradient, not an ad-hoc
penalty. This combination of a generative distance restraint with the project's best-measured
torsion channel is the classical limit of the quantum Hamiltonian in §7.

### 6.4 Optimisation and the selection rule

L-BFGS-B with the analytic Jacobian (`maxiter` 400, `maxcor` 20, `ftol` 1e-12, `gtol` 1e-10), from a
bank of **native-free starts**: the retrieval circular-mean torsions, an ideal α-helix, an extended
strand, and draws from the retrieved pool.

> **Multi-start selection uses the OBJECTIVE ONLY. RMSD is read post hoc, never to choose.**

This is stated as a rule because violating it is the single easiest way to manufacture a result in
this problem, and because a related failure — an `argmin` over a tied signal silently reading the
ORACLE sort order — invented a 1.386 Å winner in an earlier sprint. Where a signal can tie, the
outcome is averaged over the tied argmin set.

### 6.5 Restraint channels

Three native-free channels, plus an ORACLE ceiling:

- **the learned distogram** — leave-fold-out, 17 non-uniform bins, per-pair mean and sd;
- **the retrieval pool** — every one of the K = 500 retrieved windows supplies a complete Cα–Cα
  distance matrix, giving a per-pair *empirical* distribution with no training at all. This channel
  is new to this programme;
- **their product**, treating the two as independent likelihoods;
- **ORACLE**, the true distance matrix, used only as a ceiling.

### 6.6 The four-stage readout

A single final RMSD hides which stage fails, so every architecture reports four numbers and the
three gaps between them:

| | |
|---|---|
| **G** | the best structure present anywhere in the generated ensemble — **ORACLE readout** |
| **S** | what the method identifies without the native — the objective's argmin |
| **A** | what an estimator recovers by combining the ensemble — coordinate consensus |
| **F** | what is finally reported, after physics — validity gate and refinement |

`S − A` is the number that decides whether selection matters at all. This decomposition is motivated
by two measurements: coordinate aggregation is worth **3.4× the entire contribution of any
objective**, and a *perfect* ranker emits 2.474 Å through a decile but 1.219 Å through a top-100 — a
1.26 Å swing from the terminal operator alone.

### 6.7 Physical validity

A steric gate on the fraction of `|i − j| ≥ 2` Cα pairs closer than 3.6 Å. **Every excluded
structure is counted and reported.** A gate that silently improves a mean by discarding structures
is a fabrication; the exclusion count is part of the result. A gate is never permitted to empty a
target.

---

## 7. THE QUANTUM EXPERIMENT

*(pending — the run is in flight; design is fixed and pre-declared below)*

The question is whether Sprint 14's decisive negative survives a change of **objective class**.
Sprint 14 found CVaR-VQE was worse than best-of-N sampled from the same circuit **untrained**
(0 wins in 12, +0.65 to +1.32 Å), and traced the mechanism to the objective having no ordering skill
inside the tail the optimiser concentrates on. Every objective it tested was **selective** — a
physical energy scoring nativeness. The restraint objective is **generative**, and its ordering skill
is independently measured to be real.

One variable changes. The apparatus is Sprint 14's, unmodified: the same nine fully enumerated
targets at k = 4 (262,144 configurations, 18 nominal qubits), the same index algebra, the same MPS
ansatz, the same hard evaluation `Counter`, the same classical searches, the same two-axis
reporting. Only the energy vector changes.

**Arms, all at an identical hard budget of objective evaluations:** best-of-N from the **untrained**
circuit (the control), CVaR-VQE at a sweep of α, uniform random, simulated annealing, and greedy
1-opt with restarts.

**Declared before running.** A positive result requires VQE to beat `untrained_bestofN` on returned
RMSD, paired across targets, with a CI excluding zero, at matched budget. **Beating uniform random
proves nothing here** — a zero-information constant α-helix beats the random control on this
instrument. A negative result — VQE losing again on a generative objective with measured ordering
skill — is a *stronger* and more publishable claim than the original, because it removes
scoring-function pathology as the explanation.

**Both axes are reported separately and never substituted:** `objective_gap` (optimisation) and
`structural_gap` (accuracy). The entire Sprint 14 story is that these come apart.

---

## 8. COMPUTE

A single machine: 15.6 GB RAM, 8 cores that are 4 fast + 4 slow, measured at **6.43
core-equivalents** on the real workload. Utilisation is targeted at ~95% and never intentionally
exceeded past 97%. Two heavy jobs fit; three do not. BLAS thread counts are pinned to 2 per process.
Long runs checkpoint to JSON every 10 targets.

A caution for anyone reading the result artefacts: **325 of 363 result JSONs carry no completeness
flag**, and of the 38 that do, **all 12 "incomplete" flags are false alarms** — the writer takes
`len()` of a dictionary of arms. A flag with a 12/12 false-alarm rate is worse than none and must
not be cited as evidence that a run finished.

---

## 9. WHAT THIS METHODOLOGY CANNOT DO

Stated here rather than left for a reviewer to find.

- **The two fitted retrieval constants were fitted on the reported targets.** `K` and `m` were not
  cross-validated. `K` is not load-bearing for the headline; `m` is.
- **`min_sep = 2` has never been swept.** All pair statistics inherit it.
- **Legacy's weights were never fitted**, and are presented as a variance-balanced starting point.
- **The FAIL18 set is a threshold artefact** and supports no argument from membership.
- **The distogram is not bit-reproducible across BLAS thread counts**, only across runs at a fixed
  count.
- **The enumerated quantum targets are nine chains at n = 9, k = 4.** Conclusions drawn there are
  conclusions about a 262,144-configuration discrete space, not about the continuous problem, and
  the two are connected only by explicit ceiling measurements.
- **No result in this document has been checked against the protected benchmark**, by design. The
  benchmark is a single-use instrument and is spent only after the protocol freeze.
