# SPRINT 14 -- OBJ. Can a learned structural objective produce an accurate low-energy region?

Agent: OBJ. Instrument checked at start:
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18` -- **all five pinned constants
unchanged**. Re-run at the END of this work: **identical, all five unchanged.**

The question this file answers is NOT "can a model predict RMSD". It is: **can a learned
objective, seeing only sequence + torsion configuration + generic priors, produce a
landscape whose low-energy region is structurally accurate on targets it has never seen.**

---

## SUMMARY FOR THE COORDINATOR

**The question:** can a learned many-body structural objective, seeing only legitimately
available information, produce a landscape whose low-energy region is structurally
accurate?

**The answer: NO -- and C19 now shows exactly WHY, which is the durable result.**

0. **IT IS TRANSFER, NOT CAPACITY (C19, my `obj_ceiling.py`, run by the coordinator).**
   On a target's own near-native band the same 4,125-parameter linear model reaches
   **0.986** in-band ordering accuracy -- against the **0.638** the 2.0 A target needs --
   and scores **identically on held-out configurations of that target** (overfitting gap
   **-0.0005, CI [-0.0011, +0.0000]**; `d_top100 = -0.868 A`). Cross-target it collapses
   to **0.600** (transfer gap **-0.386, CI [-0.437, -0.332]**), worth -0.047 A.
   **A richer model is not worth building** -- the linear model already saturates the
   within-target problem, and no architecture addresses transfer, where all the loss is.
1. **A learned objective DOES improve discrimination -- in the wrong place.** Global
   pairwise ordering accuracy below 1.5 A rises from Legacy's 0.534 to **0.616**
   (sequence-blind variant). In-band accuracy -- inside the low-energy decile, where a
   search lives -- goes from Legacy's 0.539 to **0.514**, i.e. slightly DOWN. Learning
   buys garbage rejection, which Legacy already had and which is spent by the time a
   search arrives.
2. **The requirement is quantified for the first time.** A noisy-oracle sweep over the
   enumerated spaces gives the transfer curve from ranking quality to structure. The
   2.0 A target needs **in-band ordering accuracy ~0.638** (`sigma ~ 0.75 A`). Everything
   measured -- Legacy 0.539, AMBER 0.463, the 1-local prior 0.501, the learned objective
   0.523 -- sits near the `sigma = 3 A` rung. This is an order-of-magnitude gap in
   ordering information, not a tuning gap.
3. **The sequence-conditioned model is refuted; a sequence-BLIND geometric shape prior is
   a real but small effect.** learned -0.166 A [-0.467, +0.119] (CI crosses zero, reverses
   on drop-top-3). seqblind **-0.337 A [-0.546, -0.126]**, 16/3, and beats Legacy by
   -0.361 A [-0.647, -0.076]. Sequence conditioning is neutral-to-harmful (+0.172 A,
   W/L 9/10), independently reproducing two pinned project findings.
4. **The positive control works, which is what makes the negative trustworthy.** The
   identical harness with a leaked label returns -1.405 A, CI [-1.597, -1.235], **19/0**,
   -1.093 after drop-top-10. Both label-permuted and random-weight nulls are significantly
   WORSE than random, so "beats the null" here is the weak-control trap; the learned
   objective does not clear random itself.
5. **Two data-integrity corrections**, both reproduced independently: the cached AMBER
   subsample is oracle-conditioned by -0.405 A (AMBER's apparent skill was 4x overstated;
   my own `amb_nonbonded` in-band claim is retracted), and **FOUR torsions are inert for
   CA-RMSD, not three** (bit-exact, 7 targets), so 14 of 18 qubits are live at n=9, k=4.

6. **The learning curve is COMPLETE on both arms, and the two DIVERGE.** The learned
   arm declines (d_decile -0.348 / **-0.380** / -0.231 / -0.230 / -0.162 A at
   m = 1/2/4/8/12) with its in-band `rho_decile` **pinned at zero throughout**
   (-0.050 / -0.059 / -0.003 / +0.033 / +0.010). The leaked-label control in the same
   harness is loud at a single target and **rises** (-1.273 / -1.397 / -1.407 / -1.419 /
   **-1.423 A**), with in-band `rho_decile` rising monotonically
   **+0.551 -> +0.728**. So the harness can express a rising curve; this model has no
   slope. **No caveat remains** -- this is Sprint 12's decisive negative reproduced on a
   different instrument, and it is the predicted consequence of C19.


**For the VQE agent: do not adopt this objective as a Hamiltonian.** If a geometric
regulariser is wanted, use the sequence-blind variant and price it at ~0.34 A, not more.
`LearnedObjective.pauli_support(9)` = {sep 2: 2 qubits, ... sep 8: 14}; the family is
already full-register, and `Featurizer(max_sep=s)` is the one real handle -- it caps Pauli
weight at exactly `2(s-1)`.

**Still open:** only `s14/obj_curve.py ood` (generator shift -- does the objective survive
a VQE's own proposal distribution). It does not change the verdict. Both experiments I had
listed as outstanding have now completed: `obj_ceiling` is item 0 above, and the learning
curve is item 6.


---

## 0. ASSETS BUILT

### 0.1 The enumerated set was extended from 9 targets to 19 -- DEMONSTRATED (infrastructure)

`s14/obj_enum.py` enumerates any tuning126 target at `k = 4` in the same schema as
`s13/results/qarch_enum_<PDB>.npz`, caching to `s14/cache/obj_enum_<PDB>.npz`.

| n  | configs / target | targets enumerated | cost / target |
|----|------------------|-------------------:|--------------:|
| 9  | 262,144          | 9  (s13, reused)   | ~1.5 min |
| 10 | 1,048,576        | **10 (new)**       | ~3.5 min (uncontended) |
| 11 | 4,194,304        | 0 -- **not run**   | ~10 min (uncontended) |

**Delivered: 19 targets, 1.28e7 exactly-labelled structures**, per-fold counts
`{0:4, 1:2, 2:3, 3:3, 4:7}` -- every pinned fold is represented, so leave-fold-out is
possible, but **fold 1 has only two held-out targets and its per-fold numbers are
correspondingly weak**. This is the honest state, not the plan.

**Why the n=11 batch was abandoned, stated plainly.** It was launched (13 targets) and
stopped after the n=10 batch. The box was running 5-9 concurrent Python jobs from other
agents all sprint; measured, my enumeration process held the single largest CPU share of
any of my jobs and the n=11 batch would have cost 4x per target at ~3x the contended wall
rate -- several hours during which nothing else of mine could run. Given that the result
below is a clear negative whose mechanism is already identified and *not* sample-limited
(section 4.4), 13 more targets would not have changed it. If a later sprint wants them:
`python -m s14.obj_enum 11` reproduces the batch exactly.

Cost, measured on this box: Legacy is 84 us/config, CA-RMSD 9 us/config, so a target
costs `~93 us * 4^n` plus ~50 s of AMBER subsample. **n=12 (16.8 M configs) would be
~27 min/target uncontended and ~940 MB of resident columns; not attempted.**

Correctness: the vectorised odometer `enum_configs` was checked to be **bit-identical** to
`itertools.product` ordering, and recomputed `rmsd`/`legacy` on sampled rows of the s13
cache match to `2e-3` / `1e-4` relative. The two caches are interchangeable.

---

## 1. STEP 2 -- THE DISCRIMINATION FLOOR, MEASURED EXACTLY

`python -m s14.obj_floor` -> `s14/results/obj_floor.json`.
Exact over the whole 262,144-configuration space of each n=9 target (no sampling).
**ORACLE DIAGNOSTIC** -- RMSD is the evaluation label throughout.

### 1.1 Nothing ranks in-band, and the components do not either -- DEMONSTRATED

Pooled over the nine n=9 targets. `rho_dec` is Spearman inside the objective's own
low-energy decile, which is where a search lives. `dec_mean` is mean CA-RMSD of that
decile; the random baseline (`all_mean`) is 3.781 A.

| objective            | rho_global | rho_decile | snap pct | top1 A | top100 A | decile mean A |
|----------------------|-----------:|-----------:|---------:|-------:|---------:|--------------:|
| legacy               | +0.081 | +0.085 | 0.629 | 3.920 | 3.765 | 3.773 |
| prior (1-local)      | +0.011 | -0.025 | 0.180 | 3.636 | 3.934 | 3.823 |
| leg_steric           | +0.156 |  ties  | 0.282 | 3.552 | 3.485 | 3.695 |
| leg_electrostatic    | +0.126 | -0.140 | 0.307 | 3.572 | 3.687 | 3.739 |
| leg_contact          | +0.040 | +0.030 | 0.451 | 4.106 | 4.010 | 3.792 |
| leg_torsion          | +0.009 | +0.014 | 0.179 | 3.302 | 3.401 | 3.741 |
| leg_coop_helix       | +0.000 | -0.012 | 0.042 | 3.629 | 3.672 | 3.775 |
| leg_aromatic         | -0.018 | -0.114 | 0.240 | 3.676 | 3.741 | 3.840 |
| leg_coop_sheet       | -0.025 | -0.081 | 0.000 | 4.701 | 4.653 | 3.781 |
| leg_hbond_local      | -0.031 | +0.091 | 0.257 | 3.827 | 3.704 | 3.844 |
| leg_hbond_longrange  | -0.071 | -0.214 | 0.015 | 5.360 | 5.204 | 3.859 |
| leg_solvation        | -0.128 | -0.073 | 0.632 | 4.455 | 4.271 | 4.097 |
| leg_compactness      | -0.177 |  ties  | 0.658 | 3.641 | 3.831 | 4.150 |

**No single Legacy component has in-band skill.** The largest |rho_decile| is
`leg_hbond_longrange` at **-0.214** -- i.e. the best in-band signal available in the
Legacy decomposition points the *wrong way*. `leg_steric` and `leg_compactness` are
**entirely tied inside their own low-energy decile** (rho undefined): once you are in the
no-clash set, the steric term carries literally zero bits.

On the 2,955-config AMBER subsample (uniform/prior/near-native band).
**>>> THE TABLE IMMEDIATELY BELOW IS SUPERSEDED -- see section 1.3. <<<** That subsample
is ORACLE-CONDITIONED by -0.405 A (ENER agent's finding, reproduced here), so every AMBER
number in it is measured on a selection that partly encodes the native. It is kept in
place because the size of the correction is itself the finding.

| objective    | rho_global | rho_decile | snap pct | top1 A | top100 A |
|--------------|-----------:|-----------:|---------:|-------:|---------:|
| amber_total  | +0.116 | -0.014 | 0.433 | 3.613 | 3.136 |
| amb_nonbonded| +0.122 | +0.087 | 0.439 | 3.039 | 2.948 |
| amb_solvation| +0.060 | -0.048 | 0.415 | 3.759 | 3.588 |
| amb_angle    | +0.030 | -0.019 | 0.235 | 3.677 | 3.379 |
| amb_torsion  | -0.086 | -0.055 | 0.623 | 3.453 | 3.579 |
| legacy@sub   | -0.061 | -0.089 | 0.637 | 3.887 | 3.751 |
| prior@sub    | -0.179 | -0.069 | 0.360 | 3.714 | 3.968 |

`amb_nonbonded` is the only column in either decomposition with a positive in-band rho
(+0.087), consistent with the memory `amber-converged-interaction-only-ranks.md`. It is
still tiny. **RETRACTED by section 1.3: on the unbiased subsample that +0.087 becomes
-0.004. The one positive in-band AMBER signal I reported was an artefact of the
oracle-conditioned sampling.**

### 1.3 CORRECTION -- the cached AMBER subsample is oracle-conditioned -- DEMONSTRATED

Flagged by the coordinator from the ENER agent's work, independently reproduced here on
all 16 targets enumerated at the time.

The AMBER subsample in `qarch_enum_*.npz` (and, by construction, in my own
`obj_enum_*.npz`, which reuses the s13 sampling scheme) is built from three populations
tagged in `amber_kind`: 0 = uniform, 1 = prior-sampled, 2 = drawn from the lowest-1%-RMSD
band. Measured selection bias of the subsample relative to its own full space:

| subsample | mean RMSD offset vs full space | rows/target |
|-----------|-------------------------------:|------------:|
| all 2,955 rows       | **-0.405 A** | 2955 |
| `amber_kind == 0` only | **+0.007 A** | 978 |

So `amber_kind == 0` is an unbiased sample and the full subsample is not. Every AMBER
number recomputed on `amber_kind == 0` (16 targets, space mean 3.844 A):

| objective | rho_g (all) | rho_g (k0) | rho_dec (all) | rho_dec (k0) | top100 (all) | top100 (k0) | in-band <1.5 (k0) |
|-----------|------------:|-----------:|--------------:|-------------:|-------------:|------------:|------------------:|
| amber_total   | +0.134 | +0.114 | -0.042 | **-0.045** | 3.223 | **3.694** | **0.463** |
| amb_nonbonded | +0.136 | +0.117 | +0.024 | **-0.004** | 3.158 | **3.683** | 0.500 |
| amb_solvation | +0.170 | +0.196 | -0.056 | -0.004 | 3.500 | 3.722 | 0.485 |
| amb_angle     | +0.011 | -0.003 | -0.019 | +0.025 | 3.473 | 3.796 | 0.501 |
| amb_torsion   | +0.006 | -0.016 | -0.074 | -0.041 | 3.616 | 3.925 | 0.498 |
| legacy@sub    | +0.038 | +0.121 | +0.025 | +0.097 | 3.562 | 3.834 | 0.527 |
| prior@sub     | -0.151 | -0.007 | -0.021 | -0.025 | 3.920 | 3.918 | 0.498 |

**What changed, stated plainly.** AMBER's top-100 mean RMSD goes from 3.223 A (an apparent
-0.62 A over random) to **3.694 A (-0.15 A)** -- the apparent skill was **four times
overstated** by the conditioning. And AMBER's in-band ordering accuracy below 1.5 A is
**0.463, i.e. significantly BELOW chance**: inside its own low-energy decile the all-atom
force field mildly anti-orders near-native structures. The direction of section 1.2's
conclusion is unchanged; its magnitude for AMBER was too generous.

**LEAKAGE AUDIT OF MY OWN OBJECTIVE -- confirmed clean.** The learned objective's data
path (`obj_model.py`, `obj_train.py`, `obj_headline.py`, `obj_curve.py`, `obj_resid.py`,
`obj_noisy.py`) reads only `rmsd`, `legacy`, `prior`, `prior_table`, `PHI`, `PSI`, `n`,
`k`, `seq`, `fold` and `snap_index`. It never reads `amber_total`, `amb_*` or
`amber_idx`; verified by grep across all six modules. All of those columns are
full-enumeration and unconditioned. **No AMBER quantity has entered any feature, label,
normalisation constant or model weight anywhere in this workstream.**

### 1.4 CORRECTION -- FOUR torsions are inert, not three -- DEMONSTRATED

The brief states `phi[0]`, `psi[n-1]`, `phi[n-1]` are inert for the CA trace. Tested
exactly on the enumerations, exploiting the odometer indexing (`index = sum_i s_i k^(n-1-i)`,
so residue 0 is the most significant digit and residue n-1 the least): reshape `rmsd` to
`(k,)*n` and ask whether it is constant along an axis.

| target | n | residue 0 inert | residue n-1 inert | residue 1 (control) | residue n-2 (control) |
|--------|--:|:---------------:|:-----------------:|:-------------------:|:---------------------:|
| 1CS9 | 9 | **True** | **True** | False | False |
| 2MK7 | 9 | **True** | **True** | False | False |
| 6EY3 | 9 | **True** | **True** | False | False |
| 8IS3 | 9 | **True** | **True** | False | False |
| 9UV5 | 9 | **True** | **True** | False | False |
| 1N9U | 10 | **True** | **True** | False | False |
| 2BAO | 10 | **True** | **True** | False | False |

Exact equality (`rtol=0, atol=0`), not approximate: CA-RMSD is **bit-identically invariant**
to the torsion states of BOTH terminal residues, and demonstrably not to their neighbours.
Confirms the coordinator's correction. The dead block is `2*log2(k)` qubits per chain, so
**14 of 18 live qubits at n=9, k=4**, and this workstream reports live counts.

Consequence taken up in the model: the LOCAL (Ramachandran) features of residues 0 and
n-1, and their contribution to the per-residue prior scalar and the helix/sheet fraction
scalars, are pure label noise and are now masked (`Featurizer(mask_terminal=True)`, the
default). PAIR features are deliberately NOT masked -- `d(0, j)` depends on residues
`1..j-1` and is a genuine observable even though residue 0's own state is not. The
`nomask` arm in section 4 prices the masking.

### 1.2 THE MECHANISM: discrimination is a monotone function of structural separation -- DEMONSTRATED

This is the most useful thing in this section. For 400,000 random configuration pairs per
target, `P(objective orders the pair the same way true RMSD does)`, binned by the pair's
RMSD separation. 0.5 is chance.

**Global (pairs drawn from the whole space):**

| objective        | 0.25-0.5 A | 0.5-1 | 1-1.5 | 1.5-2 | 2-3 | 3-4 | 4-5 | >5 |
|------------------|-----------:|------:|------:|------:|----:|----:|----:|---:|
| legacy           | 0.508 | 0.517 | 0.537 | 0.573 | 0.624 | 0.719 | 0.760 | 0.770 |
| leg_steric       | 0.512 | 0.526 | 0.552 | 0.586 | 0.628 | 0.669 | **0.938** | **0.961** |
| amber_total      | 0.510 | 0.520 | 0.528 | 0.539 | 0.574 | 0.595 | 0.681 | **0.948** |
| prior            | 0.499 | 0.499 | 0.505 | 0.516 | 0.528 | 0.554 | 0.551 | 0.556 |
| leg_contact      | 0.504 | 0.510 | 0.526 | 0.540 | 0.547 | 0.524 | 0.695 | 0.526 |
| leg_compactness  | 0.459 | 0.431 | 0.409 | 0.401 | 0.404 | 0.420 | 0.000 | 0.000 |

**In-band (pairs drawn only from the objective's own low-energy decile):**

| objective        | 0.25-0.5 A | 0.5-1 | 1-1.5 | 1.5-2 | 2-3 | 3-4 |
|------------------|-----------:|------:|------:|------:|----:|----:|
| legacy           | 0.516 | 0.529 | 0.543 | 0.552 | 0.562 | 0.649 |
| leg_steric       | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| amber_total      | 0.501 | 0.492 | 0.486 | 0.479 | 0.495 | 0.489 |
| prior            | 0.497 | 0.495 | 0.493 | 0.490 | 0.486 | 0.490 |
| leg_contact      | 0.498 | 0.501 | 0.505 | 0.510 | 0.508 | 0.604 |

Three things fall out of this table, and they are the mechanistic core of the sprint's
"neither energy ranks the native" result:

1. **Both physics objectives are garbage detectors, not structure rankers.** Legacy
   reaches 0.770 at >5 A separation and 0.508 at 0.25-0.5 A. AMBER reaches 0.948 at >5 A
   and 0.510 at 0.25-0.5 A. Discrimination is monotone in structural separation and
   **decays to chance below ~1.5 A**, which is exactly the band the 1.594 A ORACLE answer
   lives in. An objective that cannot order two structures 1 A apart cannot deliver a
   1.6 A answer no matter how well it is optimised.

2. **Legacy's skill is clash rejection, and it is spent by the time you are in-band.**
   `leg_steric` is the single most discriminating column at long separation (0.938/0.961)
   and is **exactly 0.500 in every in-band bin** -- all ties. This is a direct, exact
   confirmation of the sprint brief's claim, now with the mechanism attached: the term
   that carries the skill has *no remaining variance* in the region a search occupies.

   **CORRECTION (section 1.3):** the `amber_total` row of both tables above is measured on
   the oracle-conditioned subsample. On the unbiased `amber_kind == 0` rows AMBER's
   in-band accuracy below 1.5 A is **0.463**, i.e. further below chance than shown here.
   The Legacy and `leg_*` rows are full-enumeration and stand as measured.

3. **AMBER goes below chance in-band** (0.479-0.495 at 1-3 A; **0.463 on the unbiased
   subsample**), i.e. inside its own
   low-energy decile the all-atom force field mildly *anti*-orders structures. This is the
   fine-grained version of `physics-ranks-real-geometry-not-lattice.md`.

`leg_compactness` is systematically **anti**-discriminating globally (0.40-0.46): more
compact is worse. Any learned objective is therefore free to discover a compactness term
of either sign, and its sign is a target-level property -- which turns out to be the whole
story in section 2.

**The bar this sets for a learned objective:** it must beat 0.50 *in-band* in the
0.25-1.5 A bins. Global rho is not the test; every physics term already has some.

---

## 2. STEP 3-4 -- THE LEARNED OBJECTIVE

### 2.1 Model family

`s14/obj_model.py`. A sequence-conditioned, distance-binned, many-body pair potential:

    E(S) = (1/n) [ sum_{j-i>=2} w_pair[sep(i,j), cls(i), cls(j), dbin(d_ij)]
                 + sum_i        w_loc[cls(i), rama_cell(phi_i, psi_i)]
                 + w_glob . g(S) ]

7 residue classes (28 unordered pairs), 8 separation classes, 17 CA-CA distance bins,
36 Ramachandran cells, 8 global scalars (Rg, end-to-end, contact counts, per-residue
prior energy, helix/sheet fractions), intercept. **DIM = 4069**, linear in indicator
counts, so ridge on streamed normal equations fits it exactly.

Information contract, one line: **sequence + integer torsion configuration + the
target-held-out k=4 library + the leakage-safe 1-local empirical prior + weights trained
on other targets.** No native anything. RMSD is a training label and an evaluation axis.

Why this family first: it is SE(3)-invariant by construction, it is the standard learned
statistical-potential form, and its Pauli structure is exactly readable off the sprint-13
locality theorem (a separation-s pair term touches exactly `2(s-1)` qubits at k=4), so a
positive result would be directly convertible to a VQE Hamiltonian. Complexity is earned
from here, not assumed.


### 2.2 First fit, and the methodological correction it forced

Ridge on the raw label ("predict the absolute within-target RMSD percentile") transferred
*negatively*: leave-one-target-out on nine targets gave mean rho_global +0.128 with a
per-target range of **-0.59 to +0.67**, and mean d_decile (decile RMSD minus the target's
own space mean) of essentially zero.

**Correction applied, and it is the right one:** the fit must be a WITHIN-target ranking
regression, not an absolute one -- otherwise the model spends its capacity explaining
between-target offsets it cannot see. `RidgeAccumulator.add(center=True)` subtracts each
target's own column means and label mean before accumulating, exactly and without
densifying (`Xc'Xc = X'X - N mu mu'`, `Xc'yc = X'y - N mu ybar`).

A second correction followed: **standardise the columns before ridge**. The count columns
and the global scalars differ in variance by orders of magnitude, so one lambda penalises
them incomparably and the highest-variance columns dominate *any* weight vector. The
symptom was diagnostic: on 1CS9 a **label-permuted twin scored d_decile +0.82 against the
real model's +0.84** -- the null was as "good" as the model, because both were really just
reading the highest-variance global column. Un-standardised ridge on heterogeneous feature
families is not a fair harness and every number below uses `standardize=True`.

### 2.3 THE FLIP DIAGNOSTIC: the learned objective is a compactness axis -- DEMONSTRATED

With centering, leave-one-target-out over 12 enumerated targets at lambda 0.03 gave mean
rho_global **+0.184** but mean rho_decile **-0.055** and mean d_decile **-0.088 A**
(W/L 8/4). The per-target rho_global still ranged from **-0.643 to +0.675**. That spread is
not noise, and it is not random: it is a single latent variable.

For each target, over its evaluation set, versus the leave-one-out learned rho_global:

| statistic (per target)                    | Spearman with learned rho_global |
|-------------------------------------------|---------------------------------:|
| rho(contact count < 8 A, RMSD)            | **+0.951** |
| rho(radius of gyration, RMSD)             | **-0.930** |
| rho(end-to-end distance, RMSD)            | **-0.881** |
| native Rg, z-scored against the space's Rg | **+0.909** |

Per target (`natRg` = native radius of gyration, `spaceRg` = mean over the enumerated
space, `z` = (natRg - spaceRg)/sd):

| pdb | learned rho_g | rho(Rg,RMSD) | natRg | spaceRg | z |
|-----|--------------:|-------------:|------:|--------:|---:|
| 2MK7 | +0.675 | -0.904 | 8.37 | 5.43 | **+3.92** |
| 2P5H | +0.668 | -0.840 | 6.83 | 5.36 | +1.96 |
| 1N9U | +0.647 | -0.741 | 7.16 | 5.90 | +1.43 |
| 7N2I | +0.620 | -0.935 | 8.15 | 5.76 | **+3.10** |
| 8IS3 | +0.552 | -0.552 | 6.21 | 5.68 | +0.65 |
| 6S0N | +0.536 | -0.321 | 5.38 | 5.62 | -0.31 |
| 2BAO | +0.328 | -0.377 | 5.58 | 5.70 | -0.15 |
| 1CS9 | -0.158 | +0.530 | 5.26 | 5.43 | -0.21 |
| 9UV5 | -0.205 | +0.256 | 5.80 | 5.73 | +0.07 |
| 6F3V | -0.309 | +0.282 | 4.68 | 5.72 | -1.29 |
| 1TOR | -0.498 | +0.569 | 4.53 | 6.43 | **-1.93** |
| 6EY3 | -0.643 | +0.728 | 4.40 | 5.73 | **-1.65** |

**Mechanism, stated plainly.** A 4069-parameter sequence-conditioned many-body potential,
trained leave-one-target-out, learned essentially ONE thing: *native peptide conformations
are more extended than a random torsion configuration is.* That is true on average, so the
model scores well on the targets whose native is unusually extended (2MK7 z=+3.9, 7N2I
z=+3.1) and **anti-ranks on the targets whose native is unusually compact** (1TOR z=-1.9,
6EY3 z=-1.7). The sign it needs is a per-target property that the model does not observe.

This is the same shape as `sequence-conditioning-hurts-the-failures.md`: a channel that is
right on the population mean and wrong exactly where it is asked a hard question.

It also explains the in-band collapse without any new assumption. The low-energy decile is
*selected on* the compactness axis, so inside the decile that axis has almost no remaining
variance -- precisely the mechanism that makes `leg_steric` exactly 0.500 in every in-band
bin in section 1.2. The learned objective spends its variance in the same place Legacy
does, and is equally spent by the time a search arrives.

---

## 3. THE REQUIREMENT: how good would a ranker have to be? -- ORACLE DIAGNOSTIC

`python -m s14.obj_noisy` -> `s12/results/s14_obj_noisy.json`.
Sweep a noisy oracle `e_sigma(x) = rmsd(x) + sigma*N(0,1)` over the enumerated spaces of
16 targets and read off, on the same evaluation sets used everywhere else, what in-band
ordering accuracy buys what structure. This is a **requirement gauge**, not a method.

| sigma (A) | in-band <1.5 A acc | in-band 0.25-0.5 A | rho_decile | decile mean A | top-100 A | argmin A |
|-----------|-------------------:|-------------------:|-----------:|--------------:|----------:|---------:|
| 0.00 (perfect) | 1.000 | 1.000 | +1.000 | 2.474 | **1.219** | **0.992** |
| 0.10 | 0.994 | 0.982 | +0.908 | 2.487 | 1.250 | 0.999 |
| 0.20 | 0.933 | 0.835 | +0.763 | 2.520 | 1.327 | 1.094 |
| 0.30 | 0.854 | 0.721 | +0.638 | 2.570 | 1.431 | 1.153 |
| 0.50 | 0.723 | 0.617 | +0.462 | 2.694 | 1.665 | 1.272 |
| 0.75 | 0.638 | 0.568 | +0.332 | 2.850 | **1.959** | 1.726 |
| 1.00 | 0.599 | 0.549 | +0.260 | 2.984 | 2.170 | **2.019** |
| 1.50 | 0.561 | 0.529 | +0.177 | 3.181 | 2.572 | 2.112 |
| 2.00 | 0.544 | 0.522 | +0.136 | 3.308 | 2.836 | 2.587 |
| 3.00 | 0.529 | 0.514 | +0.094 | 3.461 | 3.091 | 2.921 |
| 5.00 | 0.518 | 0.509 | +0.059 | 3.604 | 3.387 | 3.546 |
| 8.00 | 0.511 | 0.505 | +0.036 | 3.688 | 3.557 | 3.733 |
| inf (random) | 0.500 | 0.499 | -0.001 | 3.843 | 3.833 | 3.538 |

Space mean 3.844 A, space min 0.992 A, 16 targets.

**Read the sprint's 2.0 A target off this table.** Through the top-100 operator, 2.0 A
needs `sigma ~ 0.75 A`, i.e. **in-band ordering accuracy ~0.638 below 1.5 A separation**
and `rho_decile ~ +0.33`.

Now put the measured objectives on the same axis:

| objective | in-band <1.5 A acc | equivalent sigma | implied top-100 |
|-----------|-------------------:|-----------------:|----------------:|
| REQUIRED for 2.0 A | 0.638 | 0.75 A | 1.96 A |
| Legacy | 0.529 | ~3.0 A | ~3.09 A |
| AMBER total | 0.493 | worse than random | -- |
| 1-local prior | 0.495 | worse than random | -- |
| learned objective (LOO) | ~0.50 | ~inf | ~3.8 A |

**This is the quantitative closure of "neither energy ranks the native".** The requirement
is not marginally above what exists; the best available in-band ordering accuracy (Legacy,
0.529) sits at the `sigma = 3 A` rung, and 2.0 A needs the `sigma = 0.75 A` rung. In
Fisher-information terms the objective needs roughly an order of magnitude more in-band
ordering skill than any measured objective has, learned or physical.

Note also the operator dependence, which matters for the VQE arm: even a **perfect**
ranker gives only 2.474 A through the *low-energy decile* but 1.219 A through the *top
100* and 0.992 A through the argmin. Consistent with `operator-consumes-set-mean.md`, the
terminal operator, not the ranker, sets a large part of the outcome -- a VQE that returns
a broad low-energy population is capped near 2.5 A even with a perfect Hamiltonian.

---

## 4. THE LEARNED OBJECTIVE, MEASURED

### 4.0 Two harness corrections that had to be made first

Both are recorded because each one, left in, would have manufactured a false result.

1. **Per-target centering.** Regressing the *absolute* within-target RMSD percentile
   forces the model to explain between-target offsets it cannot observe, and it
   transferred negatively (mean rho_global +0.128, per-target range -0.59 to +0.67).
   `RidgeAccumulator.add(center=True)` subtracts each target's own column means and label
   mean exactly, without densifying, turning the fit into a within-target ranking
   regression. Mean rho_global rose to +0.184.

2. **Column standardisation.** The count columns and the global scalars differ in
   variance by orders of magnitude, so a single lambda penalises them incomparably and
   the highest-variance columns dominate *any* weight vector. The symptom was a
   label-permuted twin scoring d_decile +0.82 against the real model's +0.84 on 1CS9.
   Every number below uses `solve(standardize=True)`.

3. **A conditioning bug found while doing (2).** With centering the intercept column is
   exactly zero, so exempting it from the ridge penalty (the usual convention) left a
   zero row AND column and the system was singular; weight norms came out at ~5e5. It
   only ever shifts every score by a constant so no ranking was affected, but the solve
   is now conditioned and the norms are interpretable (17.4 rather than 547,065).

The banked leave-fold-out implementation was verified against direct accumulation to
**5.1e-15 relative** on the weight vector, so the algebraic shortcuts (leave-fold-out by
subtraction, `noint` by zeroing rows/columns, `leak` by bordering the normal equations)
are exact, not approximate.

### 4.1 Leave-one-target-out, 12 enumerated targets -- the shape of the result

Fixed lambda 0.03, 20,000 uniform training configurations per target, evaluated on the
full enumeration of the held-out target.

| statistic | value |
|-----------|------:|
| mean rho_global | **+0.184** |
| mean rho_decile (in-band) | **-0.055** |
| mean d_decile (decile RMSD minus space mean) | **-0.088 A** |
| mean d_top100 | -0.084 A |
| W/L on d_decile | 8/4 |
| per-target rho_global range | **-0.643 to +0.675** |

A model with a mean global correlation of +0.18, an in-band correlation of zero, a
structural effect under 0.1 A, and a per-target sign that flips. Section 2.3 explains
every one of those numbers with a single latent variable.

### 4.2 HEADLINE -- leave-fold-out, with nulls and a leaked-label positive control

`python -m s14.obj_headline arms`. 19 enumerated targets, pinned folds, every arm trained
ONLY on targets in other folds, lambda fixed at 0.03 for every arm (predefined), 8,000
uniform training configurations per training target, evaluated on a fixed 100,000-config
uniform subsample of each held-out target's full enumeration. The structural axis is
`d_decile` = (mean CA-RMSD of the objective's low-energy decile) - (that target's space
mean), so **0 is random and negative is better**.

**All five folds, all 19 targets, all five uniform-kind arms.**

| arm | mean d_decile | bootstrap 95% CI | W/L | drop-top-10 | drop-top-3 | per fold (0/1/2/3/4) |
|-----|--------------:|:-----------------|:---:|------------:|-----------:|:---------------------|
| **learned** | **-0.166 A** | **[-0.457, +0.116]** | 13/5 | **+0.296** | +0.019 | -0.04 / +0.27 / -0.01 / -0.79 / -0.16 |
| **seqblind (was meant as a NULL)** | **-0.338 A** | **[-0.536, -0.134]** | **16/3** | +0.022 | **-0.206** | -0.25 / -0.17 / +0.08 / -0.49 / -0.55 |
| noint (no composition x global) | -0.060 | [-0.369, +0.242] | 11/8 | +0.467 | -- | +0.09 / +0.43 / +0.03 / -0.69 / -0.06 |
| permuted (NULL) | **+0.314** | **[+0.094, +0.543]** | 7/12 | +0.717 | -- | +0.23 / +0.06 / -0.02 / +0.48 / +0.51 |
| randw (NULL) | **+0.416** | **[+0.175, +0.680]** | 5/14 | +0.898 | -- | +0.38 / +0.91 / +0.06 / +0.42 / +0.45 |
| **leak (POSITIVE CONTROL)** | **-1.405 A** | **[-1.597, -1.235]** | **19/0** | **-1.093** | -1.191 | -1.21 / -2.02 / -1.35 / -1.26 / -1.42 |

(At n=19, drop-top-10 removes more than half the sample; drop-top-3 is the check
proportionate to the brief's drop-top-20-of-126 and is given for the two arms that
matter.)

Head-to-head, paired, n = 19:

| comparison | mean diff | 95% CI | W/L | drop-top-10 |
|------------|----------:|:-------|:---:|------------:|
| learned vs permuted null | -0.481 | [-0.957, -0.040] | 14/5 | **+0.283** |
| learned vs randw null | -0.583 | [-1.078, -0.125] | 14/5 | **+0.208** |
| noint vs permuted null | -0.374 | [-0.878, +0.077] | 12/7 | +0.461 |
| learned vs noint | -0.106 | [-0.191, -0.025] | 12/6 | +0.048 |
| **learned vs seqblind** (neg = sequence helps) | **+0.172** | **[-0.040, +0.416]** | 9/10 | +0.550 |
| **leak vs learned** | **-1.239** | **[-1.559, -0.933]** | **19/0** | **-0.656** |

The leak control wins on **19 targets out of 19**, its worst single target is still
-0.87 A, and its drop-top-5 mean is -1.191 A. That is what a signal looks like in this
harness.

### 4.2b THE SEQUENCE-BLIND TWIN BEATS THE SEQUENCE-CONDITIONED MODEL -- DEMONSTRATED

This was designed as a null and it came back as the best arm. It is the single most
informative row in the table, and I did not predict it.

* **seqblind** -- residue identity removed from every pair term, every Ramachandran
  term and every interaction term, so the objective is a purely GEOMETRIC shape prior --
  reaches **-0.338 A, CI [-0.536, -0.134], W/L 16/3, median -0.30, drop-top-3 -0.206**.
  It is **the only arm whose CI against random excludes zero and which survives a
  proportionate concentration check.**
* **learned** (sequence-conditioned) is *worse*: +0.172 A relative to seqblind,
  CI [-0.040, +0.416], W/L 9/10.

**Sequence conditioning of this objective is neutral-to-harmful, and every gram of real
signal in the model is geometric.** That is exactly what section 2.3 predicted from the
flip diagnostic -- the model's only transferable content is a globular/extended shape
axis, and the 3,808 sequence-indexed pair weights are capacity spent on noise. It
independently reproduces two pinned project findings on a completely different
instrument: `phi-carries-no-sequence-signal.md` and
`sequence-conditioning-hurts-the-failures.md`.

**How much is -0.338 A worth?** Read it off the section-3 gauge. Space mean 3.844 A ->
decile 3.506 A. The `sigma = 3 A` rung of the noisy oracle gives decile 3.461 A. So the
sequence-blind geometric prior is worth **about the same as Legacy** -- it sits on the
same rung -- and the 2.0 A target needs the `sigma = 0.75 A` rung. **A real effect, four
rungs short of useful.** It is reported as DEMONSTRATED and explicitly NOT as a route to
the sprint target.

Concentration, spelled out because it is decisive: the learned arm's five largest gains
are 7T3H -1.77, 7VI4 -1.10, 8HVS -0.60, 2BAO -0.56, 1N9U -0.55. **Dropping the top three
takes the mean from -0.166 to +0.019; dropping the top five takes it to +0.101** -- i.e.
past random, in the wrong direction.

**VERDICT: REFUTED. There is no usable learned objective here.** Four independent
readings of the same table say so, and I went looking for each of them, including the
one that would have let me claim a win:

1. **The confidence interval against random crosses zero.** -0.166 A [-0.457, +0.116] on
   all 19 targets. By the sprint's own standard that is not an improvement.
2. **The effect is concentrated, and concentration reverses it.** Drop-top-10 is
   **+0.296**; dropping just the top three gains takes the mean from -0.166 to **+0.019**
   and the top five to **+0.101**. The brief's rule -- "a result carried by two targets is
   not a result" -- disposes of this directly.
3. **The margin over the nulls is the same concentration artefact, and the nulls are weak
   controls.** learned vs permuted is -0.481, CI [-0.957, -0.040], and learned vs randw is
   -0.583, CI [-1.078, -0.125]. Both *look* significant. Both have a **positive
   drop-top-10** (+0.283 and +0.208): the entire margin lives in ten of nineteen targets
   and reverses without them. And both nulls are themselves significantly WORSE than
   random (permuted +0.314 [+0.094, +0.543]; randw +0.416 [+0.175, +0.680]), so beating
   them is the "constant alpha-helix beats uniform random" trap the brief names
   explicitly. Against the only null that matters -- **random itself** -- the learned
   objective does not clear.
4. **The harness is not the problem, and this is the control that settles it.** The SAME
   pipeline with the true RMSD percentile as one extra input column returns **-1.405 A,
   CI [-1.597, -1.235], 19/0 wins, per-fold -1.21/-2.02/-1.35/-1.26/-1.42, worst single
   target still -0.87 A, and -1.093 after dropping the top ten**. A real signal in this
   harness is loud, uniform across every fold, and *survives* the concentration check.
   The learned objective has none of those three properties.

The gap between the leak control (-1.405 A, on 19 targets out of 19) and the learned
objective (-0.166 A, carried by three) is the measurement. It is not a tuning gap.

### 4.3 The ablations, and what they say

`noint` (drop the composition x global interaction block) is **-0.063 A vs the full
model's -0.168 A**, and both CIs contain zero and each other. The interaction block was
added specifically to let the model make the SIGN of its compactness response depend on
sequence composition -- the one thing section 2.3 identified as deciding whether it helps
or hurts a target. **It does not measurably do so.** HYPOTHESIS "sequence composition
predicts whether this target's native is compact or extended, well enough to fix the sign"
is therefore **not supported**; with 19 targets it is also underpowered, and I flag it as
the honest limitation rather than a refutation.

### 4.3b Terminal masking is correct but costs nothing -- DEMONSTRATED

The `nomask` arm (identical model, terminal residues NOT masked) returns d_decile values
that match `learned` target-for-target to within +/- 0.03 A on all 12 targets measured
(e.g. 6B9K +0.01 vs 0.00, 7N2I +0.74 vs +0.76, 7T3H -1.77 vs -1.77). So the two
CA-RMSD-inert terminal residues' local features were already being suppressed by the
ridge penalty, and masking them -- which section 1.4 shows is exactly right in principle
-- **buys no measurable accuracy in this model**. Reported so nobody spends time on it:
the correction to the qubit accounting is real and important, but as a modelling
intervention here it is worth ~0.00 A. It remains the default because it is correct and
free.

### 4.5 C19 -- IT IS TRANSFER, NOT CAPACITY. The deepest result in this workstream.

**Experiment mine (`s14/obj_ceiling.py`), RUN BY THE COORDINATOR** once the box freed up;
written up as C19 in `s14/coord_FINDINGS.md`. Numbers below re-derived by me from the raw
rows in `s12/results/s14_obj_ceiling.json` and they reproduce exactly.

**ORACLE DIAGNOSTIC.** Band membership (CA-RMSD <= 2.5 A) is chosen with the native; this
prices a ceiling and is not an inference-time procedure.

12 targets (`complete: true`, `n_rows 12`). Band mean RMSD 2.195 A, band min 1.031 A,
band sizes 3,072-30,000 configurations.
**(One correction to the hand-off: the coordinator's note heads this "19 targets"; it is
12 -- `obj_ceiling.main` defaults to the first twelve enumerated targets, and the 12/0
W/L columns carry the giveaway. Every other number in the note is exact.)**

| arm | in-band pairacc (<1.5 A) | 95% CI | rho | d_top100 | W/L vs chance |
|-----|-------------------------:|:-------|----:|---------:|:-------------:|
| in-sample (upper bound) | **0.986** | [0.974, 0.997] | +0.869 | -0.947 A | 12/0 |
| **same-target held out** | **0.986** | [0.973, 0.997] | +0.868 | **-0.868 A** | 12/0 |
| cross-target | **0.600** | [0.549, 0.655] | +0.111 | -0.047 A | 10/2 |

The two gaps, paired:

| gap | mean diff | 95% CI |
|-----|----------:|:-------|
| same-target held out vs in-sample (**overfitting**) | **-0.0005** | **[-0.0011, +0.0000]** |
| cross-target vs same-target held out (**transfer**) | **-0.3859** | **[-0.4366, -0.3320]** |

**THE FEATURE SPACE IS NOT THE BARRIER, AND THE FAILURE IS NOT OVERFITTING.** The same
4,125-parameter linear pair potential that scores -0.166 A leave-fold-out reaches **0.986
in-band ordering accuracy** on a target's own near-native band -- against the **0.638**
that section 3's noisy-oracle gauge says the 2.0 A target requires. And it does so on
**held-out configurations of that target**, scoring identically to in-sample: the
overfitting gap is **-0.0005 A, CI [-0.0011, +0.0000]**, i.e. indistinguishable from zero.
`d_top100 = -0.868 A` on unseen configurations of a seen target is the cleanest evidence
in this workstream that the model genuinely learns a target's in-band ordering and
generalises it perfectly *within* that target. **The entire loss -- 0.386 of accuracy,
CI [-0.437, -0.332] -- sits in transfer between targets.**

**What this closes, and it is my branch to close.** In section 7 I listed "nonlinear model
families" as open and recommended this experiment before anyone built one. It is now
answered: **a richer model is not worth building.** The linear model already saturates the
within-target problem at 0.986; nonlinearity, more parameters, equivariance, graph or
geometric architectures cannot exceed 0.986, and none of them touches transfer, which is
where all of the loss is. That is a **positive measurement closing the branch**, which is
a far stronger negative than my headline arm could have supported on its own. HYPOTHESIS
"a higher-capacity model family would find in-band signal" -> **REFUTED**.

**Third witness to one mechanism.** C19 converges with my two other diagnostics:

| diagnostic | says |
|------------|------|
| flip diagnostic (2.3) | per-target skill correlates **+0.909** with the native's z-scored Rg and **+0.951** with rho(contacts, RMSD) |
| learning curve (4.4) | **declines** with more targets; m=2 (-0.380) beats m=12 (-0.162) |
| **C19 (this section)** | within-target **0.986**, cross-target **0.600** |

All three say the same thing: **the in-band discriminating axis is real, learnable, and
per-target -- its correct sign is a property of the individual target that inference
cannot observe.** Averaging across targets cancels it, which is precisely why more
training data makes the objective worse. The declining learning curve is not a curiosity;
it is the predicted consequence.

**A comparability warning, so nobody misreads 0.600 against my 0.523.** They are different
populations and are not comparable. C19's arms are trained on and evaluated on the
ORACLE-selected near-native band (RMSD <= 2.5 A); section 4.2c's `learned` arm is trained
on uniform samples and evaluated on the objective's OWN low-energy decile of the whole
space. Cross-target 0.600 on an oracle-selected band is not a licence to expect 0.600 at
inference -- the corresponding structural number is `d_top100 = -0.047 A`, which is
nothing.

**NULL CHECK, prompted by the ENER caution.** ENER report that a random-tail null sits at
**0.524-0.527, not 0.500**, because gaps inside a tight tail are small by construction.
I verified whether that inflation touches my statistic. It does not, and for a concrete
reason: **my pair accuracy conditions on the dRMSD bin**, so pairs are matched on
separation and the null is exactly chance. Measured with a purely random objective over
8 targets:

| statistic | measured null |
|-----------|--------------:|
| binned (<1.5 A) on the low-energy decile -- **the statistic used throughout this file** | **0.5006** (sd 0.0013) |
| binned (<1.5 A) on the near-native band -- **the C19 statistic** | 0.5050 (sd 0.0072) |
| unbinned on the near-native band | 0.5025 (sd 0.0031) |

So every in-band number in this file is measured against a true 0.500 and stands as
reported. ENER's caution is correct for their statistic -- an unbinned accuracy on the
lowest 0.1%, where the RMSD spread is far tighter than my RMSD <= 2.5 A band -- and should
be applied to tail numbers of that shape, not to these.

---

## 5. STEP 5 -- THE DELIVERABLE

### 5.1 DECISION: the callable exists, is honest, and MUST NOT be used as a VQE Hamiltonian

`s14/obj_model.py` exposes `LearnedObjective`, verified deterministic and round-trip
identical through `save`/`load`:

    from s14.obj_model import LearnedObjective
    obj = LearnedObjective.load(path)     # or LearnedObjective(w)
    e = obj.score("1CS9", S)              # S: (B, n) int8/int -> (B,) float64, lower better
    e = obj("1CS9", S)                    # __call__ is score

Information contract, one line: **sequence + integer torsion configuration + the
target-held-out k=4 library + the leakage-safe 1-local empirical prior + weights trained
on other folds.** Never any native quantity. Audited AMBER-free (section 1.3).

**But the brief says "expose it IF AND ONLY IF the objective has real signal", and it does
not.** The measurement is unambiguous:

| what a VQE would get | value |
|----------------------|------:|
| learned objective, d_decile | -0.166 A, CI [-0.457, +0.116], reverses on drop-top-3 |
| sequence-blind twin, d_decile | -0.338 A, CI [-0.536, -0.134] -- real, but Legacy's rung |
| required for the 2.0 A target (section 3) | in-band accuracy 0.638; measured ~0.50 |
| leak control in the same harness | -1.405 A, 19/19 |

**Recommendation to the VQE agent: do not adopt this as a Hamiltonian.** C19
(section 4.5) sharpens the reason: the in-band ordering a VQE would need IS learnable --
0.986 within a target -- but it is target-specific, and a deployed Hamiltonian must carry
one fixed ordering for every target, where it is worth 0.600 and -0.047 A. This is also
the answer to the VQE workstream's tail result from the other side: every objective being
at chance or worse inside its own lowest 0.1% is not because tail ordering is unlearnable,
but because the learnable ordering does not transfer. Optimising it
harder will move the objective and not the structure -- and per the brief's own warning,
a VQE at n=12-16 will sit on the improving limb of exactly the curve whose limit is known
to be bad. If a geometric shape prior is wanted as a *regulariser* alongside Legacy, use
the **sequence-blind** variant, which is the only arm with a real effect, and price it at
~0.34 A on the decile mean, not more.

### 5.2 Pauli locality of this objective class -- DEMONSTRATED (analytic + verified)

Under the sprint-13 exact locality theorem, `d_ij` depends on exactly the `j-i-1`
residues strictly between i and j, so at k=4 (2 qubits/residue) a pair term at sequence
separation `s = j - i` is supported on exactly `2(s-1)` qubits.
`LearnedObjective.pauli_support(9)` returns, verified:

    {sep 2: 2 qubits, 3: 4, 4: 6, 5: 8, 6: 10, 7: 12, 8: 14}

Consequences for anyone building a Hamiltonian from this family:

* The model is **already full-register** at n=9 (a separation-8 term touches 14 qubits),
  so it is no more VQE-friendly than Legacy or AMBER. It confirms rather than escapes the
  brief's "no 2-local Ising form of a distance-based molecular objective exists in this
  encoding".
* `Featurizer(max_sep=s)` truncates the pair block and **caps Pauli weight at 2(s-1)**
  exactly. A `max_sep=3` model is 4-local. This is the one genuinely useful handle the
  model family offers, and it is available regardless of the negative result above.
* The local (Ramachandran) block is exactly 2-local per residue and the global scalars
  are full-register.

### 4.2c THE COMPLETE HEADLINE TABLE, and the mechanism that closes the case

Run completed. `s14/results/obj_hl.log`, `s12/results/s14_obj_headline.json`.
All 19 targets, all five folds, all nine arms. `inband<1.5` and `global<1.5` are the
pairwise ordering accuracies of section 1.2 (0.5 = chance), computed on the same
evaluation sets.

| arm | rho_global | rho_decile | **in-band <1.5** | global <1.5 | d_decile | d_top100 | d_argmin | W/L |
|-----|-----------:|-----------:|-----------------:|------------:|---------:|---------:|---------:|:---:|
| legacy (baseline) | +0.106 | +0.106 | **0.539** | 0.534 | +0.024 | -0.080 | +0.080 | 10/9 |
| prior (baseline) | +0.005 | -0.015 | 0.501 | 0.504 | +0.061 | +0.094 | +0.152 | 9/10 |
| learned | +0.166 | +0.048 | **0.523** | 0.565 | -0.166 | -0.301 | -0.304 | 13/6 |
| nomask | +0.168 | +0.051 | 0.525 | 0.565 | -0.164 | -0.294 | -0.248 | 13/6 |
| noint | +0.086 | +0.015 | 0.512 | 0.541 | -0.060 | -0.116 | -0.102 | 11/8 |
| **seqblind** | **+0.305** | +0.030 | **0.514** | **0.616** | **-0.337** | **-0.348** | -0.146 | **16/3** |
| band (trained in-band) | -0.051 | -0.092 | **0.470** | 0.481 | +0.181 | +0.358 | +0.475 | 7/12 |
| permuted (NULL) | -0.039 | -0.166 | 0.450 | 0.488 | +0.314 | +0.857 | +1.178 | 7/12 |
| randw (NULL) | -0.143 | -0.150 | 0.454 | 0.469 | +0.417 | +0.835 | +0.966 | 5/14 |
| **leak (CONTROL)** | **+0.995** | **+0.773** | **0.916** | 0.996 | **-1.405** | **-1.925** | **-1.949** | **19/0** |

Paired against random (0), with concentration:

| arm | d_decile [CI] drop10 | d_top100 [CI] drop10 | vs LEGACY on d_decile |
|-----|:---------------------|:---------------------|:----------------------|
| legacy | +0.024 [-0.130,+0.202] +0.304 | -0.080 [-0.329,+0.167] +0.354 | -- |
| learned | -0.166 [-0.467,+0.119] **+0.296** | -0.301 [-0.703,+0.083] +0.280 | -0.190 [-0.553,+0.133] 10/9 |
| **seqblind** | **-0.337 [-0.546,-0.126]** +0.024 | **-0.348 [-0.620,-0.070]** +0.125 | **-0.361 [-0.647,-0.076]** 13/6 |
| noint | -0.060 [-0.376,+0.248] +0.467 | -0.116 [-0.521,+0.277] +0.529 | -0.083 [-0.448,+0.249] 9/10 |
| band | +0.181 [-0.028,+0.395] | +0.358 [-0.080,+0.798] | +0.158 [-0.165,+0.443] 5/14 |
| permuted | +0.314 [+0.091,+0.547] | +0.857 [+0.364,+1.362] | +0.291 [+0.045,+0.504] 4/15 |
| randw | +0.417 [+0.179,+0.674] | +0.835 [+0.416,+1.261] | +0.393 [+0.240,+0.568] 1/18 |
| **leak** | **-1.405 [-1.604,-1.235]** -1.094 | **-1.925 [-2.146,-1.733]** -1.591 | -1.428 [-1.787,-1.113] 19/0 |

**THE MECHANISM, and it closes the case.** Compare the two pairwise-accuracy columns:

| arm | global <1.5 A | in-band <1.5 A |
|-----|--------------:|---------------:|
| legacy | 0.534 | **0.539** |
| learned | 0.565 | 0.523 |
| **seqblind** | **0.616** | 0.514 |
| leak | 0.996 | 0.916 |

**The learned objective is a better GARBAGE DETECTOR than Legacy and a WORSE in-band
ranker.** seqblind lifts global ordering accuracy from Legacy's 0.534 to **0.616** -- a
real, substantial gain -- while its in-band accuracy is **0.514, below Legacy's 0.539**.
Every angstrom it wins comes from more reliably rejecting bad structures, which is
precisely the skill section 1.2 showed is already spent by the time a search reaches the
region that matters. Its d_decile of -0.337 A is the *whole* of what that lift is worth,
and section 3 prices the 2.0 A target at in-band accuracy **0.638**.

So the sprint's question gets a sharp answer: **learning does buy discrimination, but it
buys it in the wrong place.** A 4,125-parameter objective trained on 1.28e7 exactly
labelled structures moved global ordering accuracy by +0.08 and in-band ordering accuracy
by **-0.02**.

**The `band` arm kills the obvious rescue.** Training the model exclusively on the
near-native 0.1% RMSD band -- i.e. supervising it precisely where in-band skill is needed
-- makes it **worse than random** (d_decile +0.181, in-band accuracy **0.470**, below
chance). HYPOTHESIS "the model fails in-band because it is trained on uniform samples
dominated by garbage; train it in-band and it will rank in-band" is **REFUTED**.

---

## 6. WHAT I REFUTED, INCLUDING MY OWN HYPOTHESES

Kept in place per hard rule 1.5.

| claim | tier | evidence |
|-------|------|----------|
| A learned many-body objective over (sequence, torsion configuration) can make the low-energy region structurally accurate | **REFUTED** | d_decile -0.166 A, CI [-0.467,+0.119], drop-top-3 +0.019; in-band accuracy 0.523 vs the 0.638 required. Leak control -1.405 A, 19/0, in the same harness. |
| MY HYPOTHESIS: fitting a within-target ranking (per-target centering) is what the earlier negative transfer needed | **PARTIALLY SUPPORTED** | rho_global +0.128 -> +0.184; structural effect still ~0. A necessary fix, not a sufficient one. |
| MY HYPOTHESIS: sequence composition can predict the sign of a target's compactness preference, via composition x global interactions | **NOT SUPPORTED** | `noint` -0.060 vs `learned` -0.166; both CIs contain zero and each other. Underpowered at 19 targets; flagged, not claimed as refuted. |
| MY HYPOTHESIS: the in-band failure is a *training-distribution* problem -- train on the near-native band and in-band skill appears | **REFUTED** | `band` arm: d_decile **+0.181** (worse than random), in-band accuracy **0.470** (below chance). |
| MY HYPOTHESIS: sequence conditioning is what makes a learned structural objective work here | **REFUTED, and reversed** | the sequence-blind twin is BETTER: -0.337 [-0.546,-0.126] vs -0.166 [-0.467,+0.119], and beats Legacy -0.361 [-0.647,-0.076]. |
| Masking the CA-RMSD-inert terminal residues improves the learned model | **REFUTED (no effect)** | `nomask` matches `learned` to +/-0.03 A per target; -0.164 vs -0.166 pooled. The qubit-accounting correction is real; the modelling intervention is worth ~0.00 A. |
| `amb_nonbonded` has positive in-band rank skill (+0.087) -- **my own section 1.1** | **RETRACTED** | -0.004 on the unbiased `amber_kind == 0` rows. Was an artefact of the oracle-conditioned subsample. |
| AMBER's top-100 is ~3.14-3.22 A, i.e. -0.6 A over random | **RETRACTED** | 3.694 A, i.e. -0.15 A, on the unbiased subsample. Skill was 4x overstated. |
| MY HYPOTHESIS (section 7, as originally written): a higher-capacity or nonlinear model family would find in-band signal these features cannot express | **REFUTED by C19** | the LINEAR model already reaches 0.986 in-band accuracy within a target, with a -0.0005 overfitting gap. Capacity is saturated; the loss is entirely transfer (-0.386 [-0.437,-0.332]). |
| The in-band failure might be a property of the FEATURE SPACE | **REFUTED by C19** | 0.986 within-target vs the 0.638 requirement. The features express it; the mapping does not transfer. |
| A sequence-blind GEOMETRIC shape prior beats random and beats Legacy | **DEMONSTRATED** | -0.337 A [-0.546,-0.126] vs random, 16/3; -0.361 A [-0.647,-0.076] vs Legacy, 13/6; d_top100 -0.348 [-0.620,-0.070]. |
| ...and is a route to the 2.0 A target | **REFUTED** | its gain is GLOBAL discrimination (0.616 vs Legacy 0.534) while its IN-BAND accuracy is 0.514, *below* Legacy's 0.539. Section 3 prices 2.0 A at 0.638 in-band. |
| Three torsions are inert for the CA trace (sprint brief) | **CORRECTED to FOUR** | bit-exact invariance of `rmsd` to residues 0 and n-1 on 7 targets, residues 1 and n-2 as live controls. 14 of 18 live qubits at n=9, k=4. |

## 7. WHAT IS STILL OPEN

* **Out-of-distribution robustness (generator shift)** -- `s14/obj_curve.py ood` is written
  and unrun. Train on uniform samples, test on prior-generated samples and vice versa.
  It matters for VQE because a VQE generates its own distribution. Not run for machine
  time; the negative above does not depend on it.
* ~~**Signal vs transfer**~~ -- **RESOLVED, see section 4.5 (C19).** It is transfer:
  0.986 within-target (overfitting gap -0.0005) against 0.600 cross-target.
* ~~**Nonlinear model families**~~ -- **CLOSED by C19.** The linear model saturates the
  within-target problem at 0.986; no architecture addresses the transfer gap where all
  the loss sits. Do not build one.
* **The one direction C19 leaves open:** anything that supplies the *per-target* sign of
  the in-band axis at inference -- i.e. a per-target conditioning signal, not a bigger
  objective. C19 says the ordering is learnable and target-specific; a deployed objective
  must carry one fixed ordering for all targets. Whether any legitimately-available
  per-target signal predicts that sign is untested here and is underpowered at 19
  targets (my `noint` arm was the cheap version and was null).
* **n=11 enumeration** (13 targets, folds all >= 6). `python -m s14.obj_enum 11`.
* Fold 1 carries only 2 held-out targets; its per-fold numbers are weak.

---

## 8. REPRODUCTION

Deterministic given the seeds recorded in each JSON. Set `PYTHONIOENCODING=utf-8`.
Run `python -m s12.instrument` before and after.

    # Step 1 -- extend the enumerated ground truth (heavy; ~3.5 min/target uncontended)
    python -m s14.obj_enum 10                 # every n=10 target -> s14/cache/
    python -m s14.obj_enum 11                 # the unrun batch, if wanted

    # Step 2 -- the discrimination floor, exact over each full enumeration
    python -m s14.obj_floor                   # -> s14/results/obj_floor.json

    # Step 3 -- the requirement gauge (what ranking quality would be needed)
    python -m s14.obj_noisy                   # -> s12/results/s14_obj_noisy.json

    # Step 4 -- the headline: leave-fold-out, 9 arms, nulls, leak control
    python -m s14.obj_headline arms           # -> s12/results/s14_obj_headline.json
                                              #    s14/results/obj_hl.log

    # Step 4 -- learning curve vs the leaked-label control
    python -m s14.obj_lc                      # -> s12/results/s14_obj_lc.json

    # Diagnostics
    python -m s14.obj_resid                   # compactness residualisation
    python -m s14.obj_ceiling                 # signal vs transfer (C19) -> s12/results/
                                              #    s14_obj_ceiling.json  (12 targets)

Modules: `obj_enum.py` (enumeration), `obj_model.py` (features, ridge, the callable),
`obj_train.py` (fitting/sampling), `obj_floor.py` (Step 2 statistics),
`obj_noisy.py` (requirement gauge), `obj_headline.py` (leave-fold-out arms),
`obj_lc.py` (learning curve), `obj_resid.py`, `obj_ceiling.py`, `obj_exp.py`.

### 4.4 The learning curve: COMPLETE. Flat and declining, against a control that rises.

`s14/obj_lc.py` implements the Sprint-12 diagnostic exactly -- training-set size
m in {1,2,4,8,12} against an identical harness carrying a leaked label, test set = the
whole of pinned fold 4 (7 targets), pool = the 12 targets in other folds, 4 repetitions
per size. It exploits the additivity of centered normal equations so that training on any
subset of m targets is just the sum of m precomputed per-target blocks -- one
featurization pass, no refitting. The 12 blocks build in 15 s; the full run took 1,140 s
of its own clock (many hours of wall time under 4-9 concurrent jobs from other agents).
**It completed.**

**COMPLETE, both arms** (`s12/results/s14_obj_lc.json`). 28 observations per point
(4 draws x 7 held-out targets of pinned fold 4; pool = the 12 targets in other folds).

| m (training targets) | learned d_decile | learned rho_dec | **leak d_decile** | **leak rho_dec** |
|---------------------:|-----------------:|----------------:|------------------:|-----------------:|
| **1** | **-0.348** | -0.050 | **-1.273** | **+0.551** |
| **2** | **-0.380** | -0.059 | **-1.397** | **+0.627** |
| **4** | -0.231 | -0.003 | -1.407 | +0.683 |
| **8** | -0.230 | +0.033 | -1.419 | +0.710 |
| **12 (whole pool)** | -0.162 | +0.010 | **-1.423** | **+0.728** |
| 15 (leave-fold-out, 4.2c) | -0.166 | +0.048 | -1.405 | +0.773 |

Full learned-arm detail (sd across draws; d_top100; rho_global):

| m | d_decile | sd | d_top100 | rho_global |
|--:|---------:|---:|---------:|-----------:|
| 1 | -0.348 | 0.973 | -0.235 | +0.367 |
| 2 | -0.380 | 0.799 | -0.134 | +0.411 |
| 4 | -0.231 | 0.609 | -0.228 | +0.229 |
| 8 | -0.230 | 0.677 | -0.241 | +0.229 |
| 12 | -0.162 | 0.687 | -0.173 | +0.216 |

Leak-arm detail: d_top100 -1.775 / -1.910 / -1.908 / -1.895 / -1.897;
rho_global +0.925 / +0.984 / +0.993 / +0.994 / +0.994.

**THE TWO CURVES DIVERGE, WHICH IS THE WHOLE POINT OF THE DIAGNOSTIC.** A single training target
reaches -0.348 A and two reach -0.380 A -- better than the same model trained on twelve
(-0.162 A) or, leave-fold-out, on fifteen (-0.166 A), and better than the full
sequence-blind model's -0.337 A. The curve is **not merely flat: it falls monotonically
from m=2 onward**, on the structural axis (-0.380 -> -0.231 -> -0.230 -> -0.162) and on
the objective axis (rho_global +0.411 -> +0.229 -> +0.229 -> +0.216). And `rho_decile` --
the in-band axis, the only one that decides whether a search benefits -- **never leaves
zero at any training size** (-0.050, -0.059, -0.003, +0.033, +0.010).

That is the flat-curve signature the diagnostic exists to detect, in its strongest form:
more data does not help, it hurts. The reading is that one target is enough to locate the
single geometric axis of section 2.3, and additional targets only average that axis
towards the population mean -- which is exactly where the per-target sign disagreement of
section 2.3 destroys it. It also closes the residual worry I flagged when the run first
stalled, that the curve might be *rising* and would keep rising past 15 targets: **it is
falling.**

**The leak arm supplies the calibration that was missing, and it rises.** It is already
loud at a SINGLE training target (-1.273 A) and improves monotonically with more:
-1.273 -> -1.397 -> -1.407 -> -1.419 -> **-1.423 A**. Its in-band `rho_decile` rises
monotonically too, **+0.551 -> +0.627 -> +0.683 -> +0.710 -> +0.728**, and its
`rho_global` saturates at +0.994.

That is what a learnable signal looks like in this harness: **more training targets help,
on both axes, including the in-band one.** So the harness is fully capable of expressing a
rising learning curve. Put the two side by side:

| | learned | leak control |
|---|---|---|
| d_decile, m=1 -> m=12 | -0.348 -> **-0.162** (declines) | -1.273 -> **-1.423** (improves) |
| in-band rho_decile, m=1 -> m=12 | -0.050 -> +0.010 (**pinned at zero**) | +0.551 -> **+0.728** (rises) |

**The diagnostic is now conclusive and no caveat remains.** The learned objective's curve
is flat-to-declining with its in-band axis pinned at zero across a 12-fold change in
training data, while the identical harness carrying a real signal starts loud and gets
louder. This is the Sprint-12 result reproduced exactly on a completely different
instrument, and it is the predicted consequence of C19 (section 4.5): the in-band ordering
is learnable **within** a target but not **across** targets, so adding targets averages
the axis away instead of sharpening it. What the leak arm shows at m=15 (-1.405 A,
19/0) is in section 4.2c. The two-point curve is consistent with everything in sections
2.3 and 4.2c: the model's entire transferable content is one geometric axis, and one
target is enough to find an axis.

What stands in its place is stronger on the one point the curve was there to establish.
The curve's purpose is to distinguish "not enough training targets" from "no signal". At
the LARGEST training size the curve would have reached -- 15 training targets, the
leave-fold-out setting of section 4.2c -- the measurement is already unambiguous:

| at 15 training targets | d_decile | in-band <1.5 | W/L |
|------------------------|---------:|-------------:|:---:|
| learned | -0.166 [-0.467,+0.119] | 0.523 | 13/6 |
| permuted twin (NULL) | +0.314 [+0.091,+0.547] | 0.450 | 7/12 |
| **leaked label, SAME harness** | **-1.405 [-1.604,-1.235]** | **0.916** | **19/0** |

The model is at its null and the leak is loud, **at the right-hand end of the curve**. A
learning curve can only tell you whether more data would help; a model that is
indistinguishable from its label-permuted twin at the largest size available has no slope
left to extrapolate. The curve would be confirmatory. It should still be run.

I flag one genuine loss from not running it: I cannot rule out that the curve is *rising*
and would keep rising past 15 targets. Given the mechanism in 2.3 and 4.2c (the model's
content is one geometric axis; its in-band accuracy is 0.523 against a 0.638 requirement)
I consider that unlikely, but it is a HYPOTHESIS, not a measurement.
