# Sprint 13 — QUANTUM-ARCHITECTURE findings

**Question.** What is the right quantum encoding and Hamiltonian for torsion-constrained
peptide conformational generation, and does the objective actually rank native structures?

Written incrementally. Every claim is tiered
DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / LITERATURE-SUPPORTED / PROPOSED and
classified PROVEN / STRONGLY SUPPORTED / SUPPORTED / OPEN / WEAK / REFUTED / UNKNOWN.

Code, in the order the argument runs:
`s13/qarch_lib.py` (shared machinery) — `qarch_enum.py` (complete k=4 enumerations) —
`qarch_locality.py` (geometric support, exact ANOVA, few-body verification, matched AMBER) —
`qarch_allatom.py` (scope of the support rule) — `qarch_sepfit.py` (sampled separability) —
`qarch_validity.py` (objective validity) — `qarch_robust.py` (AMBER's scale) —
`qarch_encoding.py` (encoding study) — `qarch_landscape.py` (neighbourhood pricing) —
`qarch_scale.py` (126 targets, k in {4,8}) — `qarch_budget.py` / `qarch_budget2.py`
(the budget curve).
Results: `s13/results/qarch_*.json` and `s13/results/qarch_enum_<PDB>.npz` (the cached
certified ground truth: 9 targets x 262,144 configurations).

---

## HEADLINE (stated first, as the brief requires)

**Does any candidate objective rank native structures? NO — none of them ranks. What Legacy
does instead is reject garbage, and its skill is entirely its clash term.**

Measured two ways that at first disagreed, and the disagreement turned out to be the result.

**(a) At a realistic optimisation budget, on the full instrument.** 126 tuning targets,
k=4, 12,001 uniform-random configurations each (`qarch_scale.py`):

| objective | mean rho(E, CA-RMSD) | best-of-12,001 | random draw | delta | W/L | native-snap percentile |
|---|---|---|---|---|---|---|
| genuine Legacy (11 terms) | +0.190 (median +0.232, sd 0.237) | 4.354 A | 4.971 A | **-0.617** | 85/41 | **0.325** |
| **Legacy with the steric term removed** | +0.134 | 5.366 A | 4.971 A | **+0.395** | 52/74 | 0.311 |
| leakage-safe 1-local torsion prior | +0.022 | 3.856 A | 4.971 A | **-1.115** | **93/33** | 0.134 |

Legacy beats a random draw — but **remove its clash term and the advantage inverts**, while
the clash term alone does nothing (+0.020 A). It is a **filter, not a ranker**: it still
places the near-native configuration at the **32.5th percentile** of its own energy, and its
best-of-12,001 at 4.354 A sits **2.76 A above the k=4 representation ceiling of 1.594 A**.
**A 1-local torsion prior beats it, on 93 of 126 targets.**

**(b) At an unlimited budget, with certified optima.** Nine n=9 targets, k=4, **every one of
the 262,144 configurations built and scored** (`qarch_enum.py`, cached):

| objective | rho | certified global argmin | random draw | verdict |
|---|---|---|---|---|
| genuine Legacy | +0.081 (sd 0.202, sign flips per target) | **3.920 A** | 3.781 A | **0.139 A worse than random**, 2W/7L |
| genuine AMBER ff14SB/GBn2 | +0.085 uniform, **-0.015 in-band** | 3.834 A | 3.789 A | +0.045 A, 3W/6L |
| every weighted combination | supremum +0.081 as lambda -> inf | 3.68-3.98 A | 3.78 A | no weight helps |
| ORACLE prior, any quality q >= 0.4 | +0.13 | **1.921 A** (= the chain-blind snap) | 3.781 A | -1.860 A, 8W/1L |

The space contains a **0.969 A** answer.

**(c) The reconciliation, and the sharpest result here: the budget curve is NON-MONOTONE.**
CA-RMSD reached by minimising Legacy, versus number of objective evaluations, on the
completely enumerated spaces:

    evaluations   10     30    100    300   1000   3000  10000  30000 100000  262144(exact)
    CA-RMSD     3.764  3.708  3.669  3.667  3.740  3.781  3.864  3.900  3.914   3.920
                                     ^best                                      ^worse than random (3.781)

**Optimising Legacy improves the structure for ~300 evaluations and then makes it steadily
worse for the next three orders of magnitude.** Sprint 12's single data point ("the certified
exact optimum emitted 0.169 A worse than its own anchor") is here a continuous curve on a
different architecture. Repeated on all 126 targets by sampling (60,000 draws each), **the
turn is controlled by the fraction of the space explored, not the evaluation count**: on the
33 targets with n <= 11 (>= 1.4% of the space seen) the curve degrades +0.267 A past its
minimum and ends worse than random; on the 93 with n >= 12 (<= 0.36% seen) it has not turned
yet. So a VQE at n = 12-16 with 10^4 evaluations sits on the monotone-improving part of a
curve whose limit is known to be bad — it will look like it is working. **The optimisation budget is a parameter of the answer, not of the
search: a better optimiser on this objective returns a worse structure**, which confounds
every optimiser comparison in the sprint unless objective quality and structural quality are
reported as two separate axes.

Four structural results support that and stand on their own:

* **Locality, exactly.** A CA-CA pair term is exactly **(|i-j|-1)-local with contiguous
  support** — agreement 1.0000, zero counterexamples over ~3,000 (pair, variable) cells, and
  a non-supporting variable moves the distance by exactly 0.000 A. **The rule does not extend
  to all-atom distances**: those have support one residue wider and populate the s=0 and s=1
  shells that CA-CA leaves exactly empty (75% of adjacent-residue atom pairs are
  non-constant vs 0% of adjacent CA-CA pairs). So AMBER's support structure genuinely
  differs from Legacy's, and the CA rule is the *Legacy* rule.
* **"AMBER is less local than Legacy" is REFUTED as stated** — the union of supports is the
  whole chain, so both are full-register n-local. What is true, on identical configurations
  in an exact 5-residue ANOVA: **AMBER's effective weight 3.054 vs Legacy's 2.013 (+1.040,
  CI [+0.709,+1.359], 11/12 cells), separable share 0.070 vs 0.388 (12/12)** — but only
  0.595 vs 0.659 after a rank transform, so ~80% of the raw gap is the vdW singularity.
* **Legacy is 95.1% interaction variance** in the full 9-variable space (Sprint 12's assembly
  Hamiltonian: 2.4%), peaking at order 3-4, and **the L2-optimal 2-body Hamiltonian explains
  only 20.4% of its variance with its minimiser at the 13th percentile of the true energy.**
  A few-body decomposition is *verified to fail*; the sampled-bitstring CVaR route is the only
  honest one. The torsion prior is the sole exception: 1-local to 2e-15.
* **Encoding: direct binary index at k=4** (25.9 qubits mean, 32 max), whole register
  feasible, no penalty. The move-size argument for one-hot is refuted (a 1-qubit binary flip
  and a 2-qubit one-hot move are the same size), but the *neighbourhood* argument is real and
  now priced: single-qubit-flip descent finds Legacy's optimum from 13.5% of starts vs 35.7%
  for residue moves, reaches 0.187 A worse structure against the oracle objective, and
  **manufactures local minima in the torsion prior, which provably has none**. Binary is still
  chosen, because 2x-4x the qubits is the difference between certified ground truth and none.

---

## 0. Instrument reproduction

Before anything else, the coordinator's representation-ceiling instrument was reproduced
through my own code path (`qarch_lib.Space`), which builds the library independently:

| quantity | coordinator (`s13/results/ceiling.json`) | reproduced here |
|---|---|---|
| 1A13, k=4, nearest-state snap CA-RMSD | 2.081447867122067 | 2.081447867122067 |

Bit-identical. `Space.PHI/PSI` therefore address the same
`torsion_lib2.library_for(seq, k, exclude_seq=seq)` state space the ceiling sweep measured,
and every number below sits on that instrument.

Timing facts measured on this box, since they set what is affordable:

| call | cost | note |
|---|---|---|
| `core.energy.components_batch` (build + all 11 terms) | **0.56 ms** / configuration | batched, n=9–14 |
| `core.amber.single_point` (ff14SB/GBn2, steps=-1) | **12.5–13.2 ms** / configuration warm | n=9, `threads=1` |
| `core.amber.builder_for(...).energy(bits)` | **≫** that | **not a single point** — the default `minimization_steps=50` means it minimises |

**Correction to the sprint's stated central fact.** The brief records AMBER at 6 ms warm.
On this box, at n=9 and `threads=1`, `single_point` is **12.5–13.2 ms**, i.e. the
AMBER : Legacy asymmetry is **~23×**, not ~20×. It is the same order and the same
conclusion, but the enumeration budgets below are sized on the measured number.
Also, `AmberHamiltonian.energy()` obtained through `builder_for` **minimises 50 steps**;
anything wanting a single point must go through `core.amber.single_point`. That trap would
silently turn a "20× costlier" experiment into a 1000× costlier one.

---

## 1. THE LOCALITY ANALYSIS — the geometric mechanism (PROVEN)

**Hypothesis (brief §4).** In torsion space the chain builds sequentially, so the CA–CA
distance between residues i and j depends on *every* torsion between them; a pairwise
distance term is |i−j|-local, not 2-local.

**Method (`s13/qarch_locality.py::geometry_locality`, DEMONSTRATED, no native read).**
64 random base configurations per target; every residue m moved to each of its k=8 library
states in turn; the median change in every CA–CA pair distance recorded. Residue m
"supports" pair (i,j) if that median |Δd_ij| exceeds 0.10 Å. Targets 1A13 (n=14),
1CS9 (n=9), 1KZ2 (n=16).

**Result — the rule is exact, not approximate.**

| target | n | mean support size | agreement with "strictly between" rule | inside-rate | outside-rate | median \|Δd\| inside |
|---|---|---|---|---|---|---|
| 1A13 | 14 | 4.00 | **1.0000** | 1.000 | 0.000 | 1.95 Å |
| 1CS9 |  9 | 2.33 | **1.0000** | 1.000 | 0.000 | 1.75 Å |
| 1KZ2 | 16 | 4.67 | **1.0000** | 1.000 | 0.000 | 1.94 Å |

The support of `d_ij` is **exactly the j−i−1 residues strictly between i and j**, contiguous,
with *zero* variables outside it on every one of the 5,000+ (pair, variable) cells tested.
Measured support size as a function of separation is `s − 1` for every separation `s`, on
every target, with no scatter:

    separation s      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    support size      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14

**Interpretation, sharper than the brief's claim.** A CA–CA pair term is not merely
|i−j|-local: it is **exactly (|i−j|−1)-local with contiguous support**, because `d_ij` is an
internal coordinate of the CA_i…CA_j sub-chain and the ideal-geometry CA–CA virtual bond
length is torsion-independent, so the two end residues drop out. Consequences:

* **s=1 pairs are constant** — adjacent-CA terms carry no torsion information at all.
* A distance term at separation s costs `(s−1)·log2(k)` qubits of support. At k=4 a
  separation-8 pair is a **14-qubit** interaction. There is no 2-local Ising form for it.
* The union of supports over all pairs is the whole chain, so **any** energy that sums over
  pairs is n-local — full-register — in the torsion variables. Legacy and AMBER both are.
* A single state flip moves an in-support pair distance by a median **1.75–1.95 Å**. The
  landscape is not smooth in the encoded variables; one qubit-block flip is a large
  geometric move.

**Classification: PROVEN** (exact agreement, no counterexample, mechanism understood).
This is the mechanism that forbids a naive 2-local QUBO mapping of any distance-based
molecular objective in torsion space, and it is the structural difference from Sprint 12's
assembly Hamiltonian, which *was* exactly 2-local because each segment's coordinates were
placed by Kabsch onto a fixed anchor and therefore did not depend on the other segments.

Result file: `s13/results/qarch_locality_geom.json`.

---

## 2. THE LOCALITY ANALYSIS — the exact functional ANOVA (PROVEN)

**Method (`qarch_locality.py::anova`, DEMONSTRATED).** The nine n=9 tuning targets at k=4
give 4^9 = 262,144 configurations, so the whole state space is enumerable
(`qarch_enum.py`, cached in `s13/results/qarch_enum_<PDB>.npz`). Under the uniform product
measure on the torsion variables every energy term has an exact functional ANOVA

    F(s) = Σ_U f_U(s_U),   Var(F) = Σ_U ‖f_U‖²   (orthogonal)

computed here by Möbius inversion over all 2⁹ subsets. The residual |Var(F) − Σ‖f_U‖²| is
reported per term and is **1e-16 to 3e-15 relative** — float noise, so these are exact
decompositions, not estimates. Definitions: `V_d` = order-d variance share; total Sobol
`ST_i`; **effective weight** `W_eff = Σ_i ST_i` (= 1 for a separable function, larger with
interaction); **participation** `W_par = (ΣST)²/Σ ST²`; interaction-graph edge when the
*pure* pair term carries ≥ 1% of Var(F).

**Result (mean over 8 of the 9 enumerated targets; the ANOVA was run before 9UV5's
enumeration finished and the memory gate has not permitted a re-run since, so section 2
is n=8 while section 3 is n=9. Residual column is the exactness check).**

| term | V₁ | V₂ | V₃ | V₄ | V₅₊ | W_eff | W_par | graph deg | residual/Var |
|---|---|---|---|---|---|---|---|---|---|
| **torsion prior** | **1.000** | −0.000 | 0.000 | −0.000 | 0.000 | **1.000** | 7.39 | 0.00 | 2.1e-15 |
| `leg_torsion` (Ramachandran) | **1.000** | −0.000 | 0.000 | −0.000 | 0.000 | **1.000** | 7.60 | 0.00 | 4.2e-16 |
| **Legacy TOTAL** | **0.049** | 0.155 | 0.311 | 0.275 | 0.210 | **3.535** | 6.16 | 1.22 | 2.4e-16 |
| `leg_steric` | 0.050 | 0.160 | 0.301 | 0.275 | 0.213 | 3.536 | 6.16 | 1.36 | 3.3e-16 |
| `leg_hbond_longrange` | 0.019 | 0.072 | 0.178 | 0.253 | **0.479** | **4.449** | 6.97 | 0.44 | 2.8e-16 |
| `leg_electrostatic` | 0.291 | 0.290 | 0.202 | 0.112 | 0.106 | 2.524 | 4.98 | 0.94 | 7.5e-17 |
| `leg_contact` | 0.273 | 0.436 | 0.183 | 0.075 | 0.033 | 2.170 | 5.56 | 1.25 | 2.3e-16 |
| `leg_compactness` | 0.381 | 0.331 | 0.192 | 0.069 | 0.026 | 2.032 | 5.16 | 1.25 | 2.6e-16 |
| `leg_solvation` | 0.477 | 0.339 | 0.120 | 0.043 | 0.021 | 1.798 | 5.91 | 1.31 | 3.1e-15 |
| **ORACLE CA-RMSD itself** | **0.164** | 0.288 | 0.308 | 0.169 | 0.070 | **2.709** | 5.13 | 1.25 | 1.6e-15 |

**Readings.**

1. **The torsion prior is 1-local to machine precision.** V₁ = 1.000000, W_eff = 1.000,
   interaction graph empty, and its exact 1-body truncation reproduces it with
   `var_explained = 1.000`, `ρ = 1.000` and the identical argmin. So the prior term is a
   genuine, exact, sum-of-single-qubit-blocks Hamiltonian — the only one here that is.
   Same for Legacy's Ramachandran term, which *is* a per-residue torsion function.

2. **Legacy is overwhelmingly many-body: 95.1% of its variance is interaction.**
   Sprint 12's assembly Hamiltonian had **2.4%**. That is a **39× larger interaction share**,
   and it is the sharpest available statement that the torsion formulation is a genuinely
   harder — and more interesting — combinatorial instance than anything this project has
   optimised before. The variance peaks at **order 3–4**, not order 2.

3. **The pure *pairwise* interaction decays fast with sequence separation**, which is the
   quantitative form of the compounding mechanism (mean share of Var(F)):

   | separation s | 1 | 2 | 3 | 4 | 5 | ≥6 |
   |---|---|---|---|---|---|---|
   | Legacy total | 0.0105 | 0.0052 | 0.0046 | 0.0012 | 0.0002 | 0.0000 |
   | ORACLE CA-RMSD | 0.0323 | 0.0034 | 0.0010 | 0.0001 | 0.0000 | 0.0000 |

   So the *order-2* structure is short-ranged, but it is only 15% of the variance; the
   long-range coupling is carried at order 3+, where a pair graph cannot see it. **Truncating
   the interaction range is not the same as truncating the interaction order, and it is the
   order that matters.**

4. **The true objective is itself 83.6% many-body** (CA-RMSD, V₁ = 0.164). This is the exact,
   quantitative form of the coordinator's snap→descent gap: any 1-local objective — any
   per-residue torsion prior, however good — can address at most **16.4%** of the variance of
   what we actually want to minimise.

**Classification: PROVEN.**

### 2b. The few-body approximation, verified — and it fails

The ANOVA truncation `F_{≤d} = Σ_{|U|≤d} f_U` is, by orthogonality, the **L2-optimal** d-body
approximation. If it fails, no hand-built d-body Hamiltonian can succeed. Verified exactly:

| term | d | var explained | ρ with full | true rank of approx argmin | approx top-1 RMSD | full top-1 | pool mean |
|---|---|---|---|---|---|---|---|
| torsion prior | 1 | **1.000** | **1.000** | **0**/262,144 | 3.575 | 3.575 | 3.814 |
| Legacy total | 1 | 0.049 | 0.269 | 53,845/262,144 | 4.020 | 3.970 | 3.814 |
| Legacy total | 2 | **0.204** | 0.297 | **33,836**/262,144 | 3.811 | 3.970 | 3.814 |
| ORACLE CA-RMSD | 1 | 0.164 | 0.383 | 40,594/262,144 | 2.349 | 0.961 | 3.814 |
| ORACLE CA-RMSD | 2 | 0.452 | 0.648 | 17,898/262,144 | 2.134 | 0.961 | 3.814 |
| `leg_contact` | 2 | 0.709 | 0.867 | 14,199/262,144 | 3.944 | 3.933 | 3.814 |
| `leg_hbond_longrange` | 2 | 0.091 | 0.146 | 436/262,144 | 4.751 | 5.432 | 3.814 |

**The best possible 2-local Ising Hamiltonian for Legacy reproduces 20.4% of its variance,
correlates ρ = +0.297 with it, and its minimiser is the 33,836-th best of 262,144
configurations (the 13th percentile).** That is not "reproduces the true objective to
numerical precision"; it is a different function. **Route (b) of the brief is therefore
closed for Legacy and for AMBER, and the sampled-bitstring CVaR route that
`run_global_cvar_vqe` already implements is the only honest option.** The single exception
is the torsion prior, which decomposes exactly at order 1 (residual 0.0).

**Classification: PROVEN (REFUTES any 2-local QUBO mapping of Legacy/AMBER in torsion space).**

Result files: `s13/results/qarch_locality.json`, `qarch_locality_geom.json`.

---

## 3. OBJECTIVE VALIDITY — the headline, and it is negative (PROVEN on this instance family)

> **PARTLY SUPERSEDED — read sections 10 and 11 with this.** Everything below is exact,
> but it is measured on nine n=9 targets with a COMPLETE enumeration, i.e. at an infinite
> optimisation budget on the weakest length bin of the instrument. On all 126 targets at a
> realistic budget Legacy's best-of-12,001 BEATS a random draw by 0.617 A (85W/41L), and the
> reconciliation is a measured non-monotone budget curve (section 11). The claims that
> survive unchanged are: the near-native configuration sits at the 32.5th percentile of
> Legacy's energy, nothing has in-band skill, the 1-local prior beats every many-body
> objective, and the certified global optimum is worse than stopping early.

**This section answers the question the coordinator called the most decisive in the sprint.
Nothing here is optimised; everything is measured over the complete state space.**

**Config.** 9 targets (all the n=9 tuning targets), k=4, library
`torsion_lib2.library_for(seq, 4, exclude_seq=seq)`, **every one of the 262,144
configurations scored** — 2,359,296 structures built and scored in total. Legacy is
`core.energy.components_batch` with `DEFAULT_WEIGHTS`; AMBER is genuine
`core.amber.single_point` (ff14SB/GBn2, no minimisation) on a labelled 2,900–3,000
configuration subsample per target (uniform / prior-sampled / ORACLE near-native band), by
cost. Seed 0. `top1_rmsd` **averages over the tied-argmin set** — the Sprint-12 trap where
`np.argmin` on a tied signal reads the enumeration order and invents a winner.

**Reference quantities (mean over the 9 targets):** the space's **best** configuration is
**0.969 A**; the chain-blind native **snap** is 1.921 A; a **random draw** is 3.781 A.
So the answer is in the space and there is 2.8 A of headroom for a working objective.

### 3a. Every native-free objective, over the complete space

| objective | mean rho(E,RMSD) | median rho | sd | argmin RMSD | random draw | **delta** | W/L | bootstrap CI | native-snap %ile | best-config %ile |
|---|---|---|---|---|---|---|---|---|---|---|
| **Legacy total** | +0.081 | +0.125 | 0.202 | 3.920 | 3.781 | **+0.139** | 2/7 | [-0.22, +0.48] | **0.629** | 0.665 |
| **torsion prior (empirical, leakage-safe)** | +0.011 | +0.008 | 0.020 | 3.636 | 3.781 | -0.145 | 6/3 | [-0.55, +0.20] | 0.180 | 0.349 |
| z(prior) + z(Legacy), derived weights | +0.054 | — | — | 3.980 | 3.781 | +0.199 | — | — | 0.181 | 0.361 |
| ORACLE prior, q = 0.4 | +0.132 | +0.132 | 0.088 | **1.921** | 3.781 | **-1.860** | 8/1 | [-2.66, -0.92] | 0.000 | 0.072 |
| ORACLE prior, q = 1.0 | +0.127 | +0.124 | 0.087 | **1.921** | 3.781 | **-1.860** | 8/1 | [-2.66, -0.92] | 0.000 | 0.057 |

**Legacy's certified global optimum over the entire 262,144-configuration space is
0.139 A WORSE than picking a configuration at random**, it loses to random on **7 of 9**
targets, and the near-native configuration sits at the **62.9th percentile** of its energy
(96.4th on 7N2I, 90.7th on 2MK7, 90.2nd on 6F3V). Its Spearman correlation with CA-RMSD has
mean +0.081 and standard deviation 0.202 — it changes sign across targets. The energy-decile
scatter is flat:

    Legacy energy decile   1     2     3     4     5     6     7     8     9    10
    mean CA-RMSD (A)      3.77  3.79  3.75  3.73  3.72  3.72  3.70  3.65  3.72  4.25

**This is Sprint 12's failure reproduced exactly, on the new architecture, with certified
optima instead of an optimiser's output.** Optimising Legacy harder in torsion space moves
the structure away from the native.

**Per-target, so no mean is hiding anything:**

| pdb | space best | snap | pool mean | Legacy argmin | rho | native %ile | best-config %ile |
|---|---|---|---|---|---|---|---|
| 1CS9 | 1.270 | 5.216 | 4.064 | 4.066 | -0.082 | 0.451 | 0.655 |
| 2MK7 | 0.675 | 0.944 | 4.869 | 5.037 | -0.029 | 0.907 | 0.936 |
| 2P5H | 1.461 | 2.135 | 3.788 | 4.222 | +0.125 | 0.754 | 0.862 |
| 6EY3 | 0.989 | 2.884 | 3.696 | 3.820 | +0.253 | 0.551 | 0.057 |
| 6F3V | 0.828 | 1.081 | 3.568 | 2.736 | +0.276 | 0.902 | 0.525 |
| 6S0N | 0.887 | 1.240 | 3.293 | 2.825 | +0.044 | 0.018 | 0.577 |
| 7N2I | 1.027 | 1.447 | 4.157 | 5.160 | -0.346 | 0.964 | 0.940 |
| 8IS3 | 0.549 | 1.009 | 3.080 | 3.893 | +0.145 | 0.834 | 0.907 |
| 9UV5 | 1.033 | 1.332 | 3.513 | 3.516 | +0.346 | 0.280 | 0.531 |

### 3b. Legacy term by term — only the 1-local term has any skill

| Legacy term | mean rho | argmin RMSD | delta vs random | W/L | bootstrap CI |
|---|---|---|---|---|---|
| `torsion` (Ramachandran, **exactly 1-local**) | +0.009 | 3.302 | **-0.479** | 6/3 | **[-0.99, -0.03]** |
| `steric` | +0.156 | 3.691 | -0.090 | 6/3 | [-0.19, +0.01] |
| `coop_helix` | +0.000 | 3.728 | -0.053 | 5/4 | [-0.28, +0.18] |
| `hbond_local` | -0.031 | 3.827 | +0.047 | 5/4 | [-0.30, +0.41] |
| `contact` | +0.040 | 3.950 | +0.169 | 5/4 | [-0.60, +1.00] |
| `aromatic` | — | 3.982 | +0.201 | 0/4 | [+0.01, +0.42] |
| `compactness` | -0.177 | 4.158 | +0.377 | 4/5 | [-0.06, +0.81] |
| `solvation` | -0.128 | 4.455 | **+0.674** | 2/7 | **[+0.03, +1.28]** |
| `coop_sheet` | — | 4.679 | +0.898 | 1/5 | [+0.23, +1.57] |
| `hbond_longrange` | -0.071 | 5.360 | **+1.579** | 2/7 | **[+0.67, +2.37]** |

**The only Legacy term whose CI excludes zero on the good side is `torsion` — the one term
that is exactly 1-local, i.e. a torsion prior in disguise.** The two terms whose CIs exclude
zero on the *bad* side (`solvation`, `hbond_longrange`) are two of the three most many-body
terms in the ANOVA table of section 2. The ordering "more many-body, worse ranking" is
visible across the whole table and is the most important structural fact in this section.

### 3c. Genuine AMBER, on its labelled populations

| population | objective | mean rho | argmin RMSD | pool mean | delta | W/L | n/target |
|---|---|---|---|---|---|---|---|
| uniform subsample | **AMBER total** | +0.085 | 3.834 | 3.789 | +0.045 | 3/6 | 1,194 |
| uniform subsample | Legacy total | +0.093 | 3.493 | 3.789 | -0.296 | 6/3 | 1,194 |
| uniform subsample | `amb_nonbonded` | +0.091 | 4.117 | 3.789 | +0.328 | 3/6 | 1,194 |
| prior-sampled | **AMBER total** | +0.123 | 3.517 | 3.844 | -0.327 | 6/3 | 1,108 |
| prior-sampled | Legacy total | +0.098 | 3.887 | 3.844 | +0.043 | 4/5 | 1,108 |
| ORACLE near-native band | **AMBER total** | -0.015 | 1.589 | 1.703 | -0.114 | 2/7 | 600 |
| ORACLE near-native band | Legacy total | -0.041 | 1.748 | 1.703 | +0.045 | 2/7 | 600 |

**AMBER does not rank either.** rho is about +0.1 on random populations and **negative in the
near-native band**, where the ranking would actually have to work; its argmin is within
0.33 A of a random draw and the win/loss record is a coin flip. On the in-band population —
the only metric this project has found to transfer — *nothing* has skill: Legacy -0.055,
AMBER -0.015, empirical prior +0.092, ORACLE prior +0.202.

One weak positive worth recording and not promoting: inside the ORACLE band, `leg_steric`
alone picks 1.171 A against a band mean of 1.693 A (-0.522, 7W/2L) although its in-band rho
is -0.045. That is an argmin effect on an oracle-defined population, so it is an
**ORACLE DIAGNOSTIC** and **OPEN**, not a deployable signal.

### 3d. AMBER's scale is not standardisable on unminimised torsion configurations (PROVEN)

Deriving the normalisation as planned (per-residue mean and sqrt(n)-scaled sd on a
native-free reference population) produced `amber_per_res_mu = 5.2e14 kcal/mol`,
`amber_per_sqrt_res_sd = 2.8e16`. The reason, measured (`s13/results/qarch_robust.json`):

| pdb | AMBER median (kcal/mol) | mean | max | fraction above 1e4 kcal/mol |
|---|---|---|---|---|
| 1CS9 | 1.3e+02 | 2.6e+11 | 1.5e+14 | 40% |
| 2P5H | 4.4e+04 | 7.4e+15 | 8.3e+18 | 60% |
| 7N2I | 3.3e+05 | 1.6e+13 | 1.9e+16 | 68% |
| 9UV5 | 2.9e+03 | 1.9e+16 | 2.3e+19 | 45% |

**40-68% of ideal-geometry discrete-torsion configurations have an AMBER single-point energy
above 1e4 kcal/mol** — van der Waals singularities from unminimised sidechains. The
distribution has no usable mean or variance, so variance standardisation is arithmetically
invalid, and the measured consequence is that `z(prior) + z(AMBER)` emits **+0.545 A worse
than random**. This is a correctness result, not a tuning detail: **any AMBER-in-the-loop
VQE must either minimise (much more than 13 ms) or use a scale-free/robust transform.**
With scale-free normalisation the conclusion does not change:

| objective (rank-normalised, uniform subsample) | mean rho | argmin RMSD | pool | delta | W/L |
|---|---|---|---|---|---|
| rank(AMBER) | +0.085 | 3.834 | 3.789 | +0.045 | 3/6 |
| rank(Legacy) | +0.093 | 3.493 | 3.789 | -0.296 | 6/3 |
| rank(prior)+rank(Legacy) | +0.069 | 3.399 | 3.789 | -0.390 | 5/4 |
| rank(prior)+rank(Legacy)+rank(AMBER) | +0.098 | 3.708 | 3.789 | -0.081 | 5/4 |
| rank(prior)+rank(AMBER) | +0.068 | 3.677 | 3.789 | -0.112 | 3/6 |

Every delta is inside the noise of a random draw, and no combination reaches even a seventh
of the 2.8 A of available headroom.

**lambda-sweep (DIAGNOSTIC — reported, never used to choose a weight).** Mean rho of
`z(prior) + lambda * z(Legacy)`: 0 -> +0.011, 0.25 -> +0.031, 0.5 -> +0.042, 1 -> +0.054,
2 -> +0.066, 4 -> +0.078, infinity -> +0.081. Monotone, and its supremum is Legacy alone at
rho = +0.081. **There is no weight at which the combination ranks; there is nothing to tune.**

### 3e. What the prior can and cannot do

The ORACLE prior parameterises the predictor agent's quality by a single knob q (probability
mass on the correct state). Its argmin is **1.921 A at every q >= 0.4** — because the argmin
of a product distribution is the per-residue mode, i.e. **exactly the chain-blind snap**. A
better predictor moves the *distribution*, not the argmin. Combined with section 2's ANOVA:

* a perfect 1-local objective returns 1.921 A; the space contains 0.969 A;
* the true objective is 83.6% many-body, so a 1-local objective can address at most 16.4%
  of its variance;
* **the 0.95 A between the snap and the space best is precisely the part of the problem that
  requires a many-body objective, and no many-body objective tested here can find it.**

That gap is the sprint's central open problem, and it is an *objective* problem, not an
encoding, ansatz or optimiser problem.

**Classification: the negative result is PROVEN on this instance family** (9 targets,
complete enumeration, certified optima, no sampling error at all for Legacy and the prior)
and **STRONGLY SUPPORTED** for AMBER (labelled subsamples, ~1,200 configurations per
population per target). **Limitation:** all nine targets are n=9 at k=4. Whether Legacy's
ranking improves at larger n or larger k is NOT established here.

Result files: `s13/results/qarch_validity.json`, `qarch_robust.json`,
`qarch_enum_<PDB>.npz` (the cached certified ground truth), `qarch_enum.json`.

---

## 4. THE ENCODING STUDY (`s13/qarch_encoding.py`, `s13/results/qarch_encoding.json`)

Five encodings of a k-state-per-residue torsion space, on 1CS9 (n=9) and 1A13 (n=14) at
k=4 and k=8. `gray_sorted` is Gray coding applied **after** reordering each residue's
library along a 1-D angular traversal of the (phi,psi) torus, so that adjacent codeword means
adjacent torsion; `hierarchical` puts a recursive 2-means basin split in the high bits.

### 4a. The book-keeping

| encoding | qubits | configurations | feasible fraction of the register | penalty | qubit flips per state change |
|---|---|---|---|---|---|
| binary / gray / gray_sorted / hierarchical | n log2 k | k^n | **1.000** | none | 1 to log2 k (mean Hamming 1.33 at k=4, 1.71 at k=8) |
| one-hot | n k | k^n | (k/2^k)^n | **required** | 2 |

Concretely: 1A13 (n=14) at k=4 is **28 qubits binary vs 56 one-hot**, feasible fraction
3.7e-09; at k=8 it is **42 vs 112 qubits**, feasible fraction 8.5e-22. Given that the
coordinator's ceiling sweep puts the deployable operating point at k=4 (25.9 qubits mean,
32 max — inside `lightning.qubit`'s 30-qubit range and inside exact statevector for the
shortest targets), **one-hot's factor of 2 to 4 in qubits is the difference between
simulable-with-certified-ground-truth and not.** Qubit count is the binding constraint.

### 4b. Move locality — the trade the brief asked to be measured, and it is not real

Flip one qubit of a random configuration; measure what moves (medians over 64 base
configurations x every residue x every bit):

| target / k | encoding | flips | torsion move (deg) | structure move (CA-RMSD, A) | median abs dE Legacy | median abs dE prior |
|---|---|---|---|---|---|---|
| 1CS9 k=4 | binary | 1 | 194.5 | 1.73 | 0.63 | 0.40 |
| 1CS9 k=4 | gray | 1 | 122.7 | 1.90 | 0.71 | 0.72 |
| 1CS9 k=4 | gray_sorted | 1 | 190.0 | 1.73 | 0.59 | 0.72 |
| 1CS9 k=4 | hierarchical | 1 | 190.0 | 1.76 | 0.59 | 0.73 |
| 1CS9 k=4 | onehot | 2 | 190.0 | 1.72 | 0.60 | 0.73 |
| 1A13 k=8 | binary | 1 | 168.8 | 2.10 | 1.49 | 0.74 |
| 1A13 k=8 | gray | 1 | 151.4 | 2.39 | 1.82 | 0.59 |
| 1A13 k=8 | gray_sorted | 1 | 146.3 | 2.11 | 1.51 | 0.39 |
| 1A13 k=8 | hierarchical | 1 | 115.7 | 2.27 | 1.64 | 0.74 |
| 1A13 k=8 | onehot | 2 | 173.1 | 2.28 | 1.50 | 0.49 |

**The brief's stated trade — "binary makes single-qubit flips move a residue between distant
states, while one-hot keeps moves local" — is REFUTED as a quantitative claim.** One-hot's
2-qubit move is 172-190 deg of torsion and 1.7-2.3 A of structure: statistically
indistinguishable from binary's 1-qubit move. `hierarchical` and `gray_sorted` do cut the
median torsion step (169 -> 146 -> 116 deg at 1A13 k=8, a real 30% reduction) but the
**structure** move barely follows (2.10 -> 2.11 -> 2.27 A), and the Legacy energy step does
not move at all.

The reason is section 1: a residue's k library states are spread over the whole Ramachandran
map, so even the geometrically nearest alternative state is a large move, and by the
locality rule that move propagates to every downstream CA. **At any k that fits on a
simulator, the state space is not geometrically local, so no bit labelling can make moves
local.** Locality would have to come from a much larger k (k=32 gives 64.8 qubits mean, i.e.
MPS-only) or from a continuous refinement stage after the discrete search.

**Classification: REFUTED (as stated); the underlying quantity is measured and is
encoding-independent.**

### 4c. Expressivity — where the encodings genuinely differ, exactly

A single RY layer with no entanglement is an independent Bernoulli per qubit, so the
per-residue distributions it can express are exactly the rank-1 (independent-bits) ones, and
the minimum KL from a target pi to that family **equals the mutual information between the
bits of the state index under that labelling** — a closed form, computed exactly, no fitting.
One-hot with post-selection on the feasible subspace is instead a full-rank family
(`p(s) ~ q_s prod_{t != s} (1 - q_t)` reaches any pi), so its minimum KL is exactly 0.

Mean min-KL per residue, in nats:

| target/k | target distribution | binary | gray | gray_sorted | hierarchical | one-hot |
|---|---|---|---|---|---|---|
| 1CS9 k=4 | empirical prior | 0.0136 | 0.0313 | 0.0791 | 0.1024 | **0.0** |
| 1CS9 k=8 | empirical prior | 0.0707 | 0.1227 | 0.1071 | 0.1369 | **0.0** |
| 1A13 k=4 | empirical prior | 0.1110 | 0.0950 | 0.0954 | 0.1095 | **0.0** |
| 1A13 k=8 | empirical prior | 0.0756 | 0.0854 | 0.0442 | 0.1364 | **0.0** |
| any k=4 | ORACLE prior q=0.8 | 0.0652 | 0.0652 | 0.0652 | 0.0652 | **0.0** |
| any k=8 | ORACLE prior q=0.8 | 0.1766 | 0.1766 | 0.1766 | 0.1766 | **0.0** |
| 1A13 k=8 | deliberately bimodal | **1.2956** | 0.0565 | 0.0565 | 0.6284 | **0.0** |
| 1CS9 k=4 | deliberately bimodal | 0.0 | 0.0 | 0.6787 | 0.4525 | **0.0** |

Three readings.

1. For a **sharp unimodal** predictor the labelling is irrelevant: the ORACLE-prior row is
   identical (0.0652 at k=4, 0.1766 at k=8) for every binary-type labelling, because a
   one-dominant-state distribution has label-invariant bit mutual information. So if the
   predictor agent delivers a confident prior, the encoding choice does not cost anything
   at one layer.
2. For a **bimodal** per-residue target the labelling matters enormously and unpredictably
   (0.0 to 1.30 nats, and the ordering of the labellings flips between targets). No fixed
   bit labelling is safe for multimodal priors.
3. **One-hot is exactly expressive at one layer, for free, always.** That is a genuine,
   provable advantage — and it costs 2x-4x the qubits plus a penalty, which at the sprint's
   operating point is the wrong trade.

### 4d. The real ansatz at depth

The genuine `core.quantum.StatevectorCircuit` (RY/CNOT, exact statevector, analytic
parameter-shift gradient, no sampling) was fitted by Adam to a deliberately **multimodal**
target: uniform over the 8 lowest-CA-RMSD configurations of a real 12-qubit sub-register
(6 residues of 1CS9 at k=4). ORACLE DIAGNOSTIC in its target only; the measured quantity is
a property of the ansatz. Achieved KL (nats, one seed, 250 Adam iterations):

| labelling | L=1 | L=2 | L=3 | L=4 | L=6 |
|---|---|---|---|---|---|
| binary | 1.386 | 0.694 | 0.693 | **0.000** | 0.316 |
| gray_sorted | 0.693 | 1.386 | 0.693 | 0.665 | 0.696 |
| hierarchical | 0.693 | 0.693 | **0.000** | 1.430 | **0.000** |

KL = ln(8/j) is exactly the value for a circuit that has put its mass on j of the 8 modes, so
1.386 = ln 4 means 2 modes captured, 0.693 = ln 2 means 4, and 0.000 means all 8.
**The ansatz IS expressive enough — depth 3-4 reaches KL = 0 exactly — but the optimisation
reaches it unreliably and non-monotonically in depth.** Expressivity is not the barrier;
trainability is (and that is the trainability agent's question).
**Caveat: one seed per cell. Classification OPEN, not SUPPORTED.** A 6-seed replication
was launched (`s13/results/qarch_expressivity_seeds.json`) and had not finished when this
was written; until it lands, the depth claim stands at OPEN and the non-monotonicity in
depth is best read as optimiser variance rather than as an expressivity statement. One-hot at this size needs
24 qubits (dim 1.7e7) and the exact fit was not affordable, so only its analytic one-layer
result is reported.

### 4e. The choice

**Direct binary index encoding, k=4, `log2 k = 2` qubits per residue, no penalty.** Justified
by: (i) it is the only encoding that keeps the operating point inside exact simulation
(25.9 qubits mean, 32 max on the 126-target set) and therefore inside certified ground truth;
(ii) the whole register is feasible, so no penalty constant has to be chosen and no
infeasible-state mass is wasted — at k=8 one-hot's feasible fraction is 8.5e-22; (iii) the
move-locality advantage that would justify one-hot or Gray **does not exist** (section 4b);
(iv) its only real deficit — one-layer product-ansatz expressivity — costs 0.065 nats per
residue against a confident prior and is removed entirely by entangling depth >= 2.
Gray and hierarchical coding buy nothing measurable and add a target-dependent failure mode.

---

## 5. HAMILTONIAN DESIGN — what survives

### 5a. The candidates, and their verdicts

| Hamiltonian | decomposes? | verified how | ranks native? | verdict |
|---|---|---|---|---|
| torsion prior only | **yes, exactly 1-local** | ANOVA V1 = 1.000000, residual 2.1e-15; 1-body truncation reproduces it with rho = 1.000 and identical argmin | no (rho +0.011); argmin = the chain-blind snap, 1.921 A at any predictor quality q >= 0.4 | **KEEP** as the generative/warm-start term. It is the only exactly-few-body piece and the only piece with a defensible argmin. |
| torsion prior + Legacy | **no** | best 2-body approximation explains 20.4% of Legacy's variance; its argmin is the 33,836-th best of 262,144 | no: rho +0.054, argmin +0.199 A worse than random | **DO NOT OPTIMISE.** |
| torsion prior + AMBER | **no** | same argument, plus AMBER has no finite variance on this population | no: rho +0.014, argmin +0.545 A worse than random with variance standardisation, -0.112 A with rank normalisation | **DO NOT OPTIMISE.** |
| any weighted combination | — | lambda sweep is monotone with supremum at Legacy alone | no at every lambda | **there is no weight to tune.** |

**Normalisation, derived not tuned.** Each channel is standardised on a native-free
reference population (4,000 configurations drawn from the leakage-safe empirical prior),
with mean and sd fitted as per-residue constants `mu = n*a`, `sd = sqrt(n)*b` pooled across
targets, and combined with unit weights so each contributes equal variance. The fitted
constants are in `qarch_validity.json::normalisation`: Legacy `a = -0.109 kcal/mol/residue`,
`b = 4.022 kcal/mol/sqrt(residue)`; prior `a = 1.146 nats/residue`, `b = 0.618`. **AMBER's
fit is degenerate** (section 3d) and its combination must use the rank/robust form instead.
No weight anywhere was chosen by inspecting a reported number.

### 5b. Route (a) or route (b)?

The brief offered two honest routes: (a) evaluate the energy on sampled bitstrings and
optimise CVaR of the sampled distribution — no decomposition needed; (b) build a genuinely
few-body Hamiltonian and verify it. **Route (b) is closed by measurement** (section 2b): the
L2-optimal 2-body Hamiltonian for Legacy is a different function, not an approximation of
it. **Route (a) is the only option**, and `core.quantum.run_global_cvar_vqe` already
implements it. The prior term is the exception and may be folded in analytically as an exact
sum of single-block operators — which is also what makes the closed-form warm start valid.

### 5c. What I would build now, and what I would not

**PROPOSED, not demonstrated.** The measured structure says the architecture should be:

1. **Encoding: binary index, k=4** (section 4e), 26 qubits mean over the 126-target set.
2. **Generative Hamiltonian: the torsion prior alone**, exactly 1-local, warm-started in
   closed form from the predictor's marginals. Its argmin is the snap (1.921 A); its *value*
   is the distribution it induces, not its minimum.
3. **The many-body energy is not a scoring function here.** Every measurement in section 3
   says Legacy and AMBER cannot rank candidates in this space, so using either as the VQE
   objective makes the answer worse. What they might still legitimately do is act as a
   *filter* (rejecting the 40-68% of configurations that are physically impossible) or as a
   *relaxation* (`refine_coords`, which changes the geometry rather than ranking it) — both
   are different uses from ranking and neither is tested here.
4. **The 0.95 A between the 1-local argmin (1.921 A) and the space best (0.969 A) is the
   prize, and it needs an objective that does not exist yet.** The ANOVA says that objective
   must be at least 3-body (Legacy's variance peaks at order 3-4; the true objective's at
   order 2-3). A trained many-body scorer over torsion configurations, fitted on training
   folds against true RMSD, is the obvious candidate and is NOT tested here.

### 5d. Small-instance ground truth, cached

`s13/results/qarch_enum_<PDB>.npz` for all nine n=9 targets, 262,144 configurations each,
each file carrying:

    rmsd (262144,)          ORACLE post-hoc CA-RMSD, float32
    legacy (262144,)        DEFAULT_WEIGHTS total
    leg_<term> (262144,)    all 11 core.energy terms
    prior (262144,)         1-local empirical-prior energy
    prior_table (9,4)       the leakage-safe per-residue state distribution
    PHI, PSI (9,4)          the library, radians
    snap_index, snap_states ORACLE nearest-state assignment
    amber_idx (~2950,)      indices of the AMBER subsample
    amber_kind              0 = uniform, 1 = prior-sampled, 2 = ORACLE near-native band
    amber_total, amb_<term> genuine ff14SB/GBn2 single points, 5 terms

Configuration index is odometer order, `np.ravel_multi_index(states, (4,)*9)`, so any
optimiser's output can be scored against the certified optimum by a single lookup. Cost:
78 s per target for the full Legacy enumeration + 37 s for the AMBER subsample. Every
optimiser claim on these nine targets can now be checked against ground truth rather than
against another optimiser.

---

## 6. LEGACY vs AMBER — the trainability mechanism

**CORRECTION, recorded rather than silently replaced.** An earlier draft of this section
classified this comparison **WEAK** on the basis of a single completed exact-ANOVA cell
(1CS9 freeze 0: Legacy W_eff 2.238 vs AMBER 3.016) plus a marginal 9-target sampled estimate,
because the AMBER run was blocked on the 1.5 GB memory gate for most of the session. After
the coordinator lifted the compute cap the run completed with **12 matched cells** and the
classification is upgraded to **STRONGLY SUPPORTED**. The full result, together with the
explicit refutation of the naive "AMBER is less local than Legacy" claim, is in
**section 8**; the sampled cross-check is `s13/qarch_sepfit.py`. Summary:
AMBER W_eff 3.054 vs Legacy 2.013 on identical configurations (+1.040, CI [+0.709, +1.359],
11/12 cells), separable share 0.070 vs 0.388 (12/12) — but only 0.595 vs 0.659 after a rank
transform, so most of the raw gap is the van der Waals singularity rather than structure.

---

## 7. SCOPE OF THE EXACT-SUPPORT RULE — it does NOT extend to all-atom distances (PROVEN)

Section 1 is now load-bearing for the Pauli-spectrum and trainability agents, so its scope
has to be pinned rather than assumed.

**Exact protocol for section 1, restated so it can be reproduced.** `qarch_locality.py::
geometry_locality`. Targets 1A13 (n=14), 1CS9 (n=9), 1KZ2 (n=16), library
`torsion_lib2.library_for(seq, 8, exclude_seq=seq)`, `k = 8`, seed 0.
64 base configurations drawn uniformly over the state space. For each residue m and each of
its k library states a, the whole batch is rebuilt with `s_m := a` and the CA trace obtained
from `s12.instrument.build_ca` (= `core.project.build_ca_exact`, the bit-exact ideal-geometry
builder the production projection uses). For each pair (i,j) the statistic is
`median_a median_batch |d_ij(s with s_m := a) − d_ij(s)|`, and residue m is said to support
(i,j) when that exceeds 0.10 Å. Agreement with the rule "support = {i+1, …, j−1}" is
**1.0000** on all three targets: inside-rate 1.000, outside-rate 0.000, over
91 + 36 + 120 = 247 pairs × n variables = 2,993 (pair, variable) cells, plus the identical
result at the earlier threshold sweep, with **zero counterexamples**. Median |Δd| for a
supporting variable is 1.75–1.95 Å; for a non-supporting variable it is **exactly 0.000**.
Support size as a function of separation s is exactly `s − 1` with no scatter at any s.

**Mechanism.** `d_ij` is an internal coordinate of the CA_i…CA_j sub-chain. Under ideal
geometry the CA–CA virtual bond length is a constant (3.80 Å) independent of every torsion,
and the sub-chain's shape is fixed by the virtual bond angles and virtual dihedrals strictly
inside it. Rotating about residue i's torsions rotates the whole downstream chain rigidly
about an axis through CA_i, and rotating about residue j's moves only atoms after CA_j, so
both end residues drop out.

**SCOPE — measured, not assumed (`s13/qarch_allatom.py`, `qarch_allatom.json`).** The same
protocol was rerun on **all-atom pair distances** over the builder's five atom types
(N, CA, C, O, CB), 2,415 atoms-pairs per target. Note this covers AMBER's real atom set for
this encoding: with `chi_bits=False` every side-chain atom is rigid in its residue's local
backbone frame, so its support equals CB's.

| quantity | CA–CA pairs | pairs involving a non-CA atom |
|---|---|---|
| mean support size | 4.00 | **4.73** |
| agreement with the `s−1` rule | **1.0000** | 0.09 |
| exact agreement with the closed interval `i..j` | — | 0.24 |
| mean support at separation s | exactly `s − 1` | `s + 0.04` (0.11, 0.97, 2.04, 3.04, 4.04, …) |
| fraction non-constant at s = 0 (intra-residue) | no such pairs | **0.107** |
| fraction non-constant at s = 1 (adjacent) | **0.000** | **0.753** |

(1A13; 1CS9 and 1KZ2 agree to within 0.01 on every cell — 0.729 and 0.681 at s = 1.)

**The rule does not extend.** All-atom pair supports are on average **one residue wider**
(`s` rather than `s − 1`), they are **not** a single clean interval rule (only 24% match the
closed interval exactly), and — the decisive part — **75% of adjacent-residue atom pairs and
11% of intra-residue atom pairs are non-constant, whereas every CA–CA pair at separation 1 is
exactly constant.** Per atom-pair type at s = 1 the support size is 2 for N–O, N–CB, CB–O,
CB–CB, 1 for most others, and exactly 0 for CA–CA, C–N, O–N, C–CA, O–CA — the covalently
rigid ones.

**So AMBER's support structure genuinely differs from Legacy's**, and in two ways at once:
each pair term reaches one residue further, and AMBER populates the s = 0 and s = 1 shells
that a CA-based objective leaves entirely empty. Cause: φ_i places C_i (and hence CB_i, built
from N_i/CA_i/C_i), ψ_i places O_i and N_{i+1}; only CA_i sits on the rotation axis. This is
a **finding, not a caveat**, and it is what the Pauli-spectrum agent's analytic prediction
must use for AMBER — the CA rule is the Legacy rule.

**Classification: PROVEN** (exact for CA–CA; the all-atom departure is measured with
agreement figures on ~7,000 pair-variable cells per target).

---

## 8. THE NAIVE HYPOTHESIS "AMBER IS LESS LOCAL THAN LEGACY" IS FALSE AS STATED (REFUTED)

Recorded explicitly, because a sprint that silently swaps a falsified hypothesis for a nearby
true one is exactly what this repository's corrections ledger exists to prevent.

**What is false.** Section 1 shows the union of pair supports is the whole chain, so **every**
pairwise molecular energy — Legacy's CA/CB terms and AMBER's ~230-atom nonbonded sum alike —
is **full-register n-local** in the torsion variables. Neither is "more local" than the other
in the support sense; both have support = all n residues. Any statement of the form "AMBER is
less local than Legacy" is a category error about this representation.

**What is true, and measured.** The difference is not *whether* they are local but *how their
variance is distributed across interaction order / Pauli weight*. On **identical
configurations**, an exact ANOVA in a 5-residue subspace (freeze n−5 residues, enumerate
4⁵ = 1,024; 12 matched cells over 1CS9, 6EY3, 7N2I; residual ≤ 6e-11 relative):

| term | V₁ | V₂ | V₃ | V₄ | V₅ | W_eff | graph degree |
|---|---|---|---|---|---|---|---|
| torsion prior | **1.000** | 0.000 | −0.000 | 0.000 | 0.000 | **1.000** | 0.00 |
| `amb_torsion` | **1.000** | 0.000 | −0.000 | 0.000 | −0.000 | **1.000** | 0.00 |
| `amb_angle` | **1.000** | −0.000 | 0.000 | −0.000 | 0.000 | **1.000** | 0.00 |
| ORACLE CA-RMSD | 0.520 | 0.315 | 0.151 | 0.014 | −0.000 | 1.660 | 1.33 |
| `leg_contact` | 0.600 | 0.307 | 0.082 | 0.010 | 0.001 | 1.504 | 1.40 |
| `amb_solvation` | 0.470 | 0.347 | 0.148 | 0.031 | 0.004 | 1.750 | 1.70 |
| **Legacy total** | **0.388** | 0.303 | 0.228 | 0.071 | 0.010 | **2.013** | 1.77 |
| `leg_steric` | 0.266 | 0.348 | 0.272 | 0.098 | 0.016 | 2.248 | 1.87 |
| **AMBER total** | **0.070** | 0.230 | 0.359 | 0.259 | 0.083 | **3.054** | 1.70 |
| `amb_nonbonded` | 0.070 | 0.230 | 0.359 | 0.259 | 0.083 | 3.054 | 1.70 |

* **AMBER's effective weight exceeds Legacy's on 11 of 12 matched cells**, mean +1.040,
  bootstrap CI **[+0.709, +1.359]**.
* **AMBER's separable share is lower on 12 of 12 cells**: V₁ = 0.070 vs Legacy's 0.388.
  AMBER's variance peaks at **order 3–4**; Legacy's at **order 1–2** *in this subspace*
  (in the full 9-variable space Legacy's V₁ falls to 0.049 and its peak moves to order 3–4 —
  freezing variables removes interaction, so subspace and full-space V₁ are not comparable
  across sections, only within them).
* **The gap is mostly, but not entirely, the van der Waals singularity.** Rank-transformed
  (which is invariant to any monotone rescaling), V₁ becomes Legacy 0.659 vs AMBER 0.595,
  AMBER lower on 10 of 12 cells. So the ordering survives a monotone transform but shrinks
  from 0.318 to 0.064 — most of the raw difference is scale, not structure. **This is worth
  stating plainly: on the rank scale the two energies are much more alike than the raw
  numbers suggest.**
* **AMBER's many-body character is entirely `nonbonded`**: its ANOVA profile is identical to
  the total's to three decimals, while `amb_bond` is constant and `amb_angle` / `amb_torsion`
  are **exactly 1-local** (V₁ = 1.000) — because with `chi_bits=False` and ideal geometry
  those terms are functions of a single residue's (φ, ψ).

Independent, no-new-compute cross-check on 9 targets × ~1,194 cached uniform configurations
(`s13/qarch_sepfit.py`): the best additive fit's cross-validated R², calibrated first against
the exact enumeration (Legacy raw 0.022 sampled vs 0.049 exact; CA-RMSD 0.136 vs 0.164, so
the estimator is mildly conservative), gives rank-transformed order-1 shares of **Legacy
0.312 vs AMBER 0.237**, paired difference −0.076, CI [−0.161, −0.002], 5W/4L. Same direction,
weak on its own; the exact subspace ANOVA is the load-bearing measurement.

**Classification: the naive claim is REFUTED; the correct claim — AMBER's variance sits at
higher interaction order than Legacy's on identical configurations — is STRONGLY SUPPORTED
(11/12 and 12/12 matched cells, CI excluding zero), with the important qualification that a
rank transform removes about 80% of the gap.**

---

## 9. PRICING THE ROUGH LANDSCAPE — what the encoding's move set costs (STRONGLY SUPPORTED)

Section 1 says a single state change moves an in-support CA–CA distance by a median
1.75–1.95 Å, so the landscape cannot be smooth in the encoded variables. `s13/qarch_landscape.py`
prices that exactly, on all nine complete 262,144-configuration spaces, by comparing the two
neighbourhoods the encoding actually determines:

* **residue move** (what a one-hot register does with two qubit flips, and what a classical
  residue-wise local search does): change one residue to *any* other state — 27 neighbours;
* **qubit flip** (what a binary index register does with one flip): 18 neighbours, and each
  residue can only reach the 2 of its 3 alternatives whose index differs in one bit;
* **gray**: same size, different reachable pair.

| move set | objective | fraction of configurations that are local minima | steepest descent reaches the global optimum | mean CA-RMSD reached | correlation length (moves) |
|---|---|---|---|---|---|
| residue | Legacy total | 0.00019 | **0.357** | 3.755 | 1.39 |
| **binary** | Legacy total | **0.00099** (5.2×) | **0.135** (2.6× worse) | 3.783 | 1.39 |
| gray | Legacy total | 0.00088 | 0.173 | 3.809 | 1.36 |
| residue | torsion prior | **0.00000** | **1.000** | 3.636 | 6.29 |
| **binary** | torsion prior | **0.00024** | **0.672** | 3.749 | 6.99 |
| gray | torsion prior | 0.00024 | 0.693 | 3.757 | 6.55 |
| residue | ORACLE CA-RMSD | 0.00471 | 0.080 | **1.704** | 1.95 |
| **binary** | ORACLE CA-RMSD | **0.01213** (2.6×) | 0.040 | **1.891** (+0.187 Å) | 1.91 |
| gray | ORACLE CA-RMSD | 0.01888 | 0.035 | 1.908 | 1.84 |

**Four things this prices.**

1. **The binary encoding manufactures local minima in a function that provably has none.**
   The torsion prior is exactly separable (section 2), so under residue moves it has
   **zero** local minima and steepest descent finds the global optimum from **100%** of
   262,144 starts. Under single-qubit flips it acquires 0.00024 local minima and the hit rate
   falls to **67.2%**. That is a pure encoding artefact, measured with certified ground truth.
2. **On Legacy the binary neighbourhood is 5.2× more rugged and 2.6× worse at finding the
   optimum** (13.5% vs 35.7% of starts). Gray coding recovers a little (17.3%) and is worse
   on RMSD; it is not a fix.
3. **The cost in structure is 0.187 Å**: descending directly on true CA-RMSD (the ORACLE
   upper bound on what any local method could reach) gets to 1.704 Å with residue moves and
   1.891 Å with qubit flips, against a space best of 0.969 Å.
4. **The landscape is extremely rugged either way.** Legacy's random-walk correlation length
   is **1.39 moves** — the objective is essentially decorrelated two moves from where you
   started — against 6.3–7.0 for the prior. Single-flip local search is a poor instrument on
   this landscape regardless of the encoding, which is a warning for the classical controls,
   not only for the VQE.

**Consequence for the encoding choice, stated honestly.** This is a *real* cost of binary
that the earlier move-size measurement (section 4b) did not see, and it partly rehabilitates
one-hot: one-hot's move set is measurably better. But (i) the cost is 0.187 Å against an
oracle objective on a landscape where **no native-free objective ranks at all** (section 3),
so it is currently unrecoverable value; and (ii) one-hot costs 2×–4× the qubits, which is the
difference between exact simulation with certified ground truth and no ground truth. **The
recommendation stays binary at k=4, but the trade is now priced rather than asserted, and if
a working objective is ever found the one-hot neighbourhood should be re-examined.**
Note also that a VQE's moves are *not* single-qubit flips — the entangling layers correlate
bits — so this penalty applies with full force to classical single-flip local search and only
by analogy to the variational optimiser.

---

## 10. SCALE, AND A CORRECTION TO MY OWN HEADLINE (this supersedes part of section 3)

The coordinator reported Legacy Spearman rho around +0.3 to +0.4 against my +0.081, and asked
that the disagreement be resolved rather than averaged. It is resolved, in my favour on the
facts and against my framing.

**Config.** All **126 tuning targets**, k in {4, 8}, **12,001 uniform-random configurations
per target** (12,000 draws plus the ORACLE native snap), genuine Legacy with `DEFAULT_WEIGHTS`,
the leakage-safe empirical prior, seed 0, 6 worker processes. `s13/qarch_scale.py`,
`s13/results/qarch_scale.json`.

| k=4, 126 targets | mean rho | median | sd | argmin RMSD | random draw | delta | W/L | native-snap %ile |
|---|---|---|---|---|---|---|---|---|
| **Legacy total** | **+0.190** | +0.232 | 0.237 | 4.354 | 4.971 | **-0.617** | **85/41** | 0.325 |
| Legacy without `steric` | +0.134 | +0.143 | 0.417 | 5.366 | 4.971 | +0.395 | 52/74 | 0.311 |
| `leg_steric` alone | +0.043 | +0.130 | 0.336 | 4.991 | 4.971 | +0.020 | 71/55 | 0.220 |
| `leg_torsion` alone | -0.073 | -0.084 | 0.084 | 4.984 | 4.971 | +0.013 | 66/60 | 0.144 |
| **torsion prior** | +0.022 | +0.002 | 0.116 | 3.856 | 4.971 | **-1.115** | **93/33** | 0.134 |

| k=8, 126 targets | mean rho | median | sd | argmin RMSD | random draw | delta | W/L | native-snap %ile |
|---|---|---|---|---|---|---|---|---|
| Legacy total | +0.169 | +0.223 | 0.187 | 4.235 | 5.013 | **-0.779** | 89/37 | 0.313 |
| torsion prior | +0.031 | +0.027 | 0.070 | 4.019 | 5.013 | **-0.995** | 90/36 | 0.183 |

**Three corrections to what section 3 said.**

1. **My rho was low because my targets were short.** Legacy's rho binned by n at k=4:
   n=9 **+0.077**, n=10 +0.129, n=11 +0.356, n=12 +0.193, n=13 +0.201, n=14 +0.194,
   n=15 +0.196, n=16 +0.148 — mean +0.190, median +0.232. The nine n=9 targets on which
   sections 2, 3, 9 and 11 are built are **the weakest bin on the whole instrument**. The
   coordinator's +0.3 to +0.4 and my +0.081 are the same phenomenon measured on different
   length distributions; my median at k=4 is +0.232, the low end of their range.
2. **"Legacy is worse than random" is FALSE at a realistic optimisation budget.** Best-of-12,001
   beats a random draw by 0.617 A at k=4 (85W/41L) and 0.779 A at k=8 (89W/37L), and the
   effect grows with length: delta by n is -0.027 (n=9), -0.312, -0.386, -0.560, -0.251,
   -0.558, -0.994, **-1.379** (n=16). Section 3's headline sentence was measured at n=9 and
   with a **complete enumeration**, and it does not generalise in that form. **Recorded as a
   correction, not quietly rewritten.**
3. **But the mechanism section 3 identified is real, and section 11 shows it directly.** The
   two facts are reconciled by the optimisation budget, not by one of them being wrong.

**What has NOT changed.** Legacy still puts the near-native configuration at the **32.5th
percentile** of its energy (k=4; 31.3rd at k=8) — a third of random garbage still scores
better than a near-native structure. Its best-of-12,001 at 4.354 A sits **2.76 A above the
k=4 representation ceiling of 1.594 A**, so what it is doing is rejecting the worst
configurations, not locating good ones — the repository's own "garbage rejection reads as
skill" hazard, in its natural habitat. And the **1-local torsion prior beats it**
(-1.115 vs -0.617 A, 93W/33L), which is the same conclusion section 3 reached: the usable
signal in this space is 1-local.

**One clean sub-result.** Removing the steric term destroys Legacy's advantage entirely
(-0.617 A becomes **+0.395 A**, 85W/41L becomes 52W/74L) while the steric term alone does
nothing (+0.020 A, rho +0.043). So Legacy's whole apparent ranking skill in torsion space is
**the steric term vetoing clashes inside an otherwise unhelpful score** — it is a filter, not
a ranker. That is the exact mirror of the coordinator's independent AMBER observation
(AMBER's rho recovers to about +0.25 when its clash term is removed) and the two point the
same way: the clash physics and the ranking physics are different jobs and should not be
carried by one number.

---

## 11. "OPTIMISE HARDER, GET WORSE" AS A MEASURED CURVE (PROVEN on the enumerated family)

The reconciliation of section 3 and section 10 is the **objective-evaluation budget**, and on
the enumerated spaces both endpoints are exact, so the whole curve can be drawn.
`s13/qarch_budget.py`: for m = 10 … 262,144 draw m configurations uniformly without
replacement, take the best by each objective (averaging over the tied-argmin set), record the
CA-RMSD; 40 independent draws per budget, and the last point is the certified global optimum.
Nine n=9 targets, k=4.

| objective-evaluation budget | 10 | 30 | 100 | 300 | 1,000 | 3,000 | 10,000 | 30,000 | 100,000 | 262,144 (exact optimum) |
|---|---|---|---|---|---|---|---|---|---|---|
| **Legacy total** | 3.764 | 3.708 | 3.669 | **3.667** | 3.740 | 3.781 | 3.864 | 3.900 | 3.914 | **3.920** |
| torsion prior | 3.859 | 3.772 | 3.922 | 3.903 | 3.941 | 3.926 | 3.881 | 3.850 | 3.822 | 3.636 |
| `leg_steric` | 3.683 | 3.703 | 3.693 | 3.694 | 3.692 | 3.690 | 3.690 | 3.691 | 3.691 | 3.691 |
| ORACLE CA-RMSD | 2.548 | 2.166 | 1.836 | 1.549 | 1.348 | 1.154 | 1.040 | 0.979 | 0.969 | **0.969** |

(random draw 3.781 A; space best 0.969 A.)

**Legacy's curve is non-monotone.** It improves to a minimum of 3.667 A at ~300 evaluations —
0.114 A better than random — and then degrades monotonically for the next three orders of
magnitude, ending at 3.920 A, **0.253 A worse than its own best-of-300 and 0.139 A worse than
doing nothing.** The true objective's curve, by contrast, falls monotonically to the space
best. `leg_steric` is flat throughout: it neither helps nor hurts, consistent with section 10's
finding that it is a veto rather than a ranker.

**This is Sprint 12's "the certified exact optimum emitted 0.169 A worse than its own anchor"
reproduced as a continuous curve rather than a single point, on a different architecture and
a different Hamiltonian.** It has three consequences that matter more than any of the
individual numbers above:

1. **The optimisation budget is a parameter of the ANSWER, not of the search.** Reporting
   "Legacy ranks" or "Legacy anti-ranks" without stating the budget is under-specified; both
   my section 3 and the coordinator's measurement are correct at their own budgets.
2. **Any comparison of optimisers on this objective is confounded.** A *better* optimiser
   finds a *worse* structure. Budget parity is not merely a fairness condition here — without
   it, the more effective method loses, and the sprint's Legacy-vs-AMBER trainability study
   must therefore compare **optimisation quality (objective value reached)** and **structural
   quality (CA-RMSD)** as two separate axes and never conflate them.
3. **A weakly-optimised bad objective is a regulariser, not a result.** The 0.114 A that
   best-of-300 Legacy buys is the value of *stopping early*, and it is smaller than what the
   1-local prior gives for free.

### 11b. The same curve on all 126 targets — and the turn is length-dependent

`s13/qarch_budget2.py`: 126 targets, k=4, **60,000 uniform draws each**, best-of-m read off
the same sample, 24 repeats per interior budget. Aggregated over all 126 the curve is
4.737 (m=10) → 4.601 (300) → **4.558 (1,000)** → 4.585 (3,000) → 4.535 (60,000): a shallow dip
and then flat. **Split by length it is not flat at all, and the splitting variable is the
FRACTION OF THE SPACE explored, not the evaluation count.**

| n | targets | fraction of 4^n seen at 60,000 | best budget | best RMSD | RMSD at 60,000 | degradation |
|---|---|---|---|---|---|---|
| 9 | 9 | 2.3e-01 | 300 | 3.688 | 3.853 | **+0.165** |
| 10 | 11 | 5.7e-02 | 1,000 | 3.999 | 4.346 | **+0.347** |
| 11 | 13 | 1.4e-02 | 30 | 4.263 | 4.676 | **+0.413** |
| 12 | 19 | 3.6e-03 | 60,000 | 3.956 | 3.956 | +0.000 |
| 13 | 23 | 8.9e-04 | 30,000 | 4.361 | 4.389 | +0.028 |
| 14 | 14 | 2.2e-04 | 60,000 | 4.861 | 4.861 | +0.000 |
| 15 | 16 | 5.6e-05 | 60,000 | 4.877 | 4.877 | +0.000 |
| 16 | 21 | 1.4e-05 | 60,000 | 5.046 | 5.046 | +0.000 |

| group | n | curve (m = 10 … 60,000) | best | at 60,000 | degradation | random draw |
|---|---|---|---|---|---|---|
| **n ≤ 11** (≥1.4% of the space seen) | 33 | 4.118 4.085 4.089 **4.075** 4.095 4.144 4.176 4.270 4.342 | 4.075 @ 300 | 4.342 | **+0.267** | 4.280 |
| **n ≥ 12** (≤0.36% seen) | 93 | 4.957 4.955 4.878 4.787 4.722 4.742 4.722 4.678 **4.604** | 4.604 @ 60,000 | 4.604 | +0.000 | 5.215 |

**Reading, and it is the most useful sentence in this document.** The degradation is not a
property of short peptides — it is a property of *how much of the space you have searched*.
It appears once roughly **1% of the configuration space** has been seen, which happens at
n ≤ 11 and does not happen at n ≥ 12 within any budget tested. On the 33 short targets the
best-of-60,000 (4.342 Å) is already **worse than a random draw** (4.280 Å); on the 93 longer
ones the curve has simply not turned yet.

**Consequence for the sprint.** A VQE at n = 12–16 with 10^4 evaluations explores ~10^-4 of
the space and will sit on the monotone-improving part of the curve. **It will look like it is
working. The enumerations prove that its objective's optimum is in the wrong place, and no
experiment confined to that budget can see it.** Any claim of the form "the VQE improved the
structure" on long targets is therefore a claim about the early part of a curve whose limit is
known to be bad, and it should be reported that way.

**Classification: PROVEN on the nine enumerated targets** (exact endpoints, 40 draws per
interior budget) and **STRONGLY SUPPORTED on the 33 targets with n ≤ 11** (sampled, 126-target
sweep, monotone degradation over three budget decades). **On n ≥ 12 the question is OPEN and
untestable at any budget reached here** — the space is too large to search far enough.

---

## 12. LEAKAGE AUDIT

* The torsion library is `torsion_lib2.library_for(seq, k, exclude_seq=seq)`, which builds
  from `peptide_db.holdout(seq)` — the target and everything above 0.6 identity to it is
  removed. The empirical prior is built from the same holdout pools.
* Native information is read in exactly three places, all labelled: `Space.nat_ca` and
  `Space.rmsd` (post-hoc evaluation only), `Space.ORACLE_native_torsions` /
  `ORACLE_snap` (diagnostics), and `ORACLE_prior(q)` (an explicitly named oracle arm
  standing in for the predictor agent's output). No native quantity enters any objective
  that is described as native-free, and no weight, threshold or state was chosen by looking
  at an RMSD.
* benchmark60 was not touched: no read of `results/benchmark_manifest.json`, no read of
  `s9/final_cache/*`. dev24 was not run. All nine targets are from the 126-target tuning
  instrument via `s12.instrument.targets()`.
* Nothing under `core/`, `s5/ s7/ s8/ s9/ s12/`, `tests/`, `verify/` was modified. No git
  writes. All new code is `s13/qarch_*.py`; all results are `s13/results/qarch_*`.
* Compute: one heavy process at a time, `OMP_NUM_THREADS=2`, `s12.instrument.free_gb()`
  checked before every AMBER stage with a 1.5 GB gate (which blocked the AMBER subspace
  ANOVA for the last part of the session, as recorded in section 6).

## 13. WHAT I DID NOT ESTABLISH

1. **No VQE was run.** By design: the brief says do not optimise an objective that does not
   rank, and none of them ranks. Budget-parity numbers, CVaR arms and optimiser comparisons
   belong to the trainability agent; the cached enumerations exist so those can be scored
   against certified optima rather than against each other. Section 11 says such a comparison
   must report objective quality and structural quality separately or it is confounded.
2. **Certified optima exist only at n=9, k=4.** Sections 2, 3, 9 and 11 all rest on the nine
   completely enumerated spaces, which section 10 shows is the *weakest* length bin of the
   instrument. The 126-target results (sections 10, and 11's sampled version) are
   sampling-based, so their "argmin" is a best-of-m, not a certified optimum. Whether the
   non-monotone budget curve keeps its shape at n=16 is being measured
   (`s13/qarch_budget2.py`) and is not yet complete at the time of writing.
3. **AMBER at scale is untested.** All AMBER results are on nine n=9 targets and ~1,200-3,000
   configurations per population. AMBER was never run on the 126-target instrument.
4. **AMBER *with* minimisation was never tested as an objective.** Only `single_point`. A
   restrained minimisation changes the geometry as well as the score and could plausibly rank
   where the single point does not.
5. **The multimodal ansatz fit is one seed per cell** (section 4d), so "depth 3-4 suffices" is
   OPEN, not SUPPORTED, and the non-monotonicity in depth is unexplained.
6. **No trained many-body objective was built.** Section 5c argues one is what is missing;
   that is a PROPOSAL with no evidence behind it.
7. **The predictor's real per-residue distribution does not exist yet**, so every prior result
   here is either the weak leakage-safe empirical prior or the ORACLE quality knob q.
8. **`leg_steric`'s in-band argmin effect** (-0.522 A, 7W/2L inside the ORACLE band, section
   3c) is unexplained and was not chased; its in-band rho is -0.045, so it is an argmin
   artefact until shown otherwise.
9. **I did not test whether early stopping is a usable strategy.** Section 11 shows
   best-of-300 Legacy is the sweet spot on nine 9-mers; that number is not a recommendation,
   it is a diagnosis, and choosing a budget by looking at RMSD would be exactly the tuning
   this brief forbids.
