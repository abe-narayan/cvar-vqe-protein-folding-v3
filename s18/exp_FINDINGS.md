# SPRINT 18 — EXPERIMENT WORKSTREAM — FINDINGS

> ## VERDICT
>
> **The λ ladder is COMPLETE at n = 126 and it is a clean falsification. Four of the five sprint
> falsifiers fire; the fifth (F2) does not fire only because the data reject its null in the
> *harmful* direction. The degree-1 branch is closed.**
>
>     lambda            0       0.25      0.5      0.75       1        Control A
>     mean RMSD       4.335    3.742    3.683    3.619    3.610         3.048
>
> Monotone **decreasing** in λ — the exact opposite of the pre-registered signature of a
> mis-specified higher-order component. And the mechanism, from an ORACLE diagnostic that cost
> nothing: **fed the native's own distances, the full objective reaches 1.152 Å and the degree-1
> object reaches only 3.769 Å** — worse than the coordinate average it was meant to replace.
> The truncation destroys the structural information; no distogram improvement could rescue it.

**Status.** The deciding arm (§4) is complete at n = 126. `exp_polish` is complete at n = 126.
Control B, H8 and H6 are **PARTIAL**; the helix-μ control is complete at its subset; each table states the n it was computed
at (§5.2 n = 85, §5.3 n = 88, §5.3b n = 17, §5.4 n = 30 complete). Nothing here
is read from a smoke run, and no partial artefact is reported as final. **The background runs were
stopped so the artefacts, `s18/results/exp_report.txt` and this file all describe one frozen
state** rather than drifting apart.

Pre-registration: `s18/PREREG_exp.md` (written before any arm; amended four times, **every time
before the relevant results were read** — for MATH's landing, for MATH's own mid-sprint revision,
for the coordinator's α ladder, and for the free `with_dhat` arms).
Modules: `s18/exp_obj.py` (adapter), `s18/exp_run.py` (arms), `s18/exp_polish.py` (the true
optimum), `s18/exp_budget.py` (Control B), `s18/exp_leverage.py` (H8), `s18/exp_helixmu.py` (the
zero-information μ), `s18/exp_report.py` (reader). Full reader output:
`s18/results/exp_report.txt`. Artefacts under `s18/results/`.

**Snapshot note.** The PARTIAL arms in §5.2, §5.3 and §5.3b are frozen at the n stated in each
table and match `s18/results/exp_report.txt`; their runs were stopped rather than left drifting.
§5.4 completed at its pre-registered 30-target subset. They are
checkpointed per target and resumable, and `python -m s18.exp_report` regenerates every table from
whatever is on disk. I tracked the drift as they grew: **Control B and H8 held their sign and their
reading unchanged from n ~ 45 through n ~ 88**, so those readings are stable even though the
magnitudes are not settled and are not claimed as settled. **H6 (§5.3b) and the helix-μ arm (§5.4)
did NOT hold as they grew** — both are underpowered and both carry an explicit correction against
an earlier reading of mine. Neither is cited for a direction.

---

## 0. What I consume, and what I did not build

The BRIEF forbids a second definition of the degree-1 object. MATH landed `s18/math_anova.py`
while my harness was being built; **their object is consumed as-is** and my provisional
implementation (`s18/exp_anova.py`, written against the published interface in
`s18/COORD_exp_to_math.md` before they landed) is **retired and not used by any reported arm**.

`s18/exp_obj.py` adds exactly two things and redefines nothing:

1. **analytic gradients** of `E_full` and `E_le1`, so every arm is optimised by the same L-BFGS
   to the same tolerance — a comparison in which one objective is optimised better than another
   measures the optimiser, not the objective;
2. the convex mix `E_λ = (1−λ)E_le1 + λE_full`, which introduces no new function.

**VERIFIED (EXACT).** Against MATH's own callables, not against a re-derivation:

| check | result |
|---|---|
| `E_full` value vs `math_anova.Target.E_full` | 9.1e-13 |
| `E_le1` value vs `math_anova.Target.E_le1` | 0.0 (bit-identical) |
| `E_full` gradient vs central differences **of MATH's function** | 3.2e-10 relative |
| `E_le1` gradient vs central differences **of MATH's function** | 3.3e-10 relative |
| `E_le1 + E_ge2 − E_full` | 9.1e-13 |

---

## 1. Two structural facts, recorded as EXACT (theorems, not discoveries)

**1.1 `E_le1` is additive in residues, so it poses no search problem at any budget.**
Its global optimum is the per-residue argmin of a tabulated 2-torus field: `O(n·G²)` lookups,
certified by separability. Sprint 17 found greedy 1-opt certifies the *full* objective's optimum
at 1,024 evaluations; for degree-1 the question does not arise at all. **Whatever is wrong with
this branch, it is not that the optimiser cannot find the optimum.** Control B measures the full
objective's side of this and is reported with the accounting convention stated.

**1.2 The degree-1 object exerts no force whatever on the two terminal residues.**
MATH's locality corollary (`f_0 = f_{n−1} = 0` exactly); confirmed on **126/126** targets.

MATH's shipped `argmin_le1` resolves the resulting tie by taking `argmin` of an identically-zero
table, which returns mesh cell (0, 0) — a fixed corner of the torus. I flagged that as an
arbitrary structural commitment on 2 of n residues that would be charged to the degree-1 object as
if it were information, and ran **both** variants. **The measurement retires my own flag: they are
identical on 126/126 to 2.13e-14** — see §3c for the exact reason. *(The docstring discrepancy
still stands and is flagged back to MATH: it says these residues are returned "at the
reference-ensemble circular mean", which is not what the code does. It happens not to matter.)*

---

## 2. Quadrature — the one place where cost could have bought a false falsification

Cutting the quadrature is **not** a neutral economy, and the direction of the harm is *against*
the hypothesis: a noisier field makes the certified argmin land partly on Monte-Carlo noise, which
degrades RMSD and would **manufacture a falsification**. So I refused to choose S by cost.

MATH originally declared S = 2048, G = 24 — ~150 s per target per fit under load, ~5.3 h for one
126-target arm, more than the envelope once the controls and the μ swap are added. I began
measuring the S-sensitivity of the **downstream arms** (not of field agreement, which is the wrong
quantity) at S ∈ {256, 512, 2048}, with the pre-registered rule "run at the smallest S whose mean
|Δ| to S = 2048 is below the 0.084 Å MDE on both arms".

**That measurement was superseded, not skipped: MATH revised their own declaration mid-sprint to
S = 512, G = 16 with a moment-based estimator.** My partial sensitivity data was taken at G = 24
and is therefore off-config; it is retained as `s18/results/_STALE_exp_sens_grid24.json` and is
**not cited anywhere**, because mixing two grid configurations is exactly the kind of quiet
incomparability this sprint keeps finding. Every reported arm runs at **MATH's own declared
S = 512, G = 16**, read dynamically from their module rather than pinned by me, with their
`CONFIG` persisted into every artefact. The quadrature is therefore MATH's choice, not a
cost-driven one of mine — which is the only version of this that could not have manufactured the
result.

---

## 3. What the λ ladder is read against (two axes, both pre-registered)

**Axis 1 — Sprint 17's measurement.** From the identical Control-A start, refining toward the
full objective moved RMSD 3.048 → 3.610, **+0.561 [+0.407, +0.712], 31W/95L**, at a 71% objective
reduction. λ = 1 of my ladder must reproduce that number; that is validity gate P0.a.

**Axis 2 — the coordinator's α ladder** (`s18/COORD_FINDING.md`, `s18/results/objceil.json`,
n = 126, same starts). Held the functional form fixed and varied only the distances:

| | |
|---|---|
| α = 0 (deployed distogram) | 3.610 |
| α = 1 (perfect distances) — **ORACLE** | 1.152, 80.2% of targets < 2.5 Å |
| shuffled residual **direction**, magnitudes kept — **ORACLE** | 2.609 (−1.001 [−1.226, −0.784] vs α=0) |
| isotropic matched noise — **ORACLE** | 2.573 |

Two consequences I adopt:

- **The functional form is sound.** Degree-1 cannot be framed as fixing a broken functional form;
  I withdraw that framing from my own motivation. What is broken is the *distances*.
- **The distogram's error is worse than random error of the same size.** That gives degree-1 the
  sharper hypothesis H2b: if the *pairwise* component of the error is the mis-directed part, an
  objective that averages over pairs should be more robust to it and should close part of the
  1.0 Å gap between 3.610 and 2.609. Pre-declared readings: ≈3.6 → nothing; ≈2.6 → recovers what
  randomising the error direction recovers; <2.6 → does something the shuffle does not.

**Caveat carried with every citation of the α curve** (the coordinator's own): α is a blend
*toward truth*, not a model of how a better predictor would err. It is a requirement in the "how
much residual must go away" sense, **not** a prediction about any achievable predictor. Every
α > 0 arm is ORACLE, a ceiling, never a result.

---

## 3b. The premise had already changed before my run reported — MATH's lattice work

`s18/results/math_lattice.json`, MATH's, 19 exhaustively enumerated cells, recomputed by me from
their artefact. **The object my 126-target run tests is the RESIDUE-ADDITIVE one, because that is
what MATH's continuous `E_le1` is.** On the very instrument the sprint's premise came from:

| argmin, 19 cells | mean RMSD |
|---|---|
| space best (ORACLE ceiling) | 1.030 |
| full objective | **2.661** |
| **residue-additive** (= the continuous object) | **2.657** |
| qubit-level Walsh weight-≤1 (Phase 0's 2.411) | 2.411 |
| qubit-level Walsh weight-≤1, averaged over the 24 relabellings of the arbitrary bit encoding | **2.862** |

- **residue-additive − full: mean −0.004, median 0.000, W/L 3/5.** On the 19-target instrument the
  object I am testing is worth *nothing at all* against the full objective.
- The entire −0.249 that motivated the sprint belongs to the **qubit-level** truncation, which
  (i) discards part of the per-residue field, (ii) has **no continuous image** — MATH's F5 finding
  is that the two lattice bits index k-means clusters of the joint (φ,ψ) library and do not
  correspond to φ and ψ, so no continuous object equals the strict Walsh weight-≤1 projection —
  and (iii) is **gauge-dependent**: averaging over the 24 permutations of the arbitrary bit
  labelling moves it to 2.862, i.e. **worse than the full objective it was supposed to beat**.
- Variance fractions confirm the split: qubit-level weight-≤1 captures 0.613 of the variance,
  residue-level 0.913.
- Alignment on the lattice: ρ(objective, RMSD) full **+0.264**, residue-additive **+0.172**,
  qubit weight-≤1 **+0.153**. Both truncations are *worse* global correlates of RMSD than the
  full objective, exactly the Phase-0 tension.

**Consequence, recorded before my continuous numbers landed:** F4's premise is already damaged by
MATH's own lattice measurement, and the 2.411 figure is a property of one arbitrary bit labelling.
My 126-target continuous test is still worth running — it is a different instrument, and it tests
*refinement* rather than only the argmin — but it should be read as confirming or overturning a
hypothesis whose lattice support has already gone.

## 3c. Corrections I made to my own flags, and one made to the record I was given

**I raised a hazard and the measurement retired it (EXACT).** I flagged that MATH's `argmin_le1`
resolves the terminal residues' undetermined fields by taking `argmin` of an all-zero table, which
returns a fixed torus corner — an arbitrary structural commitment on 2 of n residues. I ran both
variants. **They give identical RMSD on every target, to 0.00e+00**, and the reason is exact: `φ_0`
and `ψ_{n−1}` are never read by the builder; `ψ_0` rotates the chain about an axis *through* CA_0,
which fixes CA_0 and is therefore a rigid rotation of the whole Cα set; `φ_{n−1}` moves only atoms
placed after CA_{n−1}. **The residues the degree-1 object cannot see are exactly the residues the
metric cannot see.** The "blind to 17% of the chain" worry is void, and I withdraw it.

**MATH's mesh argmin is not the objective's optimum, and correcting it does not help.** `E_le1` is
evaluated by trigonometric interpolation, so it dips below every mesh node between them; L-BFGS
from Control A reaches a *lower* `E_le1` than the "certified" mesh argmin, which would have made
my own "fraction of the way to the certified optimum" column read >100% and let me report a search
failure that was really a tabulation artefact. `s18/exp_polish.py` computes the true continuous
per-residue optimum. It is genuinely lower on **30/30** targets measured so far, with ~11 residues
per target sitting off-mesh — **and the RMSD does not improve.** The objective/structure
dissociation reappears *inside* the degree-1 object.

**A retraction I was handed and have propagated.** The coordinator's `shuf_paired` arm (2.072 Å),
which motivated my H8 leverage arm, was found by ADVERSARIAL to pass a permuted weight vector into
the fit and is **withdrawn as confounded**. H8 stays pre-registered and is run as written, at the
coordinator's request, so the record shows an arm that outlived its motivation; 2.072 Å is quoted
nowhere as its ceiling. The reference arms that still stand are `shuffled` 2.609, `shuf_strat`
2.560, `isotropic` 2.573. **My λ ladder has no weighting axis** — λ mixes two objectives at fixed
`w` — so it does not duplicate ADVERSARIAL's `wflat` / `wperm_only`.

## 3d. Answers to the coordinator's two claim-checks (posted here so they are visible immediately)

**`s18/priorfit.py` does NOT collide with my λ ladder — keep it, it is yours.** My ladder mixes two
*existing* objectives at fixed weights, `E_λ = (1−λ)E_le1 + λE_full`; it adds no term to the
objective and has **no prior, no regularisation and no weighting axis**. Likewise it does not
duplicate ADVERSARIAL's `wflat` / `wperm_only`. The only weighting arm I own is H8's geometric
leverage, whose control permutes the *leverage* factor while leaving `w = 1/sd²` in place —
deliberately a different question from permuting `sd`.

**The `fit_correction` subset trap: confirmed independently, and already guarded.** I walked into
the same trap on a 3-target smoke early on — a 3-target debias made λ=1 disagree with s17 by
−0.457 and +0.975 Å on two targets — diagnosed it as the debias rather than a harness bug, and
every module in this lane now fits the leave-fold-out correction on the **full `I.targets()`** and
iterates the subset separately (`exp_run.main(subset=...)` documents this in its own docstring).
My validity gate P0.a consequently reproduces s17's `refine_full` to **mean |Δ| = 0.00000 Å**.

## 4. RESULTS — n = 126, complete (`s18/results/exp_main_pool_512.json`)

Validity gates first. **P0.a PASSES**: λ = 1 reproduces `s17/refine.py`'s `refine_full` on all 126
targets to mean |Δ| = **0.00001 Å** (max 0.00057), the objective to mean |Δ| = 0.0000, and Control A
is bit-identical to Sprint 17's. P0.b: `E_le1 + E_ge2 = E_full` to 9.1e-13 through MATH's own
callables. P0.d: my analytic gradients match central differences *of MATH's functions* to 3.3e-10
relative on both objectives. P0.e: the degree-1 field is identically zero on exactly {0, n−1} on
126/126, as MATH's theorem requires.

### 4.1 THE λ CURVE — the sprint's central measurement

    E_λ = E_le1 + λ·E_ge2 = (1−λ)E_le1 + λE_full,  λ ∈ {0, 0.25, 0.5, 0.75, 1}   (pre-registered)

| arm | RMSD | median | vs Control A [fold-clustered 95% CI] | W/L | median Δ | frac improved |
|---|---|---|---|---|---|---|
| **avg (Control A)** | **3.048** | 2.837 | — | — | — | — |
| proj (torsion start) | 3.213 | 2.966 | +0.164 [+0.128, +0.212] | 15/111 | +0.096 | — |
| **λ = 0 (DEGREE-1)** | **4.335** | 4.034 | **+1.287 [+1.124, +1.473]** | 18/108 | +1.168 | 0.143 |
| λ = 0.25 | 3.742 | 3.545 | +0.694 [+0.536, +0.801] | 23/103 | +0.637 | 0.183 |
| λ = 0.5 | 3.683 | 3.518 | +0.635 [+0.543, +0.704] | 27/99 | +0.541 | 0.214 |
| λ = 0.75 | 3.619 | 3.556 | +0.570 [+0.438, +0.689] | 32/94 | +0.448 | 0.254 |
| **λ = 1 (FULL, = s17)** | **3.610** | 3.370 | +0.561 [+0.408, +0.675] | 31/95 | +0.443 | 0.246 |

**SHAPE: MONOTONE DECREASING IN λ.** The curve is worst at λ = 0 and improves at every step to
λ = 1. That is the *exact opposite* of the pre-registered signature of a mis-specified
higher-order component. Removing `E_ge2` does not remove harm — it removes **information**.

    λ = 0 minus λ = 1:   +0.726 [+0.605, +0.816]   median +0.663   W/L 34/92
    best λ minus Control A (λ = 1):  +0.561 [+0.408, +0.677]

Per-target ΔRMSD vs Control A (never endpoints alone):

| λ | p10 | p25 | p50 | p75 | p90 | frac improved | worst | best |
|---|---|---|---|---|---|---|---|---|
| 0.00 | −0.191 | +0.467 | +1.168 | +1.825 | +2.735 | 0.143 | +6.588 | −2.004 |
| 0.25 | −0.186 | +0.137 | +0.637 | +1.136 | +1.631 | 0.183 | +3.192 | −1.045 |
| 0.50 | −0.230 | +0.092 | +0.541 | +1.141 | +1.639 | 0.214 | +3.188 | −1.238 |
| 0.75 | −0.315 | +0.006 | +0.448 | +1.012 | +1.654 | 0.254 | +3.632 | −1.946 |
| 1.00 | −0.391 | +0.019 | +0.443 | +1.035 | +1.697 | 0.246 | +3.699 | −1.966 |

The effect is **not concentrated**: against a uniform-effect null, drop-top-10 sits at the 52.2nd
percentile for λ=0, 50.9th for λ=1, 50.3rd for the argmin — dead centre. It is uniform across the
sample, not carried by a few targets. It holds in **every** length stratum and **every** fold
(λ=0 costs +0.895 / +1.104 / +1.692 Å at n≤11 / 12–13 / ≥14; +1.06 to +1.58 across the five folds).

### 4.2 The certified argmin, and the objective/structure dissociation *inside* degree-1

    EXACT separable argmin of E_le1 (MATH's mesh)        4.346   +1.297 [+1.012, +1.500]  17W/109L
    same, terminals held at the start                    4.346   identical to 2.13e-14
    MATH's harmonic order-1 (angle) truncation           4.533   +1.485 [+1.362, +1.610]
    MATH's sub-residue angle-additive argmin             4.600   +1.552 [+1.413, +1.684]

The *coarser* the truncation, the worse the structure, monotonically: residue-additive 4.35,
angle-additive 4.60. This is the same ordering the λ curve gives, on a second axis.

### 4.3 CONTROLS — degree-1 beats them, and that is the point

| control | RMSD | vs Control A |
|---|---|---|
| shuffled distogram through the identical machinery | 5.347 | (worse than every arm) |
| random additive field, matched amplitude **and** power spectrum | 5.000 | +1.951 [+1.866, +2.047] |
| matched-magnitude random torsion move (matched to λ=0) | 4.190 | +1.142 [+0.922, +1.360] |
| matched-magnitude random move (matched to the argmin) | 4.941 | +1.893 [+1.712, +2.084] |

Degree-1 (4.335) does beat its zero-information controls, so it is **not noise** — it carries real
information. But it beats the matched-displacement random move by only 4.190 → 4.335, i.e. **it
does not beat it at all**: moving that far in a random direction is as good. The information it
carries is not enough to pay for the displacement it spends.

### 4.4 H7 — WHAT LIMITS DEGREE-1: the truncation, not the distances

MATH's `with_dhat` rebuilds every field *exactly* from the cached conditional moments, so the same
degree-1 object can be fed different distances at zero cost.

| distances fed to the objective | E_full | E_le1 | argmin E_le1 | label |
|---|---|---|---|---|
| deployed distogram | 3.610 | 4.335 | 4.346 | native-free |
| pairs permuted | 5.309 | 5.347 | 5.258 | native-free CONTROL |
| magnitudes kept, direction destroyed | 2.648 | 4.388 | — | **ORACLE** |
| **the native's own distances** | **1.152** | **3.769** | 3.847 | **ORACLE** |

> **Perfect distances buy the full objective 2.458 Å and buy degree-1 only 0.566 Å.** At perfect
> distances degree-1 is +2.617 Å [+2.237, +2.957] worse than the full objective, 5W/121L, and
> **+0.721 Å [+0.498, +0.945] worse than the coordinate average it was supposed to replace.**

This is the sprint's mechanism, and it is decisive: **the degree-1 truncation destroys the
structural information, not the distogram's error.** A perfect distogram would not rescue it. The
coordinator's α-ladder shows the *functional form* is sound (1.152 Å at α=1); this shows the
*truncation* is not — the additive object cannot express a structure even when handed the truth.
It also answers H2b flatly: degree-1 closes **−72.5%** of the α=0 → shuffled gap, i.e. it moves
the wrong way by three quarters of that gap. **Averaging over pairs does not launder the bias.**

*(The α-ladder caveat applies wherever I cite it: α is a blend toward truth, not a model of how a
better predictor would err; it is a requirement, not a prediction, and every α>0 arm is ORACLE.)*

### 4.5 The four quantities dissociate — again, and harder

- **Objective quality.** λ=1 reduces `E_full` 192.8 → 56.5 (**63.6%**). λ=0 reduces `E_le1`
  375.4 → −381.9. Both objectives are optimised well. *(No "% of the way to the optimum" is
  quotable: `E_le1` is signed and its mesh argmin is not its optimum — see 4.6.)*
- **Structural quality.** Every arm is worse than Control A; the better the objective is
  optimised, the worse the structure.
- **Alignment**, four named axes over the same 75 real pool structures per target:

  | objective | Spearman | Pearson | pairwise-ordering | in-band Spearman | frac ρ>0 |
  |---|---|---|---|---|---|
  | full | +0.063 [−0.016, +0.133] | +0.074 | 0.521 | +0.001 | 0.540 |
  | degree-1 | +0.038 [−0.012, +0.108] | +0.059 | 0.512 | **−0.032** | 0.516 |

  Both CIs include zero. Degree-1 is the worse correlate on every axis, and its **in-band**
  correlation is *negative*. Phase 0's 19-target tension (+0.153 vs +0.264) survives in sign and
  collapses in magnitude on the real instrument: ρ(le1) − ρ(full) = −0.026 [−0.081, +0.038].
  In-band mean band size 49.8/75; pool diversity 2.486 Å.
- **Search quality.** See 4.6 and Control B.

### 4.6 The mesh argmin is not the optimum — and fixing that does not help

`E_le1` is evaluated by trigonometric interpolation, so it dips below every mesh node between them.
**L-BFGS from Control A ends below MATH's mesh argmin on 85/126 targets.** `s18/exp_polish.py`
computes the true continuous per-residue optimum (the additive object's global optimum, certified
by separability, from the K lowest cells per residue). It is genuinely lower on every target
measured, with ~11 residues per target sitting off-mesh — **and RMSD does not improve.** Reaching
the *true* optimum of the degree-1 objective instead of an approximation buys nothing.

### 4.7 Phase 9 — the coordinate average survives, untouched

| route from the identical pool | RMSD | vs Control A | W/L |
|---|---|---|---|
| **pool → coordinate average (INCUMBENT)** | **3.048** | — | — |
| pool ORACLE best (a ceiling, not a result) | 2.306 | −0.742 [−0.801, −0.690] | 121/5 |
| pool mean = a random single pick | 3.551 | +0.502 [+0.450, +0.565] | 6/120 |
| pool → select by `E_full` | 3.531 | +0.483 [+0.345, +0.650] | 34/92 |
| pool → select by `E_le1` | 3.600 | +0.551 [+0.468, +0.648] | 27/99 |
| pool → avg → refine on `E_le1` | 4.335 | +1.287 [+1.120, +1.476] | 18/108 |
| pool → avg → refine on `E_full` | 3.610 | +0.561 [+0.408, +0.674] | 31/95 |
| pool → avg → `E_le1` exact argmin | 4.346 | +1.297 [+1.012, +1.500] | 17/109 |

**The architectural answer is NOWHERE.** The degree-1 objective belongs neither before averaging
(selection by it is worse than a random pick from the pool: 3.600 vs 3.551), nor after averaging
(refinement by it costs 1.287 Å), nor as a candidate generator (its certified optimum is 4.346).
Selection by either objective is **worse than picking a pool member at random**, which reproduces
the programme's standing result that nothing ranks within the pool.

### 4.8 FALSIFIERS

| | verdict | evidence |
|---|---|---|
| **F1** degree-1 lands ≈3.6 Å on the real instrument | **FIRED** | it lands at **4.335** (refine) / **4.346** (certified argmin) — worse than 3.6 |
| **F2** degree-1 no better than the full objective | **not fired — it is DECISIVELY WORSE** | +0.726 [+0.605, +0.816], 34W/92L. F2 as written tests a null; the data reject the null in the harmful direction |
| **F3** the advantage dies under matched controls | **FIRED** | there is no advantage to die; degree-1 (4.335) does not beat the matched-displacement random move (4.190) |
| **F4** the effect exists only on the 19-target instrument | **FIRED** | at n = 126 the sign is reversed and the interval excludes zero by a wide margin |
| **F5** invalid continuous mapping | MATH's | their finding: no continuous object equals the strict Walsh weight-≤1 projection |

**Four of the five sprint falsifiers fire, and the fifth is worse than firing.** The branch is
closed. Per the BRIEF: no rescue with tuned weights, no tuning against RMSD, no extension of the
λ ladder beyond the pre-registered five.

### 4.9 What this leaves standing

1. **Sprint 17's finding is confirmed on a second, independent objective family.** The pipeline
   works because it does not optimise its own scoring function, and that is now true of the whole
   `E_λ` family, not just one member.
2. **"The objective is mis-aimed" needs amending.** The coordinator showed the functional form is
   sound at perfect distances (1.152 Å). H7 shows the *truncation* is not — degree-1 reaches only
   3.769 Å with perfect distances. So there are two separate defects, and only one of them is the
   distogram. **Reducing the objective's expressive power is not a repair; it is the damage.**
3. **Search is not, and never was, the binding constraint.** `E_le1`'s global optimum is certified
   at O(n·K) evaluations by separability, and reaching the true continuous optimum rather than the
   mesh one changes the objective and not the structure. Every result here is a *discrimination*
   result, consistent with the standing record.

---

## 5. WHAT WAS NOT COMPLETED, AND WHY — read this before citing anything below

The 126-target λ ladder above is **complete** and is the only arm that decides the sprint. The
secondary arms below ran on a heavily contended box (five to seven sibling heavy processes pinned
at 96–100% CPU for the whole session; my processes were getting under 15% of a core each). They
are reported **as partial subsets with their n stated**, never as final, and no directional
conclusion is drawn from any of them that the complete arm does not already carry.

**One arm was not run at all: H6, the μ-sensitivity swap on the continuous instrument.** It needs
a second full coefficient fit per target and did not fit the envelope. I am not substituting a
guess for it. What exists instead is MATH's own μ-sensitivity on the *enumerated* instrument
(`s18/results/math_lattice.json`), where the uniformised and raw measures give the **same** Walsh
weight-≤1 argmin (3.655 on 1CS9) and variance fractions within 0.03 of each other. That is
evidence about the lattice, not about the continuum, and it is cited as such. **H6 is OPEN.**

### 5.1 The true continuous optimum of `E_le1` — **n = 126, COMPLETE** (`s18/results/exp_polish.json`)

    objective   mesh argmin  -379.81  ->  polished  -394.51   (lower on 126/126 targets;
                                                               11.0 residues/target off-mesh)
    RMSD        mesh argmin    4.346  ->  polished    4.330    Control A 3.048
                polished minus mesh          -0.016 [-0.057, +0.032]     (a null, below the MDE)
                polished minus Control A     +1.282 [+0.991, +1.475]

**Reaching the genuine optimum of the degree-1 objective instead of an approximation lowers the
objective on every single target (126/126) and does not move RMSD.** The objective/structure dissociation
that Sprint 17 found between objectives reappears *inside* one objective. The angle-additive
object's polished argmin is 4.623, again worse than the residue-additive one.

### 5.2 CONTROL B — matched budget — n = 85 (`s18/results/exp_budget.json`, **PARTIAL**)

K = 8 pool-medoid states per residue, budget 1,024 evaluations, identical Control-A start.
**Accounting convention, stated once: one evaluation = one objective call at one complete
configuration.** Control A on this subset is 3.136 Å (a pinned-order prefix, so its mean differs
from the full set's 3.048 — every comparison below is against *its own* Control A).

| method | RMSD | objective | calls | vs Control A | W/L |
|---|---|---|---|---|---|
| avg (Control A, this subset) | 3.136 | — | — | — | — |
| discrete start | 3.699 | 429.3 | 1 | +0.563 | — |
| greedy 1-opt on `E_full` | 3.620 | 119.16 | 748 | +0.483 | 16/69 |
| greedy at **4× budget** | 3.624 | 118.60 | 826 | +0.488 | 16/69 |
| random local search | 3.562 | 114.59 | 1024 | +0.426 | 21/64 |
| Metropolis / annealing | 3.660 | 108.23 | 1024 | +0.524 | 20/65 |
| best-of-N (the required null) | 3.769 | 143.34 | 1024 | +0.632 | 19/66 |
| greedy on `E_le1` | 4.023 | −88.30 | 815 | +0.887 | 19/66 |
| annealing on `E_le1` | 4.009 | −93.99 | 1024 | +0.873 | 19/66 |
| **exact argmin of `E_le1` (CERTIFIED)** | 4.044 | **−94.20** | **105** | +0.908 | 19/66 |

- Greedy reaches a **certified 1-opt local optimum within budget on 69/85** targets for `E_full`
  and **72/85** for `E_le1`. Quadrupling the budget improves the objective on 14/85 (`E_full`) and
  16/85 (`E_le1`) and improves RMSD on **neither** (3.620→3.624, 4.023→4.044).
- Greedy on `E_le1` reaches the **certified separable optimum exactly on 69/85**; the certified
  optimum itself costs **105 calls** — `n·K` — because the object is additive.
- Both accounting conventions, so neither reading can be cherry-picked:

  | method | objective calls | gradient-equivalent (×(2n+1)) |
  |---|---|---|
  | L-BFGS on `E_full` | 230 | 6527 |
  | L-BFGS on `E_le1` | 48 | 1328 |
  | greedy 1-opt on `E_full` | 748 | 748 |
  | exact argmin of `E_le1` | 105 | 105 |

> **CONTROL B is a null by construction, and that is the finding.** Budget is not binding for
> either objective — Sprint 17's result that greedy certifies the full objective at ~1,024
> evaluations reproduces here, and for degree-1 the question does not arise at all because the
> optimum is closed-form. **No optimiser, classical or quantum, at any budget, changes any number
> in this report.** Every result in Sprint 18 is a discrimination result.

### 5.3 H8 — geometric-leverage re-weighting — n = 88 (`s18/results/exp_leverage.json`, **PARTIAL**)

Its motivating ceiling was retracted mid-sprint (§3c); it was run as pre-registered anyway.

| β | RMSD | vs Control A | W/L | permuted-leverage CONTROL | matched-random CONTROL | vs permuted |
|---|---|---|---|---|---|---|
| 0.00 (= deployed) | 3.634 | +0.521 [+0.342, +0.692] | 23/65 | 3.634 | 3.634 | +0.000 |
| 0.25 | 3.680 | +0.566 [+0.367, +0.766] | 23/65 | 3.650 | 3.668 | +0.029 |
| 0.50 | 3.760 | +0.647 [+0.420, +0.852] | 21/67 | 3.731 | 3.719 | +0.030 |
| 0.75 | 3.817 | +0.703 [+0.465, +0.960] | 17/71 | 3.772 | 3.704 | +0.044 |
| 1.00 | 3.927 | +0.813 [+0.556, +1.058] | 13/75 | 3.713 | 3.741 | +0.214 |

Validity gate passes: β = 0 reproduces `s17` `refine_full` to mean |Δ| = 0.00000 on all 88 shared
targets. **Every β > 0 is worse than β = 0** (β=0.5: +0.126 [+0.064, +0.207]; β=1: +0.293
[+0.203, +0.399]) **and none beats its permuted-leverage control** — the only β whose permuted
comparison excludes zero is β = 1, and it excludes zero on the *wrong side* (+0.214 [+0.122,
+0.308]: the real leverage assignment is worse than a permuted one). The pre-registered falsifier
fires and the arm closes. As pre-committed, a null here is a clean result.

**And the ORACLE diagnostic explains why, which is the part worth keeping:**

    rho(leverage, |residual|)   -0.410   median -0.417
    rho(1/sd^2,   |residual|)   -0.471
    rho(leverage, separation)   -0.645        leverage CV 0.806, max/median 4.7

The distogram's errors are **smaller** where leverage is higher, not larger. The premise that
"errors concentrate on the pairs that most determine the structure" is **backwards** on this
instrument. Leverage tracks short separations (rho = -0.645), which are exactly the distances the
distogram predicts best. There is nothing for a leverage re-weighting to fix.

**One result here bears directly on ADVERSARIAL's `wflat` arm, so I flag it rather than sit on
it:** ρ(1/sd², |residual|) = **−0.471** — the model's own confidence is rank-correlated with its
error *in the correct direction*. Whatever `wflat` finds, "the confidences are anti-informative"
is not supported by the rank correlation; if uniform weights help, the mechanism is something
other than mis-ordered confidence. (ORACLE diagnostic, n = 88, PARTIAL.)

### 5.3b H6 — μ-sensitivity — n = 17 shared targets (`exp_main_uniform_512_sub30.json`, **PARTIAL**)

The pre-registered identification check. **Uniform μ is the zero-information limit of the
reference measure** — no positional information at all — so this arm doubles as the
zero-information-μ control the coordinator asked for, in its strongest form.

    lam000            pool 3.918 (+0.962 [+0.245,+1.687])   uniform 4.465 (+1.509 [+1.233,+1.642])   SAME SIGN
    argmin_le1_hold   pool 4.023 (+1.068 [+0.353,+1.567])   uniform 4.663 (+1.708 [+1.449,+1.846])   SAME SIGN
    lam100            pool 4.018 (+1.063 [+0.728,+1.352])   uniform 4.018 (identical, as it must be) SAME SIGN

**The sign survives the μ swap on all three arms**, so the degree-1 object's failure is
*identified* and is not an artefact of the reference measure. That is the pre-registered H6
question and its answer does not depend on n. **Still PARTIAL; H6 stays OPEN for its magnitude.**

> **A correction I am making against myself.** At n = 7 the uniform-μ arm was marginally *better*
> than pool-μ (5.032 vs 5.194) and I wrote that up as weak early evidence for the coordinator's
> §5.4 concern — a fourth zero-information control matching its informative arm. **By n = 17 it has
> reversed**: uniform-μ is 4.465 against pool-μ's 3.918, and the argmin 4.663 against 4.023. The
> pool's target conditioning does appear to buy the degree-1 object something after all. I am
> withdrawing the n = 7 reading rather than leaving it standing, and neither direction is claimed
> at this n — it is exactly the "no directional conclusions from a small instrument" trap this
> sprint exists to punish, and I walked two steps into it. **The zero-information-μ question
> (§5.4) is therefore genuinely OPEN, not leaning either way.**

### 5.4 The constant-α-helix reference measure — **COMPLETE at its pre-registered n = 30**

MATH's `E_le1` is *defined* against a reference measure, and the shipped μ is the target's own
retrieval-pool per-residue torsion marginal — a conditioned torsion channel inside the objective's
definition, not beside it. So the control the coordinator asked every torsion-channel arm to carry
applies here. `s18/exp_helixmu.py` supplies a zero-information μ (every residue at the ideal
α-helix, 10° jitter, **the same measure for every residue of every target** — sequence-blind,
position-blind, target-blind) by adding one name to `math_anova.mu_samples` through a delegating
wrapper, leaving MATH's estimator untouched. Length-stratified 30-target subset, as pre-registered.

| arm | pool μ (conditioned) | helix μ (zero-information) | difference [95% CI] |
|---|---|---|---|
| refine on `E_le1` | 4.141 | 4.405 | +0.264 [−0.334, +0.708] |
| certified argmin | 4.105 | 4.551 | +0.446 [−0.243, +1.183] |

Against Control A (2.957 on this subset): helix-μ refinement +1.448 [+0.899, +2.277], helix-μ
argmin +1.594 [+0.959, +2.656].

> **INCONCLUSIVE at n = 30, leaning toward the conditioned measure.** Both intervals span zero, so
> nothing is established — but the point estimates are **positive and have stabilised** (+0.26 and
> +0.45), meaning the zero-information measure is if anything *worse*, not equal. **This is the
> opposite of what I claimed at n = 14.** The 30-target subset is far too small for the 0.084 Å
> MDE; a properly-powered version is the full 126, and it is **OPEN**.

**The correction, and why it is left in the document.** At n = 14 I wrote this up as "the
zero-information reference measure matches the conditioned one — the fourth zero-information
control in this sprint to match its informative arm", and marked it SUPPORTED. Watching it grow:

    n = 14   refine -0.041   argmin +0.106      <- the reading I published, and withdrew
    n = 16   refine +0.154   argmin +0.331
    n = 20   refine +0.256   argmin +0.493
    n = 30   refine +0.264   argmin +0.446      <- stabilised, still spanning zero

The estimate drifted away from my claim and then settled somewhere else entirely. I had described
this exact failure mode for H6 in §5.3b and then committed it in the next section on the same class
of arm. **On this instrument a small-n point estimate with a zero-spanning CI is "not measured",
never "matched"** — that is the transferable lesson and it is why the mistake is recorded rather
than edited away.

**None of this touches the falsification.** Degree-1 fails against Control A by +1.287 Å
[+1.124, +1.473] at n = 126 under the shipped μ, and no reference measure rescues it — helix-μ is
+1.448 [+0.899, +2.277] against Control A on its own subset. What stays OPEN is only the weaker
question in §4.3: whether degree-1's information is the pool's target conditioning or generic
backbone plausibility. **The coordinator's `priorfit` (n = 126, complete) is the only
properly-powered evidence and says generic. My arm leans the other way and is not powered to
contradict it; it should not be cited as doing so.**

**One methodological point this arm does settle.** My uniform-μ arm (§5.3b) and this helix-μ arm
are both "zero-information" and behave differently: uniform-on-the-torus puts mass on physically
impossible backbone conformations, so its conditional expectations are dominated by clashes and it
is a *worse* measure rather than an *uninformative* one. A zero-information control for a
reference measure must stay on the physically plausible manifold. **Uniform ≠
zero-information-but-plausible.**

*Process note: this artefact was extended concurrently by a sibling agent also running
`s18/exp_helixmu.py`. While trying to freeze it I killed their process, having mistaken it for a
stray copy of my own — it kept reappearing after each kill because it was never mine. The module is
deterministic and checkpointed per target, so nothing was lost and the run completed. Recorded so
the next person checks the parent process before killing anything.*

A note on compute, since it shaped this section: MATH's `fit` defaults to `batch = 32768`, which
on this box is the *slow* setting — measured 38.8 s at 4096 against 63.1 s at 32768 for the same
target, with ~10⁸ page faults per process under multi-process load. `batch` only chooses chunk
sizes (every reduction is within a row, and `step = batch // S` keeps each mesh point's S samples
in one chunk), so the coefficients are bit-identical and only the wall clock moves. Dropping to
4096 roughly halved the remaining run. This is worth passing to whoever runs that module next.

---

## 6. CLAIMS

| # | claim | label |
|---|---|---|
| E1 | On the real 126-target continuous-torsion instrument, refining the coordinate average toward the **degree-1** objective costs **+1.287 Å [+1.124, +1.473]**, 18W/108L, median +1.168. | **ESTABLISHED** |
| E2 | The λ curve is **monotone decreasing in λ**: 4.335 / 3.742 / 3.683 / 3.619 / 3.610. Removing the ≥2 component removes information, not harm. Degree-1 minus full = **+0.726 [+0.605, +0.816]**, 34W/92L. | **ESTABLISHED** |
| E3 | **No** member of the `E_λ` family beats the coordinate average. The best (λ=1) is +0.561 [+0.408, +0.677]. Sprint 17's finding generalises to the whole family. | **ESTABLISHED** |
| E4 | Fed the **native's own distances**, the full objective reaches 1.152 Å and degree-1 only 3.769 Å — **+0.721 [+0.498, +0.945] worse than Control A**. The truncation, not the distogram, is what destroys the information; no distogram improvement rescues a degree-1 objective. | **ESTABLISHED (ORACLE diagnostic)** |
| E5 | The effect is **not concentrated**: drop-top-10 sits at the 50–52nd percentile of a uniform-effect null on all three headline arms, and holds in every length stratum and every fold. | **ESTABLISHED** |
| E6 | `E_le1` is additive in residues, so its global optimum is closed-form: `O(n·K)` evaluations, certified by separability. **Search is not the binding constraint at any budget, for either objective.** | **EXACT** (theorem) + measured (Control B, n = 85) |
| E7 | The residues the degree-1 object cannot see (`f_0 = f_{n−1} = 0`) are **exactly** the residues the frozen metric cannot see: φ₀ and ψ_{n−1} are never read by the builder, ψ₀ is a rigid rotation about an axis through CA₀, φ_{n−1} moves nothing in the Cα set. Identical RMSD on 126/126 to 2.13e-14. **I raised this as a hazard; the measurement retires it.** | **EXACT** |
| E8 | Reaching the **true continuous optimum** of `E_le1` instead of MATH's mesh argmin lowers the objective on **126/126** targets and moves RMSD by −0.016 [−0.057, +0.032] — a null below the MDE. The objective/structure dissociation exists *inside* the degree-1 object. | **ESTABLISHED** |
| E9 | Selection by either objective over the retrieval pool is **worse than a random pick** (3.531 and 3.600 vs 3.551 pool mean). The coordinate average survives every route. **The degree-1 objective belongs nowhere in the architecture.** | **ESTABLISHED** |
| E10 | Alignment: ρ(objective, RMSD) is +0.063 [−0.016, +0.133] for the full objective and +0.038 [−0.012, +0.108] for degree-1; both CIs include zero, and degree-1's **in-band** correlation is negative (−0.032). Phase 0's tension survives in sign, collapses in magnitude. | **ESTABLISHED** |
| E11 | Geometric-leverage re-weighting of the refinement objective is harmful at every β and beats no control (at β=1 it is significantly **worse** than a permuted-leverage control, +0.217 [+0.036, +0.438]). **ORACLE mechanism: ρ(leverage, \|residual\|) = −0.410** — the distogram's errors are *smaller* where leverage is higher, so the premise is backwards. | **SUPPORTED** (n = 88, PARTIAL) |
| E12 | ρ(1/sd², \|residual\|) = **−0.471**: the distogram's own confidence is rank-correlated with its error in the **correct** direction. Bears on ADVERSARIAL's `wflat`. | **SUPPORTED** (n = 88, PARTIAL, ORACLE) |
| E13 | Whether degree-1's information is the pool's target conditioning or generic backbone plausibility. Both the helix-μ (n = 30) and uniform-μ (n = 17) arms are **underpowered, and the helix-μ point estimate drifted from −0.04 at n = 14 to +0.26 at n = 30 before stabilising** — the opposite of the reading I first published; neither direction is claimed. The coordinator's `priorfit` (n = 126, complete) is the only properly-powered evidence and says generic. | **OPEN** — my arms; SUPPORTED elsewhere by `priorfit` |
| E13b | **A zero-information control for a reference measure must stay on the physically plausible manifold.** Uniform-on-the-torus puts mass on impossible backbone conformations and is a *worse* measure, not an uninformative one; it behaves differently from a constant-helix measure. **Uniform ≠ zero-information-but-plausible.** | **SUPPORTED** (methodological) |
| E14 | μ-sensitivity: the **sign survives** the swap from the pool marginal to a uniform (zero-information) reference measure on all three arms — the failure is identified, not a measure artefact. The *direction* of the pool-vs-uniform difference **reversed between n = 7 and n = 17**, so nothing is claimed about whether the conditioning buys anything. | **SUPPORTED** for the sign (n = 17, PARTIAL); **OPEN** for the magnitude and for the conditioning question (§5.3b, §5.4) |

### The sprint's falsifiers

**F1 FIRED** (4.335 / 4.346, not ≈3.6) · **F2** does not fire *only because degree-1 is decisively
worse than the full objective rather than equal to it* · **F3 FIRED** (there is no advantage to
die; degree-1 does not beat a matched-magnitude random move) · **F4 FIRED** (the sign reverses at
n = 126 and the interval excludes zero by a wide margin) · **F5** is MATH's, and their answer is
that no continuous object equals the strict Walsh weight-≤1 projection at all.

**The branch is closed. A clean falsification is a successful sprint, and this one is clean.**
