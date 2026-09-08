# SPRINT 20 — LEDGER

---

## L1 — SPRINT 19 §5 IS DELIVERED: the quantum SAMPLER hypothesis is REFUTED, and two real positives survive (2026-09-07, WORKSTREAM D)

`s19/results/qb_main.json` was **complete on disk** (126/126 targets, 17 arms, zero missing cells)
and never analysed when Sprint 19 froze. WORKSTREAM D executed `s19/PREREG_B.md` sections 6-9
**exactly as written** — no re-running, no new structures, pure analysis of the prior lane's own
artefact. `s20/results/D_QB_CLOSE/`, COMPLETE; code `s20/d_qb_close.py`.

**Basis is clean throughout** (the L20 rule applied on its first day): `realised_ORACLE` is the
**projected built chain** for every arm including the incumbent, `avg_ORACLE` is the point cloud,
both recorded per arm. The (M,D) identity `readout^2 = M^2 - D^2` asserts at max |resid|
**5.3e-14 A^2**, PASS.

### Realised (built chain), n = 126

    pool500      3.230   <- the shipped retrieval pool, ZERO objective evaluations
    c_marg       3.326      zero-information: matched empirical torsion marginals
    c_metroH     3.407      best classical sampler
    c_helix      3.426      zero-information constant helix
    q_a1.00      3.486      BEST QUANTUM
    c_lbfgs      3.658      multi-start L-BFGS on the objective, 8193 evals
    q_untrained  3.696

**P1 (PRIMARY) FAILS.** `q_a1.00 - c_metroH = +0.079 [+0.001, +0.163]`, fold-clustered
[-0.039, +0.203], 60W/66L. **The quantum arm is worse**, and nowhere near the pre-registered -0.084
with a CI excluding zero.

**FOUR of the nine pre-registered kill rules FIRE**: loses to the zero-information marginal sampler
(+0.160 [+0.084, +0.240], 5/5 folds); loses to the incumbent pool (+0.256 [+0.131, +0.390], 5/5);
**matched by a classical thermostat** (+0.001 [-0.089, +0.094]); matched by annealing
(-0.001 [-0.097, +0.100]).

**P3 (generation) FAILS.** Best quantum generation ceiling 1.672 vs `c_cem` 1.605 = +0.067
[+0.013, +0.121]; at matched count 500, +0.041 NOT MEASURED. Against `pool500` at matched count:
**+0.404 [+0.270, +0.546], 5/5.** *The one axis the lane expected to have a chance does not.*

**P2 (eps-dominance) FAILS.** Every trained quantum point is dominated by the classical point set
with a CI excluding zero. The leave-one-out control does its job: the only non-dominated points are
`pool500`, `c_helix` and `q_untrained` — and `q_untrained` is non-dominated **only because it is the
most diverse** (D = 2.581), which Sprint 19's C3 already established is an accounting identity and
not a value.

### THE FINDING THAT MATTERS MOST FOR PRIORITY 1

> **No sampler on the board, quantum or classical, at 8192 objective evaluations, beats the
> ZERO-EVALUATION shipped retrieval pool on realised RMSD.**
>
> **And it is not a min-of-N artefact**: `pool500` selects its 75 from 500 candidates while every
> sampler selects its 75 from **8192**. `c_lbfgs`, which optimises the deployed objective hardest, is
> the **second-worst arm on the board.**

This is `search-saturates-discrimination-binds` reproduced on a **fifth instrument**, and it prices
Sprint 20's Q2 before Q2 starts: **15 independent generators, every one worse than retrieval after
the terminal operator.**

### TWO POSITIVES OUT OF THE WRECKAGE

**1. THE VQE GENUINELY TRAINS — for the first time in this programme.**
`q_a1.00 - q_untrained = -0.210 [-0.339, -0.082], 5/5 folds`; on generation **-0.299 [-0.407,
-0.195], 5/5.** Against the **mandatory** control — best-of-N from the *untrained* circuit, not an
initialisation mean.

> **This is a SCOPE CORRECTION to a project memory.** `concentration-is-wrong-when-discrimination-binds`
> recorded *"running the VQE is WORSE than not running it, 0/12 cells, +0.65 to +1.32 A"* — measured
> on the **k=4 lattice register**. On the **continuous basin-latent encoding at matched budget the
> sign REVERSES**: training beats best-of-N from the untrained circuit, decisively, 5/5 folds.
> **PRIORITY 2 is satisfied on its own terms — the pillar is genuine and it works. It simply does not
> beat anything classical.**

**2. ENTANGLEMENT BUYS NOTHING, measured directly.** `q_a1.00 - c_prod0.25` — identical pipeline with
the **CNOTs deleted**, same CVaR estimator, same seed — is **-0.013 [-0.095, +0.077], NOT MEASURED**.
The correlated basin latent, the specific mechanism the lane was built to test (*"a whole region
flipping between two conformer families"*), is worth 0.013 A with a CI spanning zero.
Nearest-neighbour mutual information in the latent is **0.045 bits** mean.

### Labels

    PREREG_B P1 / P2 / P3                REFUTED (n=126, all pre-registered)
    "the VQE trains on the continuous encoding"   ESTABLISHED (ORACLE scoring, native-free decisions)
    "entanglement contributes nothing here"       NOT MEASURED, -0.013 [-0.095, +0.077]

**The last label is the lane's own and it is the right one** — a CI spanning zero is not a refutation,
and calling it one would have been the easier and worse choice.

---

## L2 — SECOND LABEL DRIFT: the +0.164 A "AMBER tax" is the PROJECTION's, and the corrected role INVERTS (2026-09-07, WORKSTREAM D)

**Coordinator-verified against the production cache, n = 126, target-level paired bootstrap** — all
four numbers reproduce to the digit:

    point cloud 3.0483 | lam=0 3.2041 | lam=0.3 3.2148 | AMBER emission 3.2355

    projection lam=0.3  vs point cloud   +0.1664 [+0.1300, +0.2018]  16W/110L   <- THIS is the +0.164
    projection lam=0    vs point cloud   +0.1557 [+0.1229, +0.1890]  18W/108L
    AMBER emission      vs point cloud   +0.1871 [+0.1474, +0.2274]  20W/106L
    AMBER on the BUILT chain (DEPLOYED)  +0.0207 [+0.0141, +0.0276]  40W/86L    <- AMBER's ACTUAL cost

**THE DEPLOYED AMBER IS NOT THE ARM ANY OF THOSE FIGURES DESCRIBED.** `core/pipeline.Config` has
`amber_k = 10.0`, `amber_steps = 0` (minimise to convergence). **Every "AMBER k=30" figure in Sprints
16 and 19 is a different restraint constant applied to a different input** — the point cloud, not the
built chain.

**AND THE SIGN INVERTS, not merely the magnitude.** `s16/CLAIMS.md` J1 measured both repair operators
against the *same* point cloud: projection 3.2052 (+0.1554), **AMBER k=30 3.1831 (+0.1333)**. So
**AMBER at k=30 is a BETTER repair operator than the projection by −0.022** — exactly what the
project memory `averaging-space-beats-the-objective` already records. **As briefed, the standing role
read AMBER as costing 0.164 A when, against the only operator it competes with, it is the cheaper of
the two.**

**PROVENANCE — Sprint 18 was RIGHT and the drift happened downstream.** `s18/FINAL_REPORT.md` §6.2 is
titled *"The repair tax is irreducible and the AMBER pass is never spent"* and its +0.164 is
explicitly the projection's; `s18/PHYS_STATUS.md:232` calls it *"the incumbent PROJECTION's +0.163
tax"*; `s19/BRIEF.md:127` still keeps the two apart. **`s20/BRIEF.md` §3 collapsed them into one
sentence** — coordinator error #8, and the third of the "label drift across sprints" class.

**Why it mattered this sprint rather than as bookkeeping**: Q1 asks whether the AMBER landscape is
genuinely harder and explicitly warns lanes not to assume "AMBER is more physical" is the
explanation. **A lane handed "AMBER costs 0.164 A" has been handed the opposite prior to the truth.**
Corrected in the live brief before any lane acted on it.

### The cross-basis sweep, reported as a NEGATIVE

WORKSTREAM D swept Sprints 16-19 for the same class of error:

* **`s18/exp_report.py:352`** prints *"pool -> coordinate average (INCUMBENT)"* and computes all seven
  [H5] rows against the point cloud, including `pool_best`/`pool_mean`/`sel_full`/`sel_le1`, which are
  **real fragments with valid geometry**. The whole table was recomputed on the built baseline
  (`proj`, 3.2126) at n = 126: **NO SPRINT-18 CONCLUSION CHANGES.** The offset is the constant +0.164
  and every arm loses to both baselines by a far larger margin (worst +1.123 [+0.891, +1.366]).
  **Recorded so the sweep is on the record as negative rather than skipped.**
* **`s19`** files record `avg` in the row; only `a_sepfix.py:164` *labels* it "the incumbent" — same
  class, already flagged at L20.
* **`s15/expand.py:147`** already sweeps an isotropic scale on the average. **The exploit was noticed
  in Sprint 15 and never connected to the reporting basis.**

> **The drift is real and recurring, and in every instance found it is confined to the LABEL — except
> this AMBER one, where it also inverts a standing role.**

---

## L3 — WORKSTREAM B's OWN REPORT: independent agreement, three EXACT structural results, and one label the coordinator is adjudicating AGAINST it (2026-09-07)

`s19/agentB_FINDINGS.md`. Continuous torsion space, no lattice, 8,192 objective evaluations per
budgeted arm. **Arrives independently of WORKSTREAM D's analysis of the same artefact (L1) and agrees
on the closure** — two lanes, two readings, same verdict.

### The closure, sharpened

**A genuine CVaR-VQE earns no role as a structured sampler, and it loses to a ZERO-INFORMATION null.**
`c_marg` draws i.i.d. from the same per-residue von Mises basin mixtures the quantum sampler draws
from and **never reads the objective for any decision**. It has the **best realised RMSD (3.326) and
the best selection ceiling (2.491) of every budgeted arm.** Best quantum arm against it:
**+0.160 A [+0.089, +0.241], 46W/80L, folds 5/5.**

**P1** fails 6 of 7 clauses, two *significantly in the wrong direction*. **P3** — a generation **tie**
with the zero-information null (+0.009 [-0.063, +0.090]) that becomes a selection **loss** (P3b:
+0.190 vs the null, +0.396 vs the incumbent, both 5/5). At **matched candidate count**, BLOSUM
retrieval generates better structures than every sampler by **0.36-0.49 A**. **P2** — every trained
quantum arm dominated; the only undominated arms are the four that read the objective least or not at
all.

Warm starts, adaptive alpha, alpha in {0.02, 0.5}, deeper `mps3f`, CVaR+uniform and CVaR+annealing
mixtures, 4 independent seeds: **not one beats the null**, and the best is the one spending half its
budget on zero-information draws.

### Three structural results that explain WHY it is not a budget problem

* **T1 (EXACT)** — **CVaR's minimiser is a FACE.** It constrains an alpha-tail and is *indifferent*
  to the rest of the distribution. Verified with three exactly-CVaR-optimal laws whose diversity
  spans **0.00 -> 3.54 A at machine-identical objective.** This is a theorem about the objective, not
  an observation about a run.
* **T4 (EXACT)** — the deployed ansatz is **bond-dimension 4, i.e. a <=16-state HMM: classically
  samplable by construction, with no separation available at ANY budget.** And training *reduces* its
  latent correlation (1.33 bits against the untrained 1.52).
* **T3** — the discrete latent carries **~53% of Var(E) and only ~27% of Var(RMSD)**.

### The gauge arm — PARTIAL n=20, readable, and it PASSES

`qb_gauge_PARTIAL_n20.json` tests the latent's **entire** labelling freedom (`b_i <-> 1-b_i`, a
`Z2^n` gauge — not a sampled corner). Identity labelling's mid-rank percentile in its own orbit:
**mean 0.500, median 0.479** — dead centre; orbit difference +0.018 [-0.130, +0.170]; the <0.05
firing rule does not fire. **Contrast Sprint 18's `W1` at percentile 0.00, p < 0.0001.**
**Nothing in this lane's conclusion lives in the encoding.**

### COORDINATOR ADJUDICATION — the entanglement label stays NOT MEASURED

The two lanes disagree. B reports the CNOT-deletion contrast as **-0.013 [-0.092, +0.079]** and calls
it *"tight enough at n=126 to exclude the 0.084 MDE — measured worth nothing, not unmeasured"*.
D reported -0.013 [-0.095, +0.077] and labelled it **NOT MEASURED**.

> **D's label is the correct one, and the direction is what settles it.** Entanglement *helping*
> would appear as a **negative** difference. The interval's lower bound is **-0.092**, so an
> MDE-sized benefit of -0.084 lies **INSIDE** the interval and is **not excluded**. The interval
> excludes a *harm* larger than +0.079 and a *benefit* larger than +0.092 — it does not exclude the
> effect of interest.
>
> **Correct statement: entanglement's contribution is NOT MEASURED; the data exclude a benefit larger
> than ~0.09 A.** B's own T4 is the stronger argument anyway — a bond-dimension-4 ansatz is a
> <=16-state HMM and cannot offer separation at any budget, which closes the question structurally
> rather than statistically.

### THE DECORRELATION NUMBER — NOT MEASURED, and B names the reason it would have been a foregone question

B states plainly that it did not compute it, that it was not in the pre-registration, and that
**`qb_main.json` persists per-arm summary statistics only, not the emitted candidate sets**, so the
artefacts cannot answer it without a rerun.

> **AND THE DESIGN CONSTRAINT B VOLUNTEERED IS THE MOST IMPORTANT THING IN THIS REPORT FOR SPRINT 20.**
> Its samplers draw torsions from **von Mises mixtures fitted to the BLOSUM retrieval pool's own
> per-residue distribution.** The candidate space is *derived from* the pool, so it is **the least
> likely thing on the board to decorrelate from the pool's error.**
>
> **Any Q2 test whose candidate space is fitted to the retrieval pool is asking a foregone question.**
> Routed to WORKSTREAM A and WORKSTREAM B as a binding design constraint.

### A BY-PRODUCT WORTH MORE THAN THE QUANTUM ARM — the operator law is WRONG ON WIDE SETS

Over **2,142 (arm, target) cells**:

    d_out = 0.803 * d_set_mean + 0.298 * d_set_best      R^2 = 0.960

**The recorded law has 0.04 on `set_best`. On wide candidate sets the coefficient is 0.298 — seven
times larger.** So a generation-ceiling improvement converts at **0.30 A per A**, not 0.04. **This
raises the value of every generation improvement by ~7x on wide sets**, and it is measured across
2,142 cells rather than fitted on one design. Project memory corrected.

**B also records, unprompted, that its own design rationale was the BRIEF §3 conformer-family
coherence claim the coordinator retracted on 2026-09-06** — *"the construction stands on its own
terms but I am not entitled to that motivation"*. Recorded as B stated it.

---

## L4 — Q1 IS PRICED BEFORE ITS LANE SPENDS A BUDGET: every circuit-side landscape metric is a difficulty proxy, and the CVaR tail is worth nothing (2026-09-07, WORKSTREAM D)

Both results n = 126, both from `s19/results/qb_main.json`, **no new computation.**

### 1. Partial out target difficulty and every distributional metric collapses

Control is `pool500`'s realised RMSD — an ORACLE difficulty proxy; rank-based partial Spearman.

    arm       metric                raw rho   partial rho
    q_a0.05   entropy_bits           +0.144      +0.065
    q_a0.05   max_prob               -0.241      -0.099
    q_a0.05   ess                    +0.204      +0.082
    q_a0.05   nn_mi_bits (latent)    +0.035      +0.167
    q_a0.05   obj_best               +0.335      +0.234
    q_a0.05   D (diversity)          +0.499      +0.356
    q_a0.05   M (member error)       +0.982      +0.930
    q_a1.00   entropy_bits           +0.222      +0.051
    q_a1.00   max_prob               -0.268      -0.066
    q_a1.00   ess                    +0.248      +0.051
    q_a1.00   obj_best               +0.290      +0.135

The CNOT-free classical control `c_prod0.25` shows the same pattern (+0.156 -> +0.055).

> **Circuit entropy, max probability, effective sample size and latent nearest-neighbour mutual
> information are all proxies for "this target is hard" and carry essentially nothing else.** Only
> `M` survives, and `M` is member error — **nearly the outcome itself, an identity rather than a
> predictor.**

This is BRIEF §4's own warning — *"a landscape metric with no predictive relationship to structural
outcome must not be over-interpreted"* — **measured rather than asserted**, and it lands directly on
**CVaR concentration**, one of the eight candidate explanations Q1 enumerates.

> **STANDING REQUIREMENT for the Q1 lane: partial out target difficulty on every metric reported.
> On this instrument a raw rho IS a difficulty measurement.**

### 2. The CVaR tail parameter is worth nothing, and alpha = 1 is the best quantum arm

Realised RMSD (built chain), paired against `q_a1.00` — **alpha = 1.00 is the plain sample mean with
NO tail truncation at all**:

    q_a0.05  - q_a1.00   +0.012 [-0.084, +0.105]  56W/70L  folds 3/5
    q_a0.25  - q_a1.00   +0.013 [-0.066, +0.089]  56W/70L  folds 2/5
    q_anneal - q_a1.00   +0.025 [-0.057, +0.108]  52W/74L  folds 3/5

**All three NOT MEASURED, and all three point estimates are on the wrong side: every tail setting is
worse than no tail.**

**PRIORITY 2 — what "preserve genuine CVaR-VQE" actually buys, stated precisely.** The CVaR machinery
in `core/quantum.py` **is genuine and correct**: `tail_indices` reproduces the stable-sort mask
exactly with a documented tie rule; `cvar_gradient`'s "const" baseline is a true control variate at
**cos +1.000000** against the exact-expectation reference; and the recorded "tail" baseline defect is
**reproduced verbatim as a named, measurable arm rather than quietly fixed.** *The pillar is intact.
Its measured contribution to RMSD on this instrument is zero to slightly negative.*

### Theory, new since Sprint 19's literature pass

**arXiv:2605.02850 (Quantum Tilted Loss, 2026)** places CVaR in a tilted-loss family and proves:
(i) it **cannot eliminate barren plateaus**, only reshape local geometry; (ii) CVaR's **hard
alpha-quantile cutoff creates optimisation discontinuities** as parameters shift samples in and out
of the tail set, where a smooth tilted loss does not; (iii) the sample complexity of resolving the
sharpened gradient scales as **O((e^{|gamma|Delta} - 1)^2 / gamma^2 eps^2)** — **exponential in the
tilt**. *Landscape sharpening is paid for in shots, exponentially.*

On a 512-shot budget an alpha = 0.05 tail is **26 samples**. **That the alpha = 1 arm wins is exactly
the predicted regime**, and this is the citation for saying so.

### The corollary the Q1 lane must not miss

> **"CVaR concentration" is not a clean landscape property of the Hamiltonian at all.** The
> tail-membership discontinuity is a property of the **CVaR RULE**, present identically under
> `H_Legacy` and `H_AMBER`. So any difference in measured CVaR concentration between the two is a
> difference in their energy **spectra** passed through a **shared discontinuous operator**.
> **Derive the operator before interpreting its statistic** (§9).

**Compute note**: the box is at **100% CPU with 7-8 python processes**. WORKSTREAM D is holding its
Block E behind Block C rather than adding a second heavy process — correct behaviour, and the
standing rule that **timing measured under contention is not a property of the code** applies to
every lane tonight.

---

## L5 — THE SHARED-REFERENT FLOOR: the "two thirds" framing dies, the ORDERING survives, and S5 was the main result all along (2026-09-07, WORKSTREAM D)

Pre-registered as Block C with falsifier **F-D1 fixed before the run**. `s20/results/D_C_ALIGNNULL/`,
COMPLETE, n = 126; code `s20/d_alignnull.py`.

### The attack

`s19/a_source.py` reports coherent-component alignment as
`corrC(A,B) = corr(d(X_A) - d_nat, d(X_B) - d_nat)`. **Both arguments are scored against the SAME
native.** Writing `d(X) - d_nat = (d(X) - d_typ) - (d_nat - d_typ)`, the second term is **identical
for every family by construction**. Any two realisable length-n peptide structures scored against one
native must correlate. **Nobody measured how much.**

### The null — three constructions, and they agree

    N1  pool x pool      two independent windows from this target's own K=500 pool, no fit   0.507 [0.472, 0.542]
    N2  fit x fit        the SAME s15/align_lib.fit, SAME cached start, SAME 1/sd^2 weights,
                         driven to those two windows' distance fields          -- MATCHED    0.505 [0.469, 0.540]
    N3  cross-target     windows of the same length from DIFFERENT targets                   0.529 [0.489, 0.568]

Eight draws per target, N1 and N2 sharing draws so their difference isolates **the fit operator and
nothing else**. **The three agreeing within 0.024 is the internal check**: the floor is a property of
*"two realisable peptides scored against one native"* — not of the fit, not of this target's pool.

### F-D1 fires at the WEAK clause, not the strong one

Pre-registered: *if N2 >= 0.575 with a CI containing it, the zero-information references are AT the
null.* **They are not** — N2 = 0.505 and its CI excludes 0.575. **Both zero-information references
carry genuine excess**: sepprior **+0.070 [+0.028, +0.112]**, helix **+0.093 [+0.037, +0.146]**, both
5/5 folds, both CIs excluding zero. **S2 is NOT refuted.**

### What dies is the "TWO THIRDS" framing

Share of the ceiling recomputed as `(corrC - null)/(ceiling - null)`:

| pair | raw | excess [95% CI] | W/L | share AS PUBLISHED -> NULL-SUBTRACTED |
|---|---|---|---|---|
| same arch, different seed (CEILING) | 0.898 | +0.394 [+0.357, +0.432] | 123/3 | 1.00 -> **1.00** |
| same arch, different regularisation | 0.833 | +0.329 [+0.290, +0.366] | 122/4 | 0.93 -> **0.83** |
| **retrieval pool mean** | 0.784 | +0.279 [+0.244, +0.316] | 117/9 | 0.87 -> **0.71** |
| different architecture (PairNet) | 0.731 | +0.226 [+0.183, +0.270] | 101/25 | 0.81 -> **0.57** |
| ZERO-INFO separation prior | 0.575 | +0.070 [+0.028, +0.112] | 83/43 | 0.64 -> **0.18** |
| ZERO-INFO constant alpha-helix | 0.598 | +0.093 [+0.037, +0.146] | 84/42 | 0.67 -> **0.24** |

> **The zero-information references reproduce ABOUT A FIFTH of the same-architecture ceiling, not two
> thirds.**

**Consequences, all corrected in place:**

* `s20/BRIEF.md` §2 — "0.64-0.67 for sequence-blind references" now reads **"0.64-0.67 raw, 0.18-0.24
  above a 0.505 shared-target null"**.
* `s19/FINAL_REPORT.md` §0 and §2 — *"Two zero-information references each reproduce ~two thirds"*
  and the bolded **"the harmful component is largely sequence-independent"** **DO NOT SURVIVE.**
  Above the null the component is **mostly sequence-DEPENDENT.** The published artifact carried this
  too and is corrected.

### S5 IS STRENGTHENED, AND THE LANE THAT WROTE IT AS A CAVEAT WAS WRITING THE MAIN RESULT

Sprint 19 stated S5 defensively — *"sequence conditioning is not worthless; 0.833 and 0.784 sit well
above 0.598."* **After null subtraction the gap between conditioned (0.83 / 0.71 / 0.57) and
zero-information (0.24 / 0.18) is three to four times larger than it looked.**

### WHAT IS UNDAMAGED — and it is the part the mechanism needs

**S1's ORDERING survives intact**: ceiling > same-arch > pool > cross-architecture >> zero-information,
**in that order before and after**. **S3** (the incoherent component does not share) is untouched.
**S7 and W5** — the pool shares the harmful component *and* is simultaneously the estimator of it —
are untouched: **the pool's null-subtracted 0.71 is the second-highest number in the table and still
above cross-architecture.**

> **THE WALL IS WHERE AGENT A PUT IT.** WORKSTREAM D recorded in advance that it expected exactly this
> shape — the ordering surviving, the "two thirds" breaking — and that is what happened.

### Limitation volunteered by the lane

The null is built from pool windows at ~4.45 A mean member RMSD while the predictors' fitted
structures sit at ~3.5-3.7 A, and the null trends mildly **upward** with member distance (N1 0.507,
N3 0.529). **A quality-matched null would be slightly LOWER and every excess slightly LARGER** — it
moves all six rows the same way and cannot change the ordering, but **0.18/0.24 is a lower bound** on
the zero-information share. A quality-matched null is the cheap follow-up.

### A NEW METHODOLOGICAL CLASS — "the shared referent floor"

Not an identity published as a discovery, but **a statistic published without deriving its operator**.

> ### Whenever two quantities are both measured as deviations from a COMMON REFERENCE, their correlation has a FLOOR set by that reference, and the floor must be measured before the correlation is interpreted.

A third class beside **Z1** (a control matched in the wrong space) and **Z2** (a measured null beats
an analytic one). **It is cheap to check and nobody in four sprints checked it.**

---

## L6 — F-D2 FIRES ON BOTH CLAUSES: the relaxation is what makes the AMBER objective DEFINED (2026-09-07, WORKSTREAM D)

Pre-registered as Block E **before the run**. `s20/results/D_E_AMBEROP/` COMPLETE; code
`s20/d_amberop.py`. Three targets (1CS9 n=9, 1A1P n=13, 1A13 n=14), 64 random bitstrings each on the
**production k=8 torsion register**, 192 calls per arm. Identical `AmberHamiltonian` object,
identical bitstrings, **only the iteration cap differs** (OpenMM's `maxIterations=0` means
*unbounded*, so 1 is the smallest honest floor).

    pdb    n  bits  finite@50  E50-E1 mean   rho(50,1)  rho(100,50)  cap bound  collapsed@50
    1CS9   9    27     48/64      -3.2e+05      0.358       0.835       1.00         1
    1A1P  13    39     29/64      -3.5e+08      0.827       0.797       1.00         4
    1A13  14    42     35/64      -3.3e+07      0.882       0.903       1.00         1

**1. THE CAP BINDS ON 100% OF CALLS — 192/192.** Doubling to 100 iterations still moves the energy by
more than 1e-3 kcal/mol on **every single call**. By Z6a, `H_AMBER` is not "AMBER energy": it is
**"50 iterations of restrained relaxation, then evaluate"** — a different operator, **on every call
rather than on an unlucky one.** This is Z6b's operational half, and the bound bound everything.

**2. THE ORDERING CHANGES, AND THE ORDERING IS WHAT CVaR EATS.**
`Spearman(E o Relax_50, E o Relax_1)` = **0.358 / 0.827 / 0.882**, against a pre-registered 0.95
threshold. And it has **not converged at 50** either — `Spearman(Relax_100, Relax_50)` is 0.797-0.903,
so **the ordering is still moving when the cap stops it.** A CVaR tail selected under `Relax_50` is a
**different subset** from the one selected under `Relax_1` or `Relax_100`.

**3. THE DECISIVE ONE, WHICH THE LANE DID NOT PREDICT — DEFINEDNESS.**

    finite at Relax_1     112/192   (58%)
    finite at Relax_50    186/192   (97%)
    finite at Relax_100   190/192   (99%)

> **THE RELAXATION IS WHAT MAKES THE OBJECTIVE DEFINED.** On **42% of the register** AMBER's energy at
> the as-built ideal-geometry coordinates **is not a finite number**, and the 50 minimisation steps
> are what turn it into one. Mean shifts of 3e5 to 3.5e8 kcal/mol say the same thing: the unrelaxed
> structures are in **hard steric overlap**.

**So `Relax_50` is not a smoothing convenience bolted onto an energy — without it there is no
objective on most of the space.** And that **kills the clean version of the lane's own
recommendation**: "just add the `E_amber o Relax_1` arm" yields a function that is **+inf on 42% of
the register** and cannot be CVaR-optimised as-is.

**4. AND THE RELAXATION MANUFACTURES PART OF THE LEGACY-AMBER AGREEMENT.**

    Spearman(E_legacy, E_amber o Relax_50) = +0.229
    Spearman(E_legacy, E_amber o Relax_1)  = +0.140

**The two Hamiltonians barely agree on ordering at all — +0.229 — and about a third of even that
agreement appears only after the relaxation.** Any Q1 statement of the form *"Legacy and AMBER rank
configurations differently"* is measuring `E_legacy` against `E_amber o Relax_50`, **not against
`E_amber`.**

### What Q1 does instead — re-specified in the brief

The comparison the brief wants (one register, one ansatz, one optimiser, two energies) is available
only if **both energies are defined on the same domain**:

* **(b) PRIMARY — hold the relaxation fixed on BOTH sides.** Give `H_Legacy` the same `Relax_50`
  preprocessing: evaluate the knowledge-based energy **at the AMBER-relaxed coordinates**. The
  relaxation becomes a **shared operator and cancels**. *This is the cleanest "change only H"
  available and it costs one wrapper.*
* **(a) DIAGNOSTIC beside it — restrict the domain.** Run both on the 58% subset where
  `E_amber o Relax_1` is finite, report the restriction **as part of the operator**, and compare
  `E_legacy`, `E_amber o Relax_1` and `E_amber o Relax_50` as three arms. The `Relax_50 - Relax_1`
  contrast is then **the relaxation's own effect, cleanly separated from the physics, on a common
  domain.**
* **(c) Otherwise** label every landscape conclusion **RELAXATION-CONFOUNDED**.

### Labels and limits, stated by the lane

**n = 3 targets, 192 calls per arm — this is NOT an RMSD claim and none is made.** The cap-binding
rate (**192/192**) and the definedness rates (**58% -> 97%**) are properties of the operator counted
over calls and are **not marginal**. The mean Spearman 0.689 is a three-target mean spanning
**0.358-0.882** and must be quoted as that range, not a point.
**Label: ESTABLISHED as a property of the operator; the Spearman magnitude INCONCLUSIVE at this n.**

### One self-correction, volunteered

The lane's smoke run wrote a **COMPLETE flag for a 1-target, 8-draw configuration** — *precisely the
partial-under-a-clean-name hazard it had flagged in three other files two messages earlier.* It
caught this on the next read, removed the stale flag, and added a guard so the module writes COMPLETE
only for the full pre-registered configuration. **Recorded because a lane that flags a class and then
commits it should say so.**

---

## L7 — WORKSTREAM D's LANE IS COMPLETE: five live claims damaged, two recommendations (2026-09-07)

`s20/PREREG_D.md` (4 blocks, **falsifiers fixed before every run**), `s20/agentD_FINDINGS.md` (1,006
lines), code `d_qb_close.py` / `d_alignnull.py` / `d_amberop.py`, artefacts `D_QB_CLOSE/`,
`D_C_ALIGNNULL/`, `D_E_AMBEROP/` each COMPLETE. **The benchmark manifest was never opened — a
SHA-256 was recorded as an integrity certificate that needs no read**, which is the right way to
prove a seal held.

Five damaged claims, all recorded above: L20 (the 3.048 point cloud), L2 (the +0.164 is the
projection's), L1/L3 (Sprint 19 §5 delivered as a clean falsification), L5 (the shared-referent
floor), L6 (Q1 is not executable as briefed).

### Two recommendations, priors against them stated first

**R1 — measure Q2's own `rho(e_coherent, e_pool)` on the generators already on disk.** The rationale
is sharp: **Sprint 19's fifteen samplers all fitted their basins to the retrieval pool, so that study
never contained a pool-independent generator** and its decorrelation question was foregone.

> **FLAGGED DISCREPANCY, not resolved.** WORKSTREAM B stated that `qb_main.json` *"persists per-arm
> summary statistics only, not the emitted candidate sets"*, so the number **cannot** be computed
> without a rerun. WORKSTREAM D proposes computing it *"on the sixteen generators already on disk"*.
> **One of these is wrong and I am not assuming which.** Referred back rather than acted on.

**R2 — the manifold-constrained Frechet mean as the repair operator**, externally supported by the
COMBO/MCORE literature, which **independently reports the same averaging artefact.**

> **The coordinator is taking R2.** It targets the largest measured loss in the pipeline that is not
> the predictor: the projection's **+0.166 [+0.130, +0.202]** tax, 1.9x MDE, paid on every target
> because the pipeline averages in **ambient R^3** and *then* projects onto the manifold. A Frechet
> mean computed **on** the manifold cannot contract by construction. See L8.

### Self-inflicted error, volunteered

The lane's smoke wrote a `COMPLETE` flag for a 1-target run — **the exact hazard it had flagged in
three other files two messages earlier.** Caught on the next read, removed, guarded, and recorded in
its own findings. **A lane that flags a class and then commits it should say so, and it did.**

---

## L8 — Q1 ANSWERED: Legacy is a COMPACTNESS model wearing a physics vocabulary, and the two potentials disagree about ordering (2026-09-07, WORKSTREAM C)

### The headline, n = 126

    Spearman(E_Legacy, E_AMBER) = -0.0886 [-0.1239, -0.0514]   median -0.095
    0/126 targets with |rho| > 0.8   ->  F-C1(a) does NOT fire
    Pearson                     = +0.0312   -- THE OPPOSITE SIGN, surviving winsorisation

> **They agree about catastrophes and disagree about ordering.** That is the Priority-3 answer in one
> line, and it is *not* "AMBER is more physical".

### What separates them — and it is the opposite of the naive prior

Partitioned against **200 matched-random partitions of identical cell sizes**, seven native-free axes
separate, and the axis is **compactness**:

> **Legacy-preferred candidates are 0.45 A MORE COMPACT** — `-0.4477 [-0.5110, -0.3911]`, **124W/2L,
> 5/5 folds, 3/3 length terciles, flat across a 5x threshold sweep** — with tighter contacts and
> better Ramachandran. **AMBER-preferred candidates are expanded with open sterics.**
>
> **Legacy is a compactness/typicality model wearing a physics vocabulary** — the opposite of the
> naive prior, since its *steric* term carries the largest weight (4.0).

**And both ORACLE axes are the two that do NOT separate**: rebuild distance `-0.0785 [-0.2120,
+0.0464]`, pool-error alignment `-0.0328 [-0.1526, +0.1010]`.

> **Complementary about GEOMETRY, jointly blind to ACCURACY.** The lane registered this prior in
> advance and it is confirmed.

### The mandated ablation table

Point-cloud basis, matched-random **one- and two-stage** controls. **Every ordered arm loses to NO
GATE with a CI excluding zero**: `legacy` +0.0763, `amber` +0.0496, `leg->amb` +0.0535, `amb->leg`
+0.0488. **Composites never beat their own first stage.**

### 53.5% OF THE CANDIDATE POOL IS A STERIC CATASTROPHE

The coordinator's L6 (unrelaxed AMBER non-finite on 42% of the *lattice* register) **reproduces on the
real top-75 domain in a worse form** — and the lane verified it before accepting it. On 9,450 top-75
rebuilds, unrelaxed AMBER is **finite everywhere** with median **+16,062**, p99 **3.1e13**, max
**5.5e23** kcal/mol; **53.5% above 1e4**, 31.1% above 1e6.

> **Finite-but-meaningless is more dangerous than +inf — nothing throws.** And the collapse is **not a
> lattice artefact: it is a property of the ideal-geometry rebuild itself.**

### EXACT, and new

`E_AMBER = f(x(theta))` but `E_Legacy = f(seq, x(theta), theta)`. **`phi_0` moves no atom**, so
**AMBER's gradient there is exactly 0.0 and Legacy's is not.** The structure-invisible channel is
**1.4-6.7% of ||grad E_Legacy||** — *the lane states the size because the categorical claim alone
oversells it.*

### SPRINT 19's Z7 IS DOWNGRADED: PLAUSIBLE -> NOT SUPPORTED

Z7 was *"any score that prefers compact, well-formed, pool-typical geometry selects toward the pool's
own error"*. The three scores that most literally instantiate it — `rg`, `leg_compact`, `min_heavy`
— are the **three lowest rows and all point the wrong way**: pooled `-0.0144 [-0.0663, +0.0437]`
against Sprint 19's own trio at `+0.0332 [+0.0102, +0.0572]`.

**NOT SUPPORTED rather than REFUTED — the difference CI spans zero.** And the lane adds that **its own
falsifier there was mis-specified (it fires on a null), and says so rather than claiming the kill.**

### THREE INSTRUMENT DEFECTS, ALL SELF-CAUGHT, TWO OF WHICH WOULD HAVE PUBLISHED NONSENSE

1. **The pre-registered finite-difference step would have published round-off as physics.** A Hessian
   divides energy noise by `h^2`; `core.amber`'s single point has ~1e-6-1e-4 kcal/mol granularity on
   E ~ 370. AMBER's negative-curvature fraction runs **0.259 at h=1e-2 (the plateau) -> 0.778 at the
   pre-registered h=1e-4**, while **Legacy is flat to 5.6e-07 across the same decade.** This is the
   lane's own falsifier F-C2b, registered against itself, firing.
2. **A spurious sqrt(n) in ALIGN**, caught by reading `s19/agentC_kv3.py:106` rather than its prose.
   Removing it turns a 4x disagreement into an **exact cross-sprint reproduction**: helix +0.0482 vs
   +0.0476, legacy +0.0256 vs +0.0256, amber +0.0258 vs +0.0258.
3. **`s12.instrument.paired`'s CI is i.i.d., NOT fold-clustered** — its `folds` argument only adds a
   per-fold breakdown. **This is programme-wide.** Every CI regenerated with `s18.phys_lib.paired`'s
   `ci_fold`, and **it changed two verdicts**: `amber - rand` (+0.0372 [-0.0025, +0.0699]) and
   `disto - rand` move from "excludes zero" to **NOT MEASURED**.

### Duplicate-work note, handled correctly

A second lane is running the identical Q3 experiment (`s19.agentC_pareto`, same `cfg_hash`).
**WORKSTREAM C stopped its own duplicate rather than the other lane's** and will merge rows under a
hash check.

---

## L9 — PROGRAMME-WIDE: `s12.instrument.paired`'s CI is i.i.d. over targets, not fold-clustered. Scope stated precisely. (2026-09-07, WORKSTREAM C found it; coordinator verified and scoped it)

**Verified at source.** `s12/instrument.py::paired` line 5:

    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])

**That is an i.i.d. resample of targets.** The `folds` argument does **not** enter the interval at
all — it only produces a `per_fold` dict of means. **293 call sites across s12-s20.**

### What this IS and IS NOT — the distinction matters and I am not overstating it

**It is not automatically wrong.** The 126-target instrument is **cluster-disjoint by construction**,
which is precisely the condition under which an i.i.d.-over-targets bootstrap is defensible.
Fold-clustering is the more **conservative** choice, and it widens intervals.

> **The defect is one of LABELLING, plus a systematic narrowing that only bites MARGINAL claims.**
> Several places in this programme's record say "fold-clustered CI" where the interval was i.i.d.
> Large effects are unaffected; claims whose interval *barely* excludes zero can flip.

**WORKSTREAM C measured the bite**: regenerating with `s18.phys_lib.paired`'s `ci_fold` **changed two
verdicts** — `amber - rand` (+0.0372 [-0.0025, +0.0699]) and `disto - rand` moved from "excludes
zero" to **NOT MEASURED.** Both were marginal to begin with.

### Which headline results are exposed, and which are not

**Not exposed** — these are far from marginal and no plausible widening reaches zero: `signflip`
**-1.202 [-1.408, -1.004]**; the in-manifold contrast **+0.298 [+0.262, +0.331], 121/5**; Legacy's
compactness axis **-0.4477 [-0.5110, -0.3911], 124W/2L, 5/5 folds**; the projection tax **+0.166
[+0.130, +0.202], 16W/110L**; the quantum arm against the zero-information null **+0.160 [+0.089,
+0.241]**; `wflat` **+0.153 [+0.045, +0.272]**.

**Exposed** — any claim whose interval barely excludes zero, in every sprint. **These must be
re-checked with clustering before being carried forward**, and several are already labelled
NOT MEASURED for other reasons.

### Standing rule, adopted

> **Quote the CI construction with the CI.** "i.i.d. over targets" and "fold-clustered" are different
> intervals and this programme has used the second name for the first object. Where a claim's
> disposition depends on which was used, **use the clustered one** — it is the conservative choice
> and the instrument's fold structure is real even if the targets are cluster-disjoint.

**Note on the coordinator's own arms**: the `_boot` helpers in `s18/objceil.py`, `s19/distobj.py`,
`s19/softsel.py` and `s20/frechet.py` are **also i.i.d. over targets**. Their headline results sit in
the "not exposed" category above, but the label must be corrected wherever "fold-clustered" was
claimed for them.

---

## L10 — THE R1 DISCREPANCY IS RESOLVED BY INSPECTION, against the audit lane (2026-09-07, coordinator)

WORKSTREAM D's R1 proposed measuring `rho(e_coherent, e_pool)` *"on the sixteen generators already on
disk"*. WORKSTREAM B stated the artefact *"persists per-arm summary statistics only, not the emitted
candidate sets"*. **I flagged it rather than assuming, and then read the file.**

`s19/results/qb_main.json`, per target: `pdb`, `n`, `fold`, `arms`. Each of the 17 arms holds
**scalars only**:

    tag, n_emitted, n_distinct, obj_best, obj_mean, gen_*_ORACLE, cov_*_ORACLE, M, D,
    readout_MD_ORACLE, identity_resid, sel_*_ORACLE, avg_ORACLE, realised_ORACLE, iters,
    param_disp, entropy_bits, max_prob, ess, nn_mi_bits_*, evals

**No candidate sets. No torsions, no coordinates, no per-pair residuals.**

> **WORKSTREAM B is right and WORKSTREAM D is wrong: R1 is NOT EXECUTABLE from disk.** The
> decorrelation number requires a rerun that retains the emitted sets.

**AND IT IS ALSO NOT WORTH EXECUTING, for B's own reason.** B recorded that its samplers draw from
**von Mises mixtures fitted to the BLOSUM retrieval pool's per-residue distribution** — the candidate
space is *derived from* the pool. **Re-running them to measure decorrelation would spend compute
establishing a foregone answer.**

> **The live version of the question belongs to WORKSTREAM A's corpus-level sources**, which differ in
> their *corpus* rather than in a fit to the current pool. R1 is closed as proposed; its scientific
> content survives in A's lane, where it was already the primary.

*Two lanes each stated a fact about the same file and they were incompatible; the file settled it in
one read. Recorded because the alternative — picking the more senior lane, or averaging the two
claims — would have produced a wrong plan either way.*

---

## L11 — THE FRECHET MEAN RECOVERS NOTHING: the projection tax is NOT the contraction (2026-09-07, coordinator)

`s20/frechet.py`, **n = 126, `complete: true`**. Pre-registered falsifier fixed before the run.

| arm | basis | RMSD | mean bond | vs incumbent [95% CI] | W/L |
|---|---|---|---|---|---|
| avg | **POINT CLOUD** | 3.0483 | **2.961** | −0.1569 [−0.1908, −0.1250] | 108/18 |
| **proj_lam0 (INCUMBENT)** | built | **3.2052** | 3.804 | +0.0000 | — |
| proj_lam03 | built | 3.2126 | 3.804 | +0.0074 [−0.0039, +0.0191] | 59/67 |
| **frechet** | built | **3.2100** | 3.804 | **+0.0048 [−0.0117, +0.0218]** | **61/65** |
| frechet_w | built | 3.2214 | 3.804 | +0.0162 [−0.0476, +0.0778] | 61/65 |
| frechet_med *(CONTROL)* | built | 3.2163 | 3.804 | +0.0111 [−0.0080, +0.0309] | 55/71 |
| torsmean *(NULL)* | built | 4.0725 | 3.804 | +0.8673 [+0.5426, +1.2045] | 48/78 |

**THE PRE-REGISTERED FALSIFIER FIRES.** `frechet − proj_lam0 = +0.0048 [−0.0117, +0.0218]`, 61W/65L —
**a dead null, an order of magnitude below the 0.084 A MDE**, and the point estimate is on the wrong
side. **The manifold-constrained Frechet mean recovers NONE of the +0.166 A projection tax.**

**THE HYPOTHESIS IS REFUTED, AND CLEANLY.** I argued the tax was the cost of un-contracting a
shrinkage estimator, so an operator whose feasible set *is* the manifold — which **cannot contract by
construction**, verified: mean bond exactly **3.804 A** on every Frechet arm — would never incur it.
**It incurs the same tax.**

> ### The projection tax is not the price of un-contracting. It is the price of the manifold itself.
>
> **Sprint 18's "irreducible" stands, and now on much stronger evidence** than when it was recorded:
> an operator built specifically to avoid the contraction pays the identical cost.

**THE STRONGEST FORM OF THE NEGATIVE, and it is what makes this worth having.** Four *different*
manifold-constrained operators — the lam=0 projection, the lam=0.3 projection, the Frechet mean from
the projection start, the Frechet mean from the medoid start, and the score-weighted Frechet mean —
**all land within 0.016 A of one another**, while the unconstrained point cloud sits 0.157 A below
them and the torsion-space mean sits 0.867 A above.

> **The manifold constraint determines the answer; the objective used to reach it does not.** Given
> this candidate set, ~0.166 A is the information cost of requiring ideal backbone geometry — not an
> artefact of how the requirement is imposed.

**Controls behaved.** The recorded-bad operator `torsmean` is 4.07 A, so the null is doing its job and
`frechet` beats it by **−0.862 [−1.200, −0.545]** — the Frechet criterion is real, it simply buys
nothing over the projection. **The medoid trap did not trigger**: `dist(frechet, medoid)` = 1.199
against `dist(proj, medoid)` = 1.097, so the Frechet arm is if anything *further* from the medoid —
it has not quietly become medoid selection.

**Coordinator's pre-registered expectation, for the record**: *"a real but partial recovery... 3.10-3.15
would support the mechanism without closing the tax."* **The measured value is 3.2100. The expectation
was wrong and the falsifier was right to be fixed in advance.**

**What this closes**: the largest measured non-predictor loss in the pipeline is **not addressable by
changing the averaging operator.** Any future attempt to recover it must change the *candidate set* or
the *geometry requirement*, not the estimator.

---

## L12 — Q2 ANSWERED: AMBER's landscape difficulty is a STERIC SINGULARITY, and 72% of its apparent damage is MOVE SIZE (2026-09-07, WORKSTREAM C)

### The headline refuted itself, by a control the lane wrote before reading its own table

Minimising AMBER in continuous torsion space moves Ca-RMSD **+0.6198 [+0.3213, +1.0037] away from
the native, 4W/26L** — the cleanest-looking result of the sprint. `c_land_null.py`, **written while
the run was still going and before its table was read**, moves the **same per-coordinate RMS torus
distance from the same starts** along the geodesic toward another pool member — realisable,
native-free, zero information:

    matched null alone      +0.4438 [+0.3884, +0.4891]
    AMBER minus that null   +0.1761 [-0.1248, +0.5741]   NOT MEASURED

> **72% of the damage is the SIZE of the move, not its direction.** The pre-registered falsifier fires
> exactly as written.

**AMBER is not moving randomly** — it beats an isotropic move of the same size by **−0.378**. **It is
moving too far**, because it starts at `E_0 = +1.9e8 kcal/mol` and **its first steps are clash
relief.** Legacy's apparent advantage dissolves identically: its move is **5.6x smaller** and costs
+0.0045 [−0.011, +0.019]. And `Legacy − AMBER` is **−0.4657 [−1.0178, +0.1427], 24W/6L** — **NOT
MEASURED despite the lopsided win/loss**, which the lane states rather than leaning on the count.

### H-C2(i) is REFUTED on its own axis — and what replaces it is better

The lane predicted AMBER would carry **more** negative curvature. **It carries slightly less**:
−0.0346 [−0.0707, +0.0116]. What actually distinguishes the two landscapes is **where the curvature
lives**:

    participation ratio   AMBER 0.0745   vs  Legacy 0.4221    30W/0L
    near-zero modes       44%
    anisotropy            18.8           vs  5.8              0W/30L
    condition number      three orders larger

**~7% of modes carry all of AMBER's curvature.**

> ### That is what a steric singularity looks like — so Q1's collapse census and Q2's spectra are ONE finding, not two.
>
> The 53.5%-of-the-pool catastrophe (L8) and the concentrated Hessian are the same object seen from
> two instruments. **This is the answer to "what mathematical property makes AMBER harder": not more
> negative curvature — CURVATURE CONCENTRATION at a steric singularity.**

### A scale-free result that reappears across instruments

    cos(grad E_Legacy, grad E_AMBER) = -0.3512 [-0.4719, -0.2221]   27 of 30 targets negative

**Q1's rank anti-correlation (−0.0886) reappearing in the differential geometry, on a different
instrument and a different statistic.** And **Legacy's structure-invisible channel does not explain
it**: that channel is **2.84%** of the gradient norm, and removing it moves the cosine from −0.3512
to **−0.3511**. *A clean null on the lane's own most quotable claim.*

### F-C2 does not fire — but on ONE convincing metric out of sixteen

AMBER's **`spec_skew` predicts final RMSD at rho = −0.487 [−0.709, −0.167], permutation p = 0.009**:
*the targets whose start is not dominated by a steric singularity are the ones where minimising AMBER
does not destroy the fold.* Legacy's `part_ratio` (−0.386, CI touching zero) sits at the n=30 power
edge and **is not claimed**.

**And the two strongest-looking rows are refused.** `n_distinct` and `barrier_rel` are
**outcome-side** — computed from the same minimisations whose RMSD is the endpoint — so the lane
**split its table into predictive and near-circular classes and declined to count them.** Everything
else is NOT MEASURED, so §3.2's enormous raw differences are **explicitly not to be over-interpreted.**

**Barriers are reported as a median and labelled a diagnostic**: a straight line in torsion space
passes through steric overlap, so the figure is an `r^-12` upper bound, **not a barrier height.**

### Still running

Q1 re-specified through the **shared `Relax` operator** (the coordinator's L6 primary) at 3/40 — and
**already reproducing L6 on the physics lane's own domain**: `collapsed = 63-75 of 75` per target and
the 50-step cap firing **8/8**. Then the repair-path experiment, then the Q3 Pareto merge.
`s20/agentC_FINDINGS.md` carries a live status table and the exact commands to finish each block.

---

## L13 — WORKSTREAM A: corpus change does NOT decorrelate, and the carrier is the SELECTOR, not the generator (2026-09-07)

All artefacts complete at n = 126, pre-registration unedited (timestamped 00:44), benchmark manifest
never opened.

### Both pre-registered falsifiers FIRE — R1 is closed

    F-A1a   rho(e_pep, e_pool)                 0.968 (cloud) / 0.923 (built)   bar 0.85
    F-A1b   corpus partition vs random ceiling 0.963 / 0.932                   bar 0.85

**Changing the corpus does not decorrelate the harmful component.** On the built basis, best of 13
arms is `fuse(pool,pep)` at **3.182 vs the incumbent's 3.205 — −0.023, a quarter of the MDE, 63W/63L,
median +0.001.** The peptide corpus alone: −0.017 [−0.096, +0.056]. **Brief candidates 1, 2 and 3 are
closed on outcome**, and `tors` — fitted to *nothing*, hardcoded literature Ramachandran constants,
the clean pool-independent case the coordinator asked for — aligns at **0.898** and lands +0.471
worse.

> **The coordinator's R1 (peptide-first corpus) is dead**, and it died on the falsifier written for
> it rather than on a judgement call.

### The reframing, and it moves the wall

A 4x2 corpus x selector cross at n = 126, **point-cloud basis throughout** (`bond` ~2.85-2.96 on
every cell, confirming the basis):

    corpus swap, selector held    blosum.disto 3.0483  ->  univ.disto 3.1571   +0.109
    SELECTOR swap, corpus held    blosum.disto 3.0483  ->  blosum.rand 3.4354  +0.387

> ### Swapping the corpus to UNIFORMLY RANDOM WINDOWS costs +0.109. Swapping the selector to random costs +0.387 — three and a half times more.
>
> **The carrier of the shared bias is the DISTOGRAM-AS-SELECTOR, not the retrieval pool.** Alignment
> agrees: corpus swap holds rho at 0.951-0.968, selector swap drops it to **0.794**.

**And it is not "any score".** A zero-information helix gate aligns 0.819 raw but only **0.268
partialled** (+0.991 A). **WORKSTREAM A states explicitly that this does not contradict Sprint 19's
C4/Z7** — those measure gate *damage*, this measures error *direction*; the helix gate is
simultaneously the most damaging and the least pool-aligned — **rather than quoting the two against
each other.**

### COORDINATOR CORRECTION TO WORKSTREAM A's WRITE-UP — a basis mismatch, the sprint's own error class

A reported *"`univ.disto` ... landing at 3.157 A"* beside the **built** incumbent. Checked in
`a_sel.json`: that cross is entirely **point-cloud basis** (`bond|blosum.disto` = 2.9614), so the
matched comparator is **3.048, not 3.205**. Read correctly, `univ.disto` is **+0.109 WORSE** than the
incumbent point cloud, not 0.048 better. **The scientific finding is unaffected and arguably
sharpened** — the selector/corpus contrast above is computed within one basis — but the number must
not be quoted across bases. *Fourth instance of `rmsd-reporting-basis-mismatch` this sprint, and the
first committed by a lane after the rule was adopted.*

### The lane's own damage, pre-committed before the number existed

On a **stricter post-hoc statistic** — alignment with the zero-information "typical peptide" mode
partialled out — the corpus partition decorrelates **materially** more than a random one: 0.679 vs
0.870, share **0.780 against the 0.85 bar**, **−0.1917 [−0.2313, −0.1564], 115W/11L, 5/5 folds.**

> **The corpus carries a real minority of the shared direction, and it converts to nothing.**
> **A recorded this commitment in writing while the run was at partial n, specifically so it could not
> choose afterwards.**

### Three independent reasons decorrelation cannot convert here

1. **Target level** — fusion gain is explained by the partner's **quality** at Spearman **0.96-0.99**
   and by its **decorrelation** at **|rho| <= 0.29 with inconsistent sign.** The mechanism dissolves.
2. **ORACLE ceiling** — with a per-target ORACLE fusion weight the peptide corpus is
   **indistinguishable from a constant alpha-helix** (+0.023 [−0.039, +0.104]); so is `tors`. All-8
   fusion is **+0.136 worse.**
3. **Diversity axes are coupled at +0.992** across sources — **nothing bought error diversity without
   buying geometric diversity**, and everything that bought geometric diversity paid in quality.
   `pep` and `prot` are *disjoint corpora* whose members' errors align at **0.673**, against **0.697
   inside the pool itself.**

**And the realisability axis moves against decorrelation**: the repair tax rises monotonically with
population diversity, **+0.157 (incumbent, matching the coordinator's 0.156) -> +0.197 (most
diverse).**

### Gates and a new form of an old defect

Top-75 reconstruction **set-identical to `shipped_record["sub"]` on 126/126**; projected `pool`
reproduces the shipped `fit_ca` to **max |delta| = 0.0004 A**; the 22.2% contraction reproduced
independently.

**A `.COMPLETE` flag that fired on the WRONG CONDITION** — `a_proj` wrote it when rows matched the
*subset* it was called with, so a 2-target smoke left a valid-looking flag and a wait-loop reported
"finished" at 8/126. Fixed to require the full instrument. **A new form of Z6: not a gate that never
fired, but one that fired on the wrong condition.**

### Handoff

> **If Q2 is reopened, the object to attack is the DISTOGRAM-AS-SELECTOR, not the pool** — and report
> the **helix-partialled** alignment, not the raw one, since a constant helix already sits at
> 0.767-0.814 raw.


---

## L-B — Q1 ANSWERED, AND IT CANNOT MOVE PRIORITY 1 (2026-09-07, WORKSTREAM B)

Pre-registration `s20/PREREG_B.md` (falsifier first, unedited). Code `s20/qb2_*.py`.
Artefacts `s20/results/qb2_*` — landscape n=20, optimiser battery n=20, relaxation arm n=10,
CVaR α n=10, encoding n=10, Q2 n=20 — all with `_COMPLETE` flags. Findings
`s20/agentB_FINDINGS.md`.

### The result that prices the whole question

Across **36 (objective, arm) cells** — Legacy, bare AMBER, log-conditioned AMBER and the deployed
distogram objective × 9 optimisers, each against the mandatory `best_of_N` at a matched 512-eval
budget:

    cells improving objective AND Ca-RMSD together, significantly:   0
    cells significantly TRADING one for the other:                   8
    within-cell target-level rho(objective gained, RMSD gained):  +0.0034 [-0.0679, +0.0717]
                                                                  positive in 20/36 cells

**Zero for thirty-six.** `search-saturates-discrimination-binds` on a **sixth** instrument.
— **ESTABLISHED**, n=20.

**And the cleanest statement of why:** over the same K=500 retrieval windows the pipeline uses,
the RMSD of each objective's own argmin is **AMBER 4.990 · Legacy 5.487 · distogram 3.676**
against a **pool mean of 4.739**. *Picking the lowest-energy structure out of 500 real windows is
worse than picking one at random, for both physics energies.* ρ(E, RMSD): AMBER **−0.071
[−0.173, +0.042]**, Legacy +0.191, distogram **+0.497**. — **ESTABLISHED**, n=20.

### Q1 itself

**AMBER is harder on every scale-free axis, and about half of it is a units artefact.** `AMBc =
sign(E)·log1p(|E|)` is a **strictly monotone** transform of the same call — identical argmin,
identical ranking (**EXACT**, verified: the two have byte-identical rank statistics). It collapses
the dynamic range **16.39 → 1.57 decades** (Legacy 1.61) and the standardised gradient norm
**1.6e13 → 26.9**, and is worth **≈1.8 IQR of optimiser progress** (on raw AMBER *no* arm beats
`best_of_N`; Adam ends **+4.469** robust units *above* it. On `AMBc` every arm reaches ≈ −1.9).
It leaves condition number (+0.115 ns), anisotropy (+0.202 ns) and local-minimum count
(+0.017 ns) **unchanged** — and `AMBc` is still **473× worse conditioned than Legacy**, with 4×
the near-zero modes and **twice the local minima per 2π**. **Two pathologies, one remedy.**

**Eliminated:** *optimizer mismatch* — SPSA/central-FD gradient cosine is **+0.15 to +0.18 and
statistically identical on all four objectives**. **REFUTED.**
**Relocated:** *CVaR concentration* belongs to the **rule**, not the Hamiltonian — ESS is a
functional of the tail's order statistics alone (**EXACT**). On raw AMBER at α=1 the
score-function gradient norm is **9.3e14** and ESS collapses to 0.090; under the log transform it
is 1.7 at every α. **CVaR on raw AMBER is a crude discontinuous conditioner.**
**Real:** *encoding* — θ vs (cos θ, sin θ), physically identical, wins on Cα-RMSD in **11/12
cells (sign test p = 0.006)**, up to **−0.649 Å**.

**One landscape metric does predict outcome.** After partialling target difficulty on ranks
(L4's binding requirement), `frac_neg` — the negative-curvature fraction at a random plausible
start — survives on **all four objectives independently**: +0.528 / +0.546 / +0.510 / +0.511,
every CI excluding zero, against 2.6 flags expected by chance in 52 tests. **Landscape geometry
is diagnostic of the TARGET and is not actionable through the OPTIMISER.** — **SUPPORTED**, n=20.

### The coordinator's relaxation confound: CORRECT in general, and NOT applicable to this lane

`H_AMBER = E ∘ Relax_50` and the bare single point agree on only **52 % of rank variance
(ρ = 0.720 [0.656, 0.778])** — **different objectives, not two views of one landscape. F-D2 does
NOT fire.** `Relax_1` is the borderline conditioner (ρ = 0.950 [0.924, 0.973]). The relaxation
compresses the range by **nine decades** and **triples the local-minimum count (5.6 → 15.2 per
2π)** — the most multimodal objective anywhere in this study. **The collapse guard, as a measured
count: `Relax_1` returns `+inf` on 1109/3302 calls (33.6 %)**; `Relax_50` on 24/3302.
**This lane's Q1 arm is the bare `core.amber._run(steps<0)` single point**, bit-exact against
`refine_coords(k=0, steps=−1)` at **0.00e+00 over 150 comparisons**, and `_run` **never calls
`_is_collapsed`** (verified in source) — so it is un-confounded by construction.

### Q2 — CLOSED

`ρ(e_VQE, e_pool) = +0.777 [+0.672, +0.870]`, against the VQE's own same-seed ceiling of +0.852:
**91.3 % of the VQE's reproducible error structure is shared with the retrieval pool.** It is
**statistically indistinguishable from a zero-information constant helix's (−0.016 ns)**, from
best-of-N's (−0.025 ns) and from a classical thermostat's (+0.049 ns). **The brief's absolute
0.75 falsifier line is MIS-SPECIFIED for generators** — the helix null already sits at +0.793 —
declared, and replaced by the null-referenced test, which returns a clean null. **The generator
role closes.** — n=20.

### Two recovered artefacts, and a methodological number every future lane needs

`s19/results/qb_ext.json` (40/40 targets) was **complete on disk and read by no one**; recovered
in `s20/qb2_s19ext.py`. **Ansatz-seed sensitivity is 0.200 Å within-target sd (untrained 0.289)
— 2.4× the instrument's 0.084 Å MDE. Any future variational arm reporting fewer than 4 seeds is
NOT MEASURED.** It also independently replicates L4: the only extension arm with a CI excluding
zero is `x_mixunif50` (**−0.109 [−0.204, −0.021], 5/5 folds**), which *dilutes* the CVaR tail
toward the plain mean, while the deepest tail (α=0.02) is the only arm on the wrong side of zero
— the regime **arXiv:2605.02850** predicts. A **deeper** ansatz buys −0.020, 2/5 folds.

**Self-caught error, recorded:** I first read `x_warm`'s 1.5187 Å generation ceiling as beating
`pool500`'s 1.7795. It is a **min-of-N artefact** (8192 vs 500); at matched count 500 it is
1.9400 and the pool wins, and two sibling arms lose with CIs excluding zero and 5/5 folds.

### Labels

    Q1 "AMBER is harder"                                    ESTABLISHED on 4 scale-free axes (n=20)
    "conditioning explains it"  (my own H1)                 HALF REFUTED -- range yes, curvature no
    "optimizer mismatch / SPSA gradient noise"              REFUTED (n=20)
    "CVaR concentration is a Hamiltonian property"          REFUTED -- it is the rule (EXACT)
    encoding is a real lever on RMSD                        SUPPORTED (n=10, p=0.006)
    frac_neg predicts final Ca-RMSD                         SUPPORTED (n=20, 4/4 objectives)
    optimisation progress predicts structure                REFUTED (0/36 cells, rho +0.003)
    H_AMBER = E o Relax is a different objective            ESTABLISHED (rho 0.720, n=10)
    Q2 decorrelation                                        CLOSED (n=20, null-referenced)
    Priority 2: the CVaR-VQE pillar                         GENUINE, TRAINS, contributes ~0

---

## L14 — WORKSTREAM B: Q1 separated experimentally, Q2 closed, and the whole optimisation question priced at ZERO for Priority 1 (2026-09-07)

`s20/PREREG_B.md` falsifier-first and unedited; **43 artefacts, all `_COMPLETE`**.

### The lane's own H1 is half refuted, by itself

It predicted the Legacy-AMBER difference is dominated by **conditioning**. **The dynamic-range half
is true and enormous** — conditioning collapses AMBER's range **16.39 -> 1.57 decades** (Legacy 1.61).
**The curvature half is false**: condition number, anisotropy and local-minimum count are
statistically **unchanged**, and conditioned AMBER is still **473x worse conditioned than Legacy with
twice the local minima.**

> **Two pathologies, one remedy** — the remedy fixes only one of them.

### AND THE WHOLE QUESTION IS PRICED AT ZERO FOR PRIORITY 1

Across **36 (objective, arm) cells** at matched 512-eval budget against the mandatory `best_of_N`:

    0 cells improve the objective AND Ca-RMSD together
    8 cells significantly TRADE one for the other
    within-cell target-level rho(optimisation gained, structure gained) = +0.003 [-0.068, +0.072]

**Sixth independent instrument for *search saturates, discrimination binds*.** And the reason, over
the same K=500 pool the pipeline actually uses:

    argmin RMSD   AMBER 4.990 | Legacy 5.487 | distogram 3.676      pool MEAN 4.739

> **Picking the lowest-energy structure out of 500 real windows is WORSE THAN PICKING AT RANDOM — for
> both physics energies.**

### Q1, separated experimentally as the brief demanded

* **Optimizer mismatch — REFUTED.** SPSA vs central-FD gradient cosine **+0.15-0.18**, statistically
  identical on all four objectives. *SPSA is not the culprit.*
* **CVaR concentration — RELOCATED, and EXACT.** ESS is a functional of the **tail's order statistics
  alone**, so it is a property of the **rule**, not the Hamiltonian. On raw AMBER at alpha=1 the
  gradient norm is **9.3e14** and ESS collapses to **0.090**; under a log transform it is **1.7 at
  every alpha**. **CVaR on raw AMBER is a crude discontinuous conditioner.**
* **Encoding IS a real lever.** `theta` vs `(cos theta, sin theta)` — *physically identical* — wins on
  RMSD in **11/12 cells (p = 0.006), up to −0.649 A.** The coordinate-artifact hypothesis in the brief
  is **supported**. *(Within the variational setting; the variational setting itself still loses to
  retrieval — see above.)*
* **One landscape metric predicts outcome**: `frac_neg` survives the difficulty partial on **all four
  objectives independently** (+0.51 to +0.55, all CIs excluding zero) against **2.6 flags expected by
  chance in 52 tests**. **Landscape geometry is diagnostic of the TARGET, not actionable through the
  optimiser.**

### The coordinator's relaxation confound: correct in general, NOT applicable to this lane's arm

**Both lanes are right, on different objects.** WORKSTREAM D measured the **deployed
`AmberHamiltonian`**, which relaxes; WORKSTREAM B built its Q1 arm as `core.amber._run(steps<0)` — **a
true single point, bit-exact at 0.00e+00 over 150 comparisons**, and verified in source that `_run`
**never calls `_is_collapsed`**. **So B's arm is un-confounded by construction and F-D2 does not fire
on it.**

B independently confirms the confound is real where it applies: `E o Relax_50` and the bare single
point agree on only **52% of rank variance (rho = 0.720)** — *different objectives*; relaxation
compresses the range **nine decades** and **triples** the local-minimum count (5.6 -> 15.2 per 2pi);
and **the collapse guard fires 1109/3302 times (33.6%) at Relax_1.**

### Q2 — CLOSED

    rho(e_VQE, e_pool) = +0.777   against the VQE's own ceiling of +0.852

> **91.3% of the VQE's reproducible error structure is the pool's.** Indistinguishable from a
> zero-information helix (−0.016 ns), from best-of-N (−0.025 ns), and from a thermostat (+0.049 ns).

**And the lane declared the coordinator's threshold mis-specified rather than passing it.** The
brief's 0.75 line is **wrong for generators** — the helix null already sits at **+0.793** — so B
replaced it with a **null-referenced** test and reported a clean null. *Declaring a bar unfit and
replacing it, before reporting against it, is the correct move.*

### A recovered artefact and a number every future variational arm needs

`s19/results/qb_ext.json` (40/40, **unread by anyone**):

> ### Ansatz-seed sensitivity is 0.200 A within-target sd — 2.4x the MDE.
> **Any variational arm with fewer than 4 seeds is NOT MEASURED.**

It independently replicates L4: the only arm with a CI excluding zero **dilutes the CVaR tail toward
the plain mean** (−0.109, 5/5 folds), and **the deepest tail is the only arm on the wrong side of
zero.**

### Self-caught error

The lane first read `x_warm`'s 1.5187 A generation ceiling as beating `pool500`'s 1.7795 — **a
min-of-N artefact** (8192 draws vs 500). At matched count it is **1.9400 and the pool wins.**

---

## L15 — THE AMBER PARETO, CLEAN RUN: F-C2 FIRES, and a scalar validity score would have crowned a broken structure (2026-09-07)

**n = 126, complete, four artefacts with config hashes.** `steps = 0` (the deployed protocol, **no
bound anywhere**), the **full** pre-registered ladder k = 1, 3, 10, 30, 100, 300, 1000, `caonly` at
all three declared strengths, **16 arms, ~2000 genuine ff14SB/GBn2 minimisations.** All three compute
deviations withdrawn.

### Axis 1 — accuracy against the incumbent k30

| arm | dCa vs k30 | n |
|---|---|---|
| `k1` | **+0.0903 [+0.0574, +0.1293]** | 123 |
| `k10` | +0.0238 [+0.0135, +0.0372] | 123 |
| `bbo_k30` / `stage_300_30` / `caonly_k30` | +0.011 / +0.008 / −0.005 | 123 |
| `k300` | −0.0737 [−0.0913, −0.0534] | 102 |
| `caonly_k100` | −0.0776 [−0.0925, −0.0666] | 123 |
| `blend50` | −0.0922 [−0.1115, −0.0793] | 123 |
| **`caonly_k300`** | **−0.1100 [−0.1355, −0.0924]** | 123 |
| `blend75` | −0.1193 [−0.1458, −0.1018] | 123 |

*Pairwise-gated: an all-arm gate would cost 64/126 targets to one divergent arm and shift the input
mean 3.050 -> 2.435, so it is printed but not read.*

### Axis 2 kills all three MDE-beaters

`blend75` clash 2.6 A **+11.63 [+9.67, +14.51]** · `blend50` **+7.49** · `caonly_k300` **cis_frac
+0.4242 [+0.3799, +0.4732]**. **Every arm whose validity matches k30 is within ±0.024 A of it**, and
three of those sit **below the rotated-frame floor and are therefore NOT MEASURED.** *Sprint 18's
irreducible repair tax reproduces by a wholly different route.*

> ### The sharpest item is methodological
> **`caonly_k300` is 0.110 A more accurate than the incumbent with a clash count STATISTICALLY
> INDISTINGUISHABLE from it (−0.008 [−0.050, +0.032]). A scalar clash-based validity score would have
> made it the sprint's headline.** Only `cis_frac` disqualifies it — and the **dose-response confirms
> the mechanism**: pinning Ca harder runs cis **0.010 -> 0.051 -> 0.101** at k = 30 -> 100 -> 300
> **while clashes stay at zero.** *Pin only Ca and the peptide plane flips.*

**Rotated-frame null, reported with its maximum as required**: n = 15, **max |dCa| = 0.014192 A**,
mean 0.001726, max |dE| 4.769 kcal/mol. **Convergence gate**: `1D6X 2NB7 7BX2` for the **sixth
sprint**, with exclusions rising **3 -> 24 -> 64** at k = 100 -> 300 -> 1000 — *a pinned backbone
cannot relieve its own strain.* **AMBER k30 as repair reproduces**: clash 2.6 A **1.113 -> 0.000**,
bond strain **0.131 -> 0.011**.

### The soft rungs that were wrongly dropped are now measured, and they close the open item

**k = 1 is +0.090 A worse than the incumbent, and carries MORE Ramachandran outliers** (+0.0340
[+0.0170, +0.0517]). **There is no better Pareto point at the soft end.** So the unjustified deviation
**cost the sprint nothing scientifically — and it is still recorded as having been made on a false
premise.**

### Two of the lane's own claims retracted in its findings

`core.amber` has **no non-termination defect** (12-15 s on a quiet box; it was CPU starvation), and
**GC3's inertness certificates were worthless both times** because they were computed on calls the
bound never touched — *the one call it did bound moved 0.1595 A.* The misleading diagnostic log is
renamed **`_DIAG_1CEK_STARVATION_ARTEFACT_not_a_result.log`**, and the capped run is preserved as
`_SUPERSEDED_*`, quoted nowhere.

**ALL LANES ARE NOW COMPLETE.** Q1 (mechanism), Q2 (frontier) and Q3 (closed) delivered at n = 126.

