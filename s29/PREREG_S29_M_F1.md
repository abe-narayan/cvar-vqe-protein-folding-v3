# PREREG S29-M-F1 -- SWAP THE SELECTION FUNCTIONAL: L1 BAYES RISK -> LOG SCORE OF THE SAME POSTERIOR
(lane M, written 2026-09-20 00:1x, BEFORE any F1 number exists; appended-to only, never edited
after a result exists -- contract rule 7)

## 0. The assignment and the motivating finding

Coordinator, `s29/STATE.md` integration note 3, after lane D's S29-L6: lane X's pair log-score is
the first cost in the record that is not ANTI-informative on the near-native ladder (ladder rho
+0.018 [-0.069, +0.123] against the shipped cost's -0.182 [-0.308, -0.053]; paired difference
+0.200 at 1.45x MDE, 5/5 folds, power 0.98, `s29/results/s29_D_cost_audit_X_cost_nll_ca.json`).
Lane D's mechanism: the shipped score is a **Bayes risk** -- an expected L1 distance error under a
posterior that is about 2x over-confident (S25 L1/L2) -- whose minimiser is a **contracted**
structure; the log score is a **proper scoring rule** of the **same** posterior, whose minimiser is
not driven to contract. Same information, different functional.

F1 asks whether that functional difference reaches the **built chain**, by changing ONE thing in the
shipped pipeline: the functional at the selection stage. It is the first test of an item on lane M's
own convenience list (`s29/CONVENIENCE_CHOICES.md` C9, "L1 Bayes risk (rather than log-likelihood,
L2, or a proper scoring rule)").

## 1. RULE 10 DECLARATION: this question has been asked before, twice, and both answers were negative

Stated up front, with the sprint and the line, as contract rule 10 requires.

**(a) The log-likelihood form was measured in S5-S7 and rejected.** `distogram.py:10-18` (the
reference module's own docstring): "The predicted distribution over bins is used through its L1
Bayes risk ... **not** through the log-likelihood of the occupied bin. The two were measured side by
side: the log-likelihood form scores +0.02 to -0.42 ... while the Bayes-risk form scores +0.16 to
+0.82 on the same predictions." Same statement at `core/predict.py:13-17`.

**(b) And it was measured on THIS instrument, at n = 126.** `docs/FINDINGS.md` (the S7 debias arm
table, artefact `s7/debias_tune.json`), arm `pmi lam=0 (pure NLL)`:

    arm                     selected   vs base    95% CI            W/L     in-band rho
    base (shipped L1)         3.454     +0.000    --                --        +0.126
    pmi lam=0 (pure NLL)      3.547     +0.093    [-0.055, +0.241]  49/40   **+0.182**

> "The pure-NLL arm has the **best in-band Spearman of any arm (+0.182 against the baseline's
> +0.126) and a WORSE selected RMSD (3.547)** -- another instance of in-band rho and the decision
> metric disagreeing."

**That is the same shape as S29-L6**: the log functional orders the near-native band better and
selects worse. It is also the project's standing law `better-matrix-worse-ranking`.

**(c) A third, adjacent closure.** S25 L12 swept the POWER family on the same 126 targets through
the same top-75 readout: q = 2 (an L2 Bayes risk, minimiser = the posterior MEAN rather than the
median) is **+0.0176, 0.96x MDE** under nested CV, and the best parameter chosen with **complete
leakage** is the shipped identity. A functional change of a neighbouring kind is measured at ~0.

### Why the old closure may not apply (the three reasons, and they are weak)

1. **The readout is different, and the record says the readout decides.** (b) is an **argmin**
   selection (baseline 3.454 = the shipped argmin on the s8 instrument). F1's endpoint is the
   **built chain through a 75-member uniform coordinate average**. The project's own law
   `operator-consumes-set-mean` (d_out = 1.16 x set mean + 0.04 x set best, R2 0.89) says these two
   readouts consume different statistics of the retained set: an arm that is worse at picking ONE
   candidate can be better at picking a SET whose mean is better. S17 L12 is the same point measured
   on the K axis (consensus improves with K where argmin degrades).
2. **The mechanism is now named and is measurable as an intermediate.** In S7 the NLL arm was a
   ranking variant among 28; nobody measured **contraction**. F1 measures the emitted cloud's mean
   virtual CA-CA bond (production 2.9614 A against a native 3.8122, `ARCHITECTURE.md` section 0) as
   the intermediate quantity (contract rule 18), so the arm is informative about the mechanism even
   if the endpoint is null.
3. **The posterior is not the S7 posterior.** (b) predates the fold pinning and the corrected
   clustering; (a) predates the ESM features. The functional is being re-tested against the
   posterior that ships.

**Registered prior, stated honestly and before the run:** the record predicts **null to worse**,
most likely **+0.00 to +0.10 A** on the built chain. I am running it because it is cheap (the pool
and the posterior are cached), because it is the coordinator's queued assignment with a named
mechanism, and because the intermediate (contraction) has never been measured for this functional.
**If it comes out positive I will be surprised, and lane D should attack it the same hour.**

## 2. The arms. Exactly one thing changes.

Everything is held at the production setting: the shipped K = 500 BLOSUM pool
(`H.Candidates.from_universe(pdb, 500)`), the same leave-fold-out 17-bin posterior
(`I.distogram`), the same tie-safe top-M with M = 75 (`np.lexsort((tiekey, E))`), the same uniform
coordinate average in the retained set's medoid frame (`H.readout_uniform`), the same production
projection (`H.readout_projected` -> `I.project`, ramah @ 0.3, multi-start, maxiter 300), the same
RMSD (`I.ca_rmsd` against the pinned `nat_ca`). No parameter is fitted anywhere; there is nothing to
leave out of fold beyond the posterior, which is already leave-fold-out.

| arm | score of candidate w, lower is better | what it changes |
|---|---|---|
| **PROD** (the comparator) | `mean_p w_p * sum_c prob[p,c] * abs(d_p(w) - C_c)`, `w_p = 1/(sd_p+0.5)` normalised to mean 1 | -- (the shipped score; reproduces `chain_rows.jsonl :: DIS` bit-for-bit, audit check 5) |
| **F1a LOG** (primary) | `sum_p -log(prob[p, bin(d_p(w))] + 1e-4)` | the functional, exactly as metered by lane D (`s29/s29_X_config.py:728 cost_nll`, `EPS_P = 1e-4`) |
| **F1b LOGW** (secondary) | `mean_p w_p * (-log(prob[p, bin(d_p(w))] + 1e-4))` | the functional ONLY -- keeps the shipped per-pair weight, so F1a minus F1b prices the weight separately |
| **C1 L2RISK** (matched functional control) | `mean_p w_p * sum_c prob[p,c] * (d_p(w) - C_c)^2` | a DIFFERENT functional of the SAME posterior that is NOT a proper scoring rule (minimiser = posterior mean). Answers: "is it the log functional specifically, or does any re-aimed functional of this posterior move the chain?" |
| **C2 LOGPERM** (information control) | F1a with the per-pair log-probability ROWS permuted across pairs by a stable per-target RNG (`rng_for(pdb,"F1perm")`) | preserves the functional's shape and its marginal, destroys the pair correspondence. If C2 moves the chain as much as F1a, the gain is not the pairs' information. |
| **C3 MONO** (stated, not run as a measurement) | any strictly monotone transform of PROD | **inert by algebra**: the top-75 SET depends only on the order, so a monotone re-ranking of the shipped score is the identity arm. This is asserted numerically once as a harness check, never reported as a control's "null". |

**Note on the coordinator's suggested control (a).** The brief proposed "a monotone re-ranking of
the shipped score, e.g. its own rank ladder passed through the log-score's empirical quantiles".
Through a top-m selection readout that control is the identity **by construction** -- it cannot move
a single member -- so reporting it as a null would be reporting an algebraic identity as a
measurement (the trap `s25/QUANTUM.md` section 5 names: "a 100% pass rate on an algebraic identity
carries no information"). It is therefore replaced by **C1**, which changes the functional's aim
while keeping its information, and **C2**, which keeps the functional and destroys the information.
This is a declared deviation from the coordinator's wording, made before any number existed.

## 3. The endpoint, the statistics, and the falsifier

**Endpoint:** mean built-chain CA-RMSD, n = 126, paired per target against PROD, `ST.compare` with
`ST.pinned_folds`, decided on the **fold-clustered CI**. MDE = 2.8016 x SE of the paired difference,
reported per comparison. Point cloud is reported beside it as the intermediate, never as the result.

**PRIMARY FALSIFIER (pre-registered).** F1a is a result only if, on the **built chain at n = 126**:
`effect < 0` (better than PROD) AND `|effect| > 1.0x its own MDE` AND the fold-clustered CI excludes
zero AND at least 4 of 5 folds agree in sign. **Anything else is a null and will be written as one.**
If `effect > 0` with the fold CI above zero, the functional is measurably WORSE and C9 in
`s29/CONVENIENCE_CHOICES.md` is closed in the direction the record predicted.

**PROBE GATE (12 targets first, contract rule 16).** The probe set is pinned and pre-existing:
`s29/s29_X_config.py:120 PROBE_TARGETS` = 1A13, 1I6Y, 1M02, 2BFI, 2LWS, 2MP9, 2P5H, 5Z5W, 6MBM,
7JGX, 8HVS, 9KAR (`P.targets()[::11][:12]`, folds 4,2,0,2,4,0,4,3,1,0,4,3; 2 of them in FAIL18).
**The probe is never quoted as evidence for the instrument.** Gate: go to 126 unless the 12-target
built-chain mean is worse than PROD by more than +0.30 A, in which case F1 stops at the probe and is
reported as a probe-stage null (with the number).

**MECHANISM, measured beside the outcome (contract rule 18).** For every arm and every target:
mean virtual CA-CA bond of the emitted point cloud (production 2.9614 A, native 3.8122 A), the
cloud's radius of gyration against the native's, and the projection price (chain minus cloud). The
registered mechanism prediction: **if lane D's contraction story is right, F1a's cloud is
measurably LESS contracted than PROD's** (bond closer to 3.81) -- and that prediction can be
confirmed while the endpoint is null, which would be the informative outcome.

**PARALLEL-BIAS CHECK (S24 L2/L3).** Per target, after optimal superposition on the native:
`e_PROD = emitted_PROD - nat`, `e_F1 = emitted_F1a - nat`; report `cos(e_PROD, e_F1)` and the
paired norms. Near 1 means F1a emits the same answer with the same error direction (a different
route to the same place); low means it is a genuinely different answer. A random reference at
3n-6 dof (|cos| ~ 0.14) is printed beside it.

**FAIL18 / 108 SPLIT with a random-18 null.** The split is reported for the primary arm only, and
the FAIL18 effect is compared against the distribution of the same statistic over 2,000 random
18-target subsets (`rng_for("F1","fail18null")`), because an 18-target mean has a large SE and the
project has been burnt by reading one (`the median-vs-mean gap is the free early warning`).
FAIL18 is an ORACLE partition (`s29/CONVENIENCE_CHOICES.md` C8) and is labelled so wherever it is
used; it selects nothing.

**MULTIPLICITY.** F1 runs **5 endpoint comparisons at n = 126** (F1a, F1b, C1, C2 against PROD on
the built chain, plus F1a minus F1b), and the primary is **F1a vs PROD, pre-specified here**. The
other four are controls and decompositions, reported with their own MDEs; the max-over-K null is
stated in the ledger entry if any of them clears its MDE while F1a does not.

**REPLICATION.** No seed enters any arm (the scores are deterministic; the only RNG is the tie key,
which is stable per target). A positive is therefore replicated by (i) the reversed fold order in
`ST.compare`, (ii) recomputing the posterior cold rather than from the s12 cache (audit check 8:
the cached order differs from a cold one on 2/126), and (iii) lane D's independent attack.

## 4. What would make me withdraw the arm before reporting it

- If F1a's top-75 sets are identical to PROD's on most targets, the arm is not a different
  operator and its endpoint difference is noise: report the set overlap first (registered
  diagnostic: mean |F1a top-75 AND PROD top-75| / 75).
- If the log score is degenerate on this pool -- e.g. saturating at the `1e-4` floor for most
  candidates, so that the ranking is dominated by a count of "impossible" pairs rather than by the
  posterior -- report the floor-hit rate per target and treat the arm as a different object from the
  one lane D metered.
- If `cost_nll`'s discreteness produces large tie sets at the m = 75 cut, report the tie fraction;
  ties are broken by the stable random key, never by array order (contract rule 12).

## 5. Cost and compute

The pool, the channels and the posterior are cached, so each target is one `channels_for` read, four
score evaluations over 500 candidates (vectorised gathers), four coordinate averages and four
projections. Measured cost of the projection in the audit: ~3 s per target per arm. Probe (12
targets x 5 arms) ~ 3 min; full run (126 x 5) ~ 30 min, one process, peak RSS ~0.6 GB by the audit's
measurement. Both go through `python s26/jobrun.py --agent S29M --tag CPU`, checkpointed per target
to `s29/results/s29_M_F1_rows.jsonl` and resumable.

---

# ADDENDUM 1 (2026-09-20 00:2x, lane M) -- THE MECHANISM CLAIM IS WITHDRAWN BEFORE THE FULL RUN; THE ARM AND THE FALSIFIER STAND; ONE CONTROL IS ADDED

Appended, not edited (contract rule 7). Written after the 12-target probe's four arms existed and
BEFORE the 126-target run, and before the LOGPERM control was repaired. The primary falsifier in
section 3 is unchanged.

## A1.1 What the coordinator corrected

Two lanes corrected the framing this prereg inherited from `s29/STATE.md` integration note 3:

- **S29-L10 (lane D):** the shipped cost's own descent direction **EXPANDS** the structure
  (bond x1.0438, Rg x1.0248; only 18 of 126 contract). There is no contraction story at the
  gradient level.
- **S29-L12 (lane L):** the Bayes estimator of an over-confident posterior is **not contracted**:
  a symmetric width error moves neither the mean nor the median. The emitted cloud's 22 to 26%
  contraction is **Jensen's inequality on the coordinate average** (`||E X - E Y|| <= E||X - Y||`,
  an arithmetic property of averaging real windows), not a property of the posterior or of the L1
  functional.

**Consequence for this prereg.** Section 0's sentence "whose minimiser is a **contracted**
structure" is **withdrawn**. F1 is no longer "remove the contraction at the selection stage". What
F1 tests is exactly this and nothing more:

> **The two functionals order the near-native half of the pool differently (+0.200 of ladder rho,
> 5/5 folds, S29-L6). Does that ordering difference survive to the built chain?**

The arms, the endpoint, the statistics and the PRIMARY FALSIFIER in section 3 are unchanged: they
never depended on the mechanism story. The registered prior is unchanged and is now better
motivated: **null to worse**, because the terminal operator consumes the retained set's MEAN
(`operator-consumes-set-mean`: d_out = 1.16 x set mean + 0.04 x set best, R2 0.89) and a re-ordering
that keeps most of the same 75 members moves that mean very little.

## A1.2 The contraction measurement STAYS, demoted to a diagnostic -- and the probe already speaks to it

The bond/Rg measurement registered in section 3 is kept, no longer as "the mechanism" but as a
diagnostic of what the functional swap does to the emitted geometry. It is worth keeping because
the 12-target probe already shows, independently of lanes D and L, that **contraction and RMSD do
not track across functionals** (`s29/results/s29_M_F1_rows.jsonl`, n = 12, PROBE -- not evidence for
the instrument):

    arm      built chain   point cloud   mean cloud bond (native 3.8122)
    PROD       3.3816        3.2529        2.9874
    LOG        3.6161        3.4711        3.1063      contracts LESS, and is WORSE
    L2RISK     3.5632        3.3633        2.8675      contracts MORE, and is WORSE

Two functionals move the contraction in **opposite** directions and both move the endpoint the
**same** way. On this evidence the contraction is not the lever, which agrees with S29-L10 and
S29-L12 from a third direction.

## A1.3 The control the corrected framing needs, and why it cannot be the one the coordinator named

The coordinator asks that the matched control -- "a monotone re-ranking of the shipped score, which
changes the functional's SHAPE but not its information" -- become the main comparison. **That
control is the identity arm and cannot be run as a measurement.** Through a top-m selection readout
the retained set depends only on the score's ORDER, and a strictly monotone map preserves the order
exactly; the emitted structure is therefore bit-identical to PROD's by algebra, not by measurement.
Reporting its null would be reporting an algebraic identity as evidence -- the failure mode
`s25/QUANTUM.md` section 5 exists to prevent. (This was already stated as C3 in section 2 and is
asserted once numerically as a harness check.)

What the corrected framing actually needs is a control matched on the thing that DOES vary: **how
many members of the retained set the re-ordering exchanges**. The probe measures LOG's top-75
overlap with PROD at **0.687**, i.e. LOG swaps about **23 of 75** members. So:

**SWAPCTL (added here; the main comparison for the corrected framing).** Per target, let
`k_swap = |PROD_top75 \ LOG_top75|`. Remove `k_swap` members of PROD's top-75 chosen uniformly at
random, and replace them with `k_swap` candidates drawn uniformly at random from PROD's ranks 76 to
500. Stable RNG `rng_for(pdb, "F1swap<d>")`, **4 independent draws**, each through the identical
uniform average and production projection; the arm's value is the mean over draws and the per-draw
sd is recorded. This is the **zero-information re-ordering at LOG's own exchange rate**, matched in
the operator's space.

**The comparison this licenses, pre-registered now:**
- If `LOG - PROD` is statistically indistinguishable from `SWAPCTL - PROD`, then the endpoint
  difference is the **price of exchanging members at that rate**, and the log functional's ordering
  information contributes nothing at the endpoint -- whatever its ladder rho.
- If `LOG - PROD` is **better** than `SWAPCTL - PROD` by more than 1.0x the MDE of that paired
  contrast with the fold CI excluding zero, the log functional's ordering carries endpoint-relevant
  information even though it does not beat PROD.
- If `LOG - PROD` is **worse** than `SWAPCTL - PROD`, the log functional is actively
  anti-informative at the endpoint -- it swaps out better members than chance would.

`LOG - SWAPCTL` is added to the multiplicity count (section 3): F1 now runs **6 endpoint
comparisons** at n = 126 (LOG, LOGW, L2RISK, LOGPERM, SWAPCTL against PROD, plus LOG - SWAPCTL);
the primary is still **LOG vs PROD**, pre-specified before any number existed.

## A1.4 The averaging bound, recorded

Lane L (S29-L12) derives from the Krogh-Vedelsby / Ueda-Nakano (1 - 1/M) coefficient that an
**infinite** pool of the same kind returns 3.040 A against the shipped 3.0483 -- "average more
members" is worth at most about **0.008 A**, with the aggregate-to-RMSD map about 3% loose. **F1
holds M = 75 fixed in every arm**, including SWAPCTL (which exchanges members, never changes the
count), so the bound does not bind on any F1 arm. It is recorded here so that no future variant of
F1 changes m without pricing it against that ceiling, and it is added to
`s29/CONVENIENCE_CHOICES.md` C21/C23.
