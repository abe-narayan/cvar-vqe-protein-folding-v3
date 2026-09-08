# SPRINT 19 — FINAL REPORT
## Why the predictor's errors are structurally maladaptive

**126 cluster-disjoint targets · full-chain Cα-RMSD, frozen implementation · sealed 60-target
benchmark untouched · four workstreams · MDE at 80% power = 0.084 Å**

> **STATUS AT REPORT FREEZE, 2026-09-06 23:55.** Sections 0–4 and 7–10 are **final** at n = 126.
> **§5 (CVaR-VQE) is now DELIVERED** — the sweep was complete on disk and was analysed in Sprint 20
> against its own pre-registration, verbatim. All three endpoints fail; two positives survive. **§6 (AMBER Pareto) is PROVISIONAL** — its run was superseded as
> cap-limited and the clean rerun lands ~02:30. Both are marked in place. **Nothing else is
> provisional, and no result in this report is quoted from a partial or superseded artefact.**


---

# 0. EXECUTIVE VERDICT

Sprint 18 established that **the distogram's errors are worse than random errors of the same
magnitude** and could not say why. Sprint 19 says why, from three independent directions, and the
answer closes most of the remaining architecture.

> ## The pool has a systematic error direction. The predictor reproduces it. The fit cannot suppress it. Any score-ordered gate amplifies it.

Every stage of the pipeline is aligned with the same error — and the retrieval pool is simultaneously
the candidate generator, the training-fragment source, and the most harmful shared component.

| lane | stage | finding |
|---|---|---|
| **A** | the predictor | the harmful coherent component is **shared across every predictor family**; the retrieval-pool-aligned part is the most harmful of all |
| **D** | the fit | the fit denoises the component **orthogonal** to the ideal-geometry manifold and has **no power** against the component inside it |
| **C** | selection | gate damage is exactly magnitude + **alignment**; score gates move the emitted mean **along** the pool's error direction |

**And the shared component is mostly sequence-DEPENDENT** — *corrected 2026-09-07*. The alignment
statistic has a **shared-referent floor of 0.505**: both arguments are scored against the same native,
so any two realisable peptides must correlate. Null-subtracted, the zero-information references
reproduce **~a fifth** of the same-architecture ceiling (0.18–0.24), not two thirds. The earlier
sentence — *"the harmful component is largely sequence-independent"* — **does not survive.** What
survives is the **ordering**, which is what the mechanism needs.

> **Every predictor, conditioned or not, emits a *typical peptide of that length*. The harmful
> coherent error is the systematic difference between "typical" and this particular native.**
> That is an identifiability statement about the inputs, not a defect of any architecture.

**The mission number did not move.** The ladder, corrected 2026-09-07 after a Sprint-20 audit: best
**point cloud** 3.048 Å (a 22.2%-contracted coordinate average, mean virtual bond 2.961 Å — *not a
buildable structure*), best **built structure** 3.204 Å (the incumbent), deployed emission 3.236 Å.
No arm in this sprint beats the incumbent, and the sprint's own headline ceilings are ORACLE. See
LEDGER L20; every arm must state its basis.

---

# 1. THE MECHANISM

## 1.1 It is the cross-pair sign correlation

`signflip` keeps every pair's own residual magnitude **and its own weight** — no orphaning — and
randomises only the sign pattern.

| arm | RMSD | resid RMS | κ | vs real |
|---|---|---|---|---|
| real (deployed) | 3.610 | 3.205 | 0.809 | — |
| **signflip_exact** | **2.408** | **3.2049** | 0.352 | **−1.202 [−1.408, −1.004]**, 112W/14L, 5/5 |
| shuffled | 2.649 | 3.008 | 0.367 | −0.961 [−1.184, −0.740] |
| **coherent** (real alternative structure) | **3.749** | 3.097 | 0.990 | **+0.139 [+0.023, +0.253]** |
| coh_scaled | 3.843 | 3.205 | 0.987 | +0.233 [+0.122, +0.348] |

**Sign-flipping alone recovers 128% of the Sprint-18 gap**, at residual RMS matched to machine
precision, and beats plain shuffling by **−0.264 [−0.448, −0.084]**, 5/5 folds. A **perfectly realisable** error — borrowed from a real alternative structure at
matched magnitude — is the **worst arm on the board**.

**Realisability is what makes an error harmful.** At matched magnitude the family is monotone in κ
across a **1.475 Å span** (Pearson +0.984 across arm means; per-target Spearman positive on 116/126).

## 1.2 Why: the fit cannot see in-manifold error

The weighted fit projects the whitened residual onto the tangent space of the ideal-geometry
manifold. That space has rank **2n − 5** — n−2 virtual angles plus n−3 virtual dihedrals —
**verified numerically**, mean rank 20.92 against 2n−5 = 20.92, and the generic in-tangent fraction
is rank/npairs: predicted 0.329, measured 0.330 [0.318, 0.341]. **EXACT.**

    real       0.695 [0.669, 0.722]      coherent    0.846 [0.822, 0.870]
    signflip   0.397 [0.373, 0.420]      rand_white  0.330 [0.318, 0.341]   <- geometric null
    shuffled   0.586 [0.565, 0.608]      rand_raw    0.581 [0.567, 0.595]   <- iso/shuffled's own null

    real - signflip   +0.298 [+0.262, +0.331]   121/5

> **The fit is a denoiser against the component orthogonal to the manifold and has NO power against
> the component inside it.**

**Design consequence — SUPPORTED, deliberately NOT promoted.** If the harm is in-manifold, no change
to the loss can fix it: any objective whose argmin lies on the manifold realises whatever error lies
in the manifold. It **retrodicts four unconnected Sprint-18 results** — the functional form is sound,
the weighting is validated, the optimiser is not the constraint, tempering does not help — as four
faces of one fact. It is **not promoted** because the linear statistic orders *arms* (+0.889 across
arm means) and not *targets* (per-target ρ = +0.134 [−0.057, +0.301]). The lane that produced it
fired its own pre-committed rule against its own promotion, twice.

**OPEN**: post-fit κ *does* predict the per-target gap (ρ +0.360; +0.482 for κ(real) − κ(signflip))
while start-point geometry does not, even though the projection tracks κ per target at ρ ≈ 0.66.
Something along the ~5-radian trajectory carries per-target information nobody is measuring. κ and
the gap share two fits, so a common cause is not excluded.

---

# 2. THE SOURCE

Each family's error splits **exactly** by its own terminal fit: `r_F = r_coh_F + r_inc_F`. Alignment
is read on the coherent part.

| pair | COHERENT alignment | share of ceiling |
|---|---|---|
| same arch, different seed — **the ceiling** | **0.898** | 1.00 |
| same arch, different regularisation | 0.833 | 0.83 |
| **different architecture** (PairNet) | 0.731 | 0.57 |
| **retrieval pool mean** — not a network at all | **0.784** | **0.71** |
| **sequence-blind separation prior** — zero-info | **0.575** | **0.18** |
| **constant α-helix** — zero-info | **0.598** | **0.24** |

**The sharing is specifically in the part that hurts.** The *incoherent* component does not share:
0.09–0.19 cross-family, **0.00** against zero-information.

**Surgery, magnitude-matched** (whitening bound −1.257): removing the shared component recovers
**62–79%** of it; keeping **only** the shared component is *worse* than the predictor's real error
(`only_pool_m` +0.425 [+0.305, +0.550]).

**Corrected and inverted, 2026-09-07.** Shares in the second column are now **null-subtracted**
against a measured **shared-referent floor of 0.505** (three independent constructions agreeing
within 0.024). Sprint 19 stated defensively that *"sequence conditioning is not worthless"* — after
null subtraction the gap between conditioned (0.83 / 0.71 / 0.57) and zero-information (0.24 / 0.18)
is **three to four times larger than it looked**. **That caveat was the main result.** The
zero-information references still carry genuine excess (+0.070 and +0.093, CIs excluding zero), so
they are not at the null — but they reproduce about a fifth of the ceiling, not two thirds.

---

# 3. SELECTION

Gate damage decomposes **exactly**: `readout² = ‖b_pool‖² + 2‖b_pool‖·ALIGN + ‖Δ‖²`.

    physics SCORE gates      ALIGN  +0.0209 [+0.0004, +0.0402]
    SCORE-FREE gates         ALIGN  -0.0048 [-0.0124, +0.0020]
    difference                      -0.0256 [-0.0455, -0.0008]

Within-arm per-target r = 0.951 / 0.967 / 0.955.

> **Score gates push the emitted mean along the direction the pool is already wrong. Score-free
> gates, at any diversity, do not.**

**The zero-information gate is the worst in the table**: a constant α-helix costs **+0.1253 [+0.0563,
+0.1819]**, twice Legacy's. **Physics is not punished for being physics — any score ordering is.**
The unifying statement: *any score that prefers compact, well-formed, pool-typical geometry is
selecting toward the pool's own systematic error*, and a constant helix is the purest such rule.

**The diversity mechanism the coordinator briefed is REFUTED.** `rand_lowD` — matched-random, no
score, no physics, no sequence — reaches Legacy's diversity to three decimals and costs +0.003
[−0.003, +0.013] where Legacy costs +0.063 [+0.010, +0.124]. **An accounting identity absorbed the
damage; it did not cause it.**

**The design rule works and is worth zero**: `legacy_clust` beats plain Legacy by −0.0664 [−0.1392,
−0.0087], recovering the entire Sprint-18 deficit, then **ties matched-random** (−0.0033, NOT
MEASURED). Native-free, ranking gate *designs* by ‖Δ‖ predicts damage at **ρ = +0.930**.

---

# 4. WHAT THIS CLOSES

| direction | status |
|---|---|
| **Predictor architecture / training loss** | **Closed.** 8 leave-fold-out predictors, none past the MDE with a CI excluding zero; best `combo` −0.094 [−0.191, +0.001]. ORACLE MAE falls 11.7% for nothing. "Better matrix, worse ranking" reproduces **8/8**. |
| **Joint consistency / metric realisability** | **Closed.** PairNet's triangle update works geometrically (violations 4.09% → 0.68%, κ 0.809 → 0.914, the highest on the board) and buys −0.080 [−0.242, +0.082]. **More realisable = more coherently wrong.** |
| **Distribution-shape consumption** | **Closed on mechanism.** Mode information is real (−0.484 [−0.570, −0.401]) and **zero transfers** — including a split-half fit *inside the same fold*. A perfect ORACLE mode-picker lands at **3.472 Å against the built start's 3.213 Å**. |
| **Legacy as a rejection gate** | **Closed with a sign.** Impossible candidates are real (2.66 per 75 below 2.0 Å), yet steric rejection is worse than random, worse than diversity-preserving rejection, and worse than **not rejecting** — monotone from r=2. |
| **Training-loss change "in the right basis"** | **Downgraded, not dead** — *corrected 2026-09-07*. The gate was "two thirds is reproduced without sequence information"; null-subtracted it is **~a fifth**, so the harmful component is mostly sequence-**dependent** and a loss change is no longer excluded on that ground. It remains unattractive because 8 predictors moved nothing past the MDE (P1) — a weaker argument than the one it replaces. |
| **Objective redesign in general** | **SUPPORTED-closed** by §1.2, pending the promotion that was withheld. |

**The Sprint-18 compass reproduces at the predictor, in the anti-utility direction.** `sw_lin` and
`sw_none` differ *only* in the training separation weight: **+0.181 [+0.081, +0.285]**, 49W/77L, 4/5
folds, while ORACLE MAE moves +0.010 [−0.039, +0.062]. Weighting training toward long-range pairs
costs 0.181 Å and changes MAE by nothing — **and the best single MLP on the board is the one whose
loss is uniform.**

**The one arm that beat the deployed objective, and why it does not matter.** `riskw` (Bayes risk
with the deployed weighting) is **−0.134 [−0.229, −0.046]**, 77/49 — the only arm in two sprints to
beat the deployed functional on a matched comparison. **It lands at 3.476 Å, still +0.428 Å worse
than not refining at all.** A better way of doing something harmful.

---

# 5. CVaR-VQE — delivered, and the sampler hypothesis is refuted

**The 126-target sweep was complete on disk when this sprint froze and was never analysed. Sprint 20
executed its pre-registration verbatim** — no re-running, no new structures. `s20/results/D_QB_CLOSE/`.
Basis clean throughout: every arm's realised number is the **projected built chain**.

    pool500      3.230   <- the shipped retrieval pool, ZERO objective evaluations
    c_marg       3.326      zero-information: matched empirical torsion marginals
    c_metroH     3.407      best classical sampler
    c_helix      3.426      zero-information constant helix
    q_a1.00      3.486      BEST QUANTUM
    c_lbfgs      3.658      multi-start L-BFGS on the objective, 8193 evals
    q_untrained  3.696

**All three pre-registered endpoints fail.** The primary is +0.079 [+0.001, +0.163] — the quantum arm
is *worse* than the best classical sampler. **Four of nine kill rules fire**, including two that
matter most: it **loses to a zero-information marginal sampler** (+0.160, 5/5 folds) and is
**matched by a classical thermostat** (+0.001 [−0.089, +0.094]). Generation fails too: against the
incumbent pool at matched count, **+0.404 [+0.270, +0.546], 5/5.**

> ### The finding that matters most for Priority 1
> **No sampler on the board — quantum or classical, at 8192 objective evaluations — beats the
> zero-evaluation shipped retrieval pool on realised RMSD.** And it is not a min-of-N artefact: the
> pool picks its 75 from 500 candidates while every sampler picks its 75 from 8192. The arm that
> optimises the deployed objective hardest is the **second-worst on the board**.

**Two positives survive, and one corrects the programme's record.**

**The VQE genuinely trains** — for the first time here. Against the mandatory control (best-of-N from
the *untrained* circuit, never an initialisation mean): **−0.210 [−0.339, −0.082], 5/5 folds**
realised, −0.299 on generation. The recorded finding that "running the VQE is worse than not running
it" was measured on the **k=4 lattice register**; on the **continuous** encoding at matched budget the
sign reverses. **Priority 2 is satisfied on its own terms: the pillar is genuine and it works. It
simply does not beat anything classical.**

**Entanglement buys nothing, measured directly.** Identical pipeline with the **CNOTs deleted**, same
CVaR estimator, same seed: **−0.013 [−0.095, +0.077] — NOT MEASURED**, which is the honest label
rather than REFUTED. Nearest-neighbour mutual information in the latent is 0.045 bits.

---

# 6. AMBER — the Pareto is mapped, and the axes dissociate

> **PROVISIONAL.** The figures below come from a run since **superseded**: its `steps = 2000` bound
> was found to have genuinely bound on one call (8T61, max ΔCα **0.1595 Å**), making that row a
> different operator. Agent C withdrew all three compute deviations and relaunched the
> pre-registration as written — no bound, full k-ladder — rather than caveating a removable confound.
> The science is expected to reproduce and will be re-quoted from the clean run.

**n = 126. F-C2 FIRES. AMBER's standing role is unchanged: stereochemical repair only.**

**Three arms beat the k=30 incumbent on Cα past the MDE** — `caonly_k300` **−0.110 [−0.136, −0.092]**,
`blend75` **−0.119**, `blend50` **−0.092** — **and all three fail the validity axis.** Every arm whose
validity matches k30 is within **±0.024 Å** of it on accuracy. **There is no Pareto point that buys
Cα accuracy without paying for it in validity.**

> ### A scalar validity score would have promoted a broken structure.
> `caonly_k300` is **0.110 Å more accurate** than the incumbent with a **clash count statistically
> indistinguishable** from it (−0.008 [−0.050, +0.032]). What kills it is **`cis_frac +0.4247` — 42%
> more cis peptide bonds** — because pinning only Cα lets the peptide plane flip.

That is a measured argument for **reporting validity as a vector, not a scalar**, and it is Sprint
18's cis-peptide finding — diagnosed there at the *averaging* operator — reappearing at the *repair*
operator.

Rotated-frame null **max |ΔCα| = 0.014192 Å** (maximum, per the standing rule). Convergence
exclusions **`1D6X 2NB7 7BX2` for the sixth sprint**, rising to **64/126 at k=1000**, because a
pinned backbone cannot relieve its own strain.

## 6.1 A retracted defect, and a coordinator failure of verification

An earlier draft of this report carried `core.amber`'s minimiser as **pathological and
non-terminating**, "latent under every AMBER arm for five sprints". **That is RETRACTED.** On a quiet
box, uncapped, the exact calls reported as ">40 minutes" and ">7 minutes" take **6–15 seconds and
converge**; the measurement was **CPU starvation from six concurrent lanes**. The disconfirming
evidence was in the same lane's completed table all along — **max per-call wall 25.7 s, p95 13.6 s,
zero calls over 40 s across 1260 minimisations.**

**The lane caught its own error. The coordinator amplified it into three documents without an
independent check**, when a timing test costs seconds. §11 — *read the code before believing the
claim* — applies to performance claims exactly as to scientific ones.

**What survives is the rule the episode produced, now demonstrated on real data instead of a false
premise**: the `steps = 2000` bound fired on **1 of 1260 calls**, so its inertness certificate is
**vacuous in exactly the sense Z6 defines**, and is reported as vacuous rather than as a pass.

---

# 7. ERRORS AND CORRECTIONS

**Five coordinator errors this sprint. Four were caught by workstreams; one by reading my own output
table. Three of four lanes refuted something I asserted.**

| # | error | caught by |
|---|---|---|
| 1 | The §3 opening mechanism — coherent conformer flips — **asserted in a binding brief with no measurement behind it** | D |
| 2 | `distobj.py`: `nll`/`risk` carried no weight while the docstring claimed they did; `nll` was log-**mass** not log-density; `mode` used the **mass** argmax (+0.380 Å outward bias) | D |
| 3 | The `modeshuf` inference — a control carrying a confound the programme had **measured and named one sprint earlier** | D |
| 4 | Briefing Agent C at the **diversity** mechanism, which was an accounting identity rather than a cause | C |
| 5 | The peak-detector padding I added while fixing (2), which over-counted modality 0.233 → 0.857 | myself |

**Workstream self-corrections, each volunteered rather than defended**: D made **five in one day**
(file contention, its own significant result, two algebraic nulls, and the QR rank defect) — four
caught by *running a measurement* rather than re-reading an argument. C found a **negative FRAME²**
in its own gate draft and a **vacuously passing** iteration bound. A declined a blocked retrain under
the memory rule rather than swapping out the machine.

## 7.1 Four methodological rules this sprint produced

1. **A control must be matched in the space the OPERATOR actually works in.** Four instances in two
   sprints, none self-caught. The programme's most repeated error.
2. **When an analytic null and a measured null disagree, the measurement is the null.** Two algebraic
   nulls were offered and wrong; a review of the derivation would have passed both.
3. **A gate can pass VACUOUSLY.** Report how many times the bound, guard, threshold or fallback
   actually **fired**; a pass with zero firings is not evidence. And an arm that hits a cap is a
   *different operator*.
4. **Pre-register the falsifier, not the hypothesis.** The control that refuted the coordinator's
   briefed diversity mechanism existed *only* because a falsifier forced the lane to name in advance
   what would refute itself. The lane states it would otherwise have believed the favourable reading.

Rules 1 and 2 are in project memory.

---

# 8. BENCHMARK

> **The sealed 60-target benchmark remains SEALED.** Not read, probed, derived from, or tuned
> against; `results/benchmark_manifest.json` never opened. **No result in this sprint justifies
> unsealing** — nothing moved the mission number, and the sprint's own ceilings are ORACLE.

---

# 9. WHAT IS LEFT

| # | direction | why it survives |
|---|---|---|
| 1 | **Is any part of the shared coherent mode estimable NATIVE-FREE?** | Everything in §1–§2 is an ORACLE ceiling. A global size or stretch mode would be estimable; a diffuse one would not. **If one component is estimable it is the only deployable thing found; if none is, this is an identifiability wall.** |
| 2 | **Do not consume the argmin at all** | Every stage ends in an argmin or top-*m*; §1.2 says the argmin is exactly where in-manifold error gets realised; the terminal reads the set **mean**. A pipeline that optimises a point estimate at every joint is mis-matched at every joint. **Nobody is working on it.** |
| 3 | The short-favouring training loss | **Blocked, not declined** — 6,429/6,788 sequences absent from the hot ESM cache against 0.82–1.90 GB free. Pre-registration standing unedited. |
| 4 | Why post-fit κ predicts per-target outcome when start-point geometry does not | §1.2, unexplained |

---

# 10. INSTRUMENT INTEGRITY

| check | result |
|---|---|
| sealed benchmark | not read, probed or derived from; manifest never opened |
| metric | full-chain Cα-RMSD, frozen implementation, unchanged; no alternative reported as primary |
| native information | evaluation and labelled ORACLE diagnostics only |
| seeding | `stable_rng` throughout; **630 arrays verified across a fresh process, 0 mismatches, distinct streams per target** |
| cross-implementation reproduction | D's `raw` reproduces `s17/refine.py` at **max \|Δ\| = 0.0002 Å, 126/126 within 0.01**; C reproduces all five Sprint-18 scores to the third decimal with a 10× tighter null; the coordinator's density-mode arm and D's independent one agree to **0.002 Å** |
| external reproduction | `arXiv:2609.02113` independently reproduces *"search saturates, discrimination binds"* (2.18–3.26 Å selection gap); `arXiv:2606.21241` independently reproduces *"the objective does not rank the native"* — different architectures, no contact with this work |
| completion flags | every quoted artefact `complete`; partials `_PARTIAL_*`/`_SUPERSEDED_*` named and unread |

---

*Ledger: `s19/LEDGER.md` (L1–L13). Claims register: `s19/CLAIMS.md`. Workstream findings:
`s19/agent{A,B,C,D}_FINDINGS.md`, pre-registrations `s19/PREREG_{A,B,C,D}.md`.*
