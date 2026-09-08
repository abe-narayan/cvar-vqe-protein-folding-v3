# SPRINT 19 — CLAIMS REGISTER

Labels: **EXACT** (theorem/identity, not a discovery) · **ORACLE** (needs the native; never a
predictive result) · **ESTABLISHED** · **SUPPORTED** · **PLAUSIBLE** · **OPEN** · **NOT MEASURED**
(a zero-spanning CI without the power to exclude the effect of interest) · **REFUTED** ·
**RETRACTED** (published here and withdrawn here).

Instrument: 126 cluster-disjoint targets, full-chain Cα-RMSD, frozen implementation. **MDE at 80%
power = 0.084 Å.** Incumbent 3.204 Å; best structure built 3.048 Å (coordinate average of the
shipped top-75).

---

## BLOCK M — THE MECHANISM (the sprint's central result)

| # | claim | label | evidence |
|---|---|---|---|
| M1 | **The harm in the distogram's error is its CROSS-PAIR SIGN CORRELATION.** `signflip_exact` — each pair keeping its own residual magnitude *and its own weight*, residual RMS matched to **machine precision** (3.2049), only the sign pattern randomised — reaches **2.408 Å**, recovering **128%** of the Sprint-18 gap, and beats plain shuffling by **−0.264 [−0.448, −0.084]**, 83W/43L, 5/5 folds. *(The first-published −1.242 / 2.368 carried a clip-induced magnitude confound; superseded, not deleted.)* | ESTABLISHED (ORACLE) | A, `a_audit.json`, n=126, L3/L15 |
| M2 | **A perfectly realisable error is MORE harmful.** The positive control — a residual taken from a real alternative structure at matched magnitude — is the worst arm on the board: +0.139 [+0.023, +0.253], and +0.233 [+0.122, +0.348] rescaled. | ESTABLISHED (ORACLE) | A, L3 |
| M3 | At matched magnitude the family is monotone in realisability across a **1.475 Å span**; **realisability is what makes an error harmful** (Pearson +0.984 across arm means; per-target Spearman positive on 116/126). | ESTABLISHED | A + D, L3/L6 |
| M4 | **The fit is a denoiser against the component orthogonal to the ideal-geometry manifold and has NO power against the component inside it.** Real field **0.695** of whitened energy in-manifold vs **0.397** for a sign-randomised version of itself: +0.298 [+0.262, +0.331], 121/5; generic direction 0.330. | ESTABLISHED (the measurement) / **SUPPORTED** (the mechanism) | D, `D_T3_tangent2`, L6 |
| M5 | **Design consequence: no change to the loss can fix an in-manifold error**, since any objective whose argmin lies on the manifold realises whatever error lies in the manifold. Retrodicts four unconnected Sprint-18 results. | **SUPPORTED — deliberately NOT promoted** | D fired its own pre-committed rule: per-target ρ = +0.134 [−0.057, +0.301], L6 |
| M6 | rank of the whitened distance Jacobian = **2n − 5** (n−2 virtual angles + n−3 virtual dihedrals); generic in-tangent fraction = rank/npairs. Predicted 0.329, measured 0.330 [0.318, 0.341]. | **EXACT** (verified numerically, mean rank 20.92 vs 2n−5 = 20.92) | D, L6 |
| M7 | Why post-fit κ predicts the per-target gap (ρ +0.360; +0.482 for κ(real)−κ(signflip)) when start-point geometry does not (+0.134), though the projection tracks κ per target at ρ ≈ +0.66. **No explanation.** κ and the gap share two fits, so a common cause is not excluded. | **OPEN** | D, L6 |

---

## BLOCK S — SOURCE: where the harmful component comes from

| # | claim | label | evidence |
|---|---|---|---|
| S1 | **The harmful coherent component is SHARED across independent predictor families, not architectural.** Share of the same-seed ceiling (0.898): cross-architecture 0.81, retrieval pool 0.87. | ESTABLISHED (ORACLE) | A, `a_source.json`, n=126, L11 |
| S2 | **Two thirds of it is reproduced by ZERO-INFORMATION references**: sequence-blind separation prior **0.64**, constant α-helix **0.67** of the ceiling. | ESTABLISHED (ORACLE) | A, L11 |
| S3 | **The sharing is specifically in the part that hurts.** The *incoherent* component does not share: 0.09–0.19 cross-family, **0.00** against zero-information. | ESTABLISHED | A, L11 |
| S4 | **Mechanism**: every predictor, conditioned or not, emits a *typical peptide of that length*; the harmful coherent error is the systematic difference between "typical" and this native. **An identifiability / prior-mean statement.** | SUPPORTED | A, L11 |
| S5 | Sequence conditioning is **not** worthless — 0.833 (same arch) and 0.784 (retrieval) sit well above 0.598 (helix). The claim is that the *majority* survives without sequence information, not that sequence buys nothing. | ESTABLISHED — stated to bound S2/S4 | A, L11 |
| S6 | **Removing the shared component recovers 62–79% of the whitening bound**; keeping *only* it, at matched magnitude, is worse than the predictor's real error (`only_pool_m` +0.425 [+0.305, +0.550]). | ESTABLISHED (ORACLE) | A, L11 |
| S7 | The most harmful shared component is the **retrieval-pool-aligned** one — and the pool is simultaneously the candidate generator and the training-fragment source. | ESTABLISHED | A, L11 |

---

## BLOCK C — SELECTION: why score gates damage an averaged set

| # | claim | label | evidence |
|---|---|---|---|
| C1 | **Gate damage is exactly a magnitude term plus an ALIGNMENT term**: `readout² = ‖b_pool‖² + 2‖b_pool‖·ALIGN + ‖Δ‖²`. | **EXACT** | C, L12 |
| C2 | **Score gates push the emitted mean ALONG the direction the pool is already wrong; score-free gates do not.** Physics gates ALIGN +0.0209 [+0.0004, +0.0402]; score-free −0.0048 [−0.0124, +0.0020]; difference −0.0256 [−0.0455, −0.0008]. Within-arm per-target r = 0.951/0.967/0.955. | ESTABLISHED | C, n=126, L12 |
| C3 | **The diversity mechanism the coordinator briefed is REFUTED.** `rand_lowD` (no score, no physics, no sequence) reaches Legacy's D to three decimals and costs +0.003 [−0.003, +0.013] where Legacy costs +0.063 [+0.010, +0.124]. **An accounting identity absorbed the damage; it did not cause it.** | **REFUTED** (coordinator's briefed mechanism) | C, L12 |
| C4 | **Physics is not being punished for being physics — any score ordering is.** The zero-information constant-α-helix gate is the *worst* in the table: +0.1253 [+0.0563, +0.1819], twice Legacy's. | ESTABLISHED | C, L12 |
| C5 | The diversity-preserving gate design works and is worth nothing: `legacy_clust` beats plain Legacy by −0.0664 [−0.1392, −0.0087] — recovering the whole Sprint-18 deficit — then ties matched-random (−0.0033 [−0.0198, +0.0112]). | ESTABLISHED / NOT MEASURED vs random | C, L12 |
| C6 | Native-free half: ranking gate *designs* by ‖Δ‖ predicts damage at **ρ = +0.930**. | ESTABLISHED | C, L12 |
| C7 | **Legacy's last unrefuted role is closed with a sign.** Impossible candidates are real (2.66 per 75 below a 2.0 Å contact), yet steric rejection is worse than matched-random (+0.0175), worse than diversity-preserving rejection (+0.0170) and worse than *not rejecting* (+0.0204); monotone from r=2. **F-C3 fires.** | **REFUTED** | C, L12 |

---

## BLOCK X — THE CONVERGENCE

| # | claim | label | evidence |
|---|---|---|---|
| X1 | **The pool has a systematic error direction; the predictor reproduces it; the fit cannot suppress it; and any score-ordered gate amplifies it.** Three lanes, three instruments, three stages of the pipeline, one direction. | **ESTABLISHED — the sprint's central statement** | M4 + S1/S7 + C2 |

---

## BLOCK P — PREDICTOR INTERVENTIONS

| # | claim | label | evidence |
|---|---|---|---|
| P1 | **Eight leave-fold-out predictors, not one moves Cα-RMSD past the MDE with a CI excluding zero.** Best is `combo` at −0.094 [−0.191, +0.001] — CI touches zero. Meanwhile ORACLE MAE falls 2.347 → 2.073 (−11.7%) for nothing. **"Better matrix, worse ranking" reproduces 8/8 on the objective path.** | ESTABLISHED | A, `a_models.py`, n=126, L11 |
| P2 | **Joint consistency is CLOSED at the architecture level.** PairNet's triangle update works geometrically (violations 4.09% → 0.68%, EDM defect 0.286 → 0.148, κ 0.809 → 0.914 — the highest on the board) and buys −0.080 [−0.242, +0.082]. **Making the field more realisable makes it more coherently wrong.** | **REFUTED** | A, L11 |
| P3 | **The Sprint-18 compass reproduces at the PREDICTOR, anti-utility direction.** `sw_lin − sw_none` = **+0.181 [+0.081, +0.285]**, 49W/77L, 4/5 folds, while ORACLE MAE moves +0.010 [−0.039, +0.062]. Training toward long-range pairs costs 0.181 Å and changes MAE by nothing; the best single MLP is the one whose loss is **uniform**. | ESTABLISHED | A, L11 |
| P4 | Metric realisability of the predicted field is real (EDM defect 0.286, 4.09% triangle violations vs the native's 0.002/0.02%) and **is not the mechanism** — a magnitude-matched incoherent field is worse on both and lands 1.24 Å better. **Nobody should design against it.** | ESTABLISHED — labelled diagnostic only | A, L3 |
| P5 | Retraining with a **short-favouring** loss — the pre-registered positive half of P3. **Blocked, not declined**: 6,429/6,788 training sequences absent from the hot ESM cache, needing 1.5 GB against 0.82–1.90 GB free. | **OPEN** — pre-registration standing unedited | A, L11 |

---

## BLOCK O — THE OBJECTIVE (coordinator's lane)

| # | claim | label | evidence |
|---|---|---|---|
| O1 | The premise stands: the predictor emits a full 17-bin distribution, the objective path keeps only (mean, sd), and the model puts **0.634** of its mass within 1 Å of the mean it reports. | ESTABLISHED | L1/L2 |
| O2 | ~~The collapse error is coherent across correlated pairs via conformer flips.~~ **Multimodal pairs are scattered (z = +0.27), the implied correction has zero spatial coherence (+0.003), and the real residual is coherent equally on unimodal (+0.0429) and multimodal (+0.0489) pairs.** | **REFUTED** — coordinator's mechanism | D, L1 |
| O3 | **Multimodality has no measurable effect on the accuracy of the mean** once separation and spread are controlled: NN-matched +0.007 [−0.291, +0.316]. Harm REFUTED, benefit NOT SUPPORTED. | ESTABLISHED | D, L1/L4 |
| O4 | Census overstated ~2.4×: **0.097 under a prominence gate**, not 0.217. The bin-edge attack *failed* (density correction changes nothing, 0.240). | ESTABLISHED | D, L1 |
| O5 | **The branch is closed on four instruments.** Primary `mode − moment` = **+0.035 [−0.016, +0.086]** (D's independent density-mode measurement: +0.035), no concentration on multimodal targets (high +0.036, low +0.033). | **REFUTED** | L2/L10 |
| O6 | **A perfect ORACLE mode-picker lands at 3.472 Å against the built start's 3.213 Å — +0.259 Å worse than doing nothing.** The branch has no headroom whether or not the mode information is real. *(Corrected 2026-09-07: as first published this compared a built chain against the 3.048 Å **point cloud** — a basis mismatch overstating the gap by 64%. Direction unchanged. See L20.)* | ESTABLISHED (ORACLE) | D, `D_P5_headroom`, L7/L20 |
| O7 | Mode information is real (−0.484 [−0.570, −0.401] per multimodal pair, above a min-of-N null) and **zero of it transfers** to any native-free rule — including a split-half fit *inside the same fold* (+0.057 [−0.075, +0.190]). **Closed on mechanism, not merely outcome.** | ESTABLISHED | D, `D_P6_extract`, L4 |
| O8 | ~~The mode targets carry ~0.096 Å of positional information.~~ Decomposed: position with alignment held **+0.029 [−0.031, +0.088]**, orphaning tax **+0.016 [−0.064, +0.096]**. Both below MDE. **Neither the coordinator's claim nor Agent D's ~0.17 Å scaling counter-argument survives.** | **NOT MEASURED** | L5/L10 |
| O9 | **`riskw` is the only arm in two sprints to beat the deployed functional on a matched comparison**: −0.134 [−0.229, −0.046], 77/49, above MDE. **And it lands at 3.476 Å — still +0.428 Å worse than not refining.** A better way of doing something harmful. | ESTABLISHED / irrelevant to the mission | L10 |
| O10 | `Distogram._risk` is convex in x with minimiser the weighted median of `CENTRES`, so the **selection path cannot express multimodality** — voiding the coordinator's own counter-argument. | **EXACT** | D, L4 |

---

## BLOCK Z — METHOD

| # | claim | label | evidence |
|---|---|---|---|
| Z1 | **A control must be matched in the space the OPERATOR actually works in.** Four instances in two sprints, none self-caught: `objceil.py:163` (weights permuted with residuals), `modeshuf` (offsets orphaned from weights), `iso`/`shuffled` (isotropic in raw, not whitened, space), GC0's negative FRAME² (free-superposition mixed with common-frame). **The programme's most repeated error.** | ESTABLISHED — methodological; in project memory | L9/L12 |
| Z2 | **When an analytic null and a measured null disagree, the measurement is the null.** Two algebraic nulls were offered and wrong (2n−2 over-counting; QR returning 2n columns from a rank-deficient matrix) before the empirical null settled it. A review of the derivation would have passed both. | ESTABLISHED — methodological; in project memory | L6 |
| Z3 | Agent D made **five corrections to its own output in one day**, four caught by running a measurement rather than re-reading an argument. Agent C and Agent A each caught and reported their own defects. | recorded as process | L6/L9/L12 |
| Z4 | ~~`core.amber`'s minimisation is pathological / non-terminating.~~ **RETRACTED — REFUTED by Agent C's own retest on a quiet box.** Uncapped (`steps=0`, the deployed protocol) the exact calls reported as ">40 min" and ">7 min" take **6–15 s and converge**. What was measured was **CPU starvation from six other lanes**, attributed to the minimiser. Corroborated by C's own completed 1260-call table: **max per-call wall 25.7 s, p95 13.6 s, zero calls over 40 s** — including the arms called pathological. **`core.amber` has no pathology and no non-termination defect.** | **RETRACTED** | C self-refuted, L17 |
| Z4a | **Coordinator failure of verification.** The claim was amplified into the ledger, the claims register and `FINAL_REPORT.md` §6 as a shared-machinery defect "latent under every AMBER arm for five sprints" — **on one lane's report, with no independent check**, when a timing test costs seconds. §11 of the brief says read the code before believing the claim; it applies to performance claims too. | ESTABLISHED — coordinator error #6 | L17 |
| Z6 | **A gate can pass VACUOUSLY.** C's first bound (10000) certified as bit-inert — and it was inert *because it was never reached*. A distinct error class from Z1: not a control in the wrong space, but a gate that never fired. **General form: report how many times the bound, guard, threshold or fallback actually FIRED; a pass with zero firings is not evidence.** | ESTABLISHED — methodological, new class | C self-caught, L13 |
| Z6a | **C applied Z6 against its own re-certification.** GC3 at the deployed `steps = 2000` is *also* never reached on its five test targets (the incumbent k=30 converges in 3–9 s). That is what GC3 is for — the bound must not perturb the incumbent **where the incumbent converges** — but it certifies **nothing** about an arm that does hit the bound. **An arm that hits the cap is a DIFFERENT OPERATOR** — "2000 iterations of restrained relaxation", not "restrained relaxation" — and is a legitimate Pareto point only under that label. `hit_cap` counted per arm. | ESTABLISHED — methodological | C, L13 |
| Z6b | **A bound must be verified on the calls it actually BOUND, not only on the ones it did not.** GC3 certified inertness on five converging calls and was silent about the one call that hit the bound — which was perturbed by **0.1595 Å / 0.418 kcal/mol**. The operational half of Z6: when the intervention *does* fire, the firing cases are the only ones the certificate is about. | ESTABLISHED — methodological | C self-caught, L19 |
| Z8 | **The pre-registration did the work, not the judgement.** `rand_lowD` — the control that refuted the coordinator's briefed diversity mechanism — existed *only* because falsifier F-C1 forced the lane to name in advance what would refute its own hypothesis, and "a score-free gate at matched D" is the only honest answer to that question. C states it would have believed the diversity reading when §3.4 returned `S_D = 1.48`. **The strongest single argument in this programme for pre-registering the falsifier rather than the hypothesis.** | ESTABLISHED — methodological | C, L13 |
| Z7 | ~~Any score that prefers compact, well-formed, pool-typical geometry selects toward the pool's own error.~~ **DOWNGRADED 2026-09-07 to NOT SUPPORTED** (not refuted — the difference CI spans zero). The three scores that most literally instantiate it (`rg`, `leg_compact`, `min_heavy`) are the three *lowest* rows and all point the wrong way: pooled **−0.0144 [−0.0663, +0.0437]** against Sprint 19's own trio at +0.0332 [+0.0102, +0.0572]. | **NOT SUPPORTED** | Workstream C, s20 L8 |
| Z5 | External reproduction: `arXiv:2609.02113` independently reproduces **"search saturates, discrimination binds"** (2.18–3.26 Å selection gap) and `arXiv:2606.21241` independently reproduces **"the objective does not rank the native"** — both on different architectures, from groups with no contact with this work. Its headline 0.623 Å is an ORACLE minimum over 7.3M native-scored snapshots on two targets; native-free it reports 2.90 Å and 5.39 Å. | ESTABLISHED | D, `agentD_FINDINGS.md` §4.1 |

---

---

## BLOCK W — THE WALL (Agent A's closure)

| # | claim | label | evidence |
|---|---|---|---|
| W1 | **The harmful coherent mode is the per-target SEPARATION PROFILE**, owning **52.5%** of the gap (−0.525 [−0.721, −0.343], magnitude-matched). Removing the per-target **offset** or **stretch** *hurts* (−17.4%, −6.9% of the gap) — **a global size error is not the problem.** | ESTABLISHED (ORACLE) | A, `a_struct.py`, n=126, L14 |
| W2 | **Variance explained does not price damage.** The per-residue additive mode explains the **most** variance (0.430) and owns **none** of the harm (+7.3%, NOT MEASURED); the separation profile explains 0.361 and owns half. | ESTABLISHED — methodological | A, L14 |
| W3 | The coherence is **local and share-a-residue**: standardised-residual correlation +0.087 for pairs sharing a residue, −0.061 for pairs sharing none, −0.017 permutation null. | ESTABLISHED | A, L14 |
| W4 | **The harmful mode is NOT estimable native-free.** `pool_sep` clears matched-random (−0.354 [−0.493, −0.214]) and the zero-information helix (−0.366 [−0.536, −0.199]) at 5/5 folds, **and fails its own pre-registered trap**: `pool_sep − poolfull` = **+0.193 [+0.070, +0.324]**, 44W/82L. The ladder is monotone interpolation toward `poolfull` — *"move toward the retrieval pool"*, not *"estimate the mode"*. **Every arm, `poolfull` included, is worse than the built start.** *(Corrected: the 3.048 Å comparator was the point cloud, not a built chain — L20.)* | **DOWNGRADED** | A, `a_sepfix.py`, L14 |
| W5 | **THE WALL, mechanically: the estimator is made of the bias.** The best native-free estimator of the separation profile available anywhere in the project is the **retrieval pool** — and the pool **shares 0.87** of the harmful component being estimated (S1). Not a gap in effort; a structural reason. | **ESTABLISHED — the sprint's closing statement** | A, L14 + L11 |
| W6 | **The harmful error is not in a detectable subset.** Magnitude-matched, **not one of seven** native-free criteria beats a size-matched random choice; five are worse with CIs excluding zero. **ORACLE-perfect detection of the largest 25% of errors is worth only −0.205 [−0.364, −0.041].** A's own error model ranks the ORACLE error at Spearman +0.537 and converts to nothing. | **REFUTED** | A, `a_subset.py`, L14 |
| W7 | **Detecting the error is partly possible; repairing what you detect is worth nothing.** | ESTABLISHED | W6 |
| W8 | The *unmatched* column of W6 shows `long_sep` −0.264 and `nf_model` −0.163 "beating" random — **a pure magnitude confound** (resid RMS 2.14 vs 2.79). Reading that column alone would report a detectable subset that does not exist. | ESTABLISHED — flagged by A against its own table | A, L14 |
| W9 | A's `short_sep` +0.170 / `lo_sd` +0.161 and Sprint 18's G2e `fix_short` −0.294 / `fix_confident` −0.308 are **different operations, not a contradiction**: G2e *shrinks* residuals within a class keeping every pair; A's arm *zeroes* a class and rescales the complement, concentrating the survivors and **increasing** coherence. Both agree with the sign-coherence mechanism. **Must not be quoted against each other.** | ESTABLISHED | A, L14 |

*Block Q (CVaR-VQE as sampler) open — Agent B running. Block F (AMBER Pareto) reported as a
non-delivery: F-C2 UNTESTED, AMBER's standing role unchanged.*
