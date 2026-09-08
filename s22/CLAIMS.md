# SPRINT 22 — CLAIMS REGISTER (FINAL CAMPAIGN)

Every row carries **n**, effect with a **CI**, its **per-comparison MDE** (`2.8016 × SE`), and a
disposition. ORACLE quantities are ceilings, never achievable arms, and are labelled at every
appearance. Point-cloud and built-chain RMSD are never compared. Fold-clustered CIs are reported
alongside iid wherever a primary is stated.

**Instrument:** 126 cluster-disjoint targets, 5 pinned folds, full-chain Cα-RMSD after Kabsch.
**Incumbent:** 3.048 Å (point cloud, score-filter → top-75 coordinate average over the K=500 pool).
**The 60-target benchmark remained SEALED** — not inspected, probed, tuned against, or used for any
architecture decision. **Targets: <2.5 Å primary, <2.0 Å aggressive. Neither was reached.**

---

## A. THE CENTRAL RESULT — A POSITIVE CONTROL PAIRED WITH A NEGATIVE

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| **A1** | **The per-target optimal averaging width m is REAL and TRANSFERABLE.** Selecting m on one independent half of a target's pool and applying it to the disjoint other half beats the fixed incumbent m. | 126 | **−0.244** SE 0.036, MDE 0.102, iid [−0.314,−0.176], **fold [−0.290,−0.191]**, 77W/49L | **The coordinator's own pre-registered falsifier FIRED.** 2.4× its own MDE. 65% of the apparent headroom transfers. |
| **A2** | **The gain is GENUINELY PER-TARGET, not a global m correction.** The single global best m chosen on half A is **m=75 — the incumbent's own choice.** | 126 | global correction **−0.001** [−0.007,+0.004]; per-target vs *global best* **−0.243** [−0.317,−0.174] | **ESTABLISHED.** The alternative explanation is dead; m=75 confirmed correct *globally*. |
| A2a | Against a random-m control. | 126 | **−0.412** [−0.482,−0.346], 119W/7L | **ESTABLISHED** |
| **A3** | **The selection is ORACLE** — half A's RMSD to native chooses m. | — | — | **This is a CEILING WITH A TRANSFERABILITY GUARANTEE, not an achievable arm.** The distinction is the whole finding. |
| A4 | Survives independent adversarial attack: fold-clustered CI [−0.292,−0.180]; jackknife leave-one-out range [−0.244,−0.228]; **and dropping the degenerate half-pool `m=500` rung keeps 82% of the effect** (−0.196, 2.15× MDE). | 126 | — | **ESTABLISHED** — the named weakness is not carrying it |
| A4a | A statistic of the coordinator's, correctly downgraded by audit. "Median tied rungs 1.00 of 6" uses an **exact floating-point** tie test on a **continuous** outcome, so it reads ~1 by construction. | — | — | **NEAR-VACUOUS.** Not in tension with S21 L23's within-1-SE tied set of 2.0/6 — a different, much weaker question. |

---

## B. AND IT IS UNREACHABLE — BY SELECTION AND BY COMBINATION

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| B1 | Held-out length-quartile router | 126 | +0.031 [MDE 0.070], wrong-signed | **~0% captured** |
| B2 | Ridge over 16 arms, **5 arm-set configurations** | 126 | **lost on all five**; `A_plus_latent` **+0.070 [+0.005,+0.146], CI EXCLUDING ZERO** | **Significantly WORSE than doing nothing** |
| B3 | Random Forest, original feature set | 126 | in-fold **−0.276** (>half the ceiling) → **+0.079 held out, SIGN REVERSED** | **Textbook overfitting signature** |
| B4 | Ridge + RF on **new candidate-set geometry features** (`spread(m)` and its shape — mechanistically the quantity m controls) | 126 | +0.014 [−0.051,+0.079], −3.8% capture; RF again **−0.185 in-fold → +0.060 [+0.004,+0.117] held out** | **Least harmful family tried; still does not beat the incumbent** |
| **B5** | **HEDGING across m — a fixed, global, native-free operator needing no per-target decision and with nothing to overfit.** | 126 | `hedge_all` **−0.0002** [MDE 0.056]; `hedge_core` **−0.0117** against **MDE 0.027**, iid [−0.030,+0.007] | **The coordinator's pre-registered falsifier FIRED. REFUTED.** Below its own detection threshold, 0.4% of baseline. |
| B5a | And the geometry cost of the hedges that *looked* best. Physical virtual Cα–Cα is 3.805 Å; the incumbent already sits at 2.961. | 126 | `hedge_all` **2.778**, `hedge_wide` **2.628**; the only arm that *improved* geometry (`hedge_narrow`, 3.054) was the worst on RMSD | **Any Å the wide hedges bought would have been bought by making the structure less physical.** |

> **Both doors — selection and combination — are closed on a signal that is provably present.**

---

## C. WHY: THE BOUND

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| **C1** | **At the actual training-fold size, no router class rich enough to express the signal is learnable, and no class simple enough to be learnable can express it.** Finite-class/Bernstein bound using **this project's own measured σ ≈ 0.41 Å**, cross-checked from two independent artefacts. | n≈100/fold | gap at n=100: **1 global threshold 0.39** · tercile 1-feature **0.96** · linear p=5 **0.57** · linear p=15 **1.18** · tree d=3 **1.82**. n needed for 0.24 Å: **180 / 513 / 281 / 650 / 1045**. | **ESTABLISHED.** Even the *simplest* router has a gap comparable to the **whole** 0.482 Å ceiling; everything actually tried needs **2–8× the entire benchmark**. |
| C1a | And the one class simple enough to possibly work at this n is measured at exactly what the bound predicts. | 126 | global best m = **−0.001 Å** | **CONSISTENT** |
| C1b | Robustness: dropping the range term and using variance alone. | — | 0.31 at the tercile row — same order | **not an artefact of the range choice** |
| C1c | **Caveats, kept.** Worst-case and **distribution-free** — *"no proof of learnability exists at this n"*, **not** *"a router is provably impossible"*. And targets are **not i.i.d.** (fold/family structure), so the bound is **optimistic**, not pessimistic. | — | — | **binding on every restatement** |
| C2 | Forward implication | — | — | more cluster-disjoint targets, or an explicit **variance-reduction estimator** (shrinkage/hierarchical, trading bias for variance deliberately) — **not another pass through the same router class at the same n** |

---

## D. THE QUANTUM PILLAR — AN A-PRIORI CLOSURE

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| **D1** | **On any exactly-enumerable candidate register, a genuine CVaR-VQE's ARGMIN readout is IDENTICAL to a zero-training circuit's, by construction.** `R_VQE_argmin = R_untrained_argmin = R_score_argmin` **exactly**, all 16 targets, mean difference **0.000**. A full-support ansatz over an exactly-enumerable diagonal Hamiltonian makes training irrelevant to an argmin. | 16 | exact | **A-PRIORI CLOSURE, not an empirical tendency.** Every pool this project uses is K≤500, i.e. **≤9 qubits, always exactly enumerable.** **Any claim of quantum argmin-selection advantage over a discrete pool this size is closed in advance, independent of ansatz or optimiser.** |
| D2 | Gauge/label-permutation gates on the candidate-identity register. | 16 | **Gate 1 (soundness) 0/320 fired.** Gate 2: **argmin 0/256 — provably invariant**; tail-average **~26% fired**, median \|Δ\| ≈ 0.03 Å, **22.7% > 0.2 Å, 5.9% > 0.5 Å** | **Real but moderate.** Not Sprint 18's near-uniform catastrophe; not zero. **The invariant readout is the useless one and the useful-looking one is gauge-sensitive.** |
| D2a | Self-caught before results were read: an exact numerical tie in one target's top-2 made a naive argmin **gauge-dependent through tie-break order alone** — the project's oldest trap, reproduced in a new substrate. Fixed to return the tied set; the buggy pilot **preserved as a separate artefact** rather than overwritten. | — | — | **correct handling** |
| D3 | Unregularised (T=0) tail-average against classical top-α sorting. | 16 | **+0.52 Å [+0.11,+0.89]**, 13/16 | **Significantly WORSE**, and indistinguishable from a size-matched random draw |
| D4 | The audit lane's derivation that CVaR degeneracy-breaking is **provably inert** — and the gap Workstream A found in it. | — | — | **Branch 1 CORRECT** (and A had already retired that design independently, on measuring ESS = 1.0 — convergent kill from two directions). **Branch 2 has a real gap:** it closes *point search*, not *set-selection-for-averaging*, which S21's own scope notes keep separate and list as **LIVE**. |
| **D5** | **The entropy-regularised CVaR tail does NOT beat the size-matched CLASSICAL TOP-m selector** — the bar Workstream A insisted on itself, because S21 priced the averaging operator at −0.28 to −0.31 Å *independent of the energy*. | 16 | T=0 **+0.423** [MDE 0.545] · T=0.05 +0.256 · T=0.2 +0.013 · **T=0.5 −0.029 [−0.110,+0.021], MDE 0.105 — a null.** Bar 3.2998. | **ITEM 2 CLOSES NEGATIVELY AND CLEANLY.** Extends S21's nine-arm list rather than escaping it. |
| D5a | **The regulariser works exactly as designed** — ESS 1.02 → 62.08, tail 2.1 → 68.5 as T rises. The degeneracy lever is real and controllable. | 16 | — | **ESTABLISHED (mechanism)** |
| **D5b** | **AND WHY THE BAR MATTERED.** Against a size-matched **random** draw the same arm wins at **every** temperature by **0.58–0.98 Å**. | 16 | −0.976 / −0.746 / −0.609 / −0.581 | **Against the wrong control this reads as a 0.6–1.0 Å quantum win.** The right bar turned a spurious positive into an honest null — the campaign's single best piece of experimental judgement, and it was the lane's own. |
| **D8** | **THE MECHANISM PROOF — D5 and D6 ARE THE SAME FACT.** `cvar_from_probs` orders states by **ENERGY**, not by trained probability, and assigns mass to a **prefix** of that order. **The tail-support SET is therefore provably always the classical top-m set**, for whatever m training realises, **regardless of ansatz, entangler or temperature.** | 16 | Verified by **set equality**, not RMSD closeness: trained-tail indices **element-for-element identical** to `argsort(scores)[:m]` (18 cells, 7 targets); argmin **max\|diff\| = 0.0 across all 64 cells**, independent of which H trained it | **The readout "lever" is *Distance ranks better than Legacy*, surfaced through a channel training cannot influence. NOT A QUANTUM LEVER** — a statement about the two Hamiltonians. **A-PRIORI, not empirical.** |
| D8a | D7's staging effect is a **third instance** of D8. | 16 | realised tail 1.25 (single) → 3.89 (staged) | The warm start leaves stage 2 less collapsed and lands a marginally better m. **Not a staging benefit.** |
| **D8b** | **THE ONE DOOR THE THEOREM LEAVES OPEN.** Every result here uses an **unweighted** average over a tail whose **membership** is classical. A **probability-WEIGHTED tail average** is **not covered** by the set-equality argument — untried, and not the shipped operator. | — | — | **THE LAST UNTRIED OBJECT ON THIS SUBSTRATE.** Named by the lane rather than letting the closure look total. |
| D8c | Reconciliation of the two lanes' readout figures — **no discrepancy.** | 16 | coordinator's pooled −1.177/−1.508 is the per-target average of the lane's two registered comparisons: **train=Legacy fixed −1.929, significant, 14W/2L**; **train=Distance fixed +0.523, same sign, n.s., 6W/10L** | explains the 11W/5L between the two tallies |
| **D6** | **THE READOUT-H LEVER REPLICATES ON A SECOND, INDEPENDENT SUBSTRATE** *(reframed by D8: real, reproducible, and NOT quantum).* Swapping which Hamiltonian is *read out* on the candidate-identity register. | 16 | tail-average **−1.177** SE 0.420, MDE 1.177, [−2.007,−0.410], 11W/5L; argmin **−1.508** SE 0.668, MDE 1.873, [−2.820,−0.312], 11W/5L | S21 measured −0.697 (Legacy) / −1.299 (AMBER) on the **torsion basin latent**; it reproduces at −1.18/−1.51 on a **different encoding and state space**. **Directions ESTABLISHED, magnitudes NOT MEASURED** (both at or below their own MDEs at n=16). **The only architectural lever surviving two substrates.** |

| **D7** | **MULTI-STAGE VQE (§36/§37): `H_Legacy` → `H_Distance`, 60+140 iters.** Staging helps the tail readout but neither arm reaches the classical bar. | 16 | staged−single **−0.173** [−0.334,−0.033] **MDE 0.223 → NOT MEASURED**; staged−BAR **+0.346**; single−BAR **+0.519**. Tail 1.2→3.9, ESS 1.00→1.09. | **Mechanism does something; the bar is not reached.** |
| **D7a** | **The a-priori closure (D1) holds under staging.** | 16 | argmin tied to the classical argmin on **100% of cells, BOTH arms** | **Third independent confirmation**, from an experiment designed to test something else. Changing H between stages cannot move an argmin over an exactly-enumerable register. |

---

## E. PHYSICS — THE MANDATED LEGACY vs AMBER COMPARISON

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| **E1** | **On COMPACTNESS, Legacy and AMBER are ANTI-CORRELATED at matched magnitude.** By **controlled direct perturbation**, 3,780 trials — not observational correlation. | 30 | Spearman(dE, dRg): **Legacy +0.454, AMBER −0.444** | **ESTABLISHED.** *This is the fundamental reason their rankings disagree*, and it independently replicates the mechanism Sprint 21 inferred from staged preconditioning. |
| E2 | **Steric response is SHARED, not specialised.** | 30 | — | **The difference is NOT "AMBER sees sterics and Legacy does not."** It is the compactness axis, where they oppose. |
| E3 | AMBER's disproportionate sensitivity to its **own** dominant curvature direction. | 30 | **~29× [4.5–543×]** | **ESTABLISHED** — after B's own pre-registered test came back *refuted backwards* because **torsion-index concentration is not Cartesian concentration** (NeRF cascades a torsion move down the whole chain). Diagnosed, addendum registered same-day, rerun on each potential's own top Hessian eigenvector. |
| E4 | **No defensible AMBER surrogate exists.** | 30 | bonded-subset tracks true AMBER at **ρ = 0.211** against a pre-registered bar of **0.7** | **REFUTED, with a mechanism: the steric singularity lives in the nonbonded and solvation terms a bonded subset cannot see.** Directive §18 closed. |

---

## F. THE EXTENDED ABLATION (DIRECTIVE §57)

*Pool instrument, shipped top-75 tail readout, n=126 — a **different instrument** from Sprint 21's
n=12 CVaR-VQE matrix; the two are never quoted as one table. Terms rank-normalised per target
(declared before any RMSD was read).*

| # | claim | disposition |
|---|---|---|
| F1 | **Distance alone wins: 3.0483.** No combination beats it; the falsifier did not fire. Closest is `D+G` at **+0.062**, below its own MDE of 0.110 with a CI touching zero — **a tie, not a win.** | n=126, ESTABLISHED |
| F2 | **The torsion prior is the single most harmful term tested: +0.899 Å alone [+0.677,+1.085], and +0.295 even when added to Distance.** The directive named it at §25; measured, it is **worse than a matched random tail by 0.52 Å.** | **ESTABLISHED** |
| F3 | **Six cells are at or worse than a MATCHED RANDOM TAIL (+0.378):** G, D+L+T+G, D+L+T, L, L+T, T — every CI excluding zero. | **Extends S21 C10 to the torsional and geometric terms at full n.** **Adding physically-motivated terms to the structural objective makes it worse, monotonically in the number added.** |
| F4 | Geometry adds nothing on its own account (`D+G` ties `D`) | as predicted — pool members are real windows, therefore already valid chains |
| F5 | `H_AMBER` **not re-run here** (deployed form is `E ∘ Relax_50`, one OpenMM context per candidate per target; cells exist at n=12 VQE and n=126 tail). | **stated, not silently skipped** |

---

## G. EXTERNAL CORROBORATION, AND THE TENSION

| # | claim | disposition |
|---|---|---|
| **G1** | **arXiv:2606.21241**: *"the average RMSD across all solutions is, on average across all peptides, better than the optimal solution according to the cost Hamiltonian"* — **peptides ≤15 residues**, this instrument's range. Spearman(cost, error) **negative** for small peptides. Obtained on a **different energy** (Miyazawa–Jennings contact) and a **different representation** (tetrahedral lattice). | **THIRD independent group.** Upgrade adopted: from *"our pipeline has a problem"* to **"this is a property of the method class."** |
| **G2** | **THE TENSION, stated rather than dropped.** Their paper attributes the failure partly to **RESOLUTION** (more interaction shells → better correlation). **Our data cuts against that specific mechanism**: AMBER is far richer than Legacy and **does not rank better** — it is worse in most configurations, and both lose to a matched random tail. | **If "more detail fixes it" were operative, AMBER should be our best-ranking energy. It is not.** Either "interaction shells" is a different axis from "force-field detail", or there is a genuine mechanistic disagreement — **the corroboration on the headline must not quietly import their proposed mechanism.** |
| G3 | **Novelty, checked not assumed.** CVaR-VQE for peptide folding is a **KNOWN combination** — QuPepFold (PLoS ONE, Feb 2026) plus two 2024 CVaR-VQE-vs-MD peptide papers. | **Nothing here may claim it as novel.** *Caught before it could be misreported.* |
| G4 | Plausibly still open: continuous-torsion + CVaR-VQE + all-atom AMBER **jointly**; the **readout-H / training-H separation** as a measured lever; the candidate-identity register's a-priori closure (D1); and E1's compactness anti-correlation (**no specific literature match found**, classified plausible empirical/mechanistic novelty with the caveat that one search pass is evidence of absence, not proof). | — |

---

## M. METHOD — THE ERRORS, AND WHO CAUGHT THEM

| # | finding |
|---|---|
| **M1** | **The coordinator's L2 was refuted algebraically by the audit lane and withdrawn.** With `corr(inc,lat) = 0.819` and OLS slope **0.580**, `d = lat − inc` is **forced** to trend in `inc` with slope `(β−1) = −0.42` by `Cov(inc,d) = Cov(inc,lat) − Var(inc)` — an identity, not a finding. The mechanical component reproduced **81/151/202/121%** of each quartile mean; the residual was **not monotone** and two of four CIs included zero. A **placebo** (within-length-decile shuffle) produced an **even bigger** sign-flipping shape from **zero real information**. |
| M1a | **And the coordinator's own defence failed for the same reason.** Re-stratifying on `pool_oracle` — believed to escape regression to the mean because it is independent of the incumbent's *realisation* — does **not** escape: the mechanical term is `Cov(strat,lat) − Cov(strat,inc)`, negative for **any** difficulty stratifier whenever the latent is less difficulty-sensitive. The "9.6-SE durable interaction" was mechanical too. **Both withdrawn.** |
| M2 | **A bookkeeping error found by independent reconstruction.** The routing ceiling was quoted as 0.482 from a **13-arm** run while all subsequent decompositions used a **12-arm** run. Recomputed cleanly: **0.505**. *0.505 computed cleanly beats 0.482 computed sloppily.* |
| M3 | **A null the coordinator got wrong, and the correct derivation.** A permutation null shuffling targets *within* each arm returned **1.115 — better than the observed oracle** — because arms are correlated through **target difficulty** and shuffling destroys that shared factor. The audit lane's replacement: **the oracle ceiling needs no null at all** (a deterministic min of K real numbers per target, like `latent_oracle`; asking whether a minimum is significant is a category error), and **the decision-theoretic floor for an achievable router is exactly the best fixed arm, by construction.** Held-out folds are therefore the *correct construction*, not a patch. |
| M4 | **A router-*construction* search carries its own best-of-K exposure**, separate from the per-target min-of-K inside the ceiling. Five arm-set configurations and any bin-count/feature choices must be counted. **"Not yet found" is the correct disposition, not "does not exist."** |
| M5 | **Three lanes refuted their own registered hypotheses.** The coordinator's m-noise hypothesis (A1) and m-hedge hypothesis (B5); Workstream B's concentration test, which came back *refuted backwards* and was rediagnosed and rerun rather than quietly dropped (E3); Workstream C's own pre-declared favourite weighting scheme, which was numerically the worst. |
| M6 | **A process failure that is the coordinator's.** The brief mandated serialised AMBER/OpenMM access. `s22/results/` shows fresh concurrently-written artefacts from multiple lanes throughout Workstream B's run — **the rule was not enforceable by the lanes themselves.** No artefact was corrupted, but the gate did not exist. |
| M7 | Self-caught defects reported rather than silently fixed: a bootstrap RNG re-seeded inside its resampling closure (degenerate CI, B); an exact-tie-induced gauge dependence with the buggy pilot **preserved** as a separate artefact (A). |
