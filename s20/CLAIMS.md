# SPRINT 20 — CLAIMS REGISTER

Labels: **EXACT** (theorem/identity) · **ORACLE** (needs the native) · **ESTABLISHED** ·
**SUPPORTED** · **PLAUSIBLE** · **OPEN** · **INCONCLUSIVE** · **NOT MEASURED** (zero-spanning CI
without power to exclude the effect of interest) · **NOT SUPPORTED** · **REFUTED** · **RETRACTED**.

**Basis is stated on every structural row** (point cloud / built chain / repaired emission — L20).
MDE at 80% power = **0.084 Å**. Ladder: point cloud 3.048 · **built structure 3.204 (incumbent)** ·
deployed emission 3.236.

*Blocks A (decorrelated generation), B2 (landscape/optimiser) and C2 (Pareto, repair path) open —
lanes running.*

---

## BLOCK F — THE READOUT OPERATOR (coordinator)

| # | claim | label | evidence |
|---|---|---|---|
| F1 | **The manifold-constrained Fréchet mean recovers NONE of the projection tax**: +0.0048 [−0.0117, +0.0218], 61W/65L — an order of magnitude below MDE, point estimate on the wrong side. Pre-registered falsifier fires. | **REFUTED** (coordinator's hypothesis) | L11 |
| F2 | **The projection tax is the price of the manifold, not of un-contracting.** An operator that *cannot* contract by construction (mean bond exactly 3.804 Å) pays the identical cost. Sprint 18's "irreducible" stands on much stronger evidence. | ESTABLISHED | L11 |
| F3 | **The manifold constraint determines the answer; the objective used to reach it does not.** Five different manifold-constrained operators land within **0.016 Å** of one another, while the unconstrained point cloud sits 0.157 Å below and the torsion-space mean 0.867 Å above. | ESTABLISHED | L11 |
| F4 | Controls behaved: the recorded-bad `torsmean` null is 4.07 Å and Fréchet beats it by −0.862 [−1.200, −0.545]; the medoid trap did not trigger (1.199 vs 1.097 — the Fréchet arm is *further* from the medoid). | ESTABLISHED | L11 |
| F5 | **The largest measured non-predictor loss is not addressable by changing the averaging operator.** Recovering it requires changing the candidate set or the geometry requirement. | ESTABLISHED | L11 |

---

## BLOCK X — CORRECTIONS TO THE PRIOR RECORD

| # | claim | label | evidence |
|---|---|---|---|
| X1 | **"Best built 3.048 Å" is a contracted POINT CLOUD, not a structure.** Mean virtual Cα–Cα bond **2.961 Å** vs physical 3.804 (22.2% contraction, global min **0.649 Å**); no side chain can be placed on it and AMBER cannot score it. Projecting it *gives the incumbent* — **3.204 Å is the best structure the pipeline builds.** | ESTABLISHED — coordinator-verified | L20 |
| X1a | Not a metric exploit others could copy: the best **constant** global rescale of `fit_ca` buys **0.034 Å**, below MDE. The averaging operator's contraction is per-target adaptive shrinkage, which is real. | ESTABLISHED | L20 |
| X1b | Corrected consequence: the ORACLE mode-picker headroom is **3.472 vs the built start's 3.213 = +0.259**, not +0.424. Direction unchanged, magnitude overstated by 64%. | RETRACTED (the figure) | L20 |
| X2 | **The "irreducible +0.164 Å AMBER tax" is the PROJECTION's.** AMBER's own deployed cost (`k=10`, `steps=0`, on the built chain) is **+0.0207 [+0.0141, +0.0276]** — a quarter of the MDE. | RETRACTED / corrected | L2 |
| X2a | **And the sign inverts**: on the same point cloud, projection +0.1554 vs **AMBER k=30 +0.1333** — AMBER is the **cheaper** of the two repair operators by −0.022. Sprint 18 had this right; the Sprint-20 brief collapsed the two. | ESTABLISHED | L2 |
| X2b | The deployed AMBER is not the arm any published figure describes: `amber_k=10.0`, `amber_steps=0`. **Every "AMBER k=30" figure in Sprints 16 and 19 is a different constant on a different input.** | ESTABLISHED | L2 |
| X3 | **The cross-basis sweep of Sprints 16–19 is NEGATIVE.** Sprint 18's `[H5]` table used the point-cloud baseline; recomputed on the built baseline at n=126, **no conclusion changes** — the offset is constant and every arm loses by far more than it. | ESTABLISHED — reported as a negative | L2 |
| X4 | **`s12.instrument.paired`'s CI is i.i.d. over targets, not fold-clustered**; the `folds` argument only adds a per-fold breakdown. **293 call sites.** Defensible on a cluster-disjoint instrument, but the label was wrong and clustering is the conservative choice. Bite measured: **two marginal verdicts moved to NOT MEASURED.** | ESTABLISHED — labelling + marginal claims only | L9 |
| X5 | **Sprint 19's Z7** ("any score preferring compact, pool-typical geometry selects toward the pool's error") — the three scores most literally instantiating it are the three *lowest* rows and all point the wrong way: pooled **−0.0144 [−0.0663, +0.0437]** vs Sprint 19's trio at +0.0332. | **NOT SUPPORTED** (not refuted; difference CI spans zero) | L8 |
| X6 | **The shared-referent floor.** Cross-family alignment was `corr(d(X_A) − d_nat, d(X_B) − d_nat)` — both arguments share `−d_nat`, so any two realisable peptides must correlate. Measured floor **0.505**, three constructions agreeing within 0.024. | **EXACT** (the construction) / ESTABLISHED (the value) | L5 |
| X6a | Null-subtracted, zero-information references reproduce **~a fifth** of the same-architecture ceiling (0.18–0.24), **not two thirds**. *"The harmful component is largely sequence-independent"* **does not survive** — above the null it is mostly sequence-**dependent**. | RETRACTED | L5 |
| X6b | **S5 is strengthened three- to fourfold.** The conditioned/zero-information gap (0.83/0.71/0.57 vs 0.24/0.18) is far larger than it appeared. **The caveat Sprint 19 wrote defensively was the main result.** | ESTABLISHED | L5 |
| X6c | **Undamaged**: S1's ordering (ceiling > same-arch > pool > cross-arch ≫ zero-info) holds before and after; S3 untouched; **S7/W5 untouched — the pool's null-subtracted 0.71 is still second-highest and above cross-architecture. The wall is where Agent A put it.** | ESTABLISHED | L5 |

---

## BLOCK Q — THE QUANTUM PILLAR (Priority 2)

| # | claim | label | evidence |
|---|---|---|---|
| Q1 | **The CVaR-VQE sampler hypothesis is REFUTED on all three pre-registered endpoints**, n=126, 17 arms, 8192 evaluations. Primary +0.079 [+0.001, +0.163]; **4 of 9 kill rules fire.** | **REFUTED** | L1, L3 |
| Q2 | It **loses to a zero-information null** — `c_marg`, drawing i.i.d. from the same von Mises basin mixtures and never reading the objective, has the **best realised RMSD (3.326) and best selection ceiling (2.491) of every budgeted arm**. Quantum vs it: **+0.160 [+0.089, +0.241], 5/5 folds.** | ESTABLISHED | L3 |
| Q3 | **No sampler — quantum or classical, at 8192 objective evaluations — beats the ZERO-EVALUATION shipped retrieval pool (3.230) on realised RMSD.** Not a min-of-N artefact: the pool picks 75 from 500, the samplers pick 75 from 8192. The arm optimising the deployed objective hardest (`c_lbfgs` 3.658) is **second-worst on the board.** Fifth independent instrument for *search saturates, discrimination binds*. | **ESTABLISHED — the sprint's hardest structural negative** | L1 |
| Q4 | **The VQE genuinely trains on the continuous encoding** — against the mandatory control (best-of-N from the *untrained* circuit): **−0.210 [−0.339, −0.082], 5/5 folds** realised, −0.299 generation. Reverses a recorded finding measured on the **k=4 lattice**. **Priority 2 is satisfied on its own terms.** | ESTABLISHED (ORACLE scoring, native-free decisions) | L1; memory scope-corrected |
| Q5 | **Entanglement contributes nothing measurable**: CNOTs deleted, everything else identical, **−0.013 [−0.095, +0.077]**. Coordinator adjudication: **NOT MEASURED**, not refuted — a benefit of MDE size lies inside the interval. Latent nearest-neighbour MI 0.045 bits. | **NOT MEASURED** | L3 (adjudicated against Workstream B's stronger label) |
| Q6 | **CVaR's minimiser is a FACE** — it constrains an α-tail and is *indifferent* to the rest. Verified with three exactly-CVaR-optimal laws whose diversity spans **0.00 → 3.54 Å at machine-identical objective**. | **EXACT** | L3 |
| Q7 | **The deployed ansatz is bond-dimension 4 — a ≤16-state HMM, classically samplable by construction, no separation available at any budget.** Training *reduces* its latent correlation (1.33 vs 1.52 bits). | **EXACT** | L3 |
| Q8 | **The CVaR tail parameter is worth nothing**: α=0.05/0.25/anneal vs α=1 (no truncation) are +0.012/+0.013/+0.025, all NOT MEASURED, **all point estimates on the wrong side.** α=1 is the best quantum arm. | ESTABLISHED | L4 |
| Q9 | Theory: **arXiv:2605.02850** proves CVaR cannot eliminate barren plateaus, that its hard α-quantile cutoff **creates optimisation discontinuities**, and that resolving the sharpened gradient costs **exponentially in the tilt**. On 512 shots, α=0.05 is **26 samples** — α=1 winning is the predicted regime. | ESTABLISHED (external) | L4 |
| Q10 | **"CVaR concentration" is not a landscape property of the Hamiltonian.** The tail-membership discontinuity belongs to the **rule**, identical under both Hamiltonians; any measured difference is two energy **spectra** passed through a **shared discontinuous operator**. | **EXACT** — derive the operator before interpreting its statistic | L4 |
| Q11 | **Every circuit-side distributional metric is a target-difficulty proxy.** Partialling out difficulty: entropy +0.222→+0.051, max_prob −0.268→−0.066, ESS +0.248→+0.051. Only `M` survives, and `M` is member error — nearly the outcome itself. **Binding: partial out difficulty on every landscape metric.** | ESTABLISHED | L4 |
| Q12 | The CVaR machinery is **genuine and correct**: `tail_indices` reproduces the stable-sort mask with a documented tie rule; the gradient baseline is a true control variate at **cos +1.000000** against the exact-expectation reference; the known defect is preserved as a *named, measurable arm*. **The pillar is intact; its measured RMSD contribution is zero to slightly negative.** | ESTABLISHED | L4 |
| Q13 | The latent's **entire** labelling freedom (`b_i ↔ 1−b_i`, a `Z2^n` gauge) is tested: identity labelling at mid-rank percentile **0.500** (median 0.479), orbit difference +0.018 [−0.130, +0.170]. **Nothing in this lane's conclusion lives in the encoding** — contrast Sprint 18's W1 at percentile 0.00. | ESTABLISHED (PARTIAL, n=20) | L3 |

---

## BLOCK L — LEGACY vs AMBER (Priority 3)

| # | claim | label | evidence |
|---|---|---|---|
| L1 | **The two potentials disagree about ordering**: `Spearman(E_Legacy, E_AMBER) = −0.0886 [−0.1239, −0.0514]`, median −0.095, **0/126 targets with \|ρ\|>0.8**. Pearson **+0.0312 — the opposite sign**, surviving winsorisation. **They agree about catastrophes and disagree about ordering.** | ESTABLISHED | L8 |
| L2c | **What separates them is COMPACTNESS.** Against 200 matched-random partitions of identical cell sizes, Legacy-preferred candidates are **0.45 Å more compact** (−0.4477 [−0.5110, −0.3911], **124W/2L, 5/5 folds, 3/3 length terciles**), with tighter contacts and better Ramachandran; AMBER-preferred are expanded with open sterics. **Legacy is a compactness/typicality model wearing a physics vocabulary** — the opposite of the naive prior, since its steric term carries the largest weight. | ESTABLISHED | L8 |
| L3c | **Both ORACLE axes are the two that do NOT separate** (rebuild distance −0.0785 [−0.2120, +0.0464]; pool-error alignment −0.0328 [−0.1526, +0.1010]). **Complementary about geometry, jointly blind to accuracy.** | ESTABLISHED | L8 |
| L4c | Mandated ablation table: **every ordered arm loses to NO GATE** with a CI excluding zero — legacy +0.0763, amber +0.0496, leg→amb +0.0535, amb→leg +0.0488. **Composites never beat their own first stage.** | ESTABLISHED | L8 |
| L5c | **53.5% of the real candidate pool is a steric catastrophe.** On 9,450 top-75 rebuilds unrelaxed AMBER is **finite everywhere** with median **+16,062**, p99 3.1e13, max **5.5e23** kcal/mol; 53.5% above 1e4. **Finite-but-meaningless is more dangerous than +inf — nothing throws.** Not a lattice artefact: a property of the ideal-geometry rebuild. | ESTABLISHED | L8 |
| L6c | **EXACT**: `E_AMBER = f(x(θ))` but `E_Legacy = f(seq, x(θ), θ)`. `φ₀` moves no atom, so **AMBER's gradient there is exactly 0.0 and Legacy's is not.** The structure-invisible channel is **1.4–6.7% of ‖∇E_Legacy‖** — size stated because the categorical claim alone oversells it. | **EXACT** | L8 |
| L7c | **`H_AMBER = E ∘ Relax₅₀`, and the relaxation is CONSTITUTIVE.** The 50-iteration cap binds on **192/192 calls**; `Spearman(Relax₅₀, Relax₁)` = 0.358–0.882 against a 0.95 threshold; and **AMBER's energy is not finite on 42% of the register without it** (58% → 97% finite). **Without the relaxation there is no objective on most of the space.** | ESTABLISHED (operator property) / INCONCLUSIVE (Spearman magnitude at n=3) | L6 |
| L8c | **The relaxation manufactures part of the Legacy–AMBER agreement**: `Spearman(E_Legacy, E_AMBER∘Relax₅₀)` = +0.229 vs **+0.140** at Relax₁ — about a third of an already-weak agreement appears only after relaxation. | ESTABLISHED | L6 |

---

## BLOCK Z — METHOD

| # | claim | label | evidence |
|---|---|---|---|
| Z1 | **The shared referent floor** — whenever two quantities are measured as deviations from a common reference, their correlation has a floor set by that reference; **measure it before interpreting the correlation.** A third class beside "control in the operator's space" and "the measurement is the null". **Nobody checked it in four sprints.** | ESTABLISHED — in project memory | L5 |
| Z2 | **A bound must be verified on the calls it actually BOUND.** Extended by measurement: the cap bound 192/192, so the operator differed on *every* call, not an unlucky one. | ESTABLISHED | L6 |
| Z3 | **Quote the CI construction with the CI.** "i.i.d. over targets" and "fold-clustered" are different intervals and this programme used the second name for the first object. | ESTABLISHED | L9 |
| Z4 | A Hessian divides energy noise by `h²`. `core.amber`'s single point has ~1e-6–1e-4 kcal/mol granularity on E≈370, so AMBER's negative-curvature fraction runs **0.259 at the plateau → 0.778 at a pre-registered h=1e-4**, while Legacy is flat to 5.6e-07. **A finite-difference Hessian on a noisy energy publishes round-off as physics.** | ESTABLISHED | L8 |
| Z5 | Two lanes stated incompatible facts about one artefact; **the file settled it in one read.** Picking the more senior lane, or averaging the claims, would have produced a wrong plan either way. | recorded as process | L10 |
| Z6 | Self-corrections this sprint: the audit lane wrote a COMPLETE flag on a smoke *after* flagging that hazard in three other files, caught and guarded it; the physics lane caught a spurious √n, a noise-dominated FD step, and an i.i.d. CI mislabelled as clustered; the quantum lane recorded that its own design rationale rested on a retracted claim. | recorded as process | L6, L8, L3 |
