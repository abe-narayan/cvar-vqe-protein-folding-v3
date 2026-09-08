# SPRINT 16 — CLAIMS REGISTER

One row per claim the sprint makes or destroys. The columns are the questions a reader should be
able to answer without reading the modules: **is it causal or merely observed; is it native-free,
trained, or ORACLE-only; did it generalise; and was it pre-registered or discovered.**

**Status vocabulary.** *ESTABLISHED* — measured at n = 126 (or the enumerated n = 9, marked) with
its controls, interval excluding zero. *NULL* — measured with adequate power, no effect.
*REFUTED* — a claim that was previously held and is now shown false. *CEILING* — an ORACLE
quantity; never a result. *THEOREM* — true by algebra, so measuring it validates an implementation
and nothing more. *PROVISIONAL* — interval pending recomputation.

**Provenance vocabulary.** *PRE-REG* — fixed in writing before the number was seen. *DISCOVERED* —
found in the data and labelled as such. *INHERITED* — carried in from Sprint 15.

---

## A. The flagship and its successor

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| A1 | Projecting the native-free direction onto the loud subspace makes it usable | **REFUTED** | every native-free arm chose step 0 on every fold; `loud4`/`loud10`/`loud13` all +0.000 [+0.000, +0.000] vs plain fit, n = 126 | direction yes, step trained | **PRE-REG** (rule 2) |
| A2 | The loud projection raises the direction's alignment as designed | **ESTABLISHED** | \|cos(P u, e_loud)\| rises 0.318 → **0.536** at rank 4, → 0.712 at rank 2 | native-free | DISCOVERED |
| A3 | A rank-r first-order RMSD share may be read as a finite-step budget | **REFUTED** | rank-4 carries 0.875 of RMSD to first order, yet stepping the **exact** loud component buys −0.028 [−0.106, +0.051] against −3.604 for the full error | ORACLE | **PRE-REG** (rule 4) |
| A4 | RMSD is strongly convex in torsion displacement toward the truth | **ESTABLISHED** | exact-error half-step realises **7.7% [−0.0, +14.9]** of the achievable gain (median 16.1%) and makes **49/126** targets worse than not moving; non-monotonic over the first third. Numbers confirmed by VERIFY; the geodesic objection refuted (`wrap` *is* the geodesic, 1.8e−15 rad) | ORACLE construction, but a property of the map | DISCOVERED |
| A4b | …because the fit's torsion errors are *compensating* | **REFUTED (own claim)** | a zero-information control — geodesic interpolation between two arbitrary retrieval windows, matched on ‖Δθ‖ — bulges **more** (−0.024 vs +0.014 at s = 0.5). Effect tracks ‖Δθ‖ (Spearman −0.571); below ≈4 rad both are linear | — | DISCOVERED, then retracted |
| A4c | It is a **large-displacement regime**, and the fit is 87% of a random torsion assignment (90° per-angle RMS) | **ESTABLISHED** | VERIFY, `s16/verify_steer.py` | — | DISCOVERED |
| A5 | Curvature explains it: the error is outside the linear regime | **ESTABLISHED** | first-order prediction ‖J e‖/√n = 9.879 Å against a realised 3.647 Å, **ratio 2.709** | ORACLE | DISCOVERED (diagnostic added pre-verdict) |
| A6 | Coordinate-space RMSD is exactly linear in motion toward the native | **THEOREM** — and **demoted from an exhibit** | reproduced on **400 random point-cloud pairs with no protein in them**, max deviation **4.6e−07 Å**. A sum of two symmetric PSD cross-covariances is PSD, so Kabsch never leaves the identity along the segment | ORACLE | DISCOVERED |
| A7 | ~~A native-free direction with usable cosine exists in coordinate space — REFUTED on signed c~~ | **RETRACTED (own claim): the wrong functional was quoted** | the law is quadratic and two-sided, so it consumes **\|c\|**. \|c\| toward `fit` = **0.305 vs a 0.128 isotropic null, +0.177 [+0.142, +0.212], W/L 100/26** | native-free | PRE-REG falsifier, mis-read |
| A7b | The magnitude is there; the **per-target sign** is not | **ESTABLISHED** | 58% of targets need contraction, 42% expansion; the residual `avg → native` is a compactness scalar | native-free measurement, ORACLE sign | DISCOVERED |
| A7c | A **zero-information constant β-strand** carries more direction than any channel the sprint built | **ESTABLISHED** | \|c\| **0.398 (+0.262 [+0.216, +0.310])** vs the best built channel's 0.305; pure isotropic breathing reaches 0.296. ORACLE two-sided ceiling **2.491 Å** | zero-information | DISCOVERED |
| A7d | `csteer.STEPS` had **no negative rung**, so half the instrument was unreachable | **ESTABLISHED (own defect)** | fraction of targets needing a negative step: 0.45 / 0.67 / 0.52 / 0.48. The published null is partly a property of the ladder | — | DISCOVERED |
| A7e | The missing bit is obtainable native-free | **REFUTED** | Sprint 15's `expand.json`, n = 126: four native-free scale references give sign accuracy **0.357 / 0.389 / 0.413 / 0.389** against a **zero-information constant-sign baseline of 0.516** — all worse than a constant guess; correlations ±0.04; using them costs +0.080 to +0.123 Å | native-free | DISCOVERED (ledger L22) |
| A8 | The large coordinate-space "effects" are improvements | **REFUTED (own claim, caught internally)** | −0.599 and −0.164 Å arms chose step 1.0 on every fold — they are endpoint substitutions reaching a structure the pipeline already computes | — | DISCOVERED |
| A9 | The ideal-geometry projection costs accuracy | **ESTABLISHED** | `proj` 3.213 Å against `avg` 3.048 Å, n = 126 — a **0.165 Å** accuracy cost bought for validity | ORACLE evaluation of two pipeline stages | DISCOVERED |

## B. What was inherited, and what is left of it

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| B1 | Every Sprint 15 number is correctly transcribed | **ESTABLISHED** | 11 load-bearing figures reproduced exactly from source | — | INHERITED, re-audited |
| B2 | The emitted structure's quiet subspace is target-specific | **REFUTED** | zero-information α-helix 0.829, random pool window 0.833, emitted 0.827; emitted significantly *worse* at the quarter. 4 of 5 null directions are fixed coordinate axes | — | INHERITED |
| B3 | The 1.855 Å quiet-rotation ceiling is a direction effect | **REFUTED** | its own random-direction control emits 2.164 Å; **direction is worth 16.9%** of the recorded 1.822 | ORACLE | INHERITED |
| B4 | The native-free surrogate carries pool information | **REFUTED** | zero-information α-helix reaches \|cos\| 0.357; **the pool adds +0.033 [−0.007, +0.075]**. Correct matched null is **0.216, not 0.157** | native-free | INHERITED |
| B5 | Quiet directions are free because RMSD is blind to them | **HALF-REFUTED** | RMSD blind (0.0018–0.0024) but **the objective is not** (0.018–0.031); a 20° quiet step raises the objective 53%. Ratio of ratios 10–17, > 1 on 99.2% of targets | — | INHERITED |
| B6 | The fusion law is a parameter-free law validated to 0.162 Å | **REFUTED → THEOREM** | Krogh–Vedelsby ambiguity decomposition; residual in one frame with the quadratic mean is **2.65e−15**. The 0.162 Å was **54% our arithmetic-mean error, 46% our frame mismatch, 0% the identity** | — | INHERITED |
| B7 | Structural disagreement `s` is a native-free a-priori fusion screen | **REFUTED** | **AUC 0.401 against a 0.500 null**; supervised LFO threshold degenerates to "always fuse" on 5/5 folds | native-free | INHERITED |
| B8 | Native-free error-direction estimation is novel | **REFUTED (prior art)** | ATOMRefine (2023), DeepAccNet (2021). **ATOMRefine's margin — GDT-HA 69.84 → 70.04 — is the correct prior** | — | INHERITED |
| B9 | The quiet subspace is novel | **REFUTED (prior art)** | manipulability ellipsoids, torsion-space NMA, IK self-motion manifolds, KGS nullspace sampling, concerted rotation; mechanism is Hansen's discrete Picard condition | — | INHERITED |
| B10 | The *combination* of B8 and B9 is unclaimed | **ESTABLISHED** | unclaimed across four literatures — and refuted by A1/A4 | — | DISCOVERED |

## C. The integrated architecture (enumerated instrument, n = 9 targets × 3 seeds)

> **RECOMPUTED 2026-09-06.** The original intervals in this block resampled rows, not targets.
> VERIFY has recomputed with the TARGET as the unit; intervals widen by up to 2.79x and C1 does not
> survive. Rows below carry the corrected status.

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| C1 | AMBER contributes as a ranker inside the architecture | **NOT ESTABLISHED** (was claimed; withdrawn) | at target level **−0.054 [−0.148, +0.040]**, CI 2.79× wider, 4 targets worse / 5 better, destroyed by dropping any of 8 of 9 targets, seed 2 alone gives −0.0003 | native-free | DISCOVERED, then withdrawn |
| C2 | Legacy contributes | **NULL** | +0.005 [−0.034, +0.044] at m = 75; at m = 5 the point estimate is +0.068 but the target-level interval is **[−0.039, +0.176]** — not a significant harm either. (70W/90L is a ROW count.) | native-free | DISCOVERED |
| C3 | AMBER beats Legacy directly | **DIRECTION ONLY** | −0.059 at m = 75, 97W/65L, but the interval inherits the same target-level widening as C1 and is not quoted as established | native-free | DISCOVERED |
| C4 | Legacy → AMBER composition adds over AMBER alone | **NULL** | composition −0.053 vs AMBER alone −0.054; interaction −0.004 [−0.042, +0.032] | native-free | DISCOVERED |
| C5 | Running the VQE beats not running it | **REFUTED** | every CVaR arm loses to uniform random at matched budget: +0.209 to +0.301, CIs excluding zero, 5–8W/19–22L. Reconfirmed on a **new terminal operator** | native-free | INHERITED, re-tested |
| C6 | α rescues it | **NULL** | no monotone trend; plain VQE (α = 1.0) among the better arms | — | DISCOVERED |
| C7 | The coordinate readout beats the torsion readout | **ESTABLISHED** | +0.214 / +0.297 / **+0.446 [+0.193, +0.706]** at m = 5/20/75, target level, **8 of 9 targets**; **advantage grows with ensemble size**. Survived every VERIFY control | native-free | **PRE-REG** (stated as a design claim to be priced) |

## D. The mechanism, and its lever

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| D1 | The Krogh–Vedelsby identity holds through the coordinate-average operator | **THEOREM** (verified) | residual −0.0000 [−0.0000, +0.0000], max 0.0000, at m = 5/20/75 in the operator's own frame | — | DISCOVERED |
| D2 | ~~The generators differ in diversity, not in conformer quality~~ | **CORRECTED — wrong by ≈2×** | the 0.056 Å was the **free-superposition** member error; the identity consumes the **common-frame** term, spread **0.122 Å**. Per-target attribution: cvar0.25 **52%** member error, cvar1.0 51%, cvar0.5 42% | diversity native-free; member error ORACLE | DISCOVERED, then corrected |
| D3 | The VQE's readout penalty is **roughly half** concentration and half worse conformers | **ESTABLISHED** (follows from D1 + corrected D2) | the decomposition is exact, so the attribution is exact; it is 42–52% member error depending on the arm | — | DISCOVERED, corrected |
| D4 | "The terminal operator consumes the set MEAN" fully accounts for the readout | **INCOMPLETE, not refuted** | the operator consumes the set mean **and** the set diversity, in equal and opposite weights; the earlier claim that the mean plays no part here was an artefact of reading the free-superposition frame | — | INHERITED |
| D5 | Diversity-aware selection improves the real pipeline | **NULL** | n = 126: **+0.015 [−0.019, +0.049]**, median +0.000, **61W/65L**. Raising λ lifted diversity +22% and member error +4.5%; readout moved ≤ 0.008 Å — **the terms rise together and cancel** | native-free direction, trained (m, λ) | **PRE-REG** (falsifier written before the run) |
| D6 | Diversity is a usable native-free selection signal | **NULL** | within-target rank correlation with the readout only −0.145 to −0.270; it orders the six generators, not individual runs | native-free | DISCOVERED |
| D7 | Within the top-200 by score, score adds nothing over diversity for the 200 → 75 cut | **ESTABLISHED, narrow** | `maxdiv75` 3.051 vs shipped 3.048; matched random 3.090. **Caveat inseparable from the claim**: the candidate set is itself the top-200 *by score*, so this locates the score's work in the 500 → 200 cut and is **not** a claim that the score is worthless | native-free | DISCOVERED |

## E. Sprint 16's own errors, preserved

| # | error | how it was caught | consequence |
|---|---|---|---|
| E1 | the `consensus` arm averaged torsions with an **arithmetic** mean (circular quantities) | review of the n = 6 smoke, not a test | n = 6 `consensus`/`c_loud*` numbers void and never quoted; fixed to a circular mean before the n = 126 run |
| E2 | the n = 8 smoke's `quiet_half` gave **−0.106 with a CI EXCLUDING zero**; truth at n = 126 is **−0.003 [−0.011, +0.007]** | the full run | the sixth small-n reversal in the programme and the first with a small-n CI excluding zero. **A CI at n = 8 is not protection** |
| E3 | `integrate.py`'s bootstrap resamples **rows**, so 162 rows sit on 9 targets; intervals are too narrow | self-identified and declared before VERIFY reported | every §C interval marked PROVISIONAL pending recomputation with the target as the unit |
| E4 | two `csteer` arms read as −0.599 and −0.164 Å improvements | their leave-fold-out step sat at the ladder boundary | recognised as endpoint substitutions; **general rule adopted** — a boundary-pinned step selector is choosing an endpoint, not a step |
| E5 | `divselect`'s selector pinned at λ = 1.5 on 3 of 5 folds | the trap flag was written into the module **before** the run | the arm still delivered nothing out of fold; recorded as the second boundary-pinning instance |
| E6 | `s15/align_jac.py:136` seeds with `hash()`, violating the sprint's seeding rule | AUDIT | value verified correct analytically, so no number changes; assigned for repair |
| E7 | **quoted the arithmetic mean of a signed cosine as evidence for a law that is quadratic and two-sided in it** | VERIFY | turned a sign problem into a phantom absence of signal; `CSTEER_FINDINGS` §2 retracted and rewritten. **The single most consequential error of the sprint** |
| E8 | published a **tautology** as "the sprint's cleanest exhibit" | VERIFY | reproducible on random point clouds at 4.6e−07 Å; §1 demoted. The module's own docstring already contained the proof |
| E9 | a **mechanism sentence with no control** ("the errors are compensating") | VERIFY | the zero-information control bulges *more* than the fit; retracted and replaced with the large-displacement account |
| E10 | built a step ladder with **no negative rung** for a direction whose sign flips per target | VERIFY | on 45–67% of targets no arm could have helped; the published null is partly a property of the ladder |
| E11 | read the **free-superposition** member error where the identity consumes the **common-frame** one | VERIFY | the diversity attribution was overstated ≈2×; corrected to 42–52% member error |

---

## F. The one thing that did not happen

**The 60-target sealed benchmark was not touched.** No module in Sprint 16 reads it. Unlocking
requires all five pre-registered conditions — architecture frozen, control frozen, success criterion
preregistered, the expected result genuinely changing the scientific conclusion, and no cheaper
internal instrument able to answer it. **The first condition alone fails: the architecture is not
frozen, because nothing in this sprint earned a place in it.**

## G. The causal ablation on the shipped ensemble (n = 126, `s16/energy_ablate.py`)

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| G1 | AMBER relaxation improves accuracy in aggregate | **ESTABLISHED, aggregate only** | −0.0234 [−0.0375, −0.0099], median −0.0074, 66W/57L, n = 123 gated; clears the floor by ~10 SE; concentration check PASS at the 58.5th percentile | native-free | PRE-REG (success criterion fixed before any filtered arm existed) |
| G2 | …and it is resolvable per target | **REFUTED** | the workstream's own pre-declared frame null FAILS: +0.0146 Å mean, max 0.140 Å on 11 converged targets. Band not loosened | — | PRE-REG |
| G3 | Legacy improves accuracy | **NULL** | +0.0243 [−0.0120, +0.0615]; vs its matched random filter +0.0118 [−0.0242, +0.0492]; random filter alone +0.0124 | native-free | PRE-REG |
| G4 | Legacy's null is the set-mean trap | **ESTABLISHED** | the filter improves the ORACLE set mean 3.551 → 3.531 while destroying the set best 2.306 → 2.416 and contracting diversity 6% | ORACLE diagnostic | DISCOVERED |
| G5 | Legacy and AMBER interact | **NULL** | +0.0073 [−0.0073, +0.0218] | native-free | PRE-REG |
| G6 | Legacy is a good detector | **ESTABLISHED** | AUROC **0.93–0.99** on `torsion` (not tautological), 0.98–0.99 on `steric` (partly tautological); its only CI-excluding-zero effects in the sprint are +0.032 Ramachandran and +0.069 Å min separation over a matched random control | native-free | DISCOVERED |
| G7 | Legacy would work if its weights were fitted | **REFUTED** | leave-fold-out fitting moves the certified optimum **+0.7785 Å [+0.0404, +1.4975]** the wrong way (8W/11L) while gaining in sample (ρ 0.105 → 0.213). Closes the training-status escape | — | DISCOVERED (the workstream's own hypothesis, refuted) |
| G8 | The eleven-term Legacy total beats its own best component | **REFUTED** | worse on all ten axes; MJ contact at chance; `compactness` an anti-detector; two components identically zero on all 126 inputs | — | DISCOVERED |
| G9 | Legacy's certified argmin beats random | **REFUTED** | +0.082 Å worse than random on 19 enumerations (5W/14L, native at the 53.5th percentile). **A zero-information constant α-helix beats it by 0.469 Å and IS that optimum on 8/19 targets.** MJ removal changes +0.082 → +0.090, so the defect is distributed | — | INHERITED, re-tested |
| G10 | The distance objective beats both energies at the certified optimum | **ESTABLISHED** | −0.926 [−1.440, −0.431] vs random (15W/4L, native at the 29.1st percentile); beats AMBER by −0.815 [−1.357, −0.310] on the identical masked set | native-free | INHERITED, re-tested |
| G11 | The ideal-geometry projection improves stereochemistry | **REFUTED** | it *degrades* Ramachandran 0.836 → 0.703 and degrades minimum separation | — | DISCOVERED |
| G12 | AMBER repair is stereochemically clean | **REFUTED (new defect)** | it repairs Ramachandran and clashes but **introduces cis peptide bonds on 44 of 126 targets** | — | DISCOVERED |
| G13 | AMBER's Ramachandran advantage is +0.408 | **REFUTED — overstated 2.39×** | `core.project.lam_path` returns *unwrapped* torsions (to +1022.72°) and `rama_ok` compares them to fixed windows. Corrected baseline 0.4658 → 0.7029 (+0.2371 [+0.1946, +0.2796], 0W/85L, 85/126 affected). **True advantage +0.171** | — | INHERITED, corrected |
| G14 | The control arm reproduces across sprints | **ESTABLISHED** | all 126 targets reproduce `s15/results/phys_repl_canon0.json` at **max \|Δ\| = 0.000e+00 Å**, across sprints, interpreters, a `core/amber.py` change, and four worker processes | — | DISCOVERED |


---

## H. Power, stated because the report calls five things nulls on this instrument

On the 9-target enumerated instrument with 3 seeds (`s16/review_checks.py`): standard error
**0.051 Å**, minimum detectable effect at 80% power **≈ 0.16 Å**, power to detect 0.05 Å **≈ 9%**.

**Every null in block C below about 0.16 Å is uninformative, not negative.** That includes C2
(Legacy), C4 (the composition) and C6 (alpha). It does not affect C5 or C7, whose effects are 0.21
to 0.45 Å.

## I. The QPHASE ladder (19 exhaustively enumerated targets, 1.28e7 configurations)

| # | claim | status | evidence |
|---|---|---|---|
| I1 | There is an objective-quality phase boundary rho* beyond which CVaR-VQE wins | **REFUTED** | two independent quality knobs give different rho* for the same control (+0.355 vs +0.140); a pooled fit of the paired difference on rho has **R2 = 0.014** and implies crossings of −2.74 and −1.10, outside the possible range of a correlation |
| I2 | The governing axis is `RMSD(argmin E) − RMSD(control)` | **ESTABLISHED** | Spearman +0.636, R2 0.427 pooled and **0.859** on the structured ladder; +0.867 across all 95 (objective x target) cells, against **−0.038** for rho |
| I3 | The VQE-vs-classical gap depends on objective quality | **REFUTED** | flat across ten rungs: vqe−greedy **+0.098 ± 0.037**, vqe−anneal +0.064 ± 0.035, vqe−anneal(1/4 budget) −0.043 ± 0.041; never crosses zero |
| I4 | CVaR-VQE beats a real classical optimiser at matched budget | **REFUTED** | significantly **worse than greedy 1-opt at 10 of 10 rungs**; never significantly better than annealing; **indistinguishable from annealing at 1/4 the budget at 9 of 10 rungs** |
| I5 | The VQE finds the objective's optimum | **REFUTED** | greedy reaches the certified global optimum in **68–100%** of cells, the VQE in **0–32%**; the VQE lands 90% of the way there (`diff = +0.051 + 0.904 Δ`) |
| I6 | Uniform random is an adequate control | **REFUTED** | labelled "WEAK CONTROL, proves nothing" by the workstream that used the strong ones. Report section 5.5's headline control is this one; section 5.8's supersedes it |
| I7 | The deployed objective's certified optimum beats the incumbent | **ESTABLISHED (ORACLE ceiling)** | `hamil` argmin **2.661 Å** against the incumbent's 3.204 Å and a space best of 1.030 Å — and a classical 1-opt reaches it **exactly, in 100% of cells**. **The 3.204 → 2.661 gap is not a search problem** |
| I8 | The untrained-circuit contrast has a stable sign | **REFUTED** | readout-dependent: QPHASE's argmin readout gives −0.309 (VQE better, null), section 5.5's 75-member coordinate average gives +0.112 to +0.204 (VQE worse) |
| I9 | CVaR's alpha is a diversity dial peculiar to the quantum arm | **SUPERSEDED** | a **diversity-matched classical thermostat reproduces alpha's whole ensemble effect** (Pearson +0.93, +0.98) and moves the readout **1.6-1.7x further**; the recorded 150x diversity range is a property of the optimisation BUDGET, not of alpha (1.35x here). The effect is real and thermal; the mechanism is classically reproducible |
| I10 | The CVaR gradient-baseline defect makes it a weak optimiser, and fixing it would rescue the arm | **REFUTED, three ways** | the defect is null at the target unit and significantly *worse* on the ensemble readout; its recorded "de-facto step-size reduction" mechanism is **impossible under Adam** (2.6x smaller gradient gives a 7% smaller step); cutting lr 16x costs +0.156 A (argmin) and +0.676 A (ensemble). **More optimisation is better here** |
| I11 | What the VQE actually does | **CHARACTERISED** | it **samples** (E_p[E] percentile 0.50 -> 0.15-0.25) and **represents** (3.6x near-native mass enrichment; mode 3.792 -> 3.014 A); it does **not optimise to completion** (0.5% mass on the certified optimum) and does **not explore** (3,551-4,020 distinct configurations vs annealing's 7,265 at the same budget) |

## J. AMBER as a coordinate-changing repair OPERATOR (n = 126, `s16/repair.py`)

| # | claim | status | evidence |
|---|---|---|---|
| J1 | The repair operators improve accuracy | **REFUTED** | both are worse than **doing nothing**: raw all-atom coordinate average **3.0498 Å**; projection 3.2052 (+0.1554 [+0.1217, +0.1911], **27W/99L**); AMBER k=30 3.1831 (+0.1333 [+0.1031, +0.1645], **24W/99L**, gated n=123) |
| J2 | AMBER's −0.022 Å is an accuracy result | **REFUTED** | it is the gap between two ways of being 0.13–0.16 Å worse than the operator that does nothing. Mechanically it is a **magnitude** difference: orthogonal-move cost +0.1177 (AMBER) vs +0.1515 (projection); the algebra closes to **0.0044 Å**. **AMBER wins by moving 0.09 Å less, not better** |
| J3 | AMBER's displacement direction is worth something | **REFUTED, in the wrong direction** | matched-magnitude random displacement reaches 3.1606 vs AMBER's 3.1831 (**+0.0225 [−0.0138, +0.0571]**, 43W/83L). ORACLE cos(v, r) = **−0.0521 [−0.0924, −0.0129]**, positive on only 37.3%; paired vs random **−0.0491 [−0.0996, −0.0009]**, fold-clustered [−0.0876, −0.0124] |
| J4 | Some (k, depth) rescues the accuracy effect on the frame-reproducible subset | **REFUTED** | five settings: k30-full −0.0095 [−0.0242, +0.0056] 39W/48L; d400 −0.0061; d100 +0.0117; d25 +0.0109; k0 **+0.1602 [+0.0502, +0.2724]**. **There is none** |
| J5 | The frame null passes | **REFUTED at every converged setting** | gated mean at k=30 is −0.00095 (reproducing the ≈0.004 Å floor) but **max 0.1373 Å — six times the claimed effect** — and 122/122 converged targets move by >1e−6. Only depth 25 passes, and the gate discards **76 of 126** targets there |
| J6 | AMBER repairs stereochemistry | **ESTABLISHED, but only as a conjunction** | rama 0.466 → 0.874 (−0.4077 [−0.4529, −0.3628], **116W/2L**), clashes 1.397 → 0.000 (**63W/0L**) — reproducing Sprint 15 exactly. **But a zero-information constant α-helix scores rama 1.000 and 0 clashes by construction and beats AMBER 0W/65L.** The defensible claim is *reaches 0.874 and zero clashes while moving only 26% of the residual* |
| J7 | The restraint constant k is a validity knob | **REFUTED** | k = 0 reaches the same validity (rama 0.868 vs 0.874, clashes 0.000 vs 0.000, geometry 0.0170 vs 0.0192) while costing **+0.2742 Å [+0.1785, +0.3823]**. **k trades displacement magnitude against nothing else**; the k-ladder is monotone toward "do not move" |
| J8 | The physics supplies a native-free direction | **REFUTED** | it is the only channel in the programme never tested for one, and its cosine with the truth is **negative and significantly worse than random**. Closed |
| J9 | AMBER moves where the score can see | **REFUTED** | the true residual is **60.2% loud**; AMBER's displacement is **55.1% non-torsional (1.5× enriched) and 23.5% loud (depleted)**, over-occupying even the exactly-null torsion directions (0.270 vs a 0.199 null) |
| J10 | The pipeline reproduces across sprints and interpreters | **ESTABLISHED** | max \|Δ\| = **0.000e+00 Å** on both arms in a fresh interpreter on a different day, final energies to 0.000e+00 kcal/mol |
| J11 | A shallower minimisation is the frame-reproducible sweet spot | **REFUTED (the workstream's own prior)** | depth 25 is the only arm passing the frame null but fails the convergence gate on **76/126** targets and carries 1.7 clashes with rama 0.639. **No depth is both frame-reproducible and valid** |
| I12 | The one place a VQE wins is an ensemble consumed without a ranker (+0.36 to 0.57 A, Sprint 15, low-power) | **REFUTED** | on the full enumerated register **both classical searches beat the VQE on the ensemble readout too — greedy at every rung, +0.27 to +0.62 A**. This was the programme's last standing positive quantum result |
