# SPRINT 18 — CLAIMS REGISTER

Every row carries a label: **ESTABLISHED** · **SUPPORTED** · **PLAUSIBLE** · **OPEN** · **REFUTED** ·
**EXACT** (a theorem, not a discovery) · **ORACLE** (needs the native; never a predictive result) ·
**RETRACTED** (published here and withdrawn here).

Provenance is named on every row. Rows are added as workstreams report; blocks D–F are open.

---

## BLOCK A — the sprint question

| # | claim | label | evidence |
|---|---|---|---|
| A1 | The 19-target degree-1 numbers **reproduce exactly** (full 2.661, degree-1 2.411, space best 1.030, weight-1 variance 0.613, cumulative ≤2 0.930). | ESTABLISHED | `s18/PHASE0.md`, coordinator |
| A2 | The **inference** drawn from them was never statistically supported: deg1 − full = −0.249, **median +0.000**, W/L 7/5, 95% CI [−0.650, +0.093]; 7 of 19 targets share an identical argmin and the mean is carried by three targets (7VI4, 1CS9, 7T3H). | ESTABLISHED | `s18/PHASE0.md`, L1 |
| A3 | ρ(objective, RMSD) is **+0.264 for the full objective and +0.153 for degree-1** — the truncation is a *worse* global correlate of RMSD even where its argmin is better. | ESTABLISHED | `s18/PHASE0.md`, L1 |
| A4 | The residue-additive (gauge-invariant) degree-1 object is worth **−0.004 Å [−0.390, +0.318]**, W/L 3/5 with 11 ties — **5% of the instrument's 0.084 Å MDE**. | ESTABLISHED | ADVERSARIAL, L5 |
| A5 | The reported degree-1 advantage is **an artefact of the bit encoding.** W1 is not invariant under relabelling the k=4 torsion states; under a random relabelling it is **+0.191 Å** vs full; the codebase's own labelling sits at gauge percentile **0.00** (p < 0.0001); encoding luck = **−0.440 [−0.844, −0.133], folds 5/5**. | ESTABLISHED | ADVERSARIAL (gauge orbit), MATH (independent reproduction: 2.862 vs 2.852; RA invariant 19/19), L5 |
| A6 | **The only interval in the entire degree-1 story that excludes zero measures the arbitrariness of a bit encoding**, not a property of the objective. | ESTABLISHED | A4 + A5 |
| A7 | Phase 0's "better argmin, worse ρ" tension resolves **against** degree-1: 143% of the −0.249 is the within-band min-of-N draw (CI spans zero); the only significant component is **band quality at +0.108 [+0.011, +0.251], W/L 4/15, folds 5/5** — the wrong direction. Tie trap checked and clean (all 19 argmin_ties = 1). | ESTABLISHED | ADVERSARIAL, L5 |
| A8 | **Falsifiers F2, F4 and F5 all fire. The degree-1 branch is CLOSED.** | REFUTED (the hypothesis) | A4–A7 |

---

## BLOCK B — mathematics

| # | claim | label | evidence |
|---|---|---|---|
| B1 | The residue-additive ANOVA object equals the projection onto Walsh coefficients supported **inside one residue** — weight 0, weight 1, **and intra-residue weight 2**. Verified to **1e−12**. | EXACT | ADVERSARIAL + MATH independently, L5 |
| B2 | `s18/BRIEF.md` §4's claim that the first-order ANOVA object *is* the Walsh weight-≤1 projection under a uniform lattice measure is **FALSE**; the two differ by 1.4 objective sd. **Coordinator error**, carried into three workstreams; corrected in place. | RETRACTED | L5 |
| B3 | RA is invariant under relabelling the k=4 states; W1 is not. Any result quoted for W1 is a statement about an encoding convention. | EXACT | MATH, ADVERSARIAL |
| B4 | W1 retains **0.613** of the objective's variance, RA **0.913** — the strict truncation discards a third of the per-residue field. Orthogonality 5.5e-17, variance budget 5.0e-16. | EXACT | MATH |
| B5 | The objective is a function of the n−2 **interior** residues: f_0 = f_{n-1} ≡ 0. | EXACT | MATH |
| B5a | ~~No arm of this objective can place the two termini the frozen metric scores.~~ **Coordinator over-reading, WITHDRAWN.** The residues degree-1 cannot see are **exactly** the residues the metric cannot see — identical on **126/126 to 2e-14**. Those torsion coordinates move no Cα, so the metric does not score them either. **Nothing is lost.** | RETRACTED | EXPERIMENT retired the hazard by measurement; L11 |
| B6 | `hamil` is **75% an exactly-additive torsion prior** (prior term 1.0000 residue-additive, distogram term 0.404, blend 0.9495), and **ρ(RA(hamil), the pure torsion prior) = 0.986** with matching argmins. **The degree-1 object on the 19-target instrument is, to three nines, the retrieval-pool torsion prior** — and `s17/refine.py`'s objective contains no prior term, so the truncation has nothing to reduce to. | ESTABLISHED | MATH; the mechanism behind F4/F5 |
| B7 | **Truncating by degree is not truncating at random**, and this must be quoted beside the closure: a matched-random projection at the same retained variance gives 3.552, which W1 beats by **−1.140 [−1.714, −0.561]** (zero-info 4.003, ORACLE ceiling 1.030). The branch closes on **transfer and gauge**, not on "degree structure is meaningless". | ESTABLISHED | MATH |
| B8 | μ-sensitivity: RA under uniform/pool/rama = −0.004 / −0.066 / −0.021, every interval covering zero, argmin identical on 14/19. Pool and rama **do not factorise over qubits on any target**, so W1 is not an ANOVA under any informative μ. | ESTABLISHED | MATH |
| B9 | `s18/math_iface.py` reproduces `s17/refine._obj` to **2.3e-13**; `E_lambda(·,1) − E_full = 0.0`; analytic gradients ≤2.9e-9; Control D PASS at 1.0e-15. | ESTABLISHED | MATH deliverable |
| B10 | ∂θ*/∂d̂: at n=29 (smoke, stopped) there is **no sign of the predicted robustification and both point estimates go the other way** (+0.407 [+0.135, +0.684] unweighted; +0.212 [−0.040, +0.479] sd-scaled). | OPEN — recorded as a smoke, points against the coordinator | MATH |

---

## BLOCK C — the quantum question

| # | claim | label | evidence |
|---|---|---|---|
| C1 | The truncation **raises** the weight-1 variance fraction (full 0.6133 → RA 0.6706 → W1 1.0000) and drives **inter-residue coupling variance to exactly zero** (full 0.0868 → 0.0000). Mean Pauli weight 1.538 → 1.329 → 1.000. | ESTABLISHED | ADVERSARIAL, `s18/q_report anova`, L6 |
| C2 | The sprint's conditional quantum test therefore **cannot fire on this objective** — predicted by the coordinator in advance of the measurement. | ESTABLISHED | C1, L6 |
| C3 | Greedy 1-opt certifies the truncations' global optimum in **100% of cells at 36 evaluations** (one coordinate pass) vs 63.2% on the full objective; both truncations' optima are closed forms at ~38 table reads. RA−full at 36 evals = +0.368 [+0.200, +0.529], folds 5/5. | ESTABLISHED | ADVERSARIAL, L6 |
| C4 | Entangled vs CNOT-free arms: **null on the consumed readout at every α, on RA, W1 and full.** Run *after* the answer was forced, with that order pre-registered. | ESTABLISHED | ADVERSARIAL, L6 |
| C5 | **Q1/Q2/Q3 all fire. The quantum branch is CLOSED**, on a genuine CVaR-VQE that was neither replaced by a classical approximation nor asked to win. | REFUTED (the hypothesis) | C1–C4 |

---

## BLOCK G — the objective's error structure (coordinator)

| # | claim | label | evidence |
|---|---|---|---|
| G1 | Refining the coordinate average toward the deployed distance objective **makes it worse**: 3.048 → 3.610, +0.561 [+0.410, +0.717], 31W/95L, while the objective falls 192.8 → 56.5 (−71%). | ESTABLISHED | `s17/refine.py`, n=126 |
| G2 | **The functional form is sound — confirmed, and the original statement was understated.** With true distances the parameterisation floor is **0.083 Å** (projecting the native itself into the same ideal-geometry parameterisation; median 0.060, max 0.435), and an ORACLE start reaches objective 0.30 against the deployed start's 3.04 on the same objective. **Perfect distances make this functional point at the native.** | ESTABLISHED | ADVERSARIAL attack 3, L8 |
| G2a | The quoted **1.152 Å is a BASIN, not a ceiling**, and is withdrawn as a ceiling. α=1 reaches the optimum on most targets and a distant basin on a minority; the median is 0.307 and the paired median difference against an ORACLE start is only −0.094, so **the mean reports the minority**. | RETRACTED (as a ceiling) | ADVERSARIAL attack 3, L8 |
| G2b | **The α ladder is not a requirement curve.** Each rung mixes better target distances with a better-conditioned landscape and larger basin, so it prices objective quality and basin reachability jointly. Every downstream use of it as a requirement curve is withdrawn. | RETRACTED | ADVERSARIAL attack 3, L8 |
| G2c | **The coordinate average is not a privileged start.** A zero-information constant α-helix start (4.070 Å) reaches 1.028 at α=1, as well as the coordinate average: −0.124 [−0.366, +0.109]. A control the coordinator's arm did not carry. | ESTABLISHED | ADVERSARIAL attack 3, L8 |
| G2d | **"Halve the residual → 2.5 Å" is withdrawn as a point requirement.** At *identical* residual RMS (1.604), outcome spans **2.145–2.697 Å** across error models — range 0.552 Å, straddling 2.5. Residual RMS does not determine the outcome. Quote as "2.15–2.70 Å depending on which pairs improve". | RETRACTED (as a point figure) | ADVERSARIAL attack 4, L8 |
| G2e | **Where an improvement lands matters more than its size, and the direction inverts engineering instinct.** Improving **already-confident** pairs (−0.308 [−0.462, −0.150]) and **short-range** pairs (−0.294 [−0.424, −0.167]) beats a uniform improvement; improving **long-range** (+0.165 [+0.028, +0.303]) and **low-confidence** (+0.244 [+0.097, +0.386]) pairs is significantly **worse** than uniform. A predictor that reduces RMS by fixing its worst, longest, least-confident pairs lands at 2.62–2.70, not 2.45. | ESTABLISHED — a live instruction for predictor work | ADVERSARIAL attack 4, L8 |
| G2f | The coordinator's framing of attack 4 was itself wrong: `d_α = d_true + (1−α)(d̂ − d_true)` **is** exactly a uniform rescaling of the residual and already preserves its direction and full correlation structure. "Shrink the residual preserving correlation structure" is what the ladder does; the real objection is non-uniform improvement across pairs. | RETRACTED (the coordinator's framing) | ADVERSARIAL, L8 |
| G3 | **The distogram's errors are worse than random errors of the same magnitude.** `shuffled` 2.609 [−1.226, −0.784], `shuf_strat` 2.560 [−1.257, −0.860], `isotropic` 2.573 [−1.273, −0.802] — all carrying the correct weights — against the real 3.610. **Reproduced independently in another lane on a different start**: `shufr_wkeep` 2.635, −0.975 [−1.170, −0.783], 102W/24L. Has now survived separation-stratified, isotropic, weight-permuted, weight-flattened and an independent reimplementation. | ESTABLISHED — the sprint's most robust result | `s18/objceil.py` L2/L3; ADVERSARIAL L10 |
| G4 | The effect is **not** a separation confound: permuting residuals *within* separation bins makes it marginally stronger, not weaker. | REFUTED (the confound) | L3, attack 1 |
| G5 | ~~The harm is *which pairs* the errors fall on: `shuf_paired` reaches 2.072 Å.~~ **`objceil.py` line 163 passes `sd[pi]`, so that arm permutes the weights too and is not the deployed functional.** The 0.537 Å between `shuffled` and `shuf_paired` prices **permuting the weights**. | **RETRACTED** | ADVERSARIAL audit; coordinator confirmed at line 163; L4 |
| G6 | ~~The distogram's 1/sd² confidences may be anti-informative inside the fit.~~ **REFUTED, and the incumbent is vindicated.** On the real distogram the deployed 1/sd² weighting is the best of the three: uniform is **+0.153 [+0.066, +0.247] WORSE** (CI excludes zero, 50W/76L), permuted +0.091 [−0.017, +0.197]. **Deleting the distogram's confidence estimates costs 0.153 Å.** Native-free positive result for the shipped design. | **REFUTED** (the hypothesis) / ESTABLISHED (the incumbent's weighting) | ADVERSARIAL, L10 |
| G6a | The 0.537 Å between `shuffled` and `shuf_paired` prices **residual–weight alignment**, not weight quality: `shufr_wperm − shufr_wkeep = −0.467 [−0.598, −0.338]`, where the first lets each residual keep its own weight and the second orphans it. The correct statement is **"orphaning a residual from its weight is bad"**, confirmed on real data by `wperm_only` (+0.091). Supersedes L4's phrasing. | ESTABLISHED | ADVERSARIAL, L10 |
| G6b | With perfect distances the weighting is irrelevant: `a1_wperm` 1.162 vs `a1_wkeep` 1.152. A clean internal consistency check on the whole arm family. | ESTABLISHED | ADVERSARIAL, L10 |
| G6c | Tempering sweep, weight = sd^−p for p ∈ {0, 0.5, 1, 1.5, 2, 3, 4} plus a matched-random control, n=126. Extended past p=2 *because* `wflat` lost: if the response is monotone in p the deployed exponent may not be optimal. Native-free and deployable if an interior or higher p wins. | **OPEN — running** | ADVERSARIAL |
| G7 | ~~α=1 may be basin-limited; the ladder is not a realistic error model.~~ **Both resolved — see G2, G2a–G2f.** | CLOSED | ADVERSARIAL attacks 3 and 4, L8 |

---


---

## BLOCK H — the prior arm (coordinator; opened by MATH's mechanism, closed by its own controls)

| # | claim | label | evidence |
|---|---|---|---|
| H1 | The reproduction gate passes at full scale: `prior0` = **3.610**, identical to `s17` refine_full. The arm is that fit plus a term. | ESTABLISHED | `s18/priorfit.py`, n=126 |
| H2 | **Adding the missing torsion prior does NOT rescue refinement.** The best rung (λ=10) recovers 0.222 Å of the 0.561 Å penalty — 40% — and still ends **+0.339 [+0.241, +0.449]** worse than doing nothing. No λ beats the coordinate average. | **REFUTED** (the primary hypothesis) | L9 |
| H3 | **The pre-registered falsifier fires.** At every λ the retrieval-pool torsion prior is statistically indistinguishable from a constant ideal α-helix (λ=10: −0.024 [−0.138, +0.088]), and at λ=0.03 it is significantly **worse** (+0.031 [+0.005, +0.063]). The 0.222 Å recovered is **regularisation**, not the prior's positional information. | ESTABLISHED | L9 |
| H4 | **A wrong positional assignment hurts; the right one is worth nothing over making no positional claim at all.** The `shuf` control (same (μ,κ) tuples permuted across residues) is significantly worse than the true prior at large λ: −0.108 [−0.211, −0.006] at λ=10, −0.138 [−0.261, −0.016] at λ=30. The prior's information is detectable and unusable in the same measurement. | ESTABLISHED | L9 |
| H5 | **Two properly-powered zero-information controls in this sprint match their informative arms**: the α-helix *start* matches the coordinate average at α=1 (G2c, n=126); the α-helix *prior* matches the retrieval prior at every λ (H3, n=126). The programme separately records a constant α-helix beating the random control in torsion space. **Generic Ramachandran plausibility reproduces most of what the project's best conditioned torsion channels deliver through an objective.** Any future torsion-channel arm without a constant-helix control is uninterpretable. *(A third such control, helix-μ, is NOT MEASURED at n=30 — see D12.)* | ESTABLISHED — cross-cutting | L8 + L9, bounded by L13 |
| H6 | The coordinator's stated *reason* for expecting a gain is refuted by his own arm: displacement moves only 5.034 → 4.960 rad (1.5%) while 0.222 Å is recovered. Whatever the recovery is, it is not the displacement shrinkage predicted. | REFUTED (the mechanism) | L9 |
| H7 | The prior alone, with no distance term, is **3.886 [+0.567, +1.122]** — worse than every combined arm and than doing nothing. | ESTABLISHED | L9 |


---

## BLOCK D — the λ ladder (EXPERIMENT)

| # | claim | label | evidence |
|---|---|---|---|
| D1 | Validity gate: λ=1 reproduces `s17` refine_full on all 126 targets to mean \|Δ\| = **0.00001 Å**; gradients match central differences of MATH's callables to 3.3e-10 rel. EXPERIMENT consumed MATH's object and retired its own provisional implementation. | ESTABLISHED | `s18/exp_FINDINGS.md` |
| D2 | **Degree-1 is not equal to the full objective — it is decisively worse.** λ=0: **4.335 Å**, +1.287 [+1.124, +1.473] vs the coordinate average, 18W/108L. The ladder is **monotone decreasing in λ**, the exact opposite of the pre-registered signature of a mis-specified higher-order component. λ=0 − λ=1 = +0.726 [+0.605, +0.816], 34W/92L. **No member of the family beats the coordinate average.** | ESTABLISHED | L11 |
| D3 | The effect is **uniform, not concentrated** (drop-top-10 at the 50–52nd percentile of a uniform-effect null) and holds in every length stratum and every fold. | ESTABLISHED | L11 |
| D4 | **The mechanism (ORACLE, pre-registered before results): the truncation destroys structural information.** Fed the native's own distances, the full objective reaches 1.152 Å and degree-1 only **3.769 Å** — +0.721 [+0.498, +0.945] *worse than the coordinate average*. Perfect distances buy the full objective 2.458 Å and degree-1 only 0.566 Å. **No distogram improvement could rescue the truncation.** | ESTABLISHED (ORACLE diagnostic) | L11 |
| D5 | This **amends** "the objective is mis-aimed": the form is sound (G2), and reducing its expressive power is not a repair — it is the damage. | ESTABLISHED | D4 + G2 |
| D6 | **Phase 9 — the answer is *nowhere*.** Selection by either objective is worse than a random pool pick (3.600, 3.531 vs 3.551). | ESTABLISHED | L11 |
| D7 | **A dissociation inside degree-1.** MATH's mesh argmin is not the objective's optimum (L-BFGS beats it 85/126); reaching the true continuous optimum lowers the objective on **126/126** while moving RMSD −0.016 [−0.057, +0.032]. | ESTABLISHED | L11 |
| D8 | **H8 leverage is a null and its premise was backwards**: ρ(leverage, \|residual\|) = **−0.43** — the distogram's errors are *smaller* where geometric leverage is higher. The coordinator's leverage lead (issued, then retracted at L4) rested on the wrong sign. | REFUTED (n=74, PARTIAL — sign held from n≈45) | L11 |
| D9 | **ρ(1/sd², \|residual\|) = −0.476: the distogram's confidences are rank-ordered in the correct direction.** An independent, differently-derived explanation of why `wflat` lost by +0.153 (G6). | ESTABLISHED (n=74, PARTIAL) | L11 |
| D10 | Control B is a **null by construction**: greedy certifies a 1-opt optimum on 59/72 within budget, 4× budget improves RMSD on neither objective, degree-1's optimum is closed-form at 106 calls. **No optimiser at any budget changes any number here.** | ESTABLISHED (n=72, PARTIAL) | L11 |
| D11 | EXPERIMENT's falsifier verdicts: **F1 FIRED** (4.335, not ≈3.6) · **F2** does not fire only because degree-1 is decisively worse rather than equal · **F3 FIRED** (degree-1 does not beat a matched-magnitude random move, 4.190) · **F4 FIRED** (sign reverses at n=126). | ESTABLISHED | L11 |
| D12 | ~~The zero-information reference measure matches the conditioned one.~~ **WITHDRAWN — not measured.** At **n=30 (complete)**: refine helix 4.405 vs pool 4.141, **+0.264 [−0.341, +0.911]**, 13W/17L; argmin 4.551 vs 4.105, **+0.446 [−0.077, +1.053]**, 14W/16L. Both intervals span zero and both point estimates favour the *conditioned* measure; the estimate drifted monotonically −0.041 (n=14) → +0.264 (n=30). | **OPEN — NOT MEASURED** (would need n=126) | coordinator's n=30 run; EXPERIMENT withdrew it first; L13 |
| D12a | **A small-n point estimate with a zero-spanning CI is "not measured", never "matched".** A match claim requires a CI tight enough to exclude the effect size of interest, or n=126. Second distinct underpowered-null trap this sprint (the first being D13). | ESTABLISHED — methodological | L13 |
| D12b | The zero-information-control pattern rests on **two** properly-powered arms, not four: the α-helix **start** matching the coordinate average at α=1 (−0.124 [−0.366, +0.109], n=126, G2c) and the α-helix **prior** matching the retrieval prior at every λ (−0.024 [−0.138, +0.088], n=126, H3). | ESTABLISHED | L13 |
| D13 | **Uniform-on-the-torus is NOT a valid zero-information measure** — it places mass on impossible backbone conformations, making it a *worse* measure rather than an uninformative one. EXPERIMENT's earlier uniform-μ reading (n=7) reversed by n=17 (pool 3.918 vs uniform 4.465) and was withdrawn. **Uniform ≠ zero-information-but-plausible.** Any future zero-information control in this project must be plausible-but-uninformative. | ESTABLISHED — methodological, generalises beyond this sprint | EXPERIMENT self-correction, L12 |
| D14 | Partial arms are frozen at a single self-consistent state: Control B n=85, H8 n=88, H6 n=17, helix-μ n=14, each labelled at its own n. Control B and H8 held sign and reading from n≈45 to n≈88; **H6 did not, and is flagged rather than smoothed.** | ESTABLISHED | L12 |


---

## BLOCK E — Legacy / `leg_contact` (PHYSICS)

| # | claim | label | evidence |
|---|---|---|---|
| E1 | **`leg_contact` fails as a term in the objective.** Pre-registered primary test λ_c=+1 vs λ_c=0: **+0.108 Å [−0.166, +0.377]**, median +0.199 (worse), 56W/70L. Falsifier fires; H-C refuted. | **REFUTED** | n=126 complete, L15 |
| E2 | **The MJ table's amino-acid identities are worth nothing measurable.** Against a zero-information control with the MJ table rebuilt on permuted residue labels (form, pair set, switch, magnitude distribution all preserved): −0.116 [−0.455, +0.164], 67W/59L. **What `leg_contact` contributes is the shape of a pair term, not its sequence content.** | ESTABLISHED | L15 |
| E3 | The pre-registered sign was right and no negative-λ rescue existed: λ_c=−1 is **+1.314 [+0.958, +1.728]**, 27W/99L. | ESTABLISHED | L15 |
| E4 | **Adding `leg_contact` destroys the objective's global ordering monotonically**: ρ 0.351 → 0.282 → 0.124 → **−0.060**, while in-band ordering stays flat and insignificant. The combined objective inherits `leg_contact`'s anti-ranking, not its information. | ESTABLISHED | L15 |
| E5 | Fourth independent confirmation of the Sprint-17 instrument: `leg_contact` alone reproduces ρ global −0.176, ρ in-band +0.071, argmin 5.634 Å **to three decimals on a rebuilt instrument**. | ESTABLISHED | L15 |
| E6 | Fifth independent corroboration of the degree-1 closure, on PHYSICS's own module: `full` beats λ_c=0 by **−0.948 [−1.160, −0.753]**, 88W/38L. | ESTABLISHED | L15 |
| E7 | L30 reproduces: `full` is +0.394 [+0.255, +0.523] worse than the projected start (+0.558 against the coordinate average, vs L30's +0.561). **Every refined arm in the ladder is worse than the structure the pipeline already builds.** | ESTABLISHED | L15 |
| E8 | PHYSICS found its own pre-registration matched the terms on **pool sd**, so the realised contact share at the primary rung is **0.85–0.94** — the ladder never tested a *small* contact correction. It **declined to re-normalise post-hoc** and recorded the untested regime as OPEN instead. | ESTABLISHED (the defect) / **OPEN** (a correctly-scaled small contact perturbation) | L15 |
| E9 | **On the DEPLOYED objective `leg_contact` is harmful, not merely unsupported**: +0.722 Å [+0.541, +0.893], 38W/88L — 1.28 Å worse than doing nothing, losing on 109/126. | **REFUTED with a sign** | L17 |
| E10 | **The coordinator's ORACLE go/no-go returns zero.** corr(contact's descent direction, the correction the distogram's error needs) = **−0.006**, median −0.012, vs an MJ-shuffled null of +0.007 [−0.014, +0.029], n=126. **The ~1.0 Å available from not trusting the distogram's error direction is unreachable through this term.** | ESTABLISHED (ORACLE diagnostic) | L17 |
| E10a | PHYSICS corrected the statistic the coordinator specified: the energy contribution and the **force** have different supports (the force is exactly zero where the switch saturates) and disagree in sign at the native (+0.175 vs −0.149). The coordinator's formulation would have measured the wrong object. | ESTABLISHED — coordinator's specification error | L17 |
| E11 | **`leg_torsion` is closed by its own workstream.** Its near-native recall gain reproduces (+0.083 [+0.019, +0.135]) but **does not convert — it converts negatively** (+0.055 worse than no filter). | **REFUTED** | L17 |
| E12 | **The operator law `d_out = 1.16·d_set_mean + 0.04·d_set_best` fails under score-based selection.** Its premise held exactly for Legacy's gate (set mean −0.024, set best +0.255) yet the output was +0.076 *worse*; miss +0.094 [+0.063, +0.132] Legacy, +0.079 [+0.060, +0.098] `leg_torsion`. The law **holds for a random gate** (+0.006 [−0.003, +0.014]). It is descriptive over random subsets and must not be used to predict the effect of any gate that orders candidates. | ESTABLISHED — corrects a recorded law | L17; project memory corrected |
| E13 | Sixth independent confirmation of the degree-1 closure at full scale: `E_res` 4.265 vs `E_full` 3.606 = **+0.658 [+0.483, +0.828], 38W/88L at n=126**, where the 19-target instrument gave −0.004. | ESTABLISHED | L17 |

---

## BLOCK F — AMBER and physics-based selection (PHYSICS)

| # | claim | label | evidence |
|---|---|---|---|
| F1 | **Every physics filter loses to a random gate of the same size**, on identical structures: AMBER single point **+0.036 [+0.006, +0.061]**, `leg_torsion` +0.041 [+0.012, +0.068], `leg_contact` +0.053, Legacy **+0.063 [+0.008, +0.125]**. Five scores, every CI excluding zero, every W/L losing. | **ESTABLISHED — the sprint's hardest negative** | L17, Phase 8, n=126 |
| F2 | **The cost is the ordering, not the truncation.** Halving the candidate set costs +0.013 Å; halving it *by any physics score* costs **4–6×** that. Holds for a force field and a knowledge-based potential alike. | ESTABLISHED | F1 |
| F3 | **The repair tax is irreducible and the AMBER pass is never spent.** The spacing tax is a *displacement* tax: a minimum-displacement correction restores 3.80 Å spacing at *less* displacement than the projection (0.801 vs 0.840) and still costs +0.023 [+0.006, +0.037] more; a random move of the projection's magnitude costs the same as the projection (+0.009 [−0.014, +0.033]). | ESTABLISHED — refutes PHYSICS's own s17 prediction | L17, Phase 10 |
| F4 | AMBER k30 reproduces s17 exactly: **+0.133 [+0.112, +0.165], 24W/99L**, with the **same three gate exclusions (`1D6X 2NB7 7BX2`) for the fifth time**. | ESTABLISHED | L17 |
| F5 | Instrument integrity: G0 contact term = genuine Legacy to **1.78e-15**; G1 gradient to 3.5e-10; G2b start **bit-identical** to the coordinator's (0.00e+00); L30 reproduces at +0.558 [+0.400, +0.680] with an **identical 31/95** W/L. | ESTABLISHED | L17 |
| F6 | **OPEN, and cheap**: *why* score-ordering damages an averaged set beyond its mean and its best — the survivors' **error covariance** under a score gate. Diversity is consistent with it (Legacy 2.213 vs random 2.484 at the same m) but was not measured. | **OPEN** | L17 |

