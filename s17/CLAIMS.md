# SPRINT 17 — CLAIMS REGISTER

One row per claim. Labels per §63 and never blurred: **ESTABLISHED** (strong evidence and controls)
· **SUPPORTED** (good evidence, not decisive) · **PLAUSIBLE** (hypothesis, incomplete evidence) ·
**OPEN** (not adequately tested) · **REFUTED** · **EXACT** (a theorem, not a discovery) ·
**ORACLE** (needs the native; never predictive).

Provenance: **PRE-REG** (fixed in writing before the number was seen) · **DISCOVERED** ·
**INHERITED** (carried in from an earlier sprint).

Instrument note: on the 126-target instrument the minimum detectable effect at 80% power is
**0.084 Å** (AUDIT A11). Any null below that magnitude is uninformative rather than negative.

---

## A. The candidate pool and its ceilings

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| A1 | The "K = 500 pool oracle of 1.711 Å" is a truncation of a 13k–27k window universe | **ESTABLISHED** | full-universe ceiling **1.313 Å**; K=500 cut costs 0.397 Å, top-75 cut costs 0.791 Å; two independent RMSD implementations agree to 4.8e−7 Å over ~2.1M windows | DISCOVERED |
| A2 | The extra windows are genuinely better structures the ranking was hiding | **REFUTED** | matched null (random permutation of the *same* universe) gives 0.473 Å where observed is 0.397; **observed − null = −0.075 [−0.132, −0.022]**. The gain is *below* the order-statistic floor | DISCOVERED (AUDIT) |
| A3 | The deep universe is near-duplicates | **REFUTED** | d\* = 1.355 Å; distinct-conformation count grows 12.2× while K grows 36.9× — a dense continuum | DISCOVERED (AUDIT) |
| A4 | Retrieval buys much over a random 500 of the same library | **REFUTED** | **+1.5 targets** of 2.0 Å recall. BLOSUM ordering has skill, but small and decaying | DISCOVERED (AUDIT) |
| A5 | A sub-2.0 Å candidate is present for most targets | **ESTABLISHED** | **81.0%** over the full universe against 46.0% in the shipped top-75 | DISCOVERED |
| A6 | The 53 "hard targets" are a real frontier | **REFUTED as sized** | 53 → **24** over the full universe (55% were truncation artefact); >2.5 Å: 27 → 8. Worth **+0.078 Å** of the mean ceiling | DISCOVERED |
| A7 | Universe size is a meaningful per-target quantity | **REFUTED** | corr(U, chain length) = **−0.988** — it is 1/n | DISCOVERED (AUDIT) |
| A8 | `u["order"]` is native-free | **ESTABLISHED** | proved by audit; the 1.711 Å oracle robust to the 11.4% tie-break (1.708 ± 0.018 over random tie-breaks) | DISCOVERED (AUDIT) |

## B. Selection — the dominant gap, closed at four levels

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| B1 | The distance objective is a real filter | **ESTABLISHED** | beats matched random by **−0.78 to −1.19 Å** at every K, all intervals excluding zero | DISCOVERED |
| B2 | …and has in-band skill | **REFUTED** | argmin within its own top-25 is **−0.014 [−0.118, +0.094], 62W/64L** against random-in-band; independently reproduced at **+0.0587 [−0.0395, +0.1606]** | PRE-REG |
| B3 | The objective can be re-engineered into a better selector | **REFUTED** | 58 functionals, leave-fold-out **−0.009 [−0.128, +0.110]** | PRE-REG (E2 falsifier fired) |
| B4 | Some native-free signal has in-band skill | **REFUTED at 29 signals** | best in-band Spearman **+0.086** (accuracy 0.529) against a ≈0.638 requirement | PRE-REG |
| B5 | Target-level calibration recovers the gap | **REFUTED** | 27 features, 5 arms: the calibrator **loses to the constant baseline** (+0.037 to +0.083) at every shrinkage; ladder boundary-pinned at maximum shrinkage | PRE-REG (E4) |
| B6 | Selector disagreement predicts selector failure | **REFUTED** | all skill vanishes when the objective's own confidence (`dist_score`, LFO 0.598) is partialled out; partial LFO −0.148 to +0.168, no consistent sign | PRE-REG |
| B7 | Widening K helps the selector | **REFUTED** | recovered fraction **negative at every width**: −8.4% / −16.4% / −17.9% | PRE-REG |
| B8 | …because widening removes the answer from the shortlist | **ESTABLISHED** | shortlist ORACLE degrades **+0.468 [+0.304, +0.642]** while its mean improves **−0.730 [−0.906, −0.560]** | DISCOVERED |
| B9 | A score-ranked shortlist retains the answer better than random | **REFUTED** | P(has sub-2 Å) null at every M; shortlist ORACLE best **significantly worse** than random at every M (+0.195/+0.242/+0.179/+0.097) | DISCOVERED |
| B10 | The recovered-fraction number needs a composition baseline | **ESTABLISHED** | the exact random pick degrades +0.529 Å over the same widening; net of drift the selector *gains* −0.388 Å [−0.486, −0.268] (fold-clustered; plain bootstrap includes zero) | DISCOVERED (AUDIT) |

## C. The readout

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| C1 | Switching from averaging to selection improves accuracy today | **REFUTED** | coordinate average dominates the medoid in **39/39 (K, m) cells**, every CI excluding zero; beats the argmin by 0.36–0.46 Å at every K. LFO: avg 3.056 / medoid 3.282 / argmin 3.454 | PRE-REG (§11 hypothesis) |
| C2 | Averaging is flat in K while selection degrades | **ESTABLISHED** | 3.056 / 3.056 / 3.064 at K = 75/500/2000 | DISCOVERED |
| C3 | The averaging readout has hidden headroom above the best-member ceiling | **SUPPORTED** | ORACLE subset average **1.598 Å at m ≈ 6** vs best-member 1.711 Å: **−0.113 [−0.240, +0.023], 90W/36L** — interval touches zero; random-greedy control −1.597 confirms the search is real | PRE-REG |
| C4 | The useful averaging subset is small | **ESTABLISHED** | optimum at m ≈ 6; the average is *worse* than the best member beyond m ≈ 17 | DISCOVERED |
| C5 | The Krogh–Vedelsby identity governs the readout | **EXACT** | per-target parallel-axis residual **1e−16** in the operator's own frame; the visible column gap is Jensen on the aggregation, not a frame error | INHERITED, re-verified |

## D. Physics

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| D1 | Legacy ranks in-band | **REFUTED** | +0.0572 [−0.0348, +0.1547], 55W/69L against matched random; vs the distance argmin +1.036 [+0.694, +1.384], 38W/88L | DISCOVERED |
| D2 | Legacy gates preserve near-native recall | **ESTABLISHED** | gate at f = 0.50: recall 0.613 vs random 0.501, **+0.112 [+0.023, +0.198]**, 49W/21L | PRE-REG |
| D3 | …without destroying the best member | **REFUTED** | set best 1.427 vs random 1.208, **+0.219 [+0.031, +0.444]**. Mechanism: in Legacy's own order the pool's best sits at 0.450 while the *median* sub-2 Å candidate sits at 0.410 (random 0.500) | PRE-REG |
| D4 | `leg_torsion` is the exception | **SUPPORTED** | recall **+0.159 [+0.073, +0.239]** with set-best damage **null** (+0.011 [−0.079, +0.128]) — the only gate measured in this programme that buys coverage without destroying the best | DISCOVERED |
| D5 | Cα-fixed AMBER removes sterics at zero Cα cost | **ESTABLISHED, magnitude corrected** | every rung reaches **0.00** clashes below 2.0 Å at Cα costs from 0 to 1.88 Å — but the gate selects a *different subset per rung*, and against its **own gated input** (`cafix` converges on the 78 least broken targets, input 0.31 clashes not 7.37) the effect is **0.31 → 0.00 (13/0)** and 2.03 → 0.01 (31/0): real, one-sided, **an order of magnitude smaller** | PRE-REG (H1a); corrected by PHYSICS in its own draft |
| D5b | `cafix` passes the rotated-frame null | **ESTABLISHED, and vacuous** | max 0.00000 — **the first AMBER arm in the programme ever to pass it**, because at ε = 0 the Cα move by 7.1e−15 Å and it buys nothing | DISCOVERED |
| D5c | Its validity cost at full n | **REFUTED (H1b)** | Ramachandran **0.906 → 0.734** (1 better / 44 worse), cis 0.000 → 0.313, ω dev 60.2°, gate fails **48/126** against the incumbent's 3 | PRE-REG |
| D6 | …and preserves torsion/peptide-bond plausibility | **REFUTED** | Ramachandran **degrades to 0.691** against a do-nothing 0.781; cis at 0.313, ω deviation 59.2°; gate fails on 15/30. **Falsifier H1b fired** | PRE-REG |
| D7 | Torsional validity can be bought cheaply in Cα displacement | **REFUTED** | Ramachandran rises 0.647 → 0.847 and cis falls 0.574 → 0.084 **monotonically in displacement**; no rung under 0.05 Å reaches the unrestrained arm's validity. **Falsifier H1c fired** | PRE-REG |
| D8 | Restraining Cα only (freeing N and C) exploits AMBER's 55.1% non-torsional displacement | **REFUTED — a null against the lane's own premise** | at matched displacement `k30` (N/CA/C, the **incumbent** set) beats `ca30` (Cα only) on accuracy (+0.122 vs +0.140), Ramachandran (0.828 vs 0.814), geometry (0.0221 vs 0.0275) **and** cis (0.086 vs 0.201) | DISCOVERED |
| D9 | The tight-Cα regime is safe for peptide geometry | **REFUTED, non-monotonically** | cis **peaks at 0.574 at `ca300`** — worse than `cafix`'s 0.313 and far worse than `free`'s 0.087 — with ω deviation peaking at 96.4°. A Cα-only restraint at high k pins the contracted spacing and lets N/C absorb it by flipping ω | DISCOVERED |
| D10 | The brief's own falsifier for the ε = 0 rung was testable | **REFUTED** | at ε = 0 the Cα coordinates do not move, so Cα-RMSD equals do-nothing **exactly, by construction**; "accuracy cost above 0.10 Å" cannot fire. Declared in `PREREG_phys.md` **before** the run | PRE-REG (the brief's error, caught in advance) |
| D11 | Level 1 — recover ≈3.05 Å with valid all-atom structures | **NOT ACHIEVED** | the repair tax cannot be eliminated; the incumbent's restraint set is the best of ten tested. The repair stage is **validated, not improved** | — |

## E. Sprint 17's own errors, preserved

| # | error | caught by | consequence |
|---|---|---|---|
| E1 | read a recorded pairwise ordering **accuracy** (null 0.500) onto the copula's **ρ** axis — accuracy 0.600 is ρ 0.309 | SELECT | the strategic read "0.600 lands at ≈2.02 Å" retracted; fixed at source with named conversion functions |
| E2 | published a **min-of-N** routing ceiling without its min-of-N null | SELECT | retracted; at n = 126 the real minimum is −0.032 [−0.067, +0.001] against its own null — at the floor, not above it |
| E3 | wrote a directional conclusion from a **5-target smoke** that reversed at n = 126 | the full run | L17 corrected by L18; averaging headroom is −0.113, not +0.019 |
| E4 | printed the **free-superposition** member error beside the **common-frame** diversity, breaking the identity by ~18% of readout² | AUDIT | fixed at source; `readout.py` already did it correctly |
| E5 | compared the coordinate average over **all of K = 500** (3.396) against the incumbent, when the right comparison is the average over the **shipped top-75** (3.048) | AUDIT | the map made averaging look 0.19 Å *worse* when it is 0.156 Å *better*; fixed at source |
| E6 | shipped a checkpointed artefact with **no completion flag**; three workstreams read it at 60/126 rows with fold 3 under-represented 2.7× | AUDIT | completion flag added |

**Every one of the coordinator's six errors was a reading or a presentation, not a computation.**
No number failed to reproduce.

---

*Blocks F (quantum) and G (features) are added as those workstreams report.*

## F. The quantum arm

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| F1 | The programme's quantum-vs-classical comparisons were budget-matched | **REFUTED** | greedy 1-opt certifies the global optimum in **100% of 152 cells at 1,024 evaluations** and 62.5% at **36**; the VQE budget used throughout Sprints 15–17 is **8,192**. Every comparison ran ≥8× past classical saturation | DISCOVERED |
| F2 | The deployed objective poses a search problem | **REFUTED — EXACT** | Walsh/Pauli-Z spectrum: **61.3% of variance at weight 1, 93.0% at weight ≤ 2**, 95.7% of weight-2 mass intra-residue. It is nearly a separable per-residue field | DISCOVERED |
| F3 | The objective's higher-degree content helps | **REFUTED** | its **weight-≤1 truncation** (closed form, zero search) has certified argmin **2.411 Å against the full objective's 2.661 Å** — better than the objective it truncates. The ORACLE truth has mean Pauli weight 3.72 | DISCOVERED |
| F4 | CVaR-VQE occupies a Pareto point classical samplers cannot | **REFUTED** | all 10 VQE configurations ε-dominated by 21 classical points, every CI excluding zero (−0.129 to −0.190 Å) | PRE-REG (P1 fired) |
| F5 | …and is uniquely bad | **REFUTED** | the classical leave-one-out null is identical (−0.139 to −0.174). **The VQE is one more classical sampler, slightly worse** | DISCOVERED |
| F6 | Entanglement contributes anything | **REFUTED** | deleting the CNOTs gives a classical product-Bernoulli VI at identical budget/estimator/optimiser/seed: **null on member error and readout at every α**, and the entangled circuit is *significantly less diverse* | DISCOVERED (control Sprint 16 never ran) |
| F7 | Local VQE beats classical local search | **REFUTED** | greedy and Metropolis certify the local optimum in 100% of windows at ¼ the exhaustive cost; VQE 68.5–72.2%, proposals worse at set-best than a random draw | PRE-REG (P3 fired) |
| F8 | VQE-informed reweighting adds information | **REFUTED** | VQE weights indistinguishable from **a permutation of their own weights** (3.299 vs 3.287 Å) and worse than doing nothing (3.135); entropy-matched Boltzmann 2.865 | PRE-REG (P4 fired) |
| F9 | VQE contributes unique elite structures in a mixture | **REFUTED** | 2.07 unique elite structures/target against the best classical arm's 5.64 | PRE-REG (P5 fired) |
| F10 | `log q_θ` carries in-band information the classical mean field does not | **OPEN** | +0.083 [+0.001, +0.179] over the exact mean-field model at α = 0.02 — but the CI clears zero by 0.001, median +0.054 ≪ mean, W/L 5/3, fold sign 2/4, **n = 8**, and it does not convert through any readout | DISCOVERED |
| F11 | α is a temperature | **SUPPORTED, refined** | `T_eff` falls 0.668 → 0.190 monotonically in physical units — but `KL(q_θ‖Boltzmann)` is 4.55–5.99 bits, so **α sets a temperature and q_θ is not that Boltzmann law** | INHERITED, refined |
| F12 | **Local refinement of a retrieval candidate is worth ~0.7 Å** | **ESTABLISHED (enumerated instrument)** | 3.636 → **2.95 Å** over 6-residue windows, certified by exhaustive enumeration at 4,096 evals/window; greedy reaches it at 1,024. Cheap, exact, native-free | DISCOVERED |

## G. Learned representations — the last open feature class

| # | claim | status | evidence | pre-reg? |
|---|---|---|---|---|
| G1 | Some representation-derived feature has in-band skill | **REFUTED** | six classes of learned sequence representation; none beats a constant ideal α-helix **and** random-in-band with an interval excluding zero, in either deployable band | PRE-REG (lane falsifier fired) |
| G2 | The recorded "ESM buys 0.288 Å on selection" is in-band content | **REFUTED** | at matched architecture ESM is **worse** than one-hot: +0.036 [−0.063, +0.139] (top-25) and +0.116 [−0.004, +0.243], fold [+0.019, +0.174] (top-75). Window-embedding cosine has ordering accuracy **0.503** against a 0.500 null. **ESM buys a filter, not a discriminator** | DISCOVERED |
| G3 | ESM's contact head carries ESM-specific in-band content | **SUPPORTED, bounded** | in-band Spearman **+0.116 [+0.047, +0.184]** at top-75; against its own ESM-free twin **+0.050 [+0.005, +0.097]**. But most is compactness (+0.000 [−0.078, +0.080] for ESM-chosen vs distogram-chosen pairs) and **none reaches the argmin in any band** | DISCOVERED |
| G4 | **The shortlist, not the ranker inside it, is the binding constraint** | **EXACT** | at B ≤ 25 the 2.5 Å target is **unreachable at ρ = 1** — a perfect ranker inside the shipped top-25 returns **2.609 Å**. And ρ_S = 0.116 delivers 3.393 Å against a band mean of 3.502: **0.109 Å against a 2.249 Å gap** | DISCOVERED |
| G5 | The ESM channel has shortlist-construction value | **REFUTED** | its shortlist is **significantly worse than matched random**, +0.315 to +0.449 | PRE-REG (F7) |
| G6 | The shipped score's shortlist beats random at retaining the best | **REFUTED — independently reproduced** | +0.209 / +0.226 / +0.170 at B = 25/75/150, against PHYSICS's +0.195 / +0.242 / +0.179 on a separate instrument | DISCOVERED |
| G7 | A mixed-channel shortlist beats the incumbent's ceiling | **REFUTED by its own control** | both CIs excluded zero — and the **pre-registered matched-random control killed it**; random does as well or better at every size | PRE-REG |
| G8 | `cf_topd_uniform` (ESM-free) is a real in-band signal | **PLAUSIBLE** | top-75 **−0.135 [−0.252, −0.022]** vs random, beating both constants and the shipped score — but **1 of ~75 uncorrected tests** and it **flips sign on the ORACLE band** | DISCOVERED, flagged not claimed |
| G9 | The ESM contact signal is length-gated | **PLAUSIBLE, post-hoc** | ρ +0.188 / −0.271 selected at length 14–16; ρ −0.047 / +0.203 at 9–11 | DISCOVERED, flagged not claimed |

## H. The decisive physics comparison, at full scale

63,000 genuine ff14SB/GBn2 single points · all 126 targets · identical K = 500 candidates · no shortlist.

| # | claim | status | evidence |
|---|---|---|---|
| H1 | AMBER can rank peptide candidates | **REFUTED** | global Spearman **−0.027** (Legacy +0.307); in-band AUROC **0.499**; argmin **5.191 Å** against random's 4.427; incremental value over the distance features **−0.004 [−0.037, +0.030]**. `angle` and `torsion` are significantly **anti-informative** |
| H2 | Sprint 16's cis-peptide defect is a defect of AMBER | **REFUTED — mechanism found** | the coordinate average **contracts the backbone 22.4%** (Cα–Cα 2.949 Å vs an ideal 3.80 Å); Spearman(input spacing, cis) = **−0.949**, and **+0.807** conditioned on convergence. **Fix the averaging operator, not AMBER** |
| H3 | `leg_contact` (MJ) is useless | **REFUTED (Sprint 16 claim G8), with a twist** | it reaches **+0.080 [+0.032, +0.128]** against a +0.010 null — real information — **and still selects 0.136 Å worse than random**. The clearest separation in the programme between *having information* and *being able to rank with it* |
| H4 | The physics instrument reproduces across sprints | **ESTABLISHED** | an independently built OpenMM System reproduces `s16/results/repair_A.json` at **max \|ΔRMSD\| = 0.000e+00 Å**, max \|ΔE\| = 2.7e−13 kcal/mol, n = 126 |
| H5 | Quoting a gated arm against an ungated input is safe | **REFUTED** | the gate selects a different subset per rung; `cafix`'s gated input carries 0.31 clashes against the instrument-wide 5.61. **Caught by the workstream in its own draft, after the coordinator had published the uncorrected version** |

## I. The objective itself — the sprint's closing measurement

| # | claim | status | evidence |
|---|---|---|---|
| I1 | Refinement toward the distance objective improves the pipeline's best structure | **REFUTED** | from the coordinate average (3.048 Å), full refinement lands at **3.610 Å**: **+0.561 [+0.407, +0.712]**, 31W/95L, while the objective falls **192.8 → 56.5 (−71%)** |
| I2 | The move shape matters (windowed vs global) | **REFUTED** | w = 4 / 6 / 8 land at 3.645 / 3.616 / 3.611 — within 0.035 Å of the full fit. The enumerated instrument's 0.7 Å local-refinement gain **does not transfer** |
| I3 | The loss is random damage from moving | **REFUTED** | matched-magnitude random move +1.499 [+1.217, +1.786]; refinement against a **shuffled** distogram +2.346. Refinement is purposeful and competent — aimed at the wrong place |
| I4 | **The objective's optimum is worse than the structure the pipeline already builds** | **ESTABLISHED** | by 0.561 Å. **The pipeline works because it does not optimise its own scoring function** |
| I5 | The objective is weak | **REFUTED — it is wrong** | a weak objective would be uninformative; this one is informative and mis-aimed. It unifies the in-band null, the 58-functional null, and the fact that its weight-≤1 truncation has a *better* argmin than the whole |
