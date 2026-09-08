# SPRINT 14 — live research ledger

Status of every claim in play. Updated as evidence arrives. Tiers:
**DEMONSTRATED** (legitimate inference-time information, CI + null) ·
**ORACLE DIAGNOSTIC** (reads the native; prices a ceiling; never a headline) ·
**LITERATURE-SUPPORTED** (with citation) · **HYPOTHESIS** · **REFUTED** · **OPEN**.

---

## The arithmetic the sprint has to beat

| quantity | value | source |
|---|---|---|
| incumbent retrieval pipeline | 3.204 A | pinned instrument, `synthesis_fit` |
| target | < 2.0 A | sprint brief |
| sigma needed for 2.0 A, i.i.d. torsion error | **15.1 deg** | S14 C4, ORACLE |
| sigma equivalent of the incumbent, i.i.d. | **27.1 deg** | S14 C4, ORACLE (refines S13's ~29) |
| best torsion channel available | phi 33.6 / psi 59.2 deg | S14 C2, DEMONSTRATED |
| best trained sequence predictor (S13) | phi 36.1 / psi 62.4 deg | S13 |
| torsion space ceiling at k=4, random start | 1.982 A | S13, ORACLE |

**The gap is not subtle.** The best torsion information obtainable from sequence and
retrieval is roughly twice as coarse as what 2.0 A requires. No aggregation trick closes a
factor of two. Either an external information channel supplies it (chemical shifts), or the
objective does work the information alone cannot.

---

## DEMONSTRATED

| # | claim | evidence |
|---|---|---|
| C1 | The shared ladder reproduces every Sprint 13 anchor; the leakage guard passes on all arms | constant helix 4.065, class prior argmax 3.969, gap -0.096 as recorded |
| C1b | Zero-information arms beat the incumbent on the failure class | helix FAIL18 5.887, class prior 5.670, incumbent 6.026 |
| C2 | The retrieval pool is the best torsion predictor in the project, beating the trained sequence model with no training | phi 33.6 vs 36.1, psi 59.2 vs 62.4 |
| C2b | And it still loses: no native-free torsion emitter beats the pipeline consuming the same windows | best 3.514, +0.310 [+0.128,+0.492]; all six CIs exclude zero |
| C3 | Mean absolute torsion error does not order emitted accuracy | 4.9 deg worse phi, 0.56 A better structure, same window set |
| C4 | Real torsion predictors have near-i.i.d. errors, lag-1 in [-0.092, +0.066] | ten emitters, 126 targets |

## ORACLE DIAGNOSTIC

| # | claim | evidence |
|---|---|---|
| C4a | Error coherence is worth nearly a factor of two in required accuracy | at sigma 12: AR(-0.7) 1.272, i.i.d. 1.617, AR(+0.7) 2.203 — a 0.93 A spread at matched sigma |
| C4b | sigma to reach 2.0 A ranges 10.6-20.2 deg depending only on coherence mode | 126 targets, 12 reps |

## LITERATURE-SUPPORTED

| # | claim | evidence |
|---|---|---|
| C5 | The Sept 2026 continuous-torsion VQE preprint is arXiv:2609.02113 (QTF, Cumbo et al.) | primary text fetched and read in full |
| C5a | Its headline 0.623 A is oracle-selected; its predictive medians are 2.90 A chignolin / 5.39 A Trp-cage on two targets, terminal-excluded | their tables |
| C5b | It independently corroborates our central negative: AMBER energy negatively correlated with RMSD, over 7.3 M structures | their Spearman analysis, no numeric rho printed |
| C5c | Its energy is never a qubit observable and its parameter count exceeds the torsion count (132 params for 43 torsions) | their Methods; the log-qubit claim is reparameterisation, not compression |
| C5d | CVaR is entirely absent from it | zero hits in full text |

## DEMONSTRATED (continued)

| # | claim | evidence |
|---|---|---|
| C6 | The cost of a torsion error is a symmetric hump peaked at mid-chain; middle 40% carries 62.8%, outer 40% carries 15.7% | 4.0x ratio; terminal tenth ~14x cheaper; explains the S13 terminal-dropout correction mechanically |
| C7 | A native-free STRUCTURAL objective orders the space better than any molecular energy — corrected downward by its own control | in-decile **+0.161** under a uniform proposal vs Legacy +0.043; ~4x, distogram-driven |
| C7b | The retrieval conditioning carries genuine target-specific information, not just generic Ramachandran | sequence-blind twin drops in-decile 0.355 -> 0.073, decile mean 3.958 -> 4.658 |
| C8 | **Searching a good objective harder buys nothing.** Search saturates by budget 300 of 20,000 while the selection gap grows to 2.04 A | 3.742 / 3.565 / 3.572; degradation past minimum +0.007 A; objective improves monotonically throughout |
| C9 | Aggregating a low-energy SET in coordinate space recovers 0.29 A over the argmin — and still loses | argmin 3.608 -> consensus at m*=200 **3.314**, vs incumbent +0.110 [+0.004,+0.214] |

## REFUTED (including my own)

| # | claim | killed by |
|---|---|---|
| — | *My hypothesis:* real torsion predictors have coherent errors, so the i.i.d. restraint surface is mispriced | C4. They are near-i.i.d. or mildly anti-correlated. The surface stands and is mildly conservative. My C3 warning is **withdrawn**. |
| — | *My hypothesis:* torsion-error cost is N-terminally dominated | C6. Front/back ratio 0.886 — no privileged terminus. Kabsch superposition is free to rotate the whole structure, so the hinge argument replaces the lever argument. |
| — | *My headline:* the structural Hamiltonian is "eight times better ordered" | C7-CORRECTION. Measured under a proposal drawn half from the prior it scores. Under a uniform proposal the w=0 term collapses from +0.355 to +0.046. Honest figure is +0.161, roughly fourfold, and it comes from the distogram not the prior. |
| — | *Prior sprint framing:* "AMBER is a less local observable than Legacy" | S13; both are full-register. Category error. |
| — | Coherence explains the C3 MAE-versus-RMSD anomaly | C4; the two arms have near-identical coherence. |
| — | Position explains the C3 MAE-versus-RMSD anomaly | C6; position-weighting makes emitter ordering *worse* (rho +0.370 vs +0.467 plain). **The anomaly remains unexplained after two failed explanations.** |
| — | *ENER's hypothesis:* `leg_torsion` selects 2.954 A and is the best native-free selector | C12 at n=126. Worst arm tested. Was concentration on two of nine targets — and their own concentration analysis called it in advance. |
| — | *ENER's own:* term reweighting helps (in-sample +0.209) | their leave-one-target-out: **-0.031**. |
| — | *ENER's own:* "AMBER refinement rescues bad structures" | +0.014 A [-0.001,+0.028]; from 4.334 A it returns 4.337 A. |
| — | *ENER's own:* Legacy's energy-structure decoupling rides on the invisible terminal DOFs | only 3.8% of variance. |
| — | H7 as **exploitable**: Legacy and AMBER are complementary | premise correct (truth-partialled error correlation +0.096) but fusion is worth +0.032 rho under a per-target ORACLE weight and **loses 0.240 rho in-sample to LOTO**; `rank(Leg)+rank(AMB)` is 0.015 A worse than Legacy alone. Pareto refuted too: frontier discards 72% of the truly-best 1%. |
| — | *Prior sprint framing:* the log-qubit continuous-torsion encoding is a compression | LIT/C5: `P ~ 2N+6n`, 132 classical parameters for 43 torsions. It is a reparameterisation, and adopting it would breach the brief's rule 1.4. |
| — | *Mine:* "the emitted structure does not move under more search" as a GENERAL claim | C8-CORRECTION. True on the 126-target instrument, FALSE on the enumerated nine, where it falls 0.38 A monotonically. Target-length dependent; the discrimination half reproduces exactly and is unaffected. |
| — | *Mine:* AMBER cannot recover the 0.157 A projection cost | Invalid inference — the refinement experiment ran on structures already ON the ideal manifold, where AMBER has nothing to fix. |
| — | *Mine:* the 0.157 A is the price of ideal geometry | **CONFIRMED at n=126:** it is the price of re-expanding a **-25.8%** contracted structure, rho(contraction, cost) = **-0.552**; any operation restoring real bonds pays it. (The estimate was unstable in flight: -0.829 at n=6, -0.224 at n=10; only the n=126 value is quotable.) |
| — | *Prior sprint framing:* QNG is refuted for this ansatz | True and complete at depth 1 for every entangler, and **does not extend to depth >= 2**: condition number 155.2 for a torsion-aware block entangler at depth 2, 478.9 at depth 3. QNG is open again in exactly the regime a problem-inspired ansatz occupies. |
| — | *Prior sprint framing:* "coverage is the make-or-break parameter", close below 75% | SHIFT: measured gaps cost 0.014 A at sigma 12, inside seed noise. The binding quantity is the fraction of TARGETS with heteronuclear shifts (0.444), not residues within them (0.952). |
| — | The 1.486 A restraint figure as the channel's value | Superseded: under TALOS-N's real error mixture the channel lands at 2.347 A deployable, and the real built model emits 3.832 A. |
| — | *VQE agent's own:* Gray coding lowers mean Pauli weight below binary | Replication over 54 cells: delta +0.006 and +0.058, Gray lower in 16/27 and 13/27 — a coin flip with sign-changing per-target swings. |
| — | *VQE agent's own:* "alpha is effectively a learning rate" | Gradients at different alpha at the SAME point are not parallel. |
| — | *ENER agent's own:* an 86% bond contraction in the averaged backbone | A frame bug caught by their own runtime assertion against another agent's number. Corrected figure 16-19%. |

| C11 | The ablation ladder closes. Aggregation is worth 3.4x the objective, and the objective's entire contribution is 0.171 A — already banked classically | 4.072 -> 3.485 (-0.587, aggregation) -> 3.314 (-0.171, objective) vs incumbent 3.204 |
| C12 | **`leg_torsion` REFUTED at n=126** — the sprint's strongest open lead was two targets | worst arm tested: argmin 4.933, in-decile -0.073, +1.729 [+1.419,+2.052]; poisons every combination |
| V1 | No arm, classical or quantum, has a selection-gap advantage; CVaR buys no discrimination | every arm returns ~1 A worse than its best evaluated; tenfold budget does not close it |
| V2 | **A 0.68 A selection gap survives INFINITE budget** on a certified landscape | certified global optimum 2.287 A vs 1.605 A best inside its own lowest decile |
| E1 | Legacy is a clash gate, not a ranking function | 98.8% of variance is `steric`, exactly zero on 73.5% of the space, constant across its own lowest decile |
| E2 | AMBER is a genuine physical validator and an inert refiner | best clash detector measured; refinement +0.014 A [-0.001,+0.028], 13W/32L |
| E3 | **The discrimination floor**: no objective exceeds 0.511 pairwise accuracy below a 0.25 A quality gap | Legacy needs 1.31 A, AMBER 1.58 A; useful range ~2.5 A wide, so ~1 bit resolved |
| E4 | Four torsions are inert, not three — both terminal residues invisible; dead block is `2*log2(k)` | 14 of 18 live at n=9,k=4; k=4 operating point 21.9 mean live qubits, not 25.9 |
| E5 | The cached AMBER subset is 40% oracle-conditioned | 0.401 A better than its space; only `amber_kind==0` usable |
| C13 | **CERTIFIED: the structural objective's global optimum is 0.885 A BETTER than random** where Legacy's is 0.139 A worse | full 262,144-config enumeration, 9 targets; "the optimum is in the wrong place" was a property of the ENERGIES |
| C13b | The two terms invert: the distogram orders the bulk, the torsion prior places the optimum | below w=0.25 the distogram's contribution to the argmin is ZERO on 7/9 targets; it is worth 0.287 A on the decile mean |
| C13c | The selection gap reproduces across disjoint target sets and survives certification | 0.006 A agreement at budget 300; 1.876 A gap at the certified optimum |
| C14 | **The chemical-shift route is closed by ARITHMETIC**, not by accuracy | 54/126 targets runnable; ORACLE-perfect torsions on all of them leave the instrument at 2.021 A; 55 needed |
| C14b | **A bimodal shift posterior is one qubit with a physical justification** — the first in this project | ORACLE search over the top-8 support vs the argmax: -2.253 A [-2.642,-1.865], 53W/1L |
| C14c | Confidently wrong costs 2-3x what absent costs | 5% error 2.459 vs 1.890; abstain rather than guess, and let a search consume abstentions as free registers |
| V3 | **VQE is the WORST optimiser tested on its own axis** | fraction reaching the certified optimum at budget 30,000: greedy 0.71, annealing 0.71, VQE 0.00-0.12; at signal 0 it is 100% vs 0% |
| V4 | Two NEW CVaR defects, plus a closed form for the known one | sampled CVaR biased upward at non-integer alpha*N (+0.134 sd at N=13); `dCVaR/dp` identically zero **iff** `p(x*) >= alpha`, 0 counterexamples in 3,000 |
| V5 | Plain expectation-value VQE returns the best structure; low alpha buys entropy, not accuracy | alpha 1.0 -> 2.266 A; alpha 0.01 -> 2.459 A with entropy 4.07 -> 10.26 bits |
| V6 | Binary encoding, for four independent reasons | 510 Pauli terms vs one-hot's 454,463; surjective; highest gradient variance; attains the qubit bound |

## CLOSED — every question this sprint opened has an answer

| # | question | answer |
|---|---|---|
| O1 | Does a learned many-body objective have an accurate low-energy region? | **No.** -0.166 A, CI crosses zero, reverses on drop-top-3; both nulls worse than random; leaked-label control loud at -1.405 A, 19/0. Learning curve **declines** from m=2. |
| O2 | At what objective quality does VQE/CVaR beat classical search? | **Never, on either axis.** No crossing at any signal level; VQE reaches the certified optimum in 0-12% of cells against greedy's 71%; and it loses to best-of-N from its own untrained circuit, 0/12. |
| O3 | Do chemical shifts exist at usable coverage? | **54/126 targets.** ORACLE-perfect torsions on all of them still leave the instrument at 2.021 A; 55 are needed. Closed by arithmetic. |
| O4 | Each energy model's role; the decoy threshold | Legacy = clash gate (98.8% one term, constant in its own decile). AMBER = validator, refinement inert (+0.014 A). No **physical** objective exceeds 0.511 below a 0.25 A quality gap. |
| O5 | Does coordinate-space aggregation recover what one torsion vector loses? | **Yes, and it is the largest single operator in the system:** +1.024 A over torsion averaging of the same windows; worth 3.4x the objective. |
| O6 | Does a native-free structural Hamiltonian have a good low-energy region? | **Yes** — certified optimum 0.885 A better than random where Legacy's is 0.139 A worse. And it still does not help: the selection gap survives certification. |
| O7 | Does torsion-error cost fall with chain position? | **It is a symmetric mid-chain hump** (4.0x middle-vs-outer), not an N-terminal slope. Explains the S13 terminal-dropout correction mechanically. Position-weighting does NOT rescue MAE as an ordering statistic. |

## STILL OPEN — carried to the next sprint

| # | question | why it matters |
|---|---|---|
| N1 | **Resolve the per-target-skill / compactness discrepancy** on a common skill definition and a larger enumerated set | +0.909 on one definition, +0.416 (p=0.177) on another. C22's chain rests on it and cannot be quoted until it is settled. |
| N2 | **Pre-register a restraint constant** and re-measure the -0.022 A | The selection rule is degenerate; a pre-registered k could land anywhere from +0.001 to -0.022. |
| N3 | Anything supplying the **per-target sign** of the in-band ordering axis at inference | The only direction with leverage. Native-free compactness proxies reach 0.24-0.37 (all CIs exclude zero) against the oracle's 0.909. |
| N4 | A **forward chemical-shift predictor** scoring candidates against the 54 targets' measured shifts | A per-target observable available at inference, and a *discrimination* channel rather than the closed generation channel. |
| N5 | Whether QNG **helps** at depth >= 2 | Only the metric's condition number was measured (155.2 at depth 2, 478.9 at depth 3), never an optimisation using it. HYPOTHESIS. |
| N6 | The out-of-distribution generator-shift test | Written, not run. Does a learned objective survive a VQE's own proposal distribution? |

## HYPOTHESIS — stated, not yet tested

| # | claim |
|---|---|
| H-a | A VQE state should not be decoded to one configuration. Measured B times, built, and averaged in COORDINATE space, it pays the compounding cost per sample and cancels it — so the quantum state chooses a distribution and the consensus extracts the structure. |
| H-b | The correct role for Legacy and AMBER is validity filtering and refinement, never the search objective. |
| H-c | The distogram and the torsion prior are complementary (non-local versus local) and their sum has a better-ordered low-energy region than either alone. |

---

## Abandoned, and why

- **Projected set consensus** (`s14/consensus.py` with `project=True`, sizes 20/50/200).
  Launched, ran two hours under heavy CPU contention, **stopped deliberately before
  completion**. Its question — does projecting the Hamiltonian-selected consensus help or
  hurt — was superseded mid-flight by the ENER agent's strictly better experiment, which
  reconstructs a genuine all-atom average with non-ideal bonds and compares projection
  against AMBER relaxation with a geometry audit. Holding a core for a marginal duplicate
  while three agents owned the sprint's remaining open questions was the wrong trade. The
  raw (unprojected) consensus result at m*=200 is complete and is the one reported.

## Controls owed

- `s14/hamil.py` draws half its configurations from the prior it also scores. A **pure
  uniform proposal** control is required before any in-decile rank correlation from that
  file is reported, or the correlation may partly measure the proposal.
- Any argmin-over-N arm must be reported with its evaluation budget beside the incumbent's,
  since best-of-4000 and a single deterministic pipeline pass are not the same budget.
- Every new arm must clear the **constant alpha-helix at 4.065 A**, not uniform random —
  "beats random" is not evidence in this representation.
