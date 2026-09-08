# SPRINT 17 — FINAL REPORT

**Programme**: quantum-assisted peptide structure prediction, 9–16 residues, full-chain Cα-RMSD
**Instrument**: 126 cluster-disjoint targets, five folds · **Date**: 2026-09-06
**Sealed benchmark**: 60 targets, **untouched** — verified by an independent audit, no `s17/` module
reads it or derives its identities

> *Written to be checked, not to persuade. Every positive claim carries its controls; every failed
> hypothesis is here; every retraction is preserved with the evidence that forced it. ORACLE
> quantities say so in the same sentence. All six workstreams have reported and are integrated.*

---

## 0. The answer in one page

**The sprint set out to switch the architecture from an averaging readout to a selection readout,
on the strength of a 1.711 Å pool ceiling against a 3.204 Å incumbent. It found the ceiling is real
and larger than believed, that the gap to it is entirely a ranking gap, and that the ranking
required does not exist in anything the programme can currently measure.** Selection is now closed
at four independent levels. The switch would be a regression today.

**The measured picture, all at n = 126, target as the unit:**

| quantity | value |
|---|---|
| incumbent (deployed) | **3.204 Å** |
| coordinate average over the shipped top-75, before repair | **3.048 Å** — and 3.048 + the 0.155 Å repair tax = 3.203 ≈ the incumbent |
| distance-objective argmin at K = 500 | 3.454 Å |
| consensus medoid at K = 500 | 3.282 Å |
| **ORACLE best member, shipped top-75** | 2.104 Å |
| **ORACLE best member, K = 500** | 1.711 Å |
| **ORACLE best member, full universe** | **1.313 Å** |
| **ORACLE best subset average, K = 500** (m ≈ 6) | **1.598 Å** |
| sub-2.0 Å candidate present, full universe | **81.0% of targets** (top-75: 46.0%) |

**The gap is 2.249 Å at the full universe. It is entirely candidate choice — and the sprint's
central result is that the shortlist, not the ranker inside it, is what binds.**

> **At B ≤ 25 the 2.5 Å target is unreachable at ρ = 1.** A *perfect* ranker inside the shipped
> top-25 returns **2.609 Å**. The programme has been trying to build a better ranker for a
> shortlist that does not contain the answer. Three workstreams reached that conclusion
> independently, from three different directions, and the third arrival is **exact**.

### The findings that matter

1. **The pool contains the answer far more often than believed.** A sub-2.0 Å candidate exists for
   **81.0%** of targets over the full universe against 46.0% in the shipped top-75, and the "53
   hard targets" fall to **24** — more than half were a truncation artefact.

2. **But widening the candidate set makes the realized answer WORSE.** The recovered fraction of
   the extra ceiling is **negative at every width**: −8.4% at K = 500, −16.4% at 2000, −17.9% at
   the full universe.

3. **The mechanism: widening removes the answer from the shortlist.** As K grows, the ORACLE best
   *inside the score's own top-75* degrades **+0.468 [+0.304, +0.642]** while the shortlist mean
   improves **−0.730 [−0.906, −0.560]**. Extra windows look more plausible without being nearer the
   native and they displace near-native members. This is Sprint 16's set-mean trap on a new axis.

4. **And it holds at fixed K too.** A score-ranked top-M shortlist is **no better than a matched
   random shortlist at containing a sub-2 Å candidate** (null at every M) and **significantly worse
   at retaining the best one** (+0.195 / +0.242 / +0.179 / +0.097 Å at M = 25/75/150/300).

5. **The distance objective is a garbage filter and nothing more.** Global rank correlation 0.568,
   **in-band 0.131**; its argmin *within its own top-25* is **−0.014 [−0.118, +0.094], 62W/64L**
   against picking at random from that same top-25 — reproduced independently by a second
   workstream at **+0.0587 [−0.0395, +0.1606]**.

6. **Selection is closed at four independent levels.** 58 functionals of the distogram
   (LFO −0.009 [−0.128, +0.110]) · 29 native-free signals measured in-band (best Spearman +0.086
   against a ≈0.638 requirement) · 27 target-level calibration features (the calibrator **loses to
   the constant baseline** at every shrinkage, ladder boundary-pinned) · and every distance/Legacy
   disagreement feature, whose skill **vanishes entirely** once the objective's own confidence is
   partialled out.

7. **The readout triangle contradicts the sprint's central architectural move.** On identical
   shortlists the coordinate average dominates the medoid in **39 of 39 (K, m) cells**, every
   interval excluding zero, and beats the distance argmin by 0.36–0.46 Å at every K. **Switching
   from averaging to selection is a 0.23–0.46 Å regression today.** The ceiling argument stands;
   the immediate move does not.

8. **Legacy likes near-natives and likes the best one least.** Its gate raises near-native recall
   **+0.112 [+0.023, +0.198]** while making the retained set's best **+0.219 [+0.031, +0.444]
   worse**. In Legacy's own order the pool's single best candidate sits at 0.450 while the *median*
   sub-2 Å candidate sits at 0.410 (random = 0.500). That is the set-mean trap resolved to a single
   candidate.

9. **Every quantum-versus-classical comparison in three sprints ran past classical saturation.**
   Greedy 1-opt certifies the objective's global optimum in **100% of cells at 1,024 evaluations**
   and 62.5% at **36**; the CVaR-VQE budget used throughout is **8,192**. *A budget the control does
   not need is not a budget.* And the reason is exact: the objective is **93% separable** by Walsh
   spectrum while the truth sits at mean Pauli weight 3.72 — **there is no search problem here for a
   quantum device**.

10. **The incumbent was not merely left standing — two of its stages were validated.** The
    coordinate average beats both alternative readouts in 39 of 39 cells, and the shipped restraint
    set (N/CA/C, k = 30) beats all nine alternatives on accuracy, Ramachandran, geometry *and* cis
    at matched displacement. **Steric validity is free at any restraint strength; torsional validity
    is a strict monotone trade; the 0.155 Å repair tax cannot be eliminated.**

11. **Selection is closed at five independent levels**, the fifth being six classes of learned
    sequence representation. The recorded "ESM buys 0.288 Å" is **global** skill: at matched
    architecture ESM is *worse* than one-hot in-band. **ESM buys a filter, not a discriminator** —
    the same object the distance objective turned out to be.

---

## 1. The gap ledger, in brief

| gap | size | character | status |
|---|---|---|---|
| **retrieval** | 1.313 Å remains at the full-universe ceiling | 24/126 targets still have no sub-2 Å candidate; the residual hard set is longer, less helical, more polar, with weaker BLOSUM matches | worth only **+0.078 Å** of the mean ceiling — a much smaller frontier than the sprint assumed |
| **filter / shortlist** | 0.791 Å of ceiling destroyed by the top-75 cut | the shortlist is no better than random at keeping the answer and **worse** at keeping the best | **open, and mis-specified** — the filter optimises the mean and the readout consumes the best |
| **selection** | **2.249 Å** | the dominant gap; the objective has no in-band skill at all | **closed at four levels** |
| **ensemble** | averaging beats selection by 0.36–0.46 Å; ORACLE subset average 1.598 Å vs best member 1.711 Å | headroom of −0.113 [−0.240, +0.023], 90W/36L — real but modest, interval touches zero | **SUPPORTED, small** |
| **repair / validity tax** | **0.155 Å** | 3.048 + 0.155 = 3.203 ≈ the incumbent | **PENDING** — the Cα-preserving lane |
| **target calibration** | +0.294 Å of headroom over a global fit (ORACLE) | but the calibrator loses to a constant baseline | **REFUTED at this feature class** |
| **quantum** | — | — | **PENDING** |
| **physical models** | Legacy null as ranker and as gate; `leg_torsion` the sole exception | AMBER pending | partially closed |

**Which gap is worth attacking next: the shortlist.** It is the only one that is both large
(0.791 Å of ceiling) and demonstrably mis-specified rather than merely hard — the pipeline
optimises a shortlist for its mean and then reads out something that consumes its best. Selection
is closed; retrieval is nearly exhausted; the ensemble headroom is 0.11 Å.

---

## 2. The design equation, and why the requirement is out of reach

`s17/selection_theory.py` derives what the programme never had: given a candidate set whose
member-RMSD distribution is measured, and a selector whose ranking has copula correlation ρ with
the truth, what does the selector return? Gaussian copula, one dependence parameter, **both limits
verified numerically** — ρ = 1 reproduces the ORACLE best exactly, ρ = 0 reproduces random
selection exactly.

**The requirement, n = 126:**

| K | ρ for 2.5 Å | ρ for 2.2 Å | ρ for 2.0 Å | ORACLE (ρ = 1) |
|---|---|---|---|---|
| 75 | 0.80 | 0.95 | unreachable | 2.104 |
| 500 | 0.65 | 0.80 | 0.90 | 1.711 |
| 2000 | 0.60 | 0.70 | 0.80 | 1.504 |
| full | **0.60** | **0.65** | **0.75** | 1.313 |

**Widening K does not merely raise the ceiling — it lowers the ranking skill required**, from
unreachable at top-75 to ρ = 0.75 for 2.0 Å over the full universe.

**And the measured skill is nowhere near it.** The shipped objective's in-band Spearman is 0.131
(copula ρ 0.137); the best of 29 native-free signals reaches Spearman 0.086. **Three different
axes** are in play and conflating them is the programme's recurring failure — pairwise ordering
accuracy (null 0.500), Spearman, and copula ρ. The conversions are named functions in the module:
`ρ = sin(π(2·acc − 1)/2)`, so a recorded accuracy of 0.600 is **ρ 0.309, not 0.600**.

---

## 3. The architecture, as the evidence leaves it

**The incumbent is internally coherent, and that is a finding rather than a compliment.** The
distance objective improves the shortlist's *mean* (the average over the score's top-300 is 3.171
against 3.402 over a matched random 300) and the coordinate average *consumes* the mean. Filter and
readout are matched. The 0.791 Å of shortlist ceiling the score destroys costs nothing to an
averaging readout, because averaging never consumes the best member — it is only fatal to the
selection readout the sprint proposed switching to.

So the pipeline is not misassembled. It is **operating at the ceiling of the operator it uses**,
and the operator it uses is the right one given the signals available.

**But both readouts have the same enormous gap, and it is a candidate-choice gap either way.**

| readout | realized | ORACLE ceiling | gap |
|---|---|---|---|
| coordinate average, shipped top-75 | **3.048** | **1.598** (best 6-member subset, K = 500) | **1.450** |
| single-candidate selection, K = 500 | 3.454 | 1.711 | 1.743 |
| single-candidate selection, full universe | 3.562 | 1.313 | 2.249 |

Choosing the best six members to average is worth 0.113 Å more than choosing the single best
member, and no more — so the two problems are of comparable difficulty and comparable reward.
**Neither is solved, and the same absence blocks both: nothing the programme can measure
discriminates among plausible candidates.**

### What would have to be true to reach the targets

From the design equation, over the full universe: **ρ = 0.60 for 2.5 Å, ρ = 0.75 for 2.0 Å**. The
best in-band correlation measured anywhere in this sprint, across 29 native-free signals plus 58
functionals of the distogram plus 27 target-level features, is **Spearman 0.086**. That is not a
near miss; it is roughly an order of magnitude short on the axis that matters.

### The one place the architecture is demonstrably mis-specified

The **shortlist** is chosen to optimise a quantity the readout does not consume in the only regime
where it matters. For an averaging readout that is harmless. For any future selection readout it is
fatal, and it is fatal *before* ranking enters — a random shortlist retains the best candidate
better than the score's shortlist at every width. **If a discriminator is ever found, the shortlist
must be rebuilt at the same time or the discriminator will have nothing to discriminate.**

---

## 4. Methodology: six coordinator errors, none of them a computation

Sprint 16 ended by adopting a rule — *before quoting a statistic as evidence for a law, derive
which functional of it the law consumes, and quote that functional.* Sprint 17 broke that rule and
five others, and **every break was caught by a workstream other than the one that made it.**

| # | error | caught by | what it cost |
|---|---|---|---|
| E1 | read a recorded pairwise ordering **accuracy** (null 0.500) onto the copula's **ρ** axis | SELECT | the strategic read "0.600 lands at ≈2.02 Å" — retracted. Accuracy 0.600 is ρ 0.309 |
| E2 | published a **min-of-N** routing ceiling with no min-of-N null | SELECT | a whole proposed research direction (§35/§36 routing), retracted within the hour |
| E3 | wrote a directional conclusion from a **5-target smoke** that reversed at n = 126 | the full run | the averaging-headroom sign: +0.019 became −0.113 |
| E4 | printed **free-superposition** member error beside **common-frame** diversity | AUDIT | broke the governing identity by ~18% of readout² inside the sprint's shared instrument |
| E5 | compared the average over **all of K = 500** (3.396) with the incumbent instead of the average over the **shipped top-75** (3.048) | AUDIT | made averaging read 0.19 Å *worse* when it is 0.156 Å *better* — inverting the sprint's central comparison |
| E6 | shipped a checkpointed artefact with **no completion flag** | AUDIT | three workstreams read it at 60/126 rows, fold 3 under-represented 2.7× |

**No number failed to reproduce.** Two independent RMSD implementations agree to 4.8 × 10⁻⁷ Å over
~2.1 M windows; two independently written instruments (`oracle_map.py` and `sel_bench.py`) agree to
**3.2 × 10⁻⁷ Å** on 120 targets at four widths; the consensus medoid reproduced Sprint 8's
−0.172 [−0.316, −0.027] as **−0.172 [−0.322, −0.031]** on a rebuilt instrument. The failures were
all in *reading*, and E1 was committed **inside the message instructing another workstream to guard
against exactly that failure**.

### Three rules this sprint adds

> **1. An ORACLE arm needs a null of its own construction.** "Best of V things chosen with labels"
> is not a ceiling until it is priced against "best of V things chosen without them". The minimum
> of N draws falls with N whether or not any structure exists — and on one instrument the *real*
> per-target minimum **lost** to that null with an interval excluding zero.

> **2. Do not write a conclusion from a smoke run at all, even a directional one.** Record the
> smoke as a smoke and wait. Three of this sprint's reversals, and one in Sprint 16 with a
> confidence interval that *excluded zero* at n = 8, were small-n reads.

> **3. State which of three axes every rank statistic is on.** Pairwise ordering accuracy (null
> 0.500), Spearman, and copula ρ are different quantities related by
> ρ = sin(π(2·acc − 1)/2) and ρ = 2 sin(π ρ_S/6). They are now named functions rather than prose.

**And an observation about the agent structure itself**, since it is the reason the errors were
caught: every one was found by a workstream with a *different* instrument and a *different*
incentive. The adversarial auditor found three; the workstream whose direction was being redirected
found two. A single-agent version of this sprint would have shipped all six.

---

## 5. The quantum contribution (§G, §H)

### 5.1 The finding that invalidates three sprints of comparisons — including this one's framing

> **Greedy 1-opt reaches the certified global optimum of the deployed objective in 100% of 152
> (target × seed) cells at 1,024 objective evaluations, and 62.5% at 36 — one coordinate-descent
> pass. The CVaR-VQE budget used throughout Sprints 15–17 is 8,192.**

**Every "matched-budget" quantum-versus-classical comparison this programme has run was conducted
at least 8× past the point where the classical control saturates. A budget the control does not
need is not a budget.** No audit caught it, because every individual number was correct — the
failure is in what "matched" meant. It invalidates the premise of the Pareto framing this sprint
briefed, and it is the sixth instance in this sprint of a correctly measured quantity read as a
different quantity.

### 5.2 Why — and this part is exact

From the Walsh/Pauli-Z spectrum on 19 targets over full registers, the deployed objective `hamil`
puts **61.3% of its variance at Pauli weight 1** and **93.0% at weight ≤ 2**, with **95.7% of the
weight-2 mass intra-residue**. It is very nearly a **separable per-residue field**. Its weight-≤1
truncation — closed form, **zero search** — has a certified argmin of **2.411 Å against the full
objective's 2.661 Å**, i.e. *better than the objective it truncates*. The ORACLE truth sits at mean
Pauli weight **3.72**, only 31.4% at weight ≤ 2.

> **There is no search problem here for a quantum device. The objective is nearly separable and the
> answer is not — and the objective's existing higher-degree content points the wrong way.**

### 5.3 Every pre-registered falsifier fired

| test | result |
|---|---|
| **P1 Pareto** | all 10 VQE configurations ε-dominated by the 21-point classical set on the readout the pipeline consumes, every CI excluding zero (−0.129 to −0.190 Å) |
| **P1, honestly stated** | the classical leave-one-out null is identical (−0.139 to −0.174). **The VQE is one more classical sampler, slightly worse — not uniquely bad** |
| **P3 local VQE** | greedy and Metropolis certify the local optimum in 100% of windows at ¼ the exhaustive cost; VQE 68.5–72.2%, proposals worse at set-best than a random draw |
| **P4 reweighting** | VQE weights are indistinguishable from **a permutation of their own weights** (3.299 vs 3.287 Å) and worse than doing nothing (3.135); entropy-matched Boltzmann reaches 2.865 |
| **P5 mixtures** | 2.07 unique elite structures per target against the best classical arm's 5.64 |
| **entanglement control** *(new — Sprint 16 never ran it)* | deleting the CNOTs gives a classical product-Bernoulli variational model at identical budget, estimator, optimiser and seed: **null on member error and on the readout at every α**, and the entangled circuit is **significantly less diverse**. Adding those three classical points collapses every remaining non-dominated cell |

### 5.4 The one live positive, and its bounds

`log q_θ` beats both the objective it was trained on (+0.120 in-band) **and** the exact mean-field
classical model (**+0.083 [+0.001, +0.179]** at α = 0.02). It is labelled **OPEN, not SUPPORTED**:
the interval clears zero by 0.001, the median (+0.054) is far below the mean, W/L is 5/3, fold sign
2/4, **n = 8**, and **it does not convert through any readout**.

**α is a temperature — with a correction.** `T_eff` falls 0.668 → 0.190 monotonically in physical
units, but `KL(q_θ‖Boltzmann)` is **4.55–5.99 bits**: α sets a temperature, and q_θ is not that
Boltzmann law.

### 5.5 The sprint's one positive lever, and it is not selection

> **Local refinement of a retrieval candidate is worth ≈0.7 Å** — 3.636 → **2.95 Å** over 6-residue
> windows, **certified by exhaustive enumeration** at 4,096 evaluations per window, with greedy
> reaching the same optimum at 1,024. Cheap, exact, native-free.

With selection closed at four levels, this is the only measured mechanism in the sprint that moves
accuracy — and it attacks **Problem A/C** (improve the candidates) rather than Problem B (choose
among them). It is also, per §5.1, a classical mechanism.

### 5.6 Honesty about power

The Pareto arm ran n = 10 of 19 targets and the representation arm n = 8; **both were stopped by
the workstream itself when free RAM hit 0.30/0.68 GB**, which the brief forbids exceeding, and were
deliberately *not* resumed while siblings held the box. The two results carrying the most
weight — the **budget floor** and the **Walsh spectrum** — both ran the full 19 targets, as did the
entanglement control. The workstream also caught and preserved **two errors of its own**: a sign
error in the primary Pareto statistic that manufactured large quantum wins (−0.48 to −0.89 Å with
intervals excluding zero), and a conditional margin undefined exactly on the cells the VQE could
have won.

---

## 6. The physics contribution (§I, §J)

### 6.1 Legacy — null as a ranker, real as a gate, and destructive in a way now measured to a single candidate

| arm | result |
|---|---|
| Legacy in-band selection vs matched random | **+0.0572 [−0.0348, +0.1547]**, 55W/69L — null, sign of harm |
| Legacy argmin vs the distance argmin, K = 500 | **+1.036 [+0.694, +1.384]**, 38W/88L — decisively worse |
| every Legacy gate vs its own matched random gate | null (worst \|mean\| 0.046) |
| **Legacy gate at f = 0.50, near-native recall** | **0.613 vs random 0.501: +0.112 [+0.023, +0.198]**, 49W/21L |
| **the same gate, set BEST** | **1.427 vs random 1.208: +0.219 [+0.031, +0.444]** — significantly worse |

> **In Legacy's own normalised order the pool's single best candidate sits at 0.450 while the
> *median* sub-2 Å candidate sits at 0.410 (random = 0.500). Legacy likes near-natives, and likes
> the best one systematically less than the typical one.**

That is Sprint 16's set-mean trap resolved to the level of an individual candidate — the sharpest
statement of that mechanism the programme has produced.

**`leg_torsion` is the single exception and the only positive Legacy result in the sprint**: the
strongest in-band component (Spearman 0.101 against the distance score's 0.126), the only score of
seventeen with a *negative* in-band selection delta (−0.0446 [−0.1378, +0.0482], 66W/57L), and the
**only gate measured anywhere in this programme that buys near-native recall (+0.159 [+0.073,
+0.239]) without significantly damaging the set best (+0.011 [−0.079, +0.128])**.

### 6.2 AMBER as constrained repair — the frontier splits in two

Pre-registered stratified 30-target subsample, reused verbatim from Sprint 16 so no new selection
freedom is created. Do-nothing on this subsample is 3.5061 Å.

| rung | Cα disp | RMSD cost vs do-nothing | ramaFav | clash < 2.0 Å | cis | ω dev |
|---|---|---|---|---|---|---|
| none | 0.000 | — | 0.781 | 7.37 | **0.000** | 0.0 |
| `cafix` (ε = 0) | **0.000** | −0.0000 *(tautology)* | 0.691 | **0.00** | 0.313 | 59.2 |
| `ca1000` | 0.185 | **+0.0081 [−0.0033, +0.0187]** | 0.647 | **0.00** | **0.539** | 89.2 |
| `ca300` | 0.317 | +0.0305 [+0.0085, +0.0519] | 0.649 | **0.00** | **0.574** | 96.4 |
| `ca100` | 0.499 | +0.0662 [+0.0331, **+0.0996**] | 0.717 | **0.00** | 0.488 | 89.1 |
| `ca30` | 0.782 | +0.1402 [+0.0843, +0.1968] | 0.814 | **0.00** | 0.201 | 47.9 |
| `ca3` | 1.111 | +0.2682 [+0.1588, +0.3833] | **0.847** | **0.00** | **0.084** | 22.5 |
| `free` | 1.878 | +0.4586 [+0.2374, +0.6876] | 0.845 | **0.00** | 0.087 | 21.6 |
| **`k30`** (N/CA/C — the incumbent) | 0.824 | +0.1217 [+0.0688, +0.1769] | 0.828 | **0.00** | **0.086** | 25.0 |

1. **Steric validity is available at zero Cα cost — but the magnitude was overstated and is
   corrected here.** *Every* rung, `cafix` included, reaches **0.00** clashes below 2.0 Å at a Cα
   cost anywhere from 0 to 1.88 Å. **However**, the convergence gate selects a *different subset at
   every rung*, and comparing a gated arm against the ungated instrument-wide input is the failure
   mode this report is about. `cafix` converges on the **78 least broken targets**, whose input
   carries **0.31** sub-2.0 Å clashes, not 7.37. Against its own gated input the effect is
   **0.31 → 0.00 (13 better / 0 worse)** and **2.03 → 0.01 (31/0)** — real, one-sided, and **an
   order of magnitude smaller** than the uncorrected figure. *The workstream caught this in its own
   draft; the coordinator had already published the uncorrected version.*
2. **Torsional and peptide-bond validity is a STRICT TRADE.** Ramachandran rises 0.647 → 0.847 and
   cis falls 0.574 → 0.084 **monotonically in the realised displacement**. No rung buys them
   cheaply, and the pre-registered target — the unrestrained arm's validity at < 0.05 Å of Cα
   cost — is **refuted**.
3. **A null against the lane's own premise, and it validates the incumbent.** The framing was that
   freeing N and C would let the force field repair in coordinates Cα-RMSD is blind to. At matched
   displacement it does the **opposite**: **`k30`, the restraint set the pipeline already ships,
   beats the Cα-only set on accuracy, Ramachandran, geometry *and* cis.** Sprint 16's finding that
   55.1% of AMBER's displacement is non-torsional is **not** an opportunity to be harvested.
4. **The tight-Cα regime is the worst place on the ladder, non-monotonically.** cis **peaks at
   0.574 at `ca300`** — worse than `cafix` and far worse than `free` — with ω deviation peaking at
   96.4°. A Cα-only restraint at high k pins the contracted spacing and lets N and C absorb it by
   flipping ω. **Holding Cα while freeing the peptide plane is the one thing not to do.**

### 6.3 Sprint 16's cis-peptide defect is explained, and it exonerates the force field

Sprint 16 filed "AMBER introduces cis peptide bonds on 44 of 126 targets" as a defect of the repair
operator. **It is not.** The coordinate average **contracts the backbone by 22.4%** — Cα–Cα spacing
**2.949 Å against an ideal 3.80 Å** — and the minimiser is handed a chain whose peptide planes
cannot be satisfied at that spacing:

    Spearman(input Cα–Cα spacing, cis fraction)        −0.949
    the same, conditioned on convergence               +0.807

> **Fix the averaging operator, not AMBER.** The force field is responding correctly to an
> unphysical input — and this also explains the non-monotone cis peak directly: a Cα-only restraint
> at high k *pins* the contracted spacing while leaving N and C free to absorb it by flipping ω.

### 6.4 The decisive AMBER-versus-Legacy comparison, at full scale

**63,000 genuine ff14SB/GBn2 single points, all 126 targets, identical K = 500 candidate sets, no
shortlist.** AMBER lost every ranking role:

| | AMBER | Legacy |
|---|---|---|
| global Spearman | **−0.027** | +0.307 |
| in-band AUROC | **0.499** | — |
| argmin | **5.191 Å** (random: 4.427) | 4.490 Å |
| incremental value over the distance features | **−0.004 [−0.037, +0.030]** | — |

`angle` and `torsion` are **significantly anti-informative**. An all-atom force field, scored
exactly, ranks peptide candidates **worse than chance and worse than random selection.**

**What Legacy knows that AMBER does not:** garbage rejection — and one specific term.
**`leg_contact`, the Miyazawa–Jennings contact term that Sprint 16 called useless, reaches
+0.080 [+0.032, +0.128] against a +0.010 null.** It carries real information and **still selects
0.136 Å worse than random** — the cleanest statement in the programme of the gap between *having
information* and *being able to rank with it*.

**What AMBER knows that Legacy does not:** where the atoms go. That is stereochemistry, not
ranking, and it is the only role either model has earned.

### 6.5 The Level-1 deliverable is not achieved, and the incumbent is validated

The repair tax cannot be eliminated. **Level 1 — recover ≈3.05 Å with valid all-atom structures —
fails**, and the incumbent's repair stage is validated rather than improved.

### 6.6 Two pieces of method discipline worth more than the numbers

**The falsifier this report's own brief handed the physics lane could not fire, and the workstream
said so *before* the run.** At ε = 0 the Cα coordinates do not move, so Cα-RMSD equals the
do-nothing RMSD **exactly, by construction** — "accuracy cost above 0.10 Å against doing nothing"
is untestable at that rung. It was recorded in the lane's pre-registration in advance and the rung
judged on validity alone. **That is the sprint's characteristic failure mode caught in the brief
itself rather than in a result.**

**And the workstream found the radians defect in its own zero-information control**:
`core.geometry.build_backbone` takes radians; its first constant-α-helix reference passed degrees,
silently building a different conformation scoring Ramachandran **0.000** instead of 1.000. Every
affected number was recomputed; nothing built on the broken control survives.

---

## 7. Learned representations: the last open feature class, closed

Every signal the sprint had tested was geometric or energetic. The one class never tried was
**learned sequence representations**, and there was a specific reason to try it: the programme's
record holds an overturned result that **ESM buys 0.288 Å over one-hot on selection**
(p = 0.005, n = 126) — measured globally, never in-band.

**It is not in-band content.** At matched architecture, varying only the per-residue block, the ESM
arm is **worse** than one-hot: +0.036 [−0.063, +0.139] in the top-25 band and
+0.116 [−0.004, +0.243], fold [+0.019, +0.174] in the top-75. A target-conditioned pseudo-likelihood
(ESM-2-650M, masked, 126 targets) does not separate from a **composition-only** control, and
residue-aligned window-embedding cosine has counted ordering accuracy **0.503** against a 0.500
null.

> **ESM buys a better filter, not a discriminator** — the same object §0.5 showed the distance
> objective to be. **No representation-derived feature has in-band skill, and the last open feature
> class is closed.**

**One real ESM-specific signal exists, on the correlation axis only.** ESM-2's contact head reaches
in-band Spearman **+0.116 [+0.047, +0.184]** at top-75, and against its own **ESM-free twin** the
increment is **+0.050 [+0.005, +0.097]**. But the twin control shows most of it is compactness
(ESM-chosen versus distogram-chosen pairs: +0.000 [−0.078, +0.080]), and **none of it reaches the
argmin in any band**.

### 7.1 The result that reframes the sprint, and it is EXACT

Fed through the design equation's band curve:

- ρ_S = +0.116 delivers **3.393 Å at B = 25** against a band mean of 3.502 — **0.109 Å against a
  2.249 Å gap**.
- **At B ≤ 25 the 2.5 Å target is unreachable at ρ = 1.** A *perfect* ranker inside the shipped
  top-25 returns **2.609 Å**.

> **The shortlist, not the ranker inside it, is the binding constraint.**

This is §0.3 (widening removes the answer from the shortlist) and the coverage result (a ceiling no
readout can use) **arriving independently from the feature side. Three directions now agree**, and
the third one is exact rather than empirical.

### 7.2 And the shortlist channel is not fixable from here either

The ESM channel has **no shortlist-construction value**: its shortlist is **significantly worse than
matched random** (+0.315 to +0.449). It also **independently reproduced** the physics result that
the shipped score's shortlist has an ORACLE best worse than matched random — **+0.209 / +0.226 /
+0.170** at B = 25/75/150 against **+0.195 / +0.242 / +0.179** measured on a separate instrument.
Two workstreams, two instruments, the same numbers.

A mixed-channel shortlist arm **beat the incumbent's ceiling with both intervals excluding zero** —
and its **pre-registered matched-random control killed it**: random does as well or better at every
size. That arm is exactly what this programme would have shipped two sprints ago.

### 7.3 Two things flagged and not claimed

- **`cf_topd_uniform`**, an *ESM-free* control — the mean realised distance over the *n* pairs the
  distogram predicts closest — is the battery's best arm (top-75 **−0.135 [−0.252, −0.022]** against
  random, beating both constant references and the shipped score). It is **1 of ~75 uncorrected
  tests** and **flips sign on the ORACLE band**. **PLAUSIBLE**, and it belongs to the selection lane.
- The ESM contact signal is **length-gated** (ρ +0.188 at length 14–16, ρ −0.047 at 9–11).
  **PLAUSIBLE, post-hoc.**

---

*Section 6 — the physics contribution — and section 10, the deliverables table, are completed once
the last workstream reports.*

---

## 8. Benchmark readiness (§N)

**The 60-target benchmark stays sealed.** The audit verified independently that no `s17/` module
reads it, probes it, or derives its identities, and that no bare `hash()` appears anywhere in the
sprint's code.

The five unlock conditions, assessed:

| condition | status |
|---|---|
| architecture frozen | **FAILS.** Nothing in this sprint earned a place in the architecture. The proposed change — averaging → selection — is a measured regression |
| controls frozen | not reached |
| evaluation script frozen | not reached |
| success criterion pre-registered | not reached |
| no cheap internal experiment could still change the architecture | **FAILS.** Two are outstanding and both are cheap: the composition-drift control on every K-effect, and resuming the two under-powered quantum arms |

**The first condition fails outright, so the question does not arise.** The honest position is the
same as Sprint 16's: the benchmark is not being withheld out of caution, it is being withheld
because spending it would answer no question the internal instrument cannot answer more cheaply.

---

## 9. Limitations

- **Two instruments, different powers.** Selection, readout and ceiling results are at **n = 126**,
  where the minimum detectable effect at 80% power is **0.084 Å** — any null below that magnitude
  is uninformative rather than negative. The quantum results are on **19 enumerated targets**, and
  two of its arms ran at n = 10 and n = 8 after being stopped on memory.
- **Every K-effect in this sprint is confounded with a composition change that nobody controlled.**
  Widening K takes the pool from **27.3% to 20.8% peptide**, and the programme's own record says the
  peptide corpus carries **≈7× the sequence–structure channel** of the protein fragments. This is
  the audit's highest-value unrun experiment and it bears on findings 2, 3 and 7 of §0.
- **The coverage sweep was stopped at 20/126** on a saturated box. Its strategic question was
  answered independently at n = 126 (a matched *random* shortlist retains the best better than the
  score's at every width), but "which coverage construction is best" is unmeasured.
- **The `m*` and step ladders are trained on labels leave-fold-out.** Fold-honest, but trained
  quantities, and never described as native-free.
- **The pre-repair and post-repair numbers are not interchangeable.** 3.048 Å is a Cα coordinate
  average before repair; the deployed 3.204 Å includes it. No pipeline claim is made from the
  former without the 0.155 Å tax attached.
- **One band-definition discrepancy is live and is reported rather than resolved.** A constant
  α-helix reference scores +0.053 in the shipped top-25 band and +0.2913 [+0.1966, +0.3899] on the
  `pool_best + 1.5 Å` band. Those are two different bands; no number from one may be quoted beside
  a number from the other.

---

## 10. Deliverables (§62)

**A. Best architecture.** No change is recommended. The incumbent —
`BLOSUM62 retrieval (K = 500) → distogram top-75 → coordinate average → ideal-geometry projection →
restrained AMBER (N/CA/C, k = 30)` — survived every attempt to improve it, and two of its stages
were **independently validated** rather than merely left alone: the coordinate average beats both
alternative readouts in 39 of 39 cells, and the shipped restraint set beats all nine alternatives on
essentially every axis at matched displacement.

**B. Best measured RMSD on the 126-target instrument.** **3.204 Å**, unchanged. No arm in this
sprint beat it. The best *pre-repair* readout is the coordinate average at **3.048 Å**, which is
the incumbent minus its 0.155 Å repair tax, not an improvement over it.

**C. Best conservative estimate.** There is no new estimate to make. The one arm that beat the
incumbent's ceiling with both intervals excluding zero — a mixed-channel shortlist — was **killed by
its own pre-registered matched-random control**.

**D. Improvement over 3.204 Å and over 3.050 Å.** **None, on either baseline.**

**E. Oracle ceilings.** Best member: **2.104 Å** (top-75) · **1.874 Å** (K = 150) · **1.782 Å**
(K = 300) · **1.711 Å** (K = 500) · **1.313 Å** (full universe). Best *subset average* at K = 500:
**1.598 Å** at m ≈ 6.

**F. Selection gap (oracle − realized).** **1.317 Å** at top-75 · **1.743 Å** at K = 500 ·
**2.249 Å** at the full universe. It grows monotonically with K.

**G. CVaR-VQE contribution.** **None, on any axis, and the comparison that would have shown one was
mis-specified.** Greedy 1-opt certifies the objective's global optimum in **100% of cells at 1,024
evaluations**; the VQE budget used for three sprints is **8,192**. All ten VQE configurations are
ε-dominated by the classical set — *and so is a classical leave-one-out null*, so the honest
statement is that **the VQE is one more classical sampler, slightly worse**. Deleting its
entanglement changes nothing. One quantity, `log q_θ`, beats the exact mean-field model in-band by
+0.083 [+0.001, +0.179] at n = 8 and is labelled **OPEN**; it converts through no readout.

**H. Classical controls.** Greedy 1-opt: certifies the optimum at 1,024 evaluations (62.5% at 36).
Simulated annealing: never significantly beaten by the VQE, and equals it at ¼ budget. Classical
thermostat: reproduces α's entire ensemble effect and moves the readout 1.6–1.7× further. Uniform:
beats every CVaR arm at matched budget on the ensemble readout. **The classical set dominates on
every axis measured.**

**I. Legacy contribution.** *Selection*: **null** (+0.0572 [−0.0348, +0.1547]); decisively worse
than the distance objective (+1.036). *Gating*: **null against its own matched random gate**, and
actively harmful to the set best (+0.219 [+0.031, +0.444]) even while raising near-native recall
(+0.112 [+0.023, +0.198]). *Validity*: its only intervals excluding zero anywhere. *Incremental
feature value*: none — disagreement features are entirely subsumed by the objective's own
confidence. **The one exception is `leg_torsion`**, the only gate in the programme that buys
recall without destroying the best.

**J. AMBER contribution.** *Ranking*: not established. *Repair*: **steric validity is free at any
restraint strength** — every rung reaches zero clashes below 2.0 Å. *Validity*: torsional and
peptide-bond plausibility are a **strict monotone trade** against Cα displacement. *Cα accuracy
tax*: **0.155 Å, and it cannot be eliminated** — the incumbent's own restraint set is the best of
ten tested.

**K. Hard targets.** The "53 targets with pool best > 2.0 Å" fall to **24** over the full universe;
above 2.5 Å, 27 → 8; **more than half were a truncation artefact**. The residual set is longer, less
helical, more polar, with weaker BLOSUM matches — a *retrieval* deficit carrying the same selection
gap as everywhere else (+2.437 vs +2.204). It is worth only **+0.078 Å** of the mean ceiling.

**L. Mechanism — what actually creates the final accuracy.** Retrieval and averaging, and nothing
else. The distance objective contributes **as a filter** (it raises the shortlist mean, worth
≈0.23 Å) and contributes **nothing as a ranker** (in-band argmin −0.014 [−0.118, +0.094] against
random inside its own top-25). The coordinate average then converts a good shortlist *mean* into
the emitted structure. Filter and readout are matched; the pipeline is coherent.

**M. Remaining ceiling.** **1.598 Å** at K = 500 (best subset average) or **1.313 Å** over the full
universe (best member). The gap to it is **entirely candidate choice** — ranking for a selection
readout, subset choice for an averaging one — and the two are of comparable difficulty and
comparable reward (0.113 Å apart).

**N. Benchmark readiness.** **Not ready.** The first unlock condition fails outright: the
architecture is not frozen because nothing earned a place in it. Two cheap internal experiments
remain outstanding.

---

## 11. Where are the missing ångströms (§64)

| gap | size | verdict |
|---|---|---|
| **retrieval** | ~1.3 Å of ceiling remains unreached; 24/126 targets lack a sub-2 Å candidate | **nearly exhausted** — the apparatus buys only +1.5 targets of 2.0 Å recall over a random 500 of the same library, and the marginal return on deeper retrieval is *below* the order-statistic null |
| **filter / shortlist** | **0.791 Å** of ceiling destroyed by the top-75 cut | **THE ONE TO ATTACK.** Large, and demonstrably mis-specified rather than merely hard |
| **selection (ranking)** | **2.249 Å** | **closed at five levels** — 58 objective functionals, 29 native-free signals, 27 calibration features, all disagreement features, and six classes of learned representation |
| **ensemble** | 0.113 Å of averaging headroom over the best member | **SUPPORTED but small**; interval touches zero |
| **repair / validity tax** | **0.155 Å** | **closed** — sterics are free, torsion is a strict trade, the incumbent restraint set is optimal among ten |
| **quantum** | 0 | **closed**, and the comparison that would have measured it was ≥8× past classical saturation |
| **physical models** | 0 as rankers or gates | **closed**, except `leg_torsion` |
| **target calibration** | +0.294 Å of ORACLE headroom | **closed** — the calibrator loses to a constant baseline |

### The one gap worth attacking next, and why

**The shortlist.** Three workstreams reached it independently and the third arrival is exact:

- widening K makes the realized answer **worse** because the extra windows displace near-native
  members out of the shortlist (+0.468 Å of shortlist-ceiling destruction against +0.141 Å of
  realized degradation);
- at fixed K, the score's shortlist is **no better than a matched random shortlist at containing a
  sub-2 Å candidate and significantly worse at containing the best one** — reproduced on two
  independent instruments to within 0.02 Å;
- and **at B ≤ 25 the 2.5 Å target is unreachable at ρ = 1**. A *perfect* ranker inside the shipped
  top-25 returns **2.609 Å**.

> **That last line is the sprint's central result. The programme has been trying to build a better
> ranker for a shortlist that does not contain the answer. No ranker, however good, can fix that —
> and the shortlist is chosen to optimise a quantity (the set mean) that only the averaging readout
> consumes.**

For the averaging readout the pipeline currently ships, this mismatch is harmless. For any selection
architecture it is fatal, and it is fatal *before* ranking enters.

### What that implies for the next sprint

1. **Build the shortlist for the readout you intend to use.** If the readout consumes the set best,
   the shortlist must be built to retain it — and a matched random shortlist already does that
   better than the score. The measured coverage constructions raise the shortlist ceiling
   (a 6-target smoke had k-means-plus-score at −0.275 [−0.473, −0.091] against the score's own
   top-B); the full comparison was stopped for compute and is the cheapest outstanding experiment.
2. **Control the composition confound before trusting any K-effect.** Widening K takes the pool from
   27.3% to 20.8% peptide, and the peptide corpus carries ≈7× the sequence–structure channel. Every
   K-result in this report is confounded with that, and nobody has controlled it.
3. **Stop looking for a fine ranker from these feature classes.** Five independent closures is
   enough. If ranking is attempted again it needs a genuinely new observable, and the design
   equation says what it must reach: ρ ≈ 0.60 for 2.5 Å over the full universe, against a measured
   best of 0.116.

4. **Use the one target-level signal that works, and it is the trivial one.** Buried in the
   refutation of disagreement features is the strongest target-level predictor in the sprint:
   **`dist_score` — how good the objective thinks its own pick is — predicts the selector's realized
   error at leave-fold-out ρ = 0.598**, and it subsumes every engineered disagreement feature. It
   cannot *select* (it is constant across candidates), but **target difficulty estimation is exactly
   what adaptive generation needs**, and the calibration lane's 27 engineered features could not
   beat a constant baseline while this does it with no new machinery.
4. **The one positive lever measured anywhere in the sprint is local refinement** — 3.636 → 2.95 Å
   over 6-residue windows on the enumerated instrument, certified by enumeration, reached by greedy
   at 1,024 evaluations. It improves *candidates* rather than choosing among them, and it is
   classical.

