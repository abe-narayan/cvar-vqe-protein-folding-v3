# PREREG S31-C — encoding / readout: native-free sparse support (§7C), candidate-index
# information allocation (§12), and the 128→512 widening's undischarged circularity gate

**Lane C. Registered before any number in this lane exists.** Nothing below was computed; every
figure quoted is a citation of an S30 artefact, named with its file.

---

## 0. Basis, labels and discipline that bind every arm below

- **Endpoint:** mean built-chain Cα RMSD on `tuning126`, **3.2105 Å**. The CA point cloud
  (**3.0483 Å**) and the *set mean* (3.5507 Å) are different objects. Every number states its basis
  in its own sentence.
- **Cloud→chain transfer for a correction is 0.92** (S31 contract rule 2). The 1.16 of
  `operator-consumes-set-mean` maps *set mean* to output and is not used here. Anything that clears
  a falsifier on the cloud is **re-measured on the chain**, not converted.
- **MDE = 2.8016 × SE, per comparison.** < 0.7× is not a result; 0.7–1.0× is NOT MEASURED.
  Fold-clustered CIs on the pinned folds.
- **ORACLE labels travel inside the sentence carrying the number.** C3 is ORACLE in every arm and
  will never be written without that word attached.
- **No deployable support, weight, temperature, cardinality or index map may be chosen on target
  native RMSD.** Grid choices are reported with a split-half transfer arm (contract rule 11) and
  only the transferred value counts as the result.
- **A stratum defined by the outcome cannot measure the thing that defines it** (contract rule 12).
  C3 exists because of exactly this rule.

## 1. What I accept as binding, and what I attack

**Accepted as binding, and not re-litigated:**

- **T1 (S30):** for a diagonal `H` the CVaR tail is always a *prefix* of the objective's order. The
  optimised state's entire structural output is one integer.
- **S30-L11 (`s30/results/s30_Q_sparse.json`):** the sparse weighted readout is closed **by price** —
  a plain argmin over the top-2^B beats the best fully-priced sparse arm at every budget B = 3…9,
  by +0.162 to +0.727 Å. On the built chain, 2-of-75 with free continuous weights is **2.1683 Å at
  11.4+ bits** against the top-128 argmin's **2.1435 Å at 7.0 bits**. *I do not repeat that
  experiment.*
- **S30-L11 F3:** a native-free support rule (score-prefix, farthest-point diversity) does beat a
  random support with ORACLE weights by +0.22 to +0.70 Å — **and every native-free support arm is
  still worse than production** (3.198 against 3.048 at s = 2 over the pool, point cloud), *even
  with ORACLE continuous weights handed to it*.
- **The averaging-terminal algebra** (coordinator's framing, S31 brief §7C): for an averaging
  terminal `set_mean² ≈ B² + S²`, where `B` is the endpoint and `S` the spread. Improving a set
  purely by **concentration** is worth **zero by algebra**; only a change in the bias `B` pays. And
  `B` is the pool's common mode, which S30 lane L proved **non-identifiable from pool data at any
  K** (the likelihood depends on `(t, μ)` only through `t + μ`).
- **S30-L16 (lane F):** widening the **filter** does not help filter-independent tails — on the
  *averaging* operator. That is a different operator from C3's, which is the **argmin**; lane F's
  own entry says so. C3 is not a repeat of it.
- `pool-error-is-68-percent-common-mode`, `consensus-is-outlier-avoidance`,
  `in-band-ordering-is-per-target`.

**What I attack:** that the −1.9004 Å widening figure means anything about hard targets (C3); that
the deployed 7-bit score-order index is the best use of 7 bits (C2); and the charter's own claim
that §7C is still open after S30 closed the 2-of-75 form (C1).

---

# C3 — THE 128→512 WIDENING, AND ITS CIRCULARITY GATE

**Run first. It is the cheapest and it either revives or closes a whole direction.**

## C3.1 The quantity

Per target, over the shipped `K = 500` BLOSUM pool in the shipped DIS score order:

```
best1(N)  =  min over the first N of the score order of the ORACLE Cα-cloud RMSD to the native
g         =  best1(128) − best1(500)        (≤ 0, more negative = widening helps more)
```

**Every number in C3 is ORACLE and NOT DEPLOYABLE.** It says the candidate *is there*, not that
anything native-free can find it; S29 closed recognition three ways and S30 confirmed it.

S30-L11's published figures, which this test re-examines: `g` = **−0.4350 Å** over all 126,
**−1.9004 Å on FAIL18**, **−0.1907 Å on the other 108**; difference of stratum means −1.7097,
SE 0.2207, **2.77× MDE**.

## C3.2 The circularity, stated precisely, because it is nearly a tautology

`FAIL18` (`s12/instrument.py:271-278`) is *defined* as the targets whose score-top-75 retained
**zero** members within 1.5 Å of the pool optimum. That definition **directly asserts that the
near-optimal band lies outside the top-75.** The quantity C3 measures — whether the best member
lies outside the top-128 — is the same statement shifted by 53 ranks. A stratum defined by "the
good candidate is below rank 75" cannot be used as evidence for "the good candidate is below rank
128."

## C3.3 Strata, including two whose bias direction I state in advance

| stratum | definition | what it is contaminated by |
|---|---|---|
| `FAIL18` | the filter's own recall, `s12/instrument.py:271-278` | **the quantity under test**, as §C3.2 |
| `worst18_poolmean` | 18 largest mean ORACLE RMSD over all 500 members | **the clean one** — a property of the pool, no argument of `g` and no reference to the score order. **This is the headline stratum.** |
| `worst18_bestpool` | 18 largest `best1(500)` | biased **against** finding an effect: `best1(500)` is the subtrahend of `g`, so selecting for it large mechanically shrinks `|g|`. Reported, not headlined |
| `defn18` | 18 largest `best1(75) − best1(500)` — the quantity FAIL18's rule thresholds | a **definition-matched** stratum: it prices how much of FAIL18's `g` is forced by its own construction |
| `random18` | 20,000 draws of 18 of the 126 | the null. Report the full distribution of the stratum mean `g`, its 2.5th percentile, and the p-value of each stratum against it |

## C3.4 Pre-registered falsifiers (Å, ORACLE point cloud)

- **F-C3a — REVIVE.** Fires iff **`worst18_poolmean` has mean `g` ≤ −1.00 Å** *and* it falls below
  the `random18` null's 2.5th percentile.
- **F-C3b — CLOSE.** Fires iff `worst18_poolmean` has mean `g` ≥ −0.50 Å, **or** it does not clear
  the `random18` null.
- **PARTIAL.** −1.00 < mean `g` < −0.50 on `worst18_poolmean` *and* it clears the null. Verdict:
  the direction survives **at reduced size**, and the number that may be quoted afterwards is the
  filter-independent one. **−1.9004 Å is never quotable again in either case.**

**Diagnostic, not a falsifier, registered so it cannot be invented afterwards:** `defn18` mean `g`
against `FAIL18` mean `g`. If `|defn18| ≥ |FAIL18|`, the published FAIL18 figure is substantially
definitional and I will say so in those words.

**Second diagnostic:** the rank of the ORACLE-best pool member per stratum against the
score-uninformative null, under which `P(rank > 128) = 372/500 = 0.744`. This separates "the score
has no skill on this stratum" from "hard targets hide their answer deeper", which are different
mechanisms with different costs.

**Level control (contract rule 12's positive form):** `g` is mechanically bounded by `best1(128)`.
Across all 126 I regress per-target `g` on `best1(128)` and report whether stratum membership adds
anything beyond the level, and the normalised `g / best1(128)` per stratum.

## C3.5 Honest prior odds

**3:1 against F-C3a** (full replication at ≤ −1.00 Å). **~2:1 in favour of PARTIAL** — a real but
roughly-halved effect. Reasoning: the shipped score's within-pool ordering degrades on every tail
definition (S30 §4.1: ρ = +0.645 on the easy 108 against +0.107 FAIL18, +0.380
`worst18_poolmean`, +0.344 `worst18_bestpool`, all below the random-18 null's lower bound
+0.4067), so *some* deep-hiding is real and not filter-defined. But FAIL18's ρ is 3.2–3.6× lower
than the other two, and that gap is the definitional part.

---

# C2 — CANDIDATE-INDEX INFORMATION ALLOCATION

## C2.1 Why re-indexing can matter at all, given T1

Under a diagonal `H` the tail is a prefix and re-indexing changes nothing about *which* candidates
the tail holds. What re-indexing changes is **which distributions the ansatz can reach and what a
partial measurement localises.** The deployed circuit is a 3-layer RY/CNOT ansatz over 7 qubits
(`core/pipeline.py:806-852`); a hardware-efficient state's measured distribution is close to a
**product distribution over the index bits**, and a product distribution's mass is organised by
*subcubes*. **The index map decides which physical candidate sets are subcubes.** Under the
deployed score-order map a subcube is a score-contiguous block; under a structure-aware map it is
a structural cluster.

So the goal, in the charter's words, is not to represent 128 candidates but to make each qubit a
scientifically meaningful candidate distinction. That is measurable and is measured here.

## C2.2 Index maps compared (all native-free constructions; the readout scoring is ORACLE)

All over the deployed top-128 of the score order, so the *candidate set* is held fixed and only the
**assignment of basis states to candidates** varies. This is the matched-space control (rule 7).

1. `score` — **deployed**: rank r → binary(r).
2. `gray` — rank r → Gray code of r. *Registered expectation, which is analytic and not a result:
   a j-bit prefix cell under Gray coding is the same contiguous block as under binary, up to
   reflection — so Gray coding can change only the Hamming geometry, never the prefix cells.*
3. `bisect` — recursive 2-medoid bisection of the 128 under pairwise Cα RMSD; bit j = the branch
   taken at level j. Cells are structural clusters.
4. `bisect_score` — `bisect`, with the branch labelled 0 at each split being the one with the better
   mean DIS score. Coarse-to-fine / hierarchical indexing.
5. `spectral` — order by the Fiedler vector of the RMSD-similarity graph, then binary index
   (locality-preserving in one dimension).
6. `perm` — a seeded random permutation of the same 128. **The matched zero-information control**
   (rule 9: a permutation of the same candidates, not a different measure).

## C2.3 What is measured

- **(M0) Pre-check, run before the contrast (rule 22).** Spearman between index Hamming distance
  and structural Cα-RMSD distance, over all `C(128,2)` pairs, per map. *This can only remove my own
  excuse:* if the structure-aware maps do **not** raise it over `score`, the family is dead and no
  ceiling needs computing.
- **(M1) Product-state reachable ceiling, ORACLE.** `p(q) = ⊗_b Bernoulli(q_b)`, 7 parameters; the
  readout is the `p`-weighted coordinate average in the deployed common frame
  (`s27/s28_A_amp.py:Frame`, the frame `average_weighted` uses). `q` is optimised **against the
  native** — ORACLE, NOT DEPLOYABLE, an upper bound on what any product-state ansatz over that
  index could emit. Compared across maps.
- **(M2) Partial-measurement curve.** For j = 1…6, measuring the first j bits selects a cell of
  `2^(7−j)` candidates. Report (a) the ORACLE-best cell's average (ceiling, NOT DEPLOYABLE) and
  (b) the **native-free** cell chosen by each of two rules — best mean DIS score, and largest
  consensus mass (the cell whose members have the smallest mean pairwise RMSD to the whole 128).
- **(M3) Bit value.** Mutual information between each bit and the indicator "this candidate is in
  the ORACLE-best decile of the 128", per map; and the implied Å value of a bit against the
  `value-of-a-bit` law `D(R) = a + c·2^(−R/γ)` (a = 1.3312, γ = 3.1636, S30 §4.2). **That law was
  fitted on an ORACLE ladder and its currencies are not exchangeable — it is used as a direction,
  never as a constant, and never to convert a measured Å into a claimed bit price.**

## C2.4 Pre-registered falsifiers (point cloud)

- **F-C2a — the product-state ceiling.** Some native-free index map raises the **ORACLE**
  product-state readout ceiling over the `score` map by **≥ 0.10 Å**, with a fold-clustered CI
  excluding zero. *Refuted otherwise.* **This is a ceiling; clearing it deploys nothing.**
- **F-C2b — the deployable one.** Some (native-free map, native-free cell rule, j) triple beats the
  deployed uniform top-75 average (**3.0483 Å, point cloud**) by **≥ 0.10 Å**, ≥ 1.0× MDE, CI
  excluding zero, **and the gain survives a split-half transfer over the (map, rule, j) grid**
  (rule 11). *Refuted otherwise.*
- **F-C2c — the pre-check.** `bisect` and `spectral` raise Spearman(Hamming, structural distance)
  over `score` by **≥ 0.15**. *If refuted, C2's structural family is dead before M1 and M2 run and
  I report that instead of a ceiling table.* **This is near-tautological by construction and is
  registered as a pre-check, not as a result** — passing it grants nothing.

## C2.5 Honest prior odds

- **F-C2a: 3:1 against.** A product state's reachable weight vectors are restricted either way, and
  the RMSD of a weighted average is dominated by the common mode, which no permutation touches.
- **F-C2b: 6:1 against.** The averaging-terminal algebra: every cell-selection rule I can build is
  a *concentration* operator, and concentration is worth zero. The one thing that could beat it is
  that a cell's centroid genuinely moves, and S12's consensus medoid (−0.172 Å) is the whole known
  size of that effect.
- **F-C2c: 2:1 in favour**, and it is the weakest claim on the page.

---

# C1 — NATIVE-FREE SUPPORT AND WEIGHTS (§7C)

## C1.1 The constraint I am not going to walk into

The deployed terminal is an **average**. For an averaging terminal `set_mean² ≈ B² + S²`, with `B`
the endpoint. **A sparse readout that merely concentrates the set is closed before it runs** — the
gain is algebraically zero. What is not closed is a readout that changes `B`. So every arm below is
reported with its `(B, S)` decomposition, and an arm that improves only `S` is reported as a null
*whatever its RMSD does*.

## C1.2 The arm families, chosen because each one has a stated mechanism for moving `B`

| family | construction (all weight rules native-free) | mechanism for `B` |
|---|---|---|
| **W1 soft score** | `w ∝ exp(−z_DIS/T)` over top-k, T on a grid | none — pure re-concentration. Carried as the **registered concentration control**: its `B` should not move |
| **W2 typicality** | `w ∝ exp(−mean-pairwise-RMSD/T)` | outlier avoidance; `B` moves only if outliers are biased |
| **W3 cluster-select** | 2-medoid / k-medoid the top-k, emit one cluster's average; cluster chosen by score-mean or by mass | the centroid genuinely moves between clusters |
| **W4 affine extrapolation** | `C = C_top + η·(C_top − C_rest)`, `C_rest` the average of ranks k…500, η on a grid. **Non-convex: this is the only family whose emitted point can leave the candidate hull** | **directly attacks `B`**: if the score orders the pool at all, the top-minus-rest direction is an estimate of the native direction. The one family here that is not concentration |
| **W5 adaptive cardinality** | per-target `m` from a native-free statistic (pool spread, score gap), leave-fold-out | `B(m)` curve; S29's m-sweep is flat, so registered as expected-null |

## C1.3 The identifiability arm, which is what I expect to carry the finding

The ORACLE convex-hull weights over the top-75 are already on disk
(`s29/results/s29_O_cloud_rows.jsonl`, arm `hull_top75`). I fit those ORACLE weight vectors from
**native-free per-candidate features** (z-scored DIS, rank, mean pairwise RMSD to the set, distance
to the set medoid, per-candidate Rg, length) **out of fold**, and report the out-of-fold R² of the
weights and the Å emitted by applying the fitted weights. *This is the S30 lane-P form and it gives
a closure by identifiability rather than by an exhausted search* (rule 30: I will keep theorem,
empirical estimate and ORACLE quantity in separate categories and will not claim an
information-theoretic impossibility).

## C1.4 Pre-registered falsifiers (point cloud first; any pass is re-run on the chain)

- **F-C1a — the deployable one.** Some arm in W1–W5 beats the deployed uniform top-75 average
  (**3.0483 Å, point cloud**) by **≥ 0.10 Å**, ≥ 1.0× MDE, fold-clustered CI excluding zero, **and
  the grid choice transfers split-half**. *Refuted otherwise.* On a pass, the arm is re-measured on
  the **built chain** on all 126 and the chain number is the result.
- **F-C1b — the mechanism.** Every arm that lowers the emitted RMSD does so by moving `S` and not
  `B`. *Refuted if any arm moves `B` toward the native by ≥ 0.10 Å.* **This is the falsifier I most
  want to fire**, because W4 is built to fire it and a `B`-mover would reopen the whole readout.
- **F-C1c — identifiability.** The out-of-fold R² of the ORACLE hull weights on native-free
  features clears **5%**. *Refuted otherwise*, which closes native-free weighting by the same route
  S30 closed prior correction.

## C1.5 Honest prior odds

- **F-C1a: 6:1 against.** S30-L11 F3 already showed every native-free support arm sits *above*
  production even when handed ORACLE continuous weights; the weight channel is the cheap one and it
  was not the problem.
- **F-C1b: 4:1 against**, i.e. I give W4 roughly a 20% chance. It is the only arm with a real
  mechanism and the only one I would be surprised to see fail for an uninteresting reason.
- **F-C1c: 5:1 against** clearing 5%.

---

## 2. Order of work, and what makes me stop

1. **C3** — analysis only, from a cached pool table. Result messaged to `main` immediately,
   whichever way it goes.
2. **C2 M0 pre-check.** If F-C2c is refuted, C2's structural family is reported dead and M1/M2 are
   not run.
3. **C2 M1–M3**, then **C1**.
4. **Chain confirmation** only for an arm that clears its own falsifier on the cloud. No arm is
   converted from cloud to chain by a coefficient.

**What makes me stop rather than continue:** if F-C1a and F-C2b are both refuted, I do not go
looking for a sixth weight family. The finding is that the readout's native-free branch is closed
by the averaging algebra plus non-identifiability, and that is reported as the answer.

## 3. Multiplicity, counted here rather than asserted later

Registered comparisons: C3 — 5 strata × 1 contrast + 1 null + 2 diagnostics + 1 regression = **9**.
C2 — 6 maps × (1 pre-check + 1 ceiling + 6 j-levels × 3 cell rules + 7 bit-MIs) = **156**, of which
**2** are read as results (F-C2a, F-C2b) and the rest are descriptive. C1 — 5 families × grids
(T: 6, η: 13, k: 5, m: 5) ≈ **120** cells + 1 identifiability fit = **121**, of which **3** are read
as results. Every grid maximum is reported with its split-half transfer and never alone.

---

*Registered by lane C, 2026-09-20. No number produced by this lane exists in this document; none
existed when it was committed.*

---

# AMENDMENT 1 — 2026-09-21, BEFORE THE FIRST LANE-C NUMBER

**Nothing above is deleted.** This amendment is appended because a record search completed after
the registration and because the coordinator posted S31-L1 (R1) and a follow-up. Every number cited
here is an existing artefact, not a lane-C measurement. **No lane-C number existed when this was
committed either.**

## A1.1 Most of C1's registered weight families are ALREADY CLOSED in the record

The families W1, W2, W3, W5 and the one-dimensional form of W4 are closed, several **by ceiling**
rather than by a failed fit. Running them again would be decoration. Sources, all CA point cloud
unless marked:

| registered family | status in the record | source |
|---|---|---|
| **W1 soft score** | **CLOSED, twice, one of them by ceiling.** 30 deployable operator arms, best `wscore_T2.0_75` **3.0445** (-0.0038, CI straddles zero, reverses on drop-top-10); the leakage-ORACLE grid over rank-power and distance-to-medoid **selects the uniform point as optimal** — "no headroom to reach, not merely an unreached one"; score-softmax **0.12x its own MDE** in sample and evaporates under nested CV; Boltzmann(T=0.3) minus the deployed state **-0.0002, 0.00x MDE** | `s12/agg_FINDINGS.md:43-70`; `s23/LEDGER.md:99-116`; `s25/LEDGER.md:282-290` |
| **W2 typicality** | **CLOSED and measurably HARMFUL.** `wcons_T1.0_75` **+0.0364 [+0.002,+0.069], 47W/79L — CI excluding zero on the bad side**; `medoid75` +0.2339; the real outlier-strip rule is **+0.142 A worse with a matched random-removal control at a dead null** — "the outliers are not noise, they are load-bearing" | `s12/agg_FINDINGS.md:48-55`; `s23/LEDGER.md:117-136` |
| **W3 cluster-select** | **CLOSED.** 192 native-free mode arms (3 clusterings x k = 2...5 x 13 rules); best `ward5_modeavg` **+0.0131**. And "every achievable cluster-choice rule (size / score / combined / random) at every (m, k) is at or worse than the incumbent", against an ORACLE ceiling of -0.26 to -0.48 A that is **unreachable by any declared rule** | `s12/agg_FINDINGS.md:129-158`; `s23/LEDGER.md:103-106` |
| **W4, 1-D form** | **CLOSED on both registered clauses, both blind definitions, both bases.** The typicality axis: ORACLE global step is **t = -0.0**, `lfo_LIB75` **+0.0000 cloud and +0.0000 chain**, bit-identical on 126/126. PC1 of the top-75 deviations: ORACLE global **eta = +0.0 exactly**, LFO eta **+0.0071 worse**, ORACLE-best eta positive on **52%** of targets — the sign is a measured coin flip | `s29/LEDGER.md:1718-1837`, `1883-1936` |
| **W4, unconstrained form** | **CLOSED catastrophically.** lambda = 0 circuit read as signed weights **6.1366, +2.930, 4.36x MDE, 5/5 folds**; untrained circuit **~ +19.1**; unconstrained signed weights 3.4880 (75) / 3.5755 (500); random 27-dim affine subspace 3.5445. **And the convex control pins it: the simplex optimum under the same objective converges to 3.0522, i.e. to production** | `s27/LEDGER.md:1059,1070-1076,1232`; `s27/s28_A_FINDINGS.md:21,42,50` |
| **W5 adaptive cardinality** | **CLOSED.** chain slope **-0.00047 A per unit m**; ORACLE global m = 72 worth -0.0018; LFO m **+0.0079 worse**; best cell priced at 112% across-target null with **-5% split-half transfer** | `s29/LEDGER.md:5208-5265`, `4288-4290` |
| **the rescale / de-contraction idea** | **CLOSED on the ENDPOINT basis.** All four native-free scale arms are worse with every CI excluding zero (distogram +0.123, pool +0.095, both +0.080), ORACLE true scale -0.275; correlation of every native-free scale estimator with the true scale **+0.027 / -0.001 / +0.106 / +0.070, i.e. zero**. And the global scale scalar is unreachable **in principle**: `s* = 1 - <c,ebar>/||c||^2` is a function of the common-mode error and of nothing else | `s15/coord_FINDINGS.md:930-996`; `s23/LEDGER.md:143-176`, `291-321` |

**Consequence:** C1 as registered would have re-run seven closed things. **C1 is narrowed** to A1.3
below, and the closures above are reported as the answer to §7C's "can a native-free mechanism
identify a small set of members and weights" — *it has been asked seven ways and the answer is on
disk*.

## A1.2 A withdrawn number my registration leaned on

**"Coordinate averaging contracts the backbone 25.8%" is WITHDRAWN** — `s15/coord_FINDINGS.md:914-921`
records the withdrawal in its author's own words ("magnitude was overstated by a factor of seven").
The defensible figures are **3.5% against true distances**, and a *monotone separation profile*
(0.773 at the virtual bond, crossing 1.00 near |i-j| = 8, reaching **1.10 at 13**) — the cloud is
short locally and **long** at long range. `s30/LEDGER.md:129` and project memory
`averaging-space-beats-the-objective` both still carry the stale 25.8%. **No lane-C arm is built on
a de-contraction premise.**

## A1.3 What C1 becomes: the constrained-affine ladder and the identifiability arm

The coordinator's follow-up asks the right question in the right place — *not* "can sparsity buy
capacity" but **"what is the right constraint set on an affine readout that is otherwise too
expressive to be safe"** — and the record shows both endpoints of that ladder are already measured
and both are bad: unconstrained affine is +2.9 to +19 A worse, and the convex optimum under the
deployed objective **is** production. **The open rung is the bounded middle.** Registered now:

- **C1-L (the ladder), ORACLE ceilings first, candidate set HELD FIXED at the top-128** — because
  S30-L11's published pair did not hold it fixed and that is the one methodological defect in it:
  `argmin` over 128 . uniform-s-of-128 for s in {2,3,5,10,20} . convex hull of 128 (**this row does
  not exist on disk and is the missing rung**) . affine with ||w||_1 <= B for B on a grid .
  unconstrained affine. Each row reports `ess`, `frac_neg`, `neg_mass` beside its RMSD, so the
  ceiling is priced by conditioning and not just by A. **Every row is ORACLE / NOT DEPLOYABLE.**
- **F-C1d — REGISTERED FALSIFIER.** The ORACLE unconstrained affine ceiling over the top-128 is
  **vacuous by dimension counting** (3n <= 48 against 128 generic windows), so it reaches ~ 0 with
  ill-conditioned cancelling weights. *Fires (i.e. the ceiling is informative rather than
  degenerate) if the ORACLE affine solution has `ess` >= 5 and `neg_mass` <= 1.0 on a majority of
  targets.* I expect it to be refuted; `s27/s28_A_FINDINGS.md:85-90` already calls the affine hull
  a tautology and says not to quote it as a bound. **This is a pre-check that can only remove my own
  excuse** (contract rule 22): if refuted, no affine ceiling is quotable and the ladder's x-axis
  must be the constraint, not the A.
- **F-C1c stands unchanged** (identifiability of the ORACLE hull weights from native-free features,
  5% out-of-fold R^2 bar). It is the one genuinely new C1 arm and it is the one that produces a
  closure rather than another null.
- **F-C1a and F-C1b stand**, but are now evaluated against the **constrained-affine** ladder rather
  than against W1-W5, whose verdicts are inherited from the table in A1.1 and not re-measured.

## A1.4 One correction to a published S30 headline, registered before I measure it

`s30/results/s30_Q_sparse.json` contains set-matched rows (`T128_*`) that S30-L11 did not use in its
headline sentence. **I register in advance that I expect the set-matched ladder to show the sparse
weighted readout beating the argmin at a fixed candidate set**, and that S30-L11's closure survives
only in its *price across sets* form. If the set-matched rows do not show that, this paragraph is a
failed prediction and will be reported as one.

## A1.5 R1's scope, registered as a disagreement before I measure anything

S31-L1 states that the quantum stage "cannot emit anything outside the pool" and "carries at most k
bits". **That is true of the SELECTION readout `consensus_medoid(block, p)` (`core/pipeline.py:795-803`)
and false of `average_weighted` (`core/pipeline.py:880-895`)**, which is also shipped, is called at
`:1110` and `:1116`, is scored as `rmsd_q_avg` / `rmsd_q_synth` / `rmsd_u_synth` at `:1179-1182`,
and emits a p-weighted coordinate average that is generically not a pool member. Because `p` is a
measurement distribution it is non-negative, so that readout is **convex**, not affine: its
reachable set is the convex hull of the posed candidates. The affine readout (`s27/s28_A_amp.py:105-117`)
is in the measurement harness, **not** in `core/pipeline.py`. Registered here so the disagreement
is on record before, not after, my numbers.

## A1.6 Revised prior odds

C1 as a source of a deployable gain: **now 15:1 against**, up from 6:1, on the strength of A1.1.
C3 and C2 are unchanged. **The lane's expected output is a closure and a correction, not a gain**,
and that is registered here rather than discovered later.

*Amendment committed before the first lane-C number.*
