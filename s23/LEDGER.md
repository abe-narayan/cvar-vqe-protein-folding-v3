# SPRINT 23 — DECISION LEDGER

## L1 — THE GEOMETRY OF THE INCUMBENT'S OUTPUT, AND WHY KABSCH MAKES IT COST RMSD

The incumbent emits a **coordinate average** of the top-75 pool members. Measured at n=126:

    mean virtual Ca-Ca bond   2.961 A     native 3.812 A     -> 22% SHORT
    radius of gyration        6.206 A     native 6.601 A     ->  6% short

**Kabsch superposition fits rotation and translation only — it does not fit scale.** So a systematic
scale error propagates directly into Cα-RMSD. That makes a scale correction a legitimate,
native-free, one-parameter lever, and it had never been tested.

**The geometry is not a uniform shrink, and that is the key observation.** The *bond* is 22% short
while the *envelope* is only 6% short. Averaging **smooths**: it shortens local bonds while largely
preserving global extent. **A single global scale cannot fix that — it trades one error for the
other** — which is exactly what the numbers below show, and it is why the local-geometry route
(torsion-space rebuild, Cα-preserving repair) is the better-motivated one.

## L2 — THE PER-TARGET SCALE IS REAL AND UNIVERSAL. THE GLOBAL CONSTANT IS WORTH NOTHING.

`s23/gscale.py`, n=126, nested CV over the pinned 5 folds.

    incumbent, s = 1.0                    3.0483
    GLOBAL s, NESTED CV (held out)        3.0486    s = 0.9738    +0.0002  63W/63L   MDE 0.017
    global s, in-sample (leaky)           3.0424    s = 0.9725    -0.0059            not measured
    per-target s* (ORACLE)                2.9110                  -0.1374  126W/  0L  BEATS

**My pre-registered falsifier FIRED.** A single global scale fitted by nested CV is worth
**+0.0002 Å — exactly nothing, at 63W/63L.** Even the *leaky* in-sample fit reaches only −0.0059,
below its own MDE. **The global-scale route is refuted.**

**But the ORACLE is 126 wins and 0 losses.** Every single target improves under its own scale. The
reason a constant cannot capture it is the spread: **s\* has sd 0.076, 10th percentile 0.90, 90th
percentile 1.10.**

### CORRECTED BY WORKSTREAM D: the ceiling is 2.5x what I reported

> I noted the censoring below but did not fix it. **D did, in closed form.** Kabsch rotation is
> provably **scale-invariant**, so the optimal scalar is exactly `s* = Ct/Cc` -- no grid, no
> clipping. **The true ceiling is -0.3403 A, not the -0.1374 I reported** (and 2x the n=40 figure
> of -0.158). It does not rescue anything; **it makes the unreachable ceiling bigger.**

### And my grid was CENSORED, which is how the error arose

    pinned at the LOW grid endpoint (s < 0.90):  42 / 126
    pinned at the HIGH grid endpoint (s > 1.10): 21 / 126
    interior optima:                             63 / 126

**Half the panel wants a scale outside the grid, and it wants it in BOTH directions** — 42 targets
want *more* contraction, 21 want *more* expansion. **The optimal scale is bimodal.** That is a
complete explanation for two things at once: why no global constant can work, and why no single
monotone feature should be expected to predict it.

## L3 — NO NATIVE-FREE FEATURE PREDICTS THE OPTIMAL SCALE

Eight native-free candidates against s\*, n=126, Spearman:

    ensemble spread   -0.099     spread/Rg      -0.102     Rg of the average  -0.014
    predicted Rg      -0.017     Rg ratio       -0.013     chain length       -0.022
    score sd          +0.071     virtual bond   +0.090

**Every one is inside ±0.11.** The ensemble's own spread — the mechanistically obvious predictor,
since a more diverse ensemble should contract more — carries **−0.099**, nothing.

> **THIS IS THE THIRD INDEPENDENT INSTANCE OF ONE PATTERN.** A per-target optimum that is
> **real, large and near-universal**, and **invisible to every native-free feature tried**:
>
> | quantity | oracle value | W/L | best native-free predictor |
> |---|---|---|---|
> | per-target averaging width m\* | −0.243 Å | 77/49 | ~0% captured, 4 router constructions |
> | per-target arm choice | −0.482 Å ceiling | — | ~0% captured, 5 configurations |
> | **per-target scale s\*** | **−0.137 Å (censored)** | **126/0** | **all \|ρ\| ≤ 0.11** |
>
> Three different quantities, three different mechanisms, the same result. This is no longer a
> series of failed experiments; it is a **property of the instrument** — and the finite-sample bound
> from Sprint 22 says why: at n≈100 per fold, no router class rich enough to express a per-target
> quantity is learnable.

## L4 — ITERATED PROCRUSTES DOES NOT HELP, AND THE CONTRACTION IS INTRINSIC

The shipped `coordinate_average` does **one** superposition onto the medoid and stops. A generalised
Procrustes mean — superpose onto the running mean, recompute, iterate to convergence — is the
textbook Fréchet mean and is a **global** rule change requiring no per-target decision. n=45:

    shipped single-pass (medoid frame)   2.9629    virtual bond 3.040
    iterated Procrustes (converged)      2.9711    virtual bond 3.035
    delta  +0.0082  SE 0.0054  MDE 0.0152  [-0.0009,+0.0204]  16W/29L

**No improvement, and — the informative part — the contraction is unchanged (3.035 vs 3.040).**

**So the contraction is not an artefact of the reference frame.** It is intrinsic to averaging
structures that genuinely differ: the mean of points on a curved manifold lies inside it. Choosing a
better frame cannot fix it. **Only changing the space you average in (torsions) or repairing
geometry afterwards can**, which is precisely where the remaining live experiments sit.

---

## L5 — WORKSTREAM A: ALL THREE FALSIFIERS FIRE, AND ONE FAILURE EXPLAINS THE WHOLE FAMILY

126/126, bit-verified against `s21/results/poolgap.json` and the pinned incumbent constant 3.0483.

**(1) Cluster-restricted averaging.** Every achievable cluster-choice rule (size / score / combined /
random) at every (m, k) is **at or worse than the incumbent**. The **oracle ceiling is real and
substantial — −0.26 to −0.48 Å, 2.2–3.5× its own MDE, CI excluding zero — and it GROWS with m.**
Unreachable by any declared rule.

**(2) Width under clustering.** The optimal m *does* shift 75 → 150 under clustering, as hypothesised.
But `clust_m` at its own best m is **3.140 against plain `avg_75` at 3.048**, and worst-target and
failure rate do not favour clustering either.

**(3) Weighted averaging — the cleanest closure in the sprint.** For rank-power and
distance-to-medoid weighting, **even the leakage-oracle in-sample grid search selects the
UNIFORM-WEIGHT point as optimal.** There is **no headroom to reach, not merely an unreached one.**
Score-softmax shows a hint in-sample at **0.12× its own MDE** and evaporates under nested CV.

### The self-damaging detail, and it is the most valuable thing in the lane

The real `combined` rule at m=75, k=3 — which strips ~4 genuine geometric outliers from the top-75 —
is **significantly WORSE than the incumbent (+0.142, CI excluding zero)**, while a **matched
random-partition control removing the same number of members is a dead null.**

> **So removing four members at random costs nothing, and removing the four most geometrically
> deviant members costs +0.142 Å. The outliers are not noise — they are load-bearing.**
>
> **Mechanism: those members carry i.i.d. error that CANCELS in the average.** Genuine clustering
> selectively removes exactly the error-cancelling members, because "geometric outlier" and
> "carries error that cancels" are the same property viewed twice.

**This is the mechanistic explanation for why the entire consensus family fails**, and it has been
sitting unexplained across nine arms and three sprints. Coordinate averaging works *because of*
diversity, not despite it — it is a variance-reduction operator, and every method that tightens the
retained set (medoid, clustering, sharp weighting, small m) attacks the mechanism that makes it work.
It also predicts, correctly, that *uniform* weights are optimal: uniform weighting is the
minimum-variance combination of exchangeable estimators with i.i.d. errors, which is exactly the
model this result implies.

**Connects to `error-coherence-decides-correctors` and `decorrelated-errors-exist-but-are-unusable`
from project memory — same physics, arrived at from the aggregation side.**

---

## L6 — WORKSTREAM D: THE SCALE SIGNAL IS REAL, TRANSFERS ALMOST PERFECTLY, AND IS UNREACHABLE IN PRINCIPLE

**(a) My ceiling was wrong and is now closed-form.** −0.3403 Å, not −0.1374. See the correction above.

**(b) The global-scalar arm is refuted at its THEORETICAL ceiling, not merely unmeasured.** D derived
the population-optimal global scale from the exact closed form — no grid — and it captures **1.7% of
the true ceiling, best case, IN-SAMPLE** (0.6% excluding FAIL18). **Mechanism: s\* splits 53/126
wanting expansion against 73/126 wanting contraction, so the corrections cancel.** D independently
reproduced my +0.0002 nested-CV null to five decimals and found no leakage bug.
**H1's global-scalar sub-arm is CLOSED. No further tuning can find anything there.**

**(c) D refuted its own registered hypothesis, in the useful direction.** It expected partial
transfer, like Sprint 22's m-ladder at 65%. Fitting s\* on half a target's pool and applying it to
the disjoint other half transfers at **98–100%**, clearing 2.2–2.35× its own MDE on both iid and
fold-clustered CIs. **The cleanest result of this shape in the project's history: the per-target
scale is real, stable, and not a selection artefact.**

**(d) THE REFRAMING, AND IT CHANGES THE DISPOSITION FROM "NOT FOUND" TO "NOT FINDABLE".** D's first
placebo was confounded; it withdrew it, documented why, and rebuilt it. The corrected version —
cross-target, same-length real natives — shows that **a randomly mismatched same-length native
produces an equal or LARGER apparent gain than the true native.**

> **So s\* is not a property of the target's fold. It is fitting the particular
> (pool, reference-structure) pair.** A different real structure of the same length wants a
> different scale and is equally well served by it.
>
> **That makes the scale correction unreachable IN PRINCIPLE from native-free information, not
> merely unreachable at n=126.** It is the first quantity in this programme with that status. Every
> previous per-target signal was *plausibly* predictable and simply was not predicted; this one
> provably carries no native-free content, because the thing it depends on is the answer.
>
> It also resolves the tension between D's own two findings — 98–100% transfer across pool halves
> (same native) and none across natives (same length). **s\* is stable in the pool and unstable in
> the reference. It is a property of the pair.**

**(e) Concentration, flagged as the brief requires.** **FAIL18 — 18 of 126 targets, 14% of the panel —
carries 57% of the true ceiling.** A per-length constant reaches 7.6% in-sample at bucket sizes of
9–23, which the standing finite-sample bound already forbids. Not pursued.

**D's closing recommendation, adopted:** the productive question is no longer scale architecture. It
is **why FAIL18's retrieval pools are size-distorted in the first place** — a candidate-generation
diagnosis, not a selection one.

---

## L7 — WORKSTREAM C: A FIFTH ROUTER FAILS, AND THE SIGNAL THAT FAILS IS A NEW KIND

`c2_rgcond.json`, complete. Nested 5-fold CV, aggregation width routed on **rg_z** — the
compactness-disagreement signal, the only native-free quantity in this project to survive every
control (partial ρ 0.34–0.38 against in-band *ordering*; survives length, difficulty and
realised-extension controls; **not** a difficulty proxy at r = +0.126).

    routed(rg_z) vs fixed m=75:   +0.0012 A   SE 0.0196   MDE 0.0550   -> 2% of its own MDE

**Essentially exactly zero on the full instrument, with no subsampling.** Secondary `rg_gap` arm
−0.0017, also null. The **reversed-direction adversarial control self-disables** — nested CV picks a
threshold outside the data range on every fold rather than manufacture a split — which rules out a
CV-construction artefact and is the right behaviour.

**Why this null is new information rather than a repetition.** The four routers that failed in
Sprint 22 all used functionals of the objective's own score distribution, and the diagnosis was
*"the objective cannot audit itself."* **rg_z is not such a functional** — it compares the
distogram's *prediction* against the pool's *realisation*. It clears every control on the ordering
task and still transfers **nothing** to the aggregation-width lever.

> **So "the objective cannot audit itself" was too narrow a diagnosis.** A signal from outside the
> objective, with demonstrated skill on a related task, also fails to move this lever. The binding
> constraint is the one the finite-sample bound names, not the provenance of the feature.

---

## L8 — THE PROBABILITY-WEIGHTED READOUT FAILS. THE LAST QUANTUM-SIDE DOOR IS CLOSED.

`c1_probweight_analysis.json` + `c1_explore_t0_analysis.json`, n=126, 4 seeds, 9 qubits, exact
`StatevectorCircuit`, α=0.15. *(Workstream C stalled in a wait loop on completed artefacts; the
analysis below is the coordinator's, on C's data and against C's own pre-registered comparisons.)*

Sprint 22 proved by set equality that the CVaR tail's **membership** is always the classical top-m
set, and named **one** object the proof did not cover: a **probability-WEIGHTED** average over that
tail, where the amplitudes enter the output rather than only selecting a subset. This is that test.

**Gate 1 — the theorem reproduces exactly.** Trained-tail SET == classical `argsort(score)[:m]` SET
at a **100% pass rate**, for both the trained and the untrained circuit, on every target.

### Registered primary, T = 0.5

    W_trained (probability-weighted)   3.0621
    C_trained (uniform, IDENTICAL set) 3.0408
    W - C  = +0.0213   SE 0.0117   MDE 0.0327   fold[-0.0012,+0.0550]   0.65x MDE   54W/72L
    -> NULL.  The weighting does not help.

    training effect on the WEIGHTED average   -0.146 [-0.232,-0.092]   MEASURED
    training effect on the UNIFORM  average   -0.016                   much smaller

**Training does move the weighted readout — by −0.146 Å — but it moves it *toward* the classical
answer, not past it.** And neither beats the incumbent (3.048).

### Exploratory, T = 0 — the regime where the distribution is genuinely non-uniform

C flagged, against its own hypothesis, that at T=0.5 the trained state holds **8.85 of 9 bits** —
near-uniform — so a weighted average is close to a uniform one *by construction*, making the primary
a weak test of weighting specifically. It registered a dated exploratory addendum at T=0 and ran it
after locking the primary. Correct handling, and the result completes the picture:

    trained entropy 4.69 bits, realised tail m collapses to 1.73
    W_trained 3.4298   C_trained 3.3913
    W - C = +0.0386 [+0.018,+0.054]  MEASURED -- weighting is significantly WORSE
    training effect on the uniform average: +0.335 -- training HURTS badly here

> **Both regimes are now covered and both close the door.** Where the distribution is near-uniform
> the weighting cannot differ from uniform and does not. Where CVaR pressure makes it genuinely
> concentrated, the weighting is **significantly worse**, and the collapse to m ≈ 1.7 reproduces the
> pipeline's own documented failure mode (α=1, T=0.1 → 0.076 bits → the plain argmin).
>
> **The probability-weighted readout was the last object the Sprint-22 set-equality theorem did not
> cover. It is now measured, in both regimes, and it does not beat the classical top-m bar.**

**Uniform weighting winning is not an accident** — it is exactly what L5's error-cancellation
mechanism predicts. Uniform is the minimum-variance combination of exchangeable estimators with
i.i.d. errors, and every departure from it (sharper weights, clustering, smaller m, amplitude
weighting) attacks the variance reduction that makes averaging work. **Three lanes reached that
conclusion this sprint by three different routes.**

---

## L9 — THE POOL'S ERROR IS 68% SHARED. THAT ONE NUMBER EXPLAINS EVERY NULL IN THIS SPRINT.

`s23/errdecomp.py`, `results/errdecomp.json`, n=126, complete. Written and forked before running.

Write the m=75 retained members, in the medoid frame the average is actually taken in, as
`w_k = t + e_k`. The emitted cloud is `c = t + ebar`. Split each member's error about the ensemble
mean, `d_k = w_k - c`. Then **exactly** (verified to 2.7e-14 Å):

    mean_k |e_k|^2   =   |ebar|^2   +   mean_k |d_k|^2
                         COMMON         IDIOSYNCRATIC
                         160.36         63.82
                         INVISIBLE from inside the pool      fully observable

**A bias shared by every member moves `c` and every `w_k` together and leaves every within-pool
statistic unchanged.** So the pool can measure the second term exactly and the first term not at all.

    common-mode fraction f = |ebar|^2 / (|ebar|^2 + mean|d_k|^2)
      mean 0.676   median 0.674   10th 0.372   90th 0.961   quartiles [0.107 .507 .674 .855 .999]
      i.i.d. model predicts f = 1/m = 0.0133      OBSERVED / PREDICTED = 50.7x
      f > 1/m on 126 / 126 targets

**Two-thirds of the squared error of a retrieved pool is a bias its members hold in common.**

### The derived estimator, and why its failure is informative

If the errors *were* i.i.d., the invisible term would be recoverable from the visible one as
`|ebar|^2 = mean|d_k|^2/(m-1)`, giving a native-free, parameter-free, *derived* scale
`s_hat = 1 - mean|d_k|^2/((m-1)|c|^2)` — the only scale estimator this project has proposed that
comes out of a model rather than a grid.

    s* (ORACLE, closed form)      mean 0.9417   sd 0.2580   range [0.252, 1.924]
    s_hat (DERIVED, native-free)  mean 0.9983   sd 0.0015   range [0.994, 1.000]
    rho(s_hat, s*) = +0.094

    s_hat vs incumbent   -0.0010  SE 0.0005  MDE 0.001  fold[-0.0018,-0.0002]  73W/53L  at the MDE
    s*    vs incumbent   -0.3403  SE 0.0561  MDE 0.157  fold[-0.4230,-0.2478] 126W/  0L  BEATS

**The i.i.d. model underestimates the common term by a factor of 1410.** `s_hat` collapses onto
1.000 and captures **0.3%** of the ceiling. The falsifier fired — and it fired *because* the model
it was derived from is wrong by three orders of magnitude, which is the informative part.

> This independently reproduces D's closed-form ceiling **−0.3403 Å at 126W/0L** on a separate
> code path, and it extends L2's censored picture: s* actually ranges **0.252 to 1.924**, 73 want
> contraction and 53 want expansion.

### What this explains, and it is most of the sprint

    s* = <c,t>/|c|^2 = 1 - <c,ebar>/|c|^2

**The optimal scale is a function of `ebar` and of nothing else** — i.e. of exactly the component
the pool cannot see. That is a *derivation* of L6d's empirical placebo, not a restatement: the
reason a mismatched same-length native serves just as well is that `s*` never depended on the
target's fold in the first place, only on the offset between the cloud and whatever reference it is
scored against.

It also prices the operator honestly, in both directions:

    typical retained member  3.7037 A        consensus average  3.0483 A
    coordinate averaging is worth  -0.6554 A  (17.7%),  126 / 126 targets

**Averaging removes the 32% that is idiosyncratic and cannot touch the 68% that is shared.** That
is the single mechanism behind L5 (removing geometric outliers costs +0.142 Å — they carry the
cancellable part), behind uniform weights being optimal, behind five failed routers (L7), behind the
probability-weighted readout (L8), and behind the scale route (L2/L6). **Every one of those methods
operates on the within-pool distribution, and the within-pool distribution is the 32%.**

`f_common` correlates **+0.351** with the incumbent's own RMSD: the targets where the consensus
fails are exactly the ones where the shared bias dominates. `f` is an ORACLE quantity (it needs the
native) and is a diagnosis, not a predictor.

> **The binding constraint is not selection, aggregation, routing, or readout. It is that the
> retrieval library hands back candidates that are wrong together.** The remaining leverage is in
> candidate GENERATION — which is where D's independent closing recommendation (diagnose why
> FAIL18's pools are size-distorted) also points, arrived at from the other side.

---

## L10 — THE MECHANISM'S OWN FORWARD PREDICTION, TESTED. RIGHT DIRECTION, NO MEASURABLE SIZE.

`s23/diverse.py`, `results/diverse.json`, n=126, complete. Forks and falsifiers written before running.

L5 and L9 say coordinate averaging works *because of* diversity: the outliers carry the error that
cancels, and the 32% that cancels is the only part any within-pool operator can touch. **Every arm
this project has run TIGHTENS the retained set. The mechanism's own prediction is the opposite
move, and it had never been tested.** All four arms draw exactly 75 members from the same top-150
score band, so they differ in diversity and in nothing else.

    arm                                   RMSD     mean pairwise spread of the retained set
    A0 incumbent  top-75 by score       3.0483            2.486 A
    A1 DIVERSE    farthest-point        3.0587            2.925 A
    A2 RANDOM     matched null          3.0706            2.705 A
    A3 TIGHT      closest-point         3.0923            2.117 A

**The predicted ORDER holds exactly — A1 < A2 < A3, monotone in the spread of the retained set —
and every effect is far under its own MDE.**

    A1 DIVERSE - A2 RANDOM   PRIMARY    -0.0119  SE 0.0204  MDE 0.057  fold[-0.045,+0.009]  75W/51L
    A3 TIGHT   - A2 RANDOM   reversed   +0.0217  SE 0.0329  MDE 0.092  fold[-0.044,+0.075]  55W/71L
    A1 DIVERSE - A0 incumbent           +0.0104  SE 0.0326  MDE 0.091                       66W/60L

**Both falsifiers fired.** 0.21x and 0.24x of their own MDEs. **NOT DEMONSTRATED.** Buying
diversity deliberately trades away score quality at almost exactly the rate it buys cancellation,
and the two cancel to within a hundredth of an Ångström. The sign of the whole ladder agrees with
L5/L9 and the magnitude does not survive contact with the instrument — which is the honest reading,
and it is recorded that way rather than as "the mechanism is confirmed".

### Two by-products worth more than the primary

**(a) The score gate is worth ~0.02 Å between ranks 75 and 150.** A2 — seventy-five members drawn
at RANDOM from the top-150 — costs only **+0.0223 Å [-0.006,+0.054], 66W/60L** against the shipped
top-75. Inside the band, the shipped score orders candidates essentially not at all. Independent
confirmation of `in-band-is-the-only-ranking-metric` from the aggregation side.

**(b) A best-of-K trap, caught by its own null.** Across the 8 random draws the per-target sd of
RMSD is 0.0676 Å, so *which* 75 you average genuinely moves the answer, and the best of 8 draws
reaches 2.9716 — **−0.0767 Å against the incumbent at 84W/42L**, which looks like a fifth
unreachable per-target oracle.

> **It is not.** `E[min of 8 standard normals] = −1.4236`, so pure draw noise predicts a best-of-8
> gain of **−0.0962 Å** against the observed **−0.0990**. The "oracle" is **101% accounted for by
> the null**. There is no exploitable subset structure here at all.
>
> This is the distinction the project's own rule demands and it separates this result from the
> others in the sprint: s\* is a **closed-form optimum at 126W/0L**, far outside any best-of-K
> null; subset choice is a **minimum over draws** and lands exactly on it.

---

## L11 — WORKSTREAM B: H3 CLOSED BY SOURCE AUDIT; H5 REPAIR STRICTLY COSTS RMSD AT EVERY SETTING

`s23/agentB_h5_carestraint.py`, `results/agentB_h5_carestraint_sweep.json`, **complete, 30/30**.
Pre-registered in `s23/PREREG_B.md`. Sole authorised AMBER/OpenMM lane this sprint; released.

### (a) H3 — the 1.024 Å torsion-averaging loss is NOT a branch-cut bug

`s14/avgspace.py:82-83` calls `s14.retprior.circ_mean`, which is `arctan2(mean sin, mean cos)` —
a **genuine circular mean**. The contrast reproduces exactly from the persisted artefact:
**+1.0241 [+0.6984, +1.3647], n=126, 41W/85L.** Per the BRIEF's own branch rule the direction stays
closed and no rerun was registered. **Not a bug; a result.**

**Named, not chased:** `circ_mean` is unweighted and not basin-conditioned, and a circular mean of a
genuinely bimodal per-residue distribution lands *between* the modes. The audit rules out a bug; it
does not certify that basin-conditioned torsion averaging is worthless. Still open.

### (b) H5 — Cα-preserving repair: the falsifier fires at every one of 17 settings

Mechanism verified from source first: `core.amber.RESTRAINED_BACKBONE = ("N","CA","C")` is a
**module-level global read once at Hamiltonian construction**, not a per-call argument — so the
deployed restrained minimisation restrains the whole rigid peptide-plane frame and *never Cα alone*.
B built a genuinely separate Hamiltonian with that global monkey-patched to `("CA",)`, asserted the
restrained-atom index counts (1× and 3× the residue count), and swept both scopes.

n=30, the project's own native-free fold×length-stratified subsample. **Point cloud on both sides.**

    incumbent point cloud (no repair at all)                    3.4928
    k = 0   unrestrained AMBER minimisation      3.9396   <b>+0.4468</b>  [+0.319,+0.592]   6W/24L  WORSE

    scope Ca   k=1     3.7790  +0.2862  ...  k=100   3.5697  +0.0769  ...  k=3000  3.5119  +0.0191
    scope full k=1     3.7588  +0.2659  ...  k=100   3.5851  +0.0923  ...  k=3000  3.5072  +0.0143

**Every cell is at or worse than doing nothing.** Eleven of seventeen are significantly worse past
their own MDE with fold CIs excluding zero; the remaining six are "null" only because at k ≥ 1000 the
restraint has turned the repair into a **no-op** (mean Cα displacement 0.16–0.22 Å). The ladder is
monotone in exactly one direction: **the less the repair moves the structure, the less it costs.**

The RMSD cost tracks the Cα displacement across the whole ladder — cost/displacement ≈ 0.11–0.26 —
which is what happens when the repair moves atoms in directions **uncorrelated with the existing
error**: it adds error roughly in quadrature and removes none.

### (c) The genuinely new part: Cα-only restraint does not buy RMSD, it relocates the damage

    ca minus full, Ca-RMSD:  null at 6 of 8 k values; +0.057 at k=3 and +0.047 at k=10 (Ca-only WORSE)
    ca minus full, cis_frac:      +0.121 (k=30)  +0.422 (k=100)  +0.550 (k=300)  +0.638 (k=3000)  MEASURED
    ca minus full, omega_dev:     +23.3   +65.1   +74.7   +83.8 degrees                            MEASURED
    ca minus full, bond_len_dev:  -0.009  -0.035  -0.081  -0.153                                   MEASURED

> **Restraining Cα alone frees the peptide planes to flip.** At k ≥ 100 the Cα-only arm carries a
> **cis fraction of 0.51–0.64 and mean ω deviation of 92–104°** — chemically absurd structures that
> happen to have the right Cα positions. The full-backbone arm keeps ω honest and pays in bond
> lengths instead (relative deviation up to 0.235 at k=3000).
>
> **So the choice is not "how much geometry do I fix" but "where do I put the damage", and neither
> placement buys a single hundredth of an Ångström.** H5 — the BRIEF's "high-priority, low-risk" item
> — is **REFUTED**, and refuted with a mechanism rather than a shrug.

**This is the same statement as L9 from the physical side.** Repair operates on the emitted cloud,
i.e. downstream of retrieval, where only the 32% idiosyncratic component is visible. It cannot see
the 68% shared bias, so the best any repair can do is leave the structure alone — which is precisely
what the k→∞ limit of this sweep does, and precisely where the ladder bottoms out.
