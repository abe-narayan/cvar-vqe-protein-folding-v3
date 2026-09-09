# WORKSTREAM D — SPRINT 24. QUANTUM INTEGRATION AND HAMILTONIAN ROLES.

Pre-registered in `s24/PREREG_D.md`; fork lists sent to the coordinator before each run.
Point-cloud basis throughout unless a row says otherwise. ORACLE / ACHIEVABLE / PRODUCTION
labelled at every appearance. No benchmark inspection anywhere in this workstream.

---

## 0. THE RESULT IN FIVE LINES

1. **Set equality is a THEOREM ABOUT THE ALGORITHM, not an artefact of this pool.** The
   selection rule is fixed by the code and no candidate manifold can change it. The direction
   is closed and should not be funded.
2. **The integration harness is built, unit-tested 18/18, and ready** for Lanes B and C.
3. **D1 (Hamiltonian disagreement) is a NULL that closes the functional lever** alongside the
   coordinator's provenance lever — **with one MEASURED positive inside it**: `AMBER_PREFERS`
   selects along a materially non-parallel direction (bias cos **0.5693**, −0.2154 past its own
   shared-referent floor, 1.41× MDE, 25W/5L), the only source measured anywhere in this sprint
   to do so. It fails on **quality** (4.13 Å, q = 1.520), which is the whole story.
4. **D2 (the AMBER-informed SCORE) is a NULL with an ORACLE ceiling of 0.0148 Å**, closing the
   functional lever in the third and last of its available forms. Nested CV chose to
   **anti-weight** AMBER, which is itself consistent with D1's mechanism.
5. **Seven of my own errors are on the record in §6**, including a green light from a dead
   assertion, a headline prediction the measurement contradicted, an over-claim withdrawn when
   L5 was retracted, and a verdict rule in my own shipped harness that flattered its user.

---

## 1. DELIVERABLE 1 — THE SOURCE-LEVEL DETERMINATION

### 1.1 The standing statement

> **The realised CVaR tail's support is always a SUBSET of an initial prefix of the energy
> order, and equals that prefix exactly when every state in the prefix carries positive
> probability. p_theta can delete a member; it can never add one outside the classical
> top-m.**

This is the coordinator's wording, adopted verbatim after it corrected a looser first draft of
mine (§1.4). It is used unchanged in `d_setequality_proof.py`, in `d_harness.py`, and in
`s24/results/d_setequality.json`.

### 1.2 The exact lines that establish it

**`core/quantum.py:229-272`, `tail_indices(energies, alpha)`.** The signature takes **no
probability vector and no circuit parameters**. `k = ceil(alpha * n)` (line 252) depends on
`(n, alpha)` only. The mask is built either as `argsort(e, kind="stable")[:k]` (lines 256-259)
or as `e < q` plus the lowest-index ties (lines 261-271). **The returned tail is a pure
function of the energy array.** `p_theta` is not in scope, so on this path there is nothing for
a pool to change.

**`core/quantum.py:336-359`, `cvar_from_probs(energies, probs, alpha)`** — the one place
`p_theta` reaches the tail at all:

    order = np.argsort(e, kind="stable")                       # line 352   ENERGY order
    cum   = np.cumsum(p[order])                                # line 353
    take  = np.clip(alpha - (cum - p[order]), 0.0, p[order])   # line 354

`take[j] > 0` requires **both** conjuncts of that clip:

* **(i)** `alpha - (cum[j] - p[order][j]) > 0`, i.e. the mass strictly before position `j` **in
  the energy order** is below `alpha`. `cum - p[order]` is the exclusive prefix sum, which is
  non-decreasing because `p >= 0`, so it crosses `alpha` exactly once. Conjunct (i) alone
  defines an initial **prefix** of the energy order, ending at a cut `m_cut`.
* **(ii)** `p[order][j] > 0`. A zero-probability state inside that prefix is **punched out**.

So the support is a prefix **with holes**, and every hole is exactly a zero-probability state.

### 1.3 It is the ALGORITHM, and the consequence

The trained state can move exactly two things and no third: **where the prefix cuts** (a rung
`m` on the classical top-m ladder, plus which prefix members it deletes where it has exact
zeros), and **the weights inside the prefix** — which s23 L8 closed at four temperatures.
**It can never ADD a candidate the classical energy order excluded.** Membership is bounded
above by `argsort(E)[:m]` as an identity, independent of the landscape's shape, degeneracy,
multimodality and scale, and therefore independent of the candidate manifold that produced it.
**Changing the pool cannot break it.** s22's 100% Gate-1 pass rate over 2016 cells was not
evidence about the pool; it was the code being read back.

Deletion is not a hidden information channel. The deleted member is chosen by the amplitude
pattern, which is a function of the **same energies** the classical rank order already uses.
The old "sampled path" scope condition is the **same fact by a second mechanism** — exact zeros
delete on the exact path, failure to draw deletes on the sampled path, neither can add — and it
is reported as one statement rather than two.

### 1.4 CORRECTION OF RECORD, AND AN ERROR OF MINE WORTH READING

**My REV1 claimed the support is always an INITIAL PREFIX, dropping conjunct (ii), and reported
0/2916 failures.** The coordinator caught the wording. The reason my test agreed with a false
claim is worse than the claim: **every probability family in REV1 carried an epsilon floor** —
`zeros_mostly` was `np.where(rng.random(n)<0.1, 1.0, 0.0) + 1e-30`, the point masses used
`1e-18`, `random` used `+1e-6`. No exact zero was ever generated, `p > 0` held on all 2916
cells, and **the prefix-hood assertion could not fire.** It was a dead assertion reporting a
pass, which is the most dangerous kind of green light: it looks like the strongest possible
evidence and carries none.

REV2 (`revision: 2` in the artefact, with REV1's failure described in `revision_note`) adds five
exact-zero families, **asserts** the true theorem, and **measures** prefix-hood separately:

    n_cells                                                3888
    subset-hood violations                                    0     <- THE THEOREM
    holes that were NOT exactly zero-probability states        0     <- the mechanism, direct
    prefix-hood violations (REV1's over-claim)             1164     = 29.9%
    full-support cells                                     2268
    full-support cells with EXACT equality           2268 / 2268     = 100%

Adversarial by construction: 6 register sizes x 9 energy structures (all-tied, heavy
degeneracy, bimodal, 1e7 scale, near-constant, a single 1e9 outlier, discrete, monotone ramp)
x 12 probability structures x 6 alphas. A separate check varies `p` over 4,000 random
full-support draws at fixed `e`: 66 distinct realised tail sizes, and **exactly one distinct
tail set for every one of the 66**.

The coordinator's independent run (17,574 trials) found 0 subset violations and a 58.8% hole
rate; mine finds 0 and 29.9%. The rates differ because the family mixes differ; **the two
assertions agree exactly.**

**Equality is the empirical regime; subset-hood is the theorem.** A trained RY/CNOT state has
generic angles and therefore full support, which is why equality holds in every cell anyone has
actually run.

`s24/d_setequality_proof.py` -> `s24/results/d_setequality.json`. Registered as a **property
test**, not a directional experiment: it reads no RMSD and has no outcome that can favour a
hypothesis, so Rule 0's fork enumeration does not apply and was not claimed. A falsifier was
registered instead and did not fire.

---

## 2. DELIVERABLE 2 — THE INTEGRATION HARNESS

`s24/d_harness.py` + `s24/d_harness_test.py`. **18/18 tests pass.** Ready before either
generator's pools.

    from s24.d_harness import Candidates, score_shipped, run_both, aggregate
    cand = Candidates.from_arrays("1L2Y", W, PHI, PSI, source="laneB_residual_v1")
    row  = run_both(cand, score_shipped(cand), alpha=0.15, T=0.5, seed=0)
    agg  = aggregate(rows)

**Why it exists in this form.** The most repeated error in project memory is
`control-must-match-the-operators-space` — three instances in two sprints. So the control is
made *structurally* impossible to get wrong: `run_both` hands the **same** `Candidates` object,
the **same** energy vector, the **same** readout and the **same** frame to both arms, and there
is no API path that varies one arm's inputs without the other's.

**The quantum arm is genuine.** Candidate-identity register `H|i> = E_i|i>` (`s22.qcand_lib`
`Encoding`, reused not reimplemented), exact `StatevectorCircuit`, Adam on the **exact
parameter-shift gradient** via `run_cvar_vqe` — so the historical `baseline="tail"` CVaR
gradient defect (cos +0.6556 at 0.758x norm) cannot arise on this path — and the tail read
exactly off `cvar_from_probs`, not sampled. No shot noise anywhere.

**Controls shipped alongside, every cell:** classical top-m at the VQE's **realised m** (the
size-matched bar), classical top-75 (the **production** rung), and a zero-information null of
random m-subsets through the **identical** coordinate average.

**Discipline built into the code rather than left to the caller:** `aggregate` **refuses** to
combine point-cloud with built-chain rows; `write` sets the completion flag on the **full key
set**, not a row count; `paired_stats` reports SE, MDE = 2.8016 x SE per comparison,
effect/MDE, iid CI **beside** fold-clustered CI, W/L, median beside mean, and worst-target
degradation, and returns the verdict `TYPE-M ZONE` for anything in 0.7-1.3x MDE; seeds come
only from `s15.seed.stable_rng`; `AmberLock` is `os.open(O_CREAT|O_EXCL)` and announces on open
and release.

**Tests that would catch a real corruption**, not just a crash: both arms provably move
together under an energy perturbation; the gate **fires** on a deliberately desynchronised set;
the score is verified native-free by recomputing it on a container with `nat_ca=None`; the
readout is invariant to a rigid motion of the whole pool; an effect placed at exactly 1.0x its
own MDE is reported as `TYPE-M ZONE`; the register ceiling **raises** rather than silently
approximating; a gauge control (permuted candidate<->basis-index label) passes.

> **One test result must not be misread, and it is labelled in the source.**
> `t_matched_control_reproduces_the_quantum_arm` passes to **1e-9**. That is **not** a success
> of the circuit — **it is the theorem being read back.** Any future reader who sees the
> quantum arm exactly matching the classical control should read §1, not celebrate.

`gate_set_equality` asserts **subset-hood** (the theorem) and **reports** equality separately,
for the reason in §1.4. AMBER is deliberately absent from this path: `run_both` never imports
`core.amber`.

---

## 3. DELIVERABLE 3 — HAMILTONIAN-DISAGREEMENT SAMPLING (D1)

`s24/d_hamiltonians.py` -> `results/d1_hamdisagree.json`; `s24/d1_analyze.py` ->
`results/d1_analysis.json`. Both `complete: true` on the full key set with source hash and git
commit. **n = 30 targets**, seeded, fold-stratified, **declared before the run**. Point-cloud
basis, ACHIEVABLE label. LOCK_AMBER held once, 541 s, announced on open and on release.

**Why the bias cosine is THE LABEL and RMSD is descriptive.** At n=30 the RMSD arm cannot clear
a meaningful MDE, so an RMSD label would hand back a near-guaranteed "not measured" that *reads*
as a null when it is really an underpowered instrument. The cosine is estimated per target and
is far better conditioned at this n. This is stated so nobody later reads the cosine label as a
way of avoiding the endpoint.

### 3.1 VERDICT: a NULL for D1's hypothesis, one MEASURED positive inside it, and a narrowly-scoped statement about selection

Three findings from one measurement, and the third has been **rescoped downward** since I first
wrote it — see §3.3.1.

> **SCOPE CORRECTION, 2026-09-08. L5 IS RETRACTED and this file no longer claims to corroborate
> it.** My earlier wording said §3.3 "corroborates L5". It does not, and it never could have:
> both my partition arm and my floor are scored against **the same target-specific `Dhat`**, so
> nothing in my design can distinguish *target-specific* prior transfer from *generic* prior
> transfer. Lane E built the floor that can — the mean `Dhat` over all other same-length
> targets, generic and crucially **not noisy** — and against it L5's beta gap is −0.0107 and its
> Primary-1 gap is −0.1225, both NOT MEASURED. A prior that has never seen the target reproduces
> 102% of beta. **What my partitions support is the weaker and still-true statement in §3.3.1.**

### 3.2 The L2 spec — **no partition clears both bars**

A source must be ≤ ~3.9 Å standalone **and** have bias cosine ≤ ~0.65 **and** q < q* = 1.274.
Both, not either.

    partition          RMSD   <=3.9    cos    <=0.65      q    <1.274   PASSES?
    AGREE_GOOD        3.7468    YES  0.8239      no    1.215    YES     FAILS
    LEGACY_PREFERS    3.6197    YES  0.7673      no    1.199    YES     FAILS
    AMBER_PREFERS     4.1342     no  0.5693     YES    1.520     no     FAILS
    STRONG_DISAGREE   3.9177     no  0.6847      no    1.376     no     FAILS
    AGREE_BAD         3.6045    YES  0.6060     YES    1.389     no     FAILS

> **The closest thing to a positive, stated as such: AMBER genuinely points somewhere else.**
> `AMBER_PREFERS` reaches **cos 0.5693**, comfortably inside the ≤0.65 bar, and against its own
> shared-referent floor that is **−0.2154, SE 0.0546, 1.41× MDE, fold CI [−0.297, −0.143],
> 25W/5L — MEASURED.** A genuinely different physics functional *does* select along a
> materially non-parallel direction. **It cannot be used, because the candidates it prefers are
> worse**: 4.13 Å and q = 1.520, failing the quality half by a wide margin. This is the
> project's recurring shape — *diversity is available for free and costs quality* — arriving
> once more, now from the physics side.

Against the floor on Primary 1, only `AMBER_PREFERS` (1.41×) and `STRONG_DISAGREE` (−0.0992,
1.70×) are MEASURED; `AGREE_BAD` is in the **Type-M zone (0.96×) and is not a result**;
`AGREE_GOOD` and `LEGACY_PREFERS` are NULL. **Legacy alone buys no new direction at all**
(−0.0166, 0.22× MDE) — the disagreement axis is carried almost entirely by AMBER.

### 3.3 PRIMARY 2 (D1-B) — the partitions do not escape the distogram's referent

    INCUMBENT top-75   cos_distogram 0.6841   beta 0.6071
    partition           cos    floor   excess    x MDE   verdict
    AGREE_GOOD        0.5334  0.5592  -0.0258    0.32    NULL
    LEGACY_PREFERS    0.4654  0.5632  -0.0979    1.11    TYPE-M ZONE, not a result
    AMBER_PREFERS     0.4862  0.5632  -0.0771    0.49    NULL
    STRONG_DISAGREE   0.5195  0.5658  -0.0463    0.36    NULL
    AGREE_BAD         0.5606  0.5583  +0.0023    0.04    NULL

**Every partition sits at its own random same-size floor.** Not one clears it. The partitions
are all far *below* the incumbent (−0.12 to −0.22, mostly MEASURED), but that is not evidence
of a different functional — **it is evidence of not having been selected.** The floor is what
separates those two readings, and it is why the null axis was built.

#### 3.3.1 What this does and does not support

**SUPPORTED, and it survives L5's retraction because the comparison is against a floor:**

> **Selection aligns an emitted cloud with whatever prior it is scored against.** Score
> selection lifts distogram-error alignment from the random same-size floor (~0.56) to the
> incumbent's **0.684**; any subset the score did *not* select sits back down at the floor,
> whatever physics chose it.

**NOT SUPPORTED, and never was, by this design:** that the *target-specific* content of the
prior is what transfers. My partition arm and my floor are both scored against **the same
target's `Dhat`**, so the design has no contrast that could separate target-specific from
generic transfer. Lane E's generic floor — the mean `Dhat` over other same-length targets — is
the instrument that can, and against it a target-blind prior reproduces 102% of beta.

The statement above is agnostic between the two: a target-blind prior is still "whatever prior
it is scored against". My floor-based framing is what kept this result standing when L5 fell,
and that is an argument for the null axis, not for my foresight — I built the floor because
`shared-referent-floor` is in project memory, not because I anticipated the retraction.

**A process failure of mine that this exposes.** I persisted per-partition *summary statistics*
and not the AMBER energies or the partition indices. So re-running §3.3 against Lane E's generic
floor — the obvious strengthening — requires repeating the 541 s AMBER hold, because the
partitions cannot be rebuilt without the energies. **The expensive raw quantity should have been
persisted and was not.** Any future D-lane run must write the energies, not just what was
computed from them.

Score selection
lifts distogram-error alignment from the random floor (~0.56) to the incumbent's **0.684**;
any subset that the score did *not* select sits back down at ~0.56, whatever physics chose it.

### 3.4 D1-C — the quality-matched arm, and the coordinator's pre-stated expectation confirmed

The prediction, stated before I looked: once the shipped score selects *inside* the partition,
alignment should move back **up** toward the incumbent, and that would be the mechanism
reproducing under my own control rather than a failure of the arm. It does.

    pre-filter          RMSD   vs incumbent   x MDE   W/L    cos_distog   raw->filtered
    AGREE_GOOD        3.7545     +0.5475      0.79   15/15     0.5416        +0.0082
    LEGACY_PREFERS    3.4657     +0.2587      0.60   12/18     0.5493        +0.0839
    AMBER_PREFERS     3.5472     +0.3402      0.69   12/18     0.6223        +0.1361
    STRONG_DISAGREE   3.5017     +0.2947      0.68   10/20     0.5853        +0.0659
    AGREE_BAD         3.5909     +0.3839      0.67   12/18     0.5599        -0.0008

**Every `raw->filtered` move is upward** (except the partition both energies reject, which is
flat), and `AMBER_PREFERS` moves furthest, **+0.1361**, from 0.4862 to 0.6223 — most of the way
back to the incumbent's 0.6841. Letting the distogram score choose re-installs the distogram's
error, regardless of which physics functional supplied the candidate set.

**And all five pre-filters are WORSE than doing nothing**, +0.26 to +0.55 Å, 5/5 in the harmful
direction. Individually these sit at 0.60–0.79× their own MDE at n=30, so **no single one is a
clean measured harm and I do not claim one**; the consistency of the sign across all five is the
reportable part.

> **Physics as a pre-filter ahead of the deployed scorer costs 0.26–0.55 Å and buys nothing.**

### 3.5 A clean independent reproduction, worth recording

    rho(Legacy, AMBER) on the RETRIEVAL manifold   -0.0829   median -0.0760   SE 0.0397   19/30 negative
    s20's value on the CONTINUOUS TORSION manifold -0.0886

Two different manifolds, two different instruments, the same number. **The two genuine
Hamiltonians are mildly ANTI-correlated in rank**, which is why the disagreement partition is
well populated at all — and it is also why disagreement still fails to buy anything usable.
Descriptively, `rho(Legacy, distogram) = +0.3875` while `rho(AMBER, distogram) = -0.0270`:
Legacy is partly aligned with the deployed scorer, AMBER is orthogonal to it.

### 3.6 What D1 settles

**The functional lever closes alongside the provenance lever.** Combined with §1 (the selection
rule is fixed by the algorithm), the coordinator's L3 (provenance is worth nothing once the
score selects) and `conf.py` (reweighting the same functional moves beta 0.520 -> 0.490 and the
endpoint not at all), the joint statement is:

> The selection RULE is fixed by the algorithm; the CANDIDATES are measured out; reweighting the
> functional moves a mechanism but not the endpoint; and a genuinely DIFFERENT physics
> functional either points the same way (Legacy) or points a different way at a quality that
> cannot pay for it (AMBER).

**"The prior is the binding constraint" is NOT asserted here.** I wrote that sentence in an
earlier draft and it is withdrawn as premature. After L5's retraction that claim rests on the
oracle-distance counterfactual (ORACLE distances give 0.36 Å pool / 0.98 Å selected) and on
Lane B's B-3 attribution ladder, which is still running. **It is the coordinator's to settle
centrally, and nothing in Lane D establishes it.** What Lane D establishes is where the
constraint is *not*: not in the selector, not in the provenance, not in the choice between the
two available physics functionals.

**My registered pessimistic prior is confirmed**, though not exactly: I predicted cosines in the
0.85–0.95 band, i.e. partitions differing in quality but not direction. `AMBER_PREFERS` at
0.5693 is a *genuine* directional difference and my prediction was wrong about that specific
number. The conclusion I registered — that this would not pay — is right, but for the reason I
listed second (quality) rather than the reason I listed first (parallel bias).

---

## 3B. D2 — AN AMBER-INFORMED **SCORE**. THE LAST LEVER, AND IT IS CLOSED.

Pre-registered in `s24/PREREG_D.md` D2; authorised by the coordinator; run as registered.
`s24/d2_amberscore.py` -> `results/d2_amberscore.json` (`complete: true`, provenance-stamped)
and the persisted energy cache `s24/cache_amber/` (126 files). **n = 126, the full dev
instrument. THE LABEL is dev-set mean RMSD**, flipped from D1's cosine because a score change
acts on all 126 targets and the full instrument has the power the n=30 AMBER panel did not.

    s(w) = zrank(s_distogram) + w * zrank(E_AMBER)          w = 0 IS the incumbent, exactly

**Validation that the parameterisation is honest:** at `w = 0` the arm returns **3.0483 Å**,
the incumbent to four decimals, because `zrank` is monotone and the ordering is therefore
bit-identical.

### 3B.1 PRIMARY — nested LOFO CV: a clean NULL

    blended - incumbent   +0.0041   SE 0.0195   MDE 0.0547   0.08x MDE
    CI iid  [-0.0359, +0.0417]      CI fold [-0.0334, +0.0564]      53W / 69L / 4 ties
    median +0.0036      worst-target degradation +0.7557      folds same sign 2 of 5
    VERDICT: NULL

At **0.08× its own MDE**, with the fold CI spanning zero, W/L slightly *losing*, and the
per-fold deltas disagreeing in sign (2 of 5), this is an absence of signal rather than an
underpowered non-result.

### 3B.2 The part I must not gloss: **nested CV did NOT pick w = 0**

My pre-registration said *"if nested CV picks w = 0, that is the result and stop."* **It did
not.** It picked **w = −0.5 on four folds and −0.03 on the fifth; 0 of 5 folds chose zero.**

**The selected weight is NEGATIVE — the fitted blend ANTI-weights AMBER**, preferring candidates
AMBER ranks *badly*. That is not noise-shaped: it is consistent with two independent D1 results,
`rho(AMBER, distogram) = -0.0270` and the `AGREE_BAD` partition (both energies reject it) coming
back at a *better* RMSD (3.6045) than `AMBER_PREFERS` (4.1342). **On this candidate manifold
AMBER's ordering is mildly anti-informative for accuracy.**

It still buys nothing. The whole `w` curve is flat to within a few hundredths on the negative
side (−0.0006 to −0.0148) and monotonically harmful on the positive side (+0.0026 at w=+0.015
rising to +0.0809 at w=+1.0). The negative-`w` improvement does not survive nested CV.

### 3B.3 The ceiling, which closes the lever by its own upper bound

    best single w on the FULL dev set (FULL LEAKAGE)   w = -0.500   3.0336   -0.0148 A
    nested CV recovers                                              3.0524   +0.0041 A

> **Even an ORACLE `w`, chosen with complete leakage on the very targets it is scored on, is
> worth 0.0148 Å.** The lever is closed not merely by a failed fit but by its own upper bound:
> there is no value of `w` anywhere on the grid worth having. Nested CV recovers none of the
> 0.0148 and the gap is the optimism.

### 3B.4 The three nulls

    w = 0 exact               3.0483 -- reproduces the incumbent bit-for-bit
    sign-flipped w            +0.0252  SE 0.0229  0.39x MDE  59W/66L   UNDERPOWERED, not a result
    rank-permuted AMBER       -0.0010  SE 0.0045  0.08x MDE  64W/53L   NULL

The **rank-permuted** control is the informative one: preserving AMBER's marginal exactly while
destroying its correspondence to candidates gives −0.0010 at 0.08× MDE. **Whatever tiny movement
the real arm shows is not distinguishable from a control that has had all of AMBER's information
removed.** Its SE (0.0045) is four times tighter than the real arm's (0.0195), which is itself
diagnostic: `w = -0.5` is a large perturbation that swings individual targets by up to ±0.83 Å
and averages to nothing (|Δ| > 0.1 Å on 44 of 126 targets).

The **sign-flip** control is reported but not leaned on: mean +0.0252 against a **median of
+0.0004** — the project's median-vs-mean warning firing, a handful of targets carrying the whole
mean — and at 0.39× MDE it is not a result.

### 3B.5 Verdict

**D2 is a NULL and my registered prior is confirmed**, this time including its mechanism: I
predicted a null at small |w|, harm at large |w|, and nested CV selecting w ≈ 0. Two of three
correct. The one I got wrong — nested CV chose w = −0.5, not 0 — makes the finding *stronger*
rather than weaker, because the direction it chose is the anti-AMBER one and its ceiling is
still 0.0148 Å.

> **The orthogonality of AMBER to the distogram (`rho = -0.0270`) is real, is the only
> structural reason anyone found to expect a different functional to help, and does not convert
> into accuracy through the score any more than it did through a filter or a partition. The
> functional lever is closed in all three of its available forms.**

### 3B.6 The durable asset

`s24/cache_amber/` holds **126 × 500 = 63,000 genuine ff14SB/GBn2 single points**, with the
distogram score and pool indices alongside, verified bit-exact (`amber_verify_max_rel = 0.0` on
every target). Any future AMBER question on this pool — including strengthening §3.3 against
Lane E's generic floor — can now be answered **without touching OpenMM**. The coordinator judged
this sufficient justification for the hold independently of the arm's outcome, and given the
result that judgement is the one that carries.

---

## 4. DELIVERABLE 4 — AMBER'S ROLE, STRICTLY BOUNDED

**No refinement arm was run and none is proposed.** s23 L11 measured 17 of 17 restrained and
unrestrained repair settings at or worse than no repair, with the ladder bottoming out at the
no-op. AMBER appears in this workstream **only** as:

* an **energy measurement** — genuine ff14SB/GBn2 **single points**, `steps=-1`,
  `k_restraint=0`, via `s20.qb2_lib.AmberSP`, whose `verify()` asserts bit-exactness against
  `core.amber.refine_coords` and returned max relative difference **0.0** on every target;
* a **constrained diagnostic** — the conditioning census in §3 and the landscape table in §5.

It is never an unconstrained refinement, never a default pipeline stage, and never a selector.

**A conditioning fact that any future AMBER consumer needs.** Raw ff14SB/GBn2 single points on
**unrelaxed retrieval windows** are clash-dominated. Across the 30-target panel a large share of
every pool sits above 1e4 kcal/mol and single candidates reach 1e18, so one clash sets the
standard deviation and **98-99.7% of every pool lands inside |z_raw| < 0.1**. A moment z-score
of raw AMBER therefore measures *which candidate has the worst steric clash*, not *which
candidate AMBER prefers*. This is `pauli-spectrum-delta-spike-artefact` reproduced on a new
instrument: **condition monotonically; 99th-percentile winsorisation is not enough.** Rank
standardisation fixes it (6.0% inside |z| < 0.1) and is strictly monotone, so it changes no
ordering, no argmin and no level set. **It is not cosmetic:** the `AMBER_PREFERS` partition
under the two normalisations overlaps by only 0.56-0.68.

---

## 5. DELIVERABLE 5 — THE LEGACY-vs-AMBER LANDSCAPE TABLE

**Reused and cited, not re-measured**, per the brief. Sprint 20 measured this on the continuous
torsion instrument at n=30 and its numbers stand; a changed candidate manifold cannot plausibly
move a Hessian spectrum or a gradient cosine, so nothing here was re-run. Only the rows marked
**[S24]** are new, and they are new because they are properties of the *candidate pool*, which
is exactly what Sprint 24 changes.

| axis | Legacy | AMBER | source |
|---|---|---|---|
| energy scale | ~1.6 decades of dynamic range | **16.39 decades** raw; 1.57 after log conditioning | s20 L (conditioning) |
| gradient scale | move 5.6x smaller | first steps are clash relief from `E_0 = +1.9e8 kcal/mol` | s20 |
| gradient direction (joint) | **cos(grad Legacy, grad AMBER) = -0.3512 [-0.4719, -0.2221]** | 27/30 targets NEGATIVE | s20 |
| Hessian condition | baseline | **three orders larger**; log-conditioned AMBER still **473x worse** | s20 |
| anisotropy | 5.8 | **18.8**, 0W/30L | s20 |
| participation ratio | 0.4221 | **0.0745**, 30W/0L — ~7% of modes carry all the curvature | s20 |
| near-zero spectrum | baseline | **44%** near-zero modes; 4x Legacy after conditioning | s20 |
| ruggedness | baseline | **2x** the local minima per 2π, unchanged by conditioning | s20 |
| negative curvature | baseline | **slightly LESS** than Legacy, -0.0346 [-0.0707,+0.0116] — s20's own prediction REFUTED | s20 |
| compactness response | **prefers candidates 0.45 Å MORE COMPACT**, -0.4477 [-0.5110, -0.3911], 124W/2L, 5/5 folds | prefers expanded, open sterics | s20 |
| steric sensitivity | steric term carries the largest weight (4.0) yet behaves as compactness/typicality | **curvature concentration at a steric singularity**; 53.5% of the pool in hard overlap | s20 |
| rank agreement (joint) | ρ = **-0.0886** torsion instrument | ρ = **-0.0829** retrieval pool **[S24]**, independent reproduction | s20, §3.5 |
| structural correlation (joint) | complementary about GEOMETRY | **jointly blind to ACCURACY**: rebuild distance -0.0785 [-0.2120,+0.0464], pool-error alignment -0.0328 [-0.1526,+0.1010] | s20 |
| CVaR behaviour | well conditioned | CVaR on raw AMBER is a **crude discontinuous conditioner**, ratio 1.7 at every α | s20 |
| **[S24] conditioning on the retrieval pool** | well conditioned | clash-dominated: 98-99.7% of the pool inside \|z_raw\|<0.1 | §4 |
| **[S24] alignment with distogram error** | filtered 0.5493 | filtered **0.6223** (incumbent 0.6841) — neither escapes | §3.3-3.4 |
| **role in the architecture** | a compactness/typicality prior in physics vocabulary | a feasibility FILTER and an energy MEASUREMENT | §1, §3 |

**Two standing corrections carried forward, so they are not re-derived wrongly.**
*"AMBER is less local"* is a **category error** — both energies are **full-register** in torsion
space (`torsion-space-locality-theorem`; `d_ij` depends on exactly the residues between `i` and
`j`, agreement 1.0000, and all-atom supports are one residue wider, but both are full-register).
And *"Legacy is the less physical model"* inverts the measured finding: **Legacy is a
compactness/typicality model wearing a physics vocabulary**, which is the opposite of the naive
prior given that its steric term carries the largest weight.

---

## 6. WHAT DAMAGED MY OWN EXPECTATIONS

Reported plainly, per the brief.

1. **My REV1 assertion was dead and reported a pass** (§1.4). An epsilon floor in every
   probability family meant the test could not fail. I found it only after the coordinator
   challenged the wording — the coordinator caught the claim, I caught why the test agreed
   with it. A green light from an assertion that cannot fire is worse than a red one.
2. **My pre-registered normalisation was unusable** (§4) and I did not anticipate it, despite
   `pauli-spectrum-delta-spike-artefact` recording exactly this failure. I caught it from a
   native-free diagnostic before computing any outcome, but I should have expected it from
   memory rather than discovered it from data.
3. **The smoke test exposed a confound in my own design** (§3, D1-C): the raw partitions are
   not score-selected, so comparing them to the incumbent conflates *"a different functional"*
   with *"no selection at all"* — a confound worth 1.4-2.4 Å. Had I reported the raw partitions
   alone I would have published a selection artefact as a mechanism.
4. **My registered prior was right for the wrong reason.** I predicted cosines of 0.85-0.95 —
   partitions differing in quality but not direction. `AMBER_PREFERS` came back at **0.5693**, a
   genuine directional difference and a MEASURED one against its floor. My headline prediction
   was wrong. The conclusion survived only because of the second condition I listed (quality),
   not the first (parallelism). Had the L2 spec been a single bar rather than a conjunction, I
   would have called this a win.
5. **The normalisation choice was even more substantive than my smoke test suggested.** On the
   full panel the rank-vs-raw partition overlap is **0.348**, against the 0.56-0.68 I quoted
   from two targets. Two-thirds of the candidates studied depend on that one amendment.
6. **My own harness's verdict rule over-claimed, and D2's control caught it.** `paired_stats`
   labelled an effect at **0.39× its own MDE** "MEASURED" purely because a 5-fold cluster
   bootstrap CI excluded zero. With five clusters that CI is unstable, and an effect below its
   own MDE is by definition one the design could not reliably detect. The rule now returns
   `UNDERPOWERED (<0.7x MDE)`. **This had already bitten D1**: three D1-C rows at 0.60–0.69×
   were machine-labelled MEASURED. My *prose* was right — I wrote "no single one is a clean
   measured harm and I do not claim one" — but the artefact's own labels disagreed with my
   text, and a later reader trusting the JSON over the prose would have been misled. Both are
   now consistent. A harness that flatters its user in exactly the direction the user wants is
   the same class of error as §1.4's dead assertion.
7. **I claimed to corroborate L5, and had to withdraw it when L5 was retracted** (§3.1, §3.3.1).
   My design could never have supported that claim: both my arm and my floor are scored against
   the same target-specific `Dhat`, so nothing in it separates target-specific from generic
   prior transfer. The measurement stands; the *interpretation* I attached to it did not, and I
   attached it in the same message in which I criticised a different over-claim of my own.
   **And I persisted summaries rather than the AMBER energies**, so re-testing §3.3 against Lane
   E's generic floor now costs a repeat of the 541 s AMBER hold. Persist the expensive raw
   quantity, not what you computed from it.

---

## 7. WHAT REMAINS LIMITING

Not a manufactured win, per BRIEF §6. The closest defensible improvement in this workstream is
**nothing**: no arm I ran moves the endpoint, and I do not propose one. What is established is
*where* the constraint sits.

* The **selector** cannot be the lever: its selection is provably the classical top-m set (§1).
* The **candidate provenance** is not the lever (coordinator L3), and D1 shows the **functional**
  is not either (§3.6).
* Both physics functionals are, on the ORACLE axes, **jointly blind to accuracy** (s20, §5) —
  they are complementary about *geometry* and agree about *catastrophes*.
* **Selection aligns the output with whatever prior it is scored against** (§3.3.1) — including,
  per Lane E, a *target-blind* one. This is narrower than the retracted L5 and it is what my
  floor-based design actually supports.

**What is NOT established by this lane, and is left conditional:** that *the prior* is the
binding constraint. That claim now rests on the oracle-distance counterfactual and on Lane B's
B-3 ladder, and it is the coordinator's to settle centrally.

The honest architectural statement for the final report: **CVaR-VQE is a genuine, exactly
simulated selector whose selection is provably the classical top-m set of whatever functional it
is given; Legacy is a compactness/typicality prior wearing a physics vocabulary; AMBER is a
feasibility filter, an energy measurement, and the one source measured in this sprint that
selects along a materially non-parallel direction — at a quality that cannot pay for it. The
leverage is not in the circuit, not in the provenance, and not in the choice between the two
available physics functionals.**
