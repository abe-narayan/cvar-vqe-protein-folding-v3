# SPRINT 21 — WORKSTREAM A PRE-REGISTRATION
## The mandatory Hamiltonian matrix, with CVaR-VQE as the selector

Written before any Sprint-21 number was produced by this lane. **Falsifiers first.**
Nothing below is edited after data are seen; corrections are appended as dated addenda with the
original text left intact.

Basis vocabulary, per BRIEF §1: **point cloud** (coordinate average, not buildable) /
**built chain** (an ideal-geometry rebuild or a real retrieved window) / **repaired emission**
(after projection and/or restrained AMBER). Every row below states its basis.

MDE at 80% power = **0.084 Å**. Anything smaller is NOT MEASURED, never "matched".

---

## STAGE 0 — THE TAIL-PRICING BOUND (`s21/tailprice.py`)

The full pre-registration is in the module docstring and was written by the coordinator before
this lane touched the file; it is reproduced here so the falsifier is in one place.

**Why it comes first.** Sprint 20 measured *argmin-by-energy* over the shipped K=500 pool:
AMBER 4.990, Legacy 5.487, distogram 3.676, against a pool mean of 4.739. Picking the
lowest-energy structure is worse than picking at random for both physics energies. **But argmin
is not CVaR.** CVaR selects an α-tail and — EXACT, Sprint 20 Q6 — its minimiser is a *face*: it
constrains the tail and is indifferent outside it. A tail **mean** regresses toward the pool
where a **minimum** chases the worst-ranked outlier. These are different operators on the same
spectrum and the difference is measurable in minutes.

    Hypothesis    a tail selector is materially better than argmin on the physics Hamiltonians.
    Prediction    tail(0.15) beats argmin by > 0.3 A on legacy and amber; the gap shrinks as
                  alpha -> 0 and vanishes on disto, which already ranks well.
    Primary       tail(0.15) - argmin per Hamiltonian, paired over 126 targets, bootstrap CI.
    FALSIFIER     if a Hamiltonian's tail arms do NOT beat their MATCHED-COUNT RANDOM subsets,
                  tail selection carries no information on that Hamiltonian, and a CVaR-VQE whose
                  selection is over this pool is bounded at chance on it no matter how well it
                  optimises.  That Hamiltonian is then priced OUT of the matrix as a selector.
    Null          random{m} at every matched tail size.  THE operative control.
    Budget        n = 126, no optimisation, pool-restricted.
    Promotion     none.  This is a bound, not a method.

**THE TRAP, pre-named.** A tail arm that beats argmin has proved nothing: as α → 1 the tail
becomes the whole-pool average, a zero-information operator already known to be decent. The
matched-count random control at the same count is the only comparison that means anything.

**BASIS (declared).** Selection ranks pool members; the emission is the **coordinate average of
the selected members' real retrieved Cα windows** — the deployed operator, on the point-cloud
basis. A **rebuild-basis** twin (coordinate average of each member's ideal-geometry rebuild) is
recorded on every arm, because Legacy and AMBER are only defined on the rebuild and it would be a
confound to score them through a geometry they never saw.

**Verified before the run, not assumed:** the shipped top-75 filter is the distogram score on the
**real windows W**, not on the rebuild — 75/75 index agreement on 1A13, 54/75 on the rebuild. The
disto arm therefore scores W; Legacy and AMBER score the rebuild. This asymmetry is declared, not
discovered.

---

## STAGE 1 — NORMALISATION, DECLARED BEFORE RMSD

Legacy and AMBER live on wildly different scales and AMBER has a heavy upper tail (Sprint 20 L5c:
median +16,062, p99 3.1e13, max 5.5e23 on 9,450 real rebuilds; 53.5% above 1e4). A raw sum is an
AMBER-outlier detector wearing a hybrid's name.

**Declared form: within-target rank-to-normal.** Each component is rank-transformed to a standard
normal *within that target's own candidate set* before any sum. Properties, stated as the
justification and not as a result: scale-free, invariant to any monotone reparameterisation of
either energy (so it cannot be gamed by a barrier transform), outlier-robust by construction,
carries no fitted constant, and needs no per-target calibration data.

**Cost, stated up front:** it discards magnitude, so it cannot express "AMBER says this one is a
catastrophe by 20 orders of magnitude". That is deliberate — Sprint 20 L5c shows the magnitude
channel is a steric singularity, not an accuracy signal.

**Audit arms, computed alongside and reported whatever they say:** raw-sum combinations, and the
per-target z-score. If a raw arm beats the rank arm the declared choice is recorded as *wrong*,
not quietly swapped.

**Measured before combining** (Stage 1 deliverable): mean, variance, skew, kurtosis, the empirical
quantile ladder, the finite fraction, and the fraction above 1e4, for each energy, per target and
pooled.

---

## STAGE 2 — THE MANDATORY MATRIX

    Legacy | AMBER | Legacy+AMBER | Distance | Distance+Legacy | Distance+AMBER | Distance+Legacy+AMBER

Same candidate problem, same representation, same budget, same seeds, **genuine CVaR-VQE as the
selector**, all seven Hamiltonians separably evaluable at all times.

    Primary endpoint   final full-chain Ca-RMSD, MEAN and MEDIAN, of the emitted structure,
                       stated on a named basis, n = 126 tuning targets.
    Reported in order  mean, median, target-level paired difference vs the declared reference,
                       bootstrap CI (i.i.d. AND fold-clustered, both quoted), W/L, per-fold means
                       -- BEFORE objective convergence, gradient/Hessian geometry, CVaR behaviour,
                       diversity or the validity vector.
    Reference arm      the structural-only baseline (Distance), which is the shipped selector.
    Seeds              4 minimum on every variational arm (Sprint 20: ansatz-seed sd 0.200 A
                       within-target = 2.4x MDE; fewer than 4 is NOT MEASURED).

**FALSIFIER F-A1.** If no Hamiltonian containing Legacy or AMBER beats the Distance-only arm by
more than the MDE with a CI excluding zero, **the physics Hamiltonians contribute nothing to
selection**, and that is this lane's result.

**FALSIFIER F-A2.** If the CVaR-VQE arm does not beat **best-of-N from the untrained circuit** at
matched N on the same Hamiltonian, the variational machinery contributes nothing over its own
initialisation and the arm is reported as NOT MEASURED or NOT SUPPORTED accordingly. An
initialisation *mean* is not this control and will not be used as one.

**FALSIFIER F-A3.** If a matched-budget classical selector (greedy, simulated annealing, exact
enumeration where the tail is small enough, and a matched-count random draw) reaches the same
RMSD, the quantum arm is reported as *reproduced by a simpler classical baseline* — kill rule,
BRIEF §12.

**Pre-committed kill rules**, any one fires: loses to the matched random control · loses to a
plausible zero-information null · CI spans the MDE · dissolves under target-level analysis ·
requires native information · rests on an arbitrary encoding that fails a gauge check ·
improves the objective while worsening RMSD.

---

## STAGE 3 — CONTROLS THAT MUST ACCOMPANY EVERY QUANTUM CLAIM

1. best-of-N from the **untrained** circuit, matched N — never an initialisation mean.
2. matched-count **random** subset at the same tail size.
3. a **zero-information** reference in the operator's own space (whole-pool average; matched
   empirical marginals), per BRIEF §7.4 — uniform is not one.
4. a matched-budget **classical** optimiser over the identical Hamiltonian.
5. the **shared-referent floor** (0.505) measured before any alignment statistic is interpreted.

Gate reporting per BRIEF §7.3: every gate, guard and cap states **how many times it FIRED**. A
gate that never fired passed vacuously and is reported as such.

---

## STAGE 4 — WHAT WOULD MAKE ME ABANDON THE LANE

If Stage 0 shows that neither physics tail beats matched random, then §4's matrix is, for six of
its seven cells, a measurement of how much damage a zero-skill selector does to a good one. That
is still worth reporting — it is the price of the mandatory matrix — but it is reported as a
**bound**, not as a search for a winner, and no amount of variational effort is offered as a
remedy. **A clean falsification is a successful result.**

---

## ADDENDA

*(appended after data; original text above unedited)*

**ADD-1 — 2026-09-07, coordinator-issued, logged before any table was read.** Stage 0's primary
`tail(0.15) − argmin` is **confounded**: `argmin` returns one structure, `tail{α}` returns the
coordinate average of α·m structures, so the contrast contains the averaging operator (~1.0 Å,
Sprint 19) and would have "confirmed" the >0.3 Å prediction on *every* Hamiltonian including a
random one. **New primary: `tail{α} − random{α}` at matched count** — already present as the
operative control. Three single-structure readouts (`tail_min`, `tail_member`, `tail_medoid`) are
added beside the averaged one, each with its own matched-count random control, so the tail
operator and the averaging operator are separable. The old primary is retained and relabelled.

**ADD-2 — 2026-09-07, coordinator-issued.** A Hamiltonian's sign is **a function of the readout**:
Legacy beats its matched control at `member` (−0.41 [−0.56, −0.26], n=126) and loses at `medoid`
(+0.30) and `avg` (+0.33). **Every row must state its readout**, and "which Hamiltonian is best"
is reported per readout, never unqualified. Mechanism is **error coherence**, not diversity — the
diversity account was refuted by its own sign control (maximising diversity is dead; *minimising*
it helps, −0.214 [−0.397, −0.040]), so **no diversity term enters any arm**.

**ADD-3 — 2026-09-07, coordinator-issued.** **MDE is per comparison**, `2.8016 × SE`. The pooled
0.084 Å constant in the header of this document is **retracted as a threshold** — it is wrong by
up to 84× in both directions on an individual contrast. Every reported contrast now carries its
own SE and its own MDE. *The original text above is left unedited; this addendum supersedes its
MDE line.*

**ADD-4 — 2026-09-07, my own error, recorded.** I reported Legacy's **marginal** rank skill
(+0.328) as its contribution. It is not: Legacy is strongly rank-correlated with the deployed
distogram (+0.478), and the **partial given the distogram is −0.0076 fold[−0.0616, +0.0681], 25/42
negative**. Every Hamiltonian's contribution is reported as a **partial**, never a marginal. The
same artefact that produced the wrong reading already contained the refutation.

**ADD-5 — compute, declared.** Stage 2's target count was reduced from 20 to **12** and Stage 0
was paused at **42 of 126**, because the shared box measured **52 ms per AMBER single point
against a 6 ms nominal** under four-way contention and the mandatory matrix is the deliverable.
**This is a compute cut, not a scientific one**: per-cell configuration is intact, only power is
reduced, and the power warning prints on every table. The untested regime — Stage 0 at the full
n=126 — is recorded as **OPEN**, and `tailprice.json` remains resumable.
