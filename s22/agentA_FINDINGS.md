# WORKSTREAM A — Sprint 22 findings (CVaR-VQE selector)

Pre-registered in `s22/PREREG_A.md` (A1/A2/A3), with two dated addenda appended there (never
edited into the frozen body): ADDENDUM 1 replaces A2's β-weighted structural term with an
entropy-regularised free energy after a pilot showed pure CVaR-over-p collapses to a single
candidate regardless of α; ADDENDUM 2 records a tie-breaking hazard caught and fixed before any
A1/A2 number was read. All work is complete: **A1, A2a, A2b, A3, and the readout-H/training-H
separation experiment have all run to completion.** No AMBER/OpenMM anywhere in this workstream.
Benchmark seal verified unchanged at start:
`a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d` (8002 bytes, byte-identical to
s21's certificate).

**Coordinator's independent read of `a2b_entropy.json` and `a_readouth.json` (ledger L13/L14) is
checked below against my own analysis and found consistent** — the apparent numerical
differences are two valid ways of slicing the same data (pooled-over-entangler vs split;
fixed-m=75 bar vs realised-m-matched bar), reconciled in §3 and §4. No error found in either
direction. What follows adds a mechanism the coordinator's numbers are also consistent with but
did not have from the raw artefact alone: **both the entropy result and the readout-H result
reduce to a single provable architectural identity**, stated in §0 and proved in §3–§4.

---

## 0. THE THROUGH-LINE, stated once so the rest of this file is read correctly

**For an order-based readout (argmin, or a mass-support tail-average) over an exactly-enumerable
diagonal candidate Hamiltonian, the trained circuit can never select a different SET of
candidates than a classical score-sort would. Training, entanglement, temperature and
cross-Hamiltonian staging can only ever change the SIZE of the selected set, never its identity.**

This is not an empirical tendency; it follows directly from the definition of
`core.quantum.cvar_from_probs` (used throughout this workstream for the exact face computation):
it sorts states by ENERGY, not by the trained probability, and assigns positive mass to a PREFIX
of that energy order until cumulative probability reaches α. So `tail = {x : mass(x) > 0}` is
always exactly "the m lowest-energy states" for whatever `m` the trained distribution's
concentration produces — regardless of the ansatz, the entangler, the temperature, or which
Hamiltonian trained it. An UNWEIGHTED coordinate average over that set (the readout this
programme uses throughout, matching the shipped operator's own convention) is therefore
**provably identical** to a classical top-m average at that same m. **Verified by direct set
comparison, not inferred**: `n_tail_cands` candidates recovered from a trained circuit are
element-for-element identical to `argsort(scores)[:m]` in every spot-check run (18 cells across
7 targets, §3/§4), and the argmin readout is separately verified exactly equal to the classical
argmin of whichever Hamiltonian scores it, in **all 64 argmin cells** of the readout-H experiment
(max|diff| = 0.0, four training/readout combinations × 16 targets).

**Consequence for how to read §2–§4.** A1's closure (training is irrelevant to the argmin
readout) and this identity (training is irrelevant to the tail-average readout's MEMBERSHIP) are
the SAME fact applied to two readouts. Every subsequent "effect" measured in this workstream — the
entropy-T sweep (§3), the readout-Hamiltonian swap (§4), and A3's staging result (§5) — is
therefore a statement about **which classical set size m gets realised**, mapped onto the
already-fully-characterised classical m-ladder (s21 G1–G5: m* interior, flat for m∈[20,150],
m=75 already correct). None of it is evidence of a distinct quantum selection mechanism. That
does not make these results uninteresting — a provable, general architectural closure is a
**better** contribution than an unexplained empirical one, and it is stated plainly here rather
than left implicit.

---

## 1. SELF-CAUGHT HAZARD, reported first

**The tie-breaking trap, reproduced in a new substrate.** A pilot run found target `5H1H`'s
top-two distogram scores exactly tied (gap 0.0). A naive `np.argmin` broke the tie by array
position, which differs between the classical control's index space (candidate-canonical order)
and the trained circuit's (label-permuted bit-index order) — manufacturing an apparent gauge
failure that was really tie arbitration. Gate 1 fired 10/20 times on that one target alone before
the fix, and 0/300 everywhere else. **Fixed**: both readouts now return the full tied candidate
SET and average ORACLE RMSD over it (gauge-invariant by construction, since it is a set of
candidate indices). The buggy pilot is preserved, not deleted, at
`s22/results/a1_gauge_PILOT_tiebug.json`.

---

## 2. A1 — CANDIDATE-STATE ENCODING: THE A-PRIORI CLOSURE

**`H|i⟩ = E_i|i⟩`, `n_qubits = ⌈log₂K⌉` (K=500 → 9 qubits, 512 states, 12 padding). Exactly
enumerable** — every quantity below is computed by full enumeration, never sampled.

### GATE 1 — soundness (exact global-optimum SET agreement): **0/320 fired**, corrected run.
(Pilot, tie-naive: 10/320, entirely the one tied target — see §1.)

### GATE 2 — achieved gauge robustness: **argmin 0/256 fired (provably always zero — see the
theorem below); tail-average/CVaR-value ~26% fired** against the canonical labelling's own
seed-noise spread. Magnitude, computed directly (excluding the tie-corrupted target): tail-average
RMSD deviation under relabelling has **mean −0.017 Å, sd 0.27 Å, median |deviation| 0.028 Å**,
with **42.6% of (target, permutation) cells exceeding 0.05 Å, 22.7% exceeding 0.2 Å, and 5.9%
exceeding 0.5 Å.**

**The honest asymmetry, worth a sentence on its own, as instructed.** The readout that is
*provably gauge-invariant* (argmin) is also the one that is *useless* — it never differs from an
untrained circuit (below). The readout that shows real, if moderate, gauge sensitivity
(tail-average) is the only one that carries any information beyond a single point. A candidate
encoding cannot have a readout that is simultaneously informative and free.

### THE A-PRIORI CLOSURE (stated as a theorem, as requested)

**`R_VQE_argmin = R_untrained_argmin = R_score_argmin` exactly, on all 16 targets** (mean
difference 0.0, all three quantities bit-identical per target). This is forced: a full-support
ansatz (any generic RY-based circuit, trained or not) assigns strictly positive probability to
every one of the `2^n` basis states, and `argmin_readout` scans the full support — so for any
register small enough to enumerate exactly, the argmin readout is the true global optimum
*independent of whether any training happened at all*.

> **THEOREM.** For a pool of size K ≤ 2^30 (i.e. every pool this project uses, K ≤ 500 → ≤ 9
> qubits, always exactly enumerable at the scale in question), a genuine CVaR-VQE's argmin
> readout over the candidate-identity register is identical to a zero-training circuit's argmin,
> and identical to the classical exact argmin. **Any claim of quantum argmin-selection advantage
> over a discrete pool this size is closed a priori, independent of ansatz, entangler, optimiser,
> or training duration.** This is a mathematical fact about the architecture class (full support +
> exact enumeration + order-based readout), not an empirical tendency, and it is expected to
> outlive the specific numbers in this file.

### RMSD decomposition, mean over 16 targets, point-cloud basis

| quantity | value | note |
|---|---|---|
| R_pool (ORACLE ceiling) | 1.780 | never achievable, consumes the native |
| R_score_argmin (classical exact) | 3.892 | = R_VQE_argmin = R_untrained_argmin, exactly |
| R_VQE_argmin (trained, canonical label, any T) | 3.892 | see theorem above |
| R_untrained_argmin (best_of_N, untrained circuit) | 3.892 | mandatory control, per BRIEF §7 |
| R_score_tailavg (classical exact top-α avg, α=0.15) | 3.300 | the incumbent-analogue for this pool |
| R_VQE_tailavg (trained, canonical, **T=0**) | 3.819 | **significantly worse** than R_score_tailavg |
| R_rand_tailavg (matched-count random) | 3.833 | zero-information null |
| R_untrained_tailavg (best_of_N, untrained circuit) | 3.345 | *below* R_VQE_tailavg(T=0) — see §3 |

**VQE_tailavg(T=0) vs classical exact tailavg**: mean **+0.519 [+0.109, +0.886]**, MDE 0.576,
13W/3L — significantly *worse*. **VQE_tailavg(T=0) vs matched random**: mean −0.015
[−0.734, +0.627], NOT MEASURED — statistically tied with random. **Classical exact tailavg vs
matched random**: mean −0.534 [−1.112, −0.018], marginally excludes zero — the plain classical
sort carries real information a random subset does not (replicating the project's most reproduced
fact, "the distogram orders in-band," on a tenth independent instrument).

**Matched controls delivered per BRIEF's hard requirements**: random (above), simulated annealing
(found the true global optimum in 69% of targets at budget = K = 500 evaluations — weaker than
exhaustive sort or the exact-enumeration argmin, which both get it 100% of the time by
construction; a genuine, if unsurprising, search-completeness gap for a stochastic method at
matched budget), best_of_N from the untrained circuit (above).

---

## 3. A2 — CVaR DEGENERACY: CLOSES NEGATIVELY, AND THE CONTROL IS WHAT MAKES IT CLEAN

### 2a — collapse dynamics (T=0, all six α schedules, checkpointed at 5/20/50/100/150/200 iters)

**Every α setting — including α=1 (the undegenerate full mean) — collapses to ESS≈1.0–1.1 and
tail_size≈1–2 by iteration 200**, and tail-average RMSD gets monotonically *worse* with more
training at every α (e.g. α=0.05: 3.635→3.849; α=0.15: 3.491→3.819; α=0.5: 3.409→3.611).

**Why α=1 also collapses, stated as the underlying fact rather than left as a puzzle**: for FIXED
energies E, `CVaR_α(p)` restricted to any candidate tail set is a LINEAR functional of `p` (a
probability-weighted sum of that set's energies). Minimising a linear functional over the
probability simplex is optimised at a vertex — a delta — for *every* α, including α=1 (the mean).
Only a genuinely non-linear regulariser (entropy) can produce an interior optimum. This is the
formal version of `core.quantum.free_energy`'s own stated empirical fact, and it is why ADDENDUM 1
replaced the β-structural-term design: there was no multi-candidate "face" to break at T=0 in the
first place, for continuous non-tied energies.

### 2b — entropy-regularised degeneracy breaking: the coordinator's numbers, checked and adopted

Pooled over both entanglers (cnot, none — I confirm this reproduces the coordinator's L13 table to
the numbers quoted, e.g. T=0.5 VQE-tail 3.271 vs my own cnot/none split of 3.264/3.278 mean-pooled
= 3.271; ESS 62.08 = mean of my cnot 52.63 and none 71.54; all figures match), against the FIXED
classical top-α=0.15 (m≈75) bar (**3.2998**, = A1's R_score_tailavg = the incumbent-analogue):

| T | VQE tail-avg | ESS | tail size | vs FIXED bar (m=75) | vs matched-random |
|---|---|---|---|---|---|
| 0.00 | 3.723 | 1.02 | 2.1 | +0.423, MDE 0.545 | −0.976 |
| 0.05 | 3.556 | 5.50 | 10.9 | +0.256, MDE 0.496 | −0.746 |
| 0.20 | 3.313 | 33.70 | 44.5 | +0.013, MDE 0.220 | −0.609 |
| 0.50 | 3.271 | 62.08 | 68.5 | **−0.029, MDE 0.105** | −0.581 |

**This is the null, and it is a clean one.** The entropy regulariser works exactly as designed
(ESS climbs 1.02→62.08, tail widens 2.1→68.5 — the degeneracy lever is real and controllable) and
beats a size-matched random draw at every T (0.58–0.98 Å) — but it **never beats the fixed
classical top-m=75 bar**, and its best cell (T=0.5) sits at −0.029 against its own MDE of 0.105,
decisively inside the null.

**My own, complementary framing (the size-MATCHED classical control I pre-registered), and why
it is not a second, different result but the mechanism for the one above.** Comparing VQE
tail-average to classical top-m at the **SAME REALISED m** (not the fixed m=75) gives **exact
equality to floating-point precision** (differences ~1e-15) at every T and every target — because,
per §0's theorem, the tail-support SET is provably always the classical top-m set for whatever m
training realises. Verified by direct set comparison (not just RMSD equality): the 69
tail-candidate indices recovered from a T=0.5-trained circuit on `1CS9` are element-for-element
identical to `argsort(scores)[:69]`, spot-checked on 12 further (target, seed) cells across 6
targets with zero mismatches.

**So the two framings say the same thing at two different resolutions.** The fixed-bar table
shows *how close* entropy regularisation gets to the incumbent's chosen m=75 as T grows (it
converges toward null from above, because the classical m-ladder is nearly flat for m∈[20,150] and
T=0.5's realised m≈68.5 lands inside that flat region). The matched-m identity explains *why*:
there is no mechanism by which the trained tail could ever be anything other than classical
top-m, so "converges to the fixed bar" is just "T pushes the realised m into the region where the
already-known-flat classical ladder no longer distinguishes it." **Both readings close item 2
negatively, and for the same underlying reason.**

**Disposition: A2 CLOSES NEGATIVELY.** Not "no native-free arm beats the incumbent" (G2's
finding, restated) but something stronger and cleaner: **an entropy-regularised CVaR-VQE's tail
readout cannot beat classical top-m even in principle**, because it is provably the same set. The
extension to a tenth (quantum-flavoured) entry in G2/G2a's list is exact, not approximate.

---

## 4. READOUT-H / TRAINING-H SEPARATION: replicates on a second substrate, and reduces to a
   classical fact once traced to its mechanism

Trained once per (target, seed) on `H_TRAIN ∈ {DIST, LEG}`, T=0, α=0.15, 150 iterations; then read
out the SAME trained distribution against **both** energies with zero further training (mirrors
s21 L16/B6's "at zero extra budget on the identical evaluated set").

| training \ readout | readout = DIST | readout = LEG |
|---|---|---|
| **DIST** | 3.819 (same-H) | 4.342 (cross-H) |
| **LEG** | 3.417 (cross-H) | 5.346 (same-H) |

**Coordinator's pooled reading (L14), checked and confirmed**: readout swap DIST-vs-LEG,
tail-average **−1.177 [−2.007,−0.410]**, SE 0.420, MDE 1.177; argmin **−1.508 [−2.820,−0.312]**,
SE 0.668, MDE 1.873. I reproduce this exactly as a pooling of my own two separately-registered
comparisons (holding training fixed at LEG: DIST-readout−LEG-readout = −1.929
[−3.106,−0.796], 14W/2L, independently significant; holding training fixed at DIST:
LEG-readout−DIST-readout = +0.523 [−0.308,+1.394], 6W/10L, same sign, not independently
significant) — averaged per target before bootstrapping, which is why the pooled W/L (11/5) sits
between my two separate tallies. **No discrepancy found; both directions agree in sign, one
independently significant, the pooled estimate clears its own CI.** Per BRIEF's Type-M rule
(|effect|/MDE ≈ 1.0 for both quantities here): **direction ESTABLISHED, magnitude NOT MEASURED**
— replicates s21's own disposition on its own matrix (B6: −0.697, B11: −1.299 on the torsion
basin latent) for the same statistical reason.

**The mechanism, verified by direct proof rather than left as "an architectural lever."** Both
numbers above are **exactly** reproduced by two zero-training classical facts, checked directly:

- **Argmin swap is exactly the classical argmin-under-the-readout-Hamiltonian, always.** Verified
  in all 64 cells (4 training/readout combinations × 16 targets): max|difference| = 0.0 between
  the trained circuit's argmin readout and the classical exact argmin of whichever Hamiltonian
  scores it — *regardless of which Hamiltonian trained the circuit*. The −1.508 Å argmin effect is
  therefore not a training or entanglement effect at all: it is the pre-existing classical fact
  that DIST's global-best pool candidate is a better structure than Legacy's, visible through a
  readout that (§2's theorem) training can never move away from either answer.
- **Tail-average swap is exactly classical top-m under the readout Hamiltonian, at whatever m the
  cross-trained circuit's concentration happens to realise there.** Verified directly: a
  Legacy-trained circuit, read out against DIST's own energy order, realises m=176 real
  candidates on `1CS9` (versus DIST's own T=0 self-collapse to m≈1–2) — and its tail-average
  RMSD (3.722) is bit-identical to classical top-176-by-DIST-score (3.722). The −1.177 Å
  tail-average effect is the SAME m-realisation mechanism as §3, not a second, independent one:
  training on the "wrong" Hamiltonian happens to leave the distribution far less collapsed when
  re-sorted by the "right" one, landing a much larger, better-positioned m on the classical
  ladder than DIST's own unregularised self-collapse does.

**Reframing, stated plainly because it revises how L14 should be read.** The effect is real, the
CIs are honest, and it replicates s21's sign and rough order of magnitude on a genuinely different
encoding — that convergence is worth keeping. But on THIS substrate it is not evidence of a
distinct quantum-selector mechanism the way it might be on the continuous torsion basin latent
(where the argmin/mode readout picks a basin CONFIGURATION that must still be decoded into a
structure, and nothing forces that decoding to coincide with a classical sort). Here, because
training is provably inert for both readout types (§0, §2), "the readout Hamiltonian is the
lever" is the *same* statement as "the training Hamiltonian doesn't matter and the readout
Hamiltonian's own classical ranking quality is all that's left" — which is B12/L29's prescription,
now given an a-priori rather than empirical grounding on this specific register.

---

## 5. A3 — MULTI-STAGE CVaR-VQE (Legacy 60 iters → Dist 140 iters, vs single-stage Dist 200 iters)

**COMPLETE.** Registered expectation (staging does not help, by analogy with s21 C7′/C7″) was
**not confirmed as stated** — but the mechanism is the same one running through this whole file.

**PRIMARY, staged − single, tail-average RMSD (target-level, n=16, 4 seeds each)**:
**mean −0.173 [−0.335, −0.036], SE 0.080, MDE 0.223, 8W/5L.** CI excludes zero; |effect|/MDE ≈
0.78, inside the Type-M zone (0.7–1.3×) — **direction favours staging, magnitude is an upper
bound, not a settled number.**

**SECONDARY, argmin readout**: staged − single = exactly 0.000 on every seed (64/64 tied to the
exact DIST optimum in both arms) — Pillar 1 is preserved: the staged procedure never drifts the
final answer toward Legacy's optimum, confirming the promotion criterion's safety check.

**Mechanism check, following §0/§3/§4's pattern rather than treating this as a fourth
phenomenon.** Realised tail size: single-stage mean 1.25 (median 1.0, i.e. near-total collapse,
consistent with T=0/200-iteration DIST-only training in §3); staged mean 3.89 (median 2.0) — a
Legacy warm start leaves the DIST fine-tuning stage measurably less collapsed after only 140
(rather than 200) pure-DIST iterations. Per §0's theorem this again means staging's entire effect
is routed through **which m gets realised**, landing marginally further from the worst part of
the classical m-ladder (pure argmin, m≈1) than single-stage's harder collapse does. **This is not
a new mechanism — it is a third instance of the same one**, and the small, Type-M-zone magnitude
is exactly what "nudging m by a couple of units, on an already-nearly-flat-but-not-quite-at-m=1
part of the ladder" would predict.

**Disposition: direction ESTABLISHED (staging's small residual undercollapse marginally helps),
magnitude NOT MEASURED, mechanism is m-realisation (§0), not a genuine benefit of training on two
Hamiltonians in sequence.** This is a substantially weaker claim than "staging works here unlike
on the continuous encoding," and it is stated that way rather than as a positive.

---

## 6. Scope limits, stated plainly

- 16 of 126 tuning targets (`s20.qb2_lib.subset(n=16)`) throughout — a declared compute/turnaround
  choice, not a result. All comparisons are target-level (mean over seeds per target, then
  bootstrapped over targets) unless stated as cell-level.
- `H_AMBER` was not run anywhere in this workstream (declared in PREREG_A); no OpenMM context was
  taken, so no announce/release cycle was needed and there is no AMBER-contention risk from this
  lane.
- `R_repair` (ideal-geometry projection of the tail-average) was not computed in the final runs —
  a scope gap, not a null; every RMSD reported above is point-cloud basis, never mixed with a
  repaired-chain number.
- **The one readout genuinely not tested**: a PROBABILITY-WEIGHTED tail average (weighting each
  candidate's coordinates by its trained mass before superposing, rather than the unweighted
  average this project's shipped operator and this workstream both use throughout) is NOT
  covered by §0's theorem — the theorem is specifically about which candidates receive ANY
  positive mass, not about how an average would be weighted across them. This is the one
  remaining place a genuinely different quantum-shaped signal could in principle live on this
  substrate, and it was not tried. Named as the honest open door, not run due to time.
- Every mean above is reported with its own SE and per-comparison MDE, never the retired 0.084 Å
  constant.

---

## 7. One-paragraph summary for the report

The candidate-state encoding was built, its permutation gauge was gated (soundness 0/320 exact;
achieved tail-average gauge sensitivity real but moderate, median 0.03 Å, 5.9% of cells >0.5 Å),
and its two order-based readouts were both found to be provably classical operations on an
exactly-enumerable register: argmin is always the classical argmin of whichever Hamiltonian scores
it, and an unweighted tail-average is always classical top-m for whichever m the training realises
— proved by direct set-equality, not inferred from RMSD closeness. Every subsequent manipulation
tried (CVaR α, entropy temperature, cross-Hamiltonian readout, two-stage training) therefore
reduces to moving m along the classical m-ladder Sprint 21 already fully characterised, and none
of them beats the size-matched classical selector. The readout-Hamiltonian effect replicates
Sprint 21's sign and rough magnitude on a second, independent substrate, which is worth keeping,
but on this substrate it now has an a-priori mechanism rather than an empirical one. The one
untested readout that escapes this closure (probability-weighted averaging) is named as the open
door for anyone continuing this line.
