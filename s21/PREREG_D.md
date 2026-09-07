# PRE-REGISTRATION — WORKSTREAM D (theory, literature, adversarial audit)
## Sprint 21

Written before the corresponding results were inspected. **Not edited after data.** Where a
pre-registration turns out to have been mis-specified I say so in `agentD_FINDINGS.md`, leave this
file unedited, and record the untested regime as OPEN (`s21/BRIEF.md` §8).

Labels per `s21/BRIEF.md` §9. **TARGET is the unit.** MDE at 80% power = **0.084 Å**.
Every structural number states its basis (point cloud / built chain / repaired emission).

**Benchmark seal, verified by hash and not by read, before anything else was done:**

    results/benchmark_manifest.json
      sha256  a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d
      size    8002 bytes      mtime 2026-09-01 18:52:50

Byte-identical to the certificate Sprint 20 recorded. `results/benchmark_manifest.json` was not
opened, probed, or derived from.

---

## D1 — THE ENCODING LEVER IS CONFOUNDED WITH STEP COUNT

**Target of attack.** `s21/BRIEF.md` §5 item 3 and Sprint 20 §8: *θ versus (cos θ, sin θ),
physically identical, wins on RMSD in 11/12 cells, p = 0.006, up to −0.649 Å.* Sprint 18's headline
was an encoding artefact; two physically identical parameterisations producing different RMSD is
that shape again.

**Facts established from the artefacts before this pre-registration was written** (they are audit,
not experiment, and are reported as such):

* `AMB|nelder` and `AMBc|nelder` are **bit-identical on all 10 targets × 2 seeds in both
  encodings**. `AMBc = sign(E)·log1p(|E|)` is strictly monotone and Nelder–Mead is comparison-only,
  so this is an **identity**. Two of the "twelve cells" are one measurement. Deduplicated: 11 cells,
  10 negative, sign-test p = 0.0117, and two — not three — distinct cells have CIs excluding zero.
* Recomputed at the brief's mandated unit (mean Δ per target over the 11 distinct cells, paired
  bootstrap over 10 targets): **−0.2167 [−0.3239, −0.1095], 9W/1L, sign p = 0.0215.** The
  measurement survives. Its interpretation is what is at issue.
* **The embedding takes fewer optimisation steps in every cell.** From the artefacts' own `iters`:
  `adam_fd` 5.1 vs 10.5 (ratio 0.49); `spsa` 159.2 vs 207.6 (0.77). `lbfgs_fd`, the arm least
  sensitive to step count, has the three smallest effects and the only positive one.
* SPSA's deficit is **not** the declared 4n-vs-2n dimension effect. It is `arm_spsa`'s
  `track_quality` diagnostic: three central-FD probes at `2·d` units each, and `d` doubles in the
  embedding, so the diagnostic eats 288/512 budget units embedded against 144/512 in θ.

**Hypothesis (D1).** The encoding effect on RMSD is a **step-size / step-count** effect, not a
property of the coordinate system. Mechanism: on this instrument optimising the deployed objectives
harder makes the *structure* worse (Q3; `search-saturates-discrimination-binds`;
`concentration-is-wrong-when-discrimination-binds`), the start is retrieval-derived and good, and
the embedding moves less far from it. Two channels: (a) fewer steps at equal budget, (b) the
embedding's radius `‖u‖` is a **pure gauge direction** — `arctan2` is scale-invariant so the true
gradient's radial component is exactly 0 — into which a finite-difference optimiser injects noise;
`adam_fd` in particular normalises per coordinate by `√v` and therefore takes a **full-size step
along a direction whose gradient is pure noise**. If `‖u‖` diffuses upward the effective angular
step decays as `1/‖u‖`: an unintended annealing schedule.

**Prediction.** (i) Realised angular displacement `‖wrap(z_best − z0)‖` is **smaller** in the
embedding in a majority of cells. (ii) `‖u‖` **grows** along the embedded trajectory. (iii) On a
budget grid, RMSD is a **increasing** function of realised step count within each encoding, and
(iv) at **matched step count** the encoding gap shrinks to inside the MDE.

**Primary endpoint.** The encoding difference at **matched realised step count**, per (kind, arm),
paired over targets, bootstrap CI over targets, with W/L, median and per-fold. Matching is done in
two directions (θ down to the embedding's step count; the embedding up to θ's) so the answer cannot
be an artefact of which arm was handicapped.

**Falsifier (pre-registered, and it falsifies MY hypothesis).** If at matched step count the
embedding still wins by more than the MDE with a CI excluding zero in **both** matching directions,
the encoding lever is **not** a step-count effect. I withdraw D1 and record the lever as SUPPORTED.

**Null / matched control.** The step-count-matched θ arm is the operative control — matched in the
space the operator works in (§7 rule 1). A budget grid `{256, 512, 1024}` in **both** encodings, so
the step-count response curve is measured rather than assumed, and the matched comparison is read
off it rather than from a single hand-picked budget.

**Second control, against my own mechanism.** `track_quality` is switched **off** in every arm here,
which removes the asymmetric diagnostic cost. If the encoding gap survives that but dies under
step-matching, the diagnostic was not the mechanism but step count still is; if it dies under
`track_quality=False` alone, the diagnostic was the whole thing. Both are reported.

**Soundness gate before any number is read.** My arm implementations are verbatim copies of
`s20.qb2_opt.arm_spsa` / `arm_adam` with recording added. With `track_quality=True` and the same
seed, start and budget they must reproduce Sprint 20's stored rows **bit-for-bit**. If the gate
fails, no D1 number is reported. The gate's firing count is reported (a gate that never fires is
not evidence, §7 rule 3).

**Budget.** 10 targets (the same 10 Sprint 20 used, so the comparison is paired to its own
instrument), kinds `{LEG, AMBc}`, arms `{spsa, adam_fd}`, encodings `{θ, emb}`, budgets
`{256, 512, 1024}`, **4 seeds** (§6 minimum; Sprint 20 used 2). One heavy process, BLAS capped at 1.

**Promotion rule.** None. This is an audit of a claim, not a method. Its only output is a label on
somebody else's lever.

**Declared limitation, in advance.** n = 10 targets is Sprint 20's panel, not the n = 126
instrument, and 10 targets cannot resolve 0.084 Å. This experiment is powered to discriminate
between "≈ −0.25 Å" and "≈ 0", which is the question, and is **not** powered to certify a residual
encoding effect of MDE size. A residual inside the MDE will be labelled **NOT MEASURED**, not
"matched".

---

## D2 — CVaR-VQE AND ARGMIN ARE THE SAME SELECTOR ON A POOL, AND THE "TAIL" ARM MEASURES AVERAGING

**Target of attack.** `s21/BRIEF.md` §2 and `s21/tailprice.py`: *"argmin is not CVaR … a tail MEAN
and a MINIMUM are different operators on the same energy … whether that difference helps is the
sprint's first measurable question."*

**The operator, derived before its statistic is interpreted (§9).** From source,
`core/quantum.py:1600–1626`, a VQE run emits `vqe_bitstring` (argmin over the final distribution's
samples), `vqe_modal_bitstring` (the mode) and `best_seen_bitstring` (argmin over everything seen).
**None of them is a tail mean.** CVaR is the *training objective*; the *readout* is an argmin or a
mode. `s21/tailprice.py`'s `tail{a}` arm is the **coordinate average of the α-tail's members**,
which is neither.

**Claim D2a (to be stated as EXACT or withdrawn).** Let the selector be restricted to a finite pool
`P`, let the training objective be `CVaR_α(H)` and let the readout be `argmin_{x ∈ supp(law)} H(x)`.
Q6 (EXACT) says every CVaR-optimal law is supported inside the α-tail `T_α(H) ⊆ P`. The α-tail of
`H` contains `argmin_P H`. Therefore the best attainable readout under a CVaR-trained law equals
`argmin_P H`, and **CVaR-VQE and argmin have the same optimum on a pool when the training objective
and the readout energy are the same `H`.** Scope conditions, stated in advance because the claim is
worthless without them: it needs (i) pool restriction, (ii) same `H` for training and readout,
(iii) an argmin (or any order-based) readout. It **fails** if the readout is an average, if the
readout energy differs from the training energy, or if the law is not converged.

**Claim D2b (measurement).** The `tail{a} − argmin` primary endpoint is confounded by the
**averaging operator**, which `averaging-space-beats-the-objective` prices at ~1.0 Å against
torsion averaging and which S19 shows dominates the objective by ~6×. Prediction: **a
zero-information energy** — a per-target random permutation of the pool — reproduces most of the
`tail(0.15) − argmin` gap. If it does, the endpoint cannot distinguish a good Hamiltonian from
noise and must be replaced.

**Primary endpoint.** `tail(α) − argmin` for a **random** energy, beside the same quantity for
`disto` and `legacy`. Paired over targets, bootstrap CI, n = 126.

**Falsifier.** If the random energy's `tail(0.15) − argmin` gap is **less than half** the physics
Hamiltonians' gaps, the endpoint is not dominated by averaging and D2b is **REFUTED**.

**The operator decomposition, which is the constructive half.** Four readouts on the same α-tail,
so that "tail vs min" and "average vs single structure" are separated instead of confounded:

    argmin(H)        single structure, lowest H                 Sprint 20's operator
    tail_member(H)   a UNIFORMLY RANDOM member of the α-tail    what a converged CVaR law + a
                                                                one-shot readout actually gives
    tail_medoid(H)   the α-tail's medoid                        single structure, set-informed
    tail_avg(H)      the α-tail's coordinate average            `tailprice.py`'s arm

each against its **matched-count random control** (random subset of the same size, same readout).
`tail_member` is the arm the sprint is missing and it is the one that prices a real CVaR readout.

**Budget.** n = 126, no optimisation, pool-restricted, distogram + Legacy + a random control.
Genuine AMBER over 126 × 500 continuous single points is deferred to a second pass and is scoped
explicitly; the operator decomposition does not depend on which `H` is used and is demonstrated on
the ones that are free.

**Promotion rule.** None. Output is a corrected endpoint for Workstream A, and a label.

---

## D3 — THE MATCHED-DISPLACEMENT NULL FOR THE λ-CONTINUATION

**Target of attack.** `s21/BRIEF.md` §5 item 1, `H(λ) = (1−λ)H_L + λH_A`, flagged as potentially the
highest-value idea.

**Hypothesis.** A continuation schedule that ends up taking **smaller total displacement** in
torsion space reproduces a "gain" that is a step-size control. Sprint 20 (Workstream C, L12)
measured **72% of AMBER's damage as move size, not direction**: any realisable displacement of
0.577 rad/coordinate costs +0.4438 [+0.3884, +0.4891], while a 0.077-rad displacement costs
+0.0045 [−0.0111, +0.0191].

**Prediction.** `Δ RMSD` across the λ-path is predicted by `‖wrap(θ_λ − θ_0)‖` alone, with the
λ-specific residual inside the MDE.

**Primary endpoint.** The λ-path arm minus a **displacement-matched isotropic control**: a random
realisable move of the *same* per-coordinate magnitude, same start, same projection and same
repair. Paired over targets, CI over targets.

**Falsifier.** If the λ-path beats its displacement-matched control by more than the MDE with a CI
excluding zero, the continuation carries direction information and D3 is **REFUTED**.

**Status at pre-registration.** Workstream C has produced no λ-path artefact this sprint. This
entry is registered **now**, before their data exists, precisely so the null cannot be chosen after
seeing it. If no artefact appears, D3 is recorded **NOT RUN**, not "matched".

---

## D4 — LITERATURE: IS A QUANTUM-RESOURCE CLAIM REACHABLE HERE AT ALL?

Not an experiment; no falsifier is claimed for a literature search. Recorded here so the questions
are fixed before the answers are read:

(a) At what bond dimension does classical simulability actually break for this circuit class?
(b) Does any ansatz that escapes it remain trainable at this qubit count (n ≈ 9–16 residues)?
(c) Has anyone demonstrated a genuine separation on a comparable continuous-variable encoding?

Sprint 20 established as **EXACT** that the deployed ansatz is bond-dimension 4 — a ≤16-state HMM,
classically samplable by construction. The output is a **reachability verdict**, and I will state it
as one, including if the verdict is that no quantum-resource claim is available at this scale.

---

## STANDING RULES I HOLD MYSELF TO

* Every artefact under `s21/results/` carries a completion flag that requires the **full**
  configuration — all targets, all kinds, all arms, all seeds, all budgets, zero skipped rows — and
  not the subset it happened to be called with. Sprint 21's first artefact
  (`s21/results/tailprice.json`) was stamped `complete: true` with 3/3 rows skipped; that is the
  hazard this rule exists for.
* CI construction is quoted with the CI. `s12.instrument.paired` is i.i.d. over targets, not
  fold-clustered, across 293 call sites; where a disposition could depend on it I use a
  fold-clustered interval and say so.
* Any alignment or correlation between two quantities measured against a common reference gets its
  **shared referent floor** measured first (0.505 here).
* A gate's firing count is reported. A bound is verified on the calls it actually bound.
* A zero-spanning CI without the power to exclude the effect of interest is **NOT MEASURED**.

---

## D8 — ADDENDUM, registered 2026-09-07 before the n = 14–16 enumeration was run

**This section was appended after the original pre-registration and before any of its own data
existed.** Nothing above it has been edited. The coordinator asked for my expectation on the
record before I look, and this is it.

**The task.** `s21/latentfull.py` enumerated all `2^n` latent configurations on the **75** targets
with `n ≤ 13`. The remaining **51** targets (`n = 14, 15, 16`; 16,384 / 32,768 / 65,536
configurations each, 2.13 M total) complete the instrument at n = 126. Those 51 are the **only**
targets on which the deployed VQE's 8192-evaluation budget was smaller than the latent space — the
only ones where it was genuinely *searching* rather than resampling (§3.1a).

**The coordinator's registered expectation:** the discrimination gap is *the same* on the 51,
because nothing in the mechanism is length-dependent.

**Mine differs, and I am saying so before looking.**

* **P1 (my primary prediction).** The discrimination gap `latent_argmin − latent_oracle` will be
  **LARGER** on `n = 14–16` than on `n ≤ 13`. Mechanism: `latent_oracle` is a **minimum over `2^n`
  structures**, so it improves with `n` by pure min-of-N, while `latent_argmin` follows the
  objective's ordering and has no such mechanism. Guessed magnitude: **+0.2 to +0.6 Å** larger.
* **P2 (where I agree with the coordinator).** The **qualitative** conclusion is unchanged: the gap
  stays strongly positive with **0 or near-0 of 51** targets going the other way, and per-fold signs
  stay uniform.
* **P3.** The matched primary `latent_argmin_bayes − pool_argmin_rb` will **not** become
  significantly negative on the 51. Exhaustive search will still not beat the pool, *even on the
  only targets where the sampler was genuinely searching.*

**Falsifier for the interesting reading.** If the gap is **smaller** on `n = 14–16` with a CI
excluding zero on the *difference between the two n-ranges*, then budget-vs-latent-size **was** the
operative variable, my min-of-N mechanism is wrong, and P1 is REFUTED. If it is larger, P1 is
supported but **see the confound below before believing it means anything.**

**THE CONFOUND, DECLARED IN ADVANCE, AND IT LIMITS WHAT P1 CAN MEAN.** Comparing the gap across
`n` compares a **min over `2^13` = 8,192** structures with a **min over `2^16` = 65,536** — an 8×
larger set on the ORACLE arm alone. So a larger gap at large `n` is *partly an artefact of the
ORACLE arm's set size*, not evidence about discrimination. The clean control is a **matched-set
ORACLE**: subsample each large latent to 8,192 configurations and re-take the minimum. That cannot
be recovered post hoc from the artefact's summary columns, and I am **not** modifying the
coordinator's schema mid-run and breaking their merge. **Consequence, accepted in advance: the
cross-`n` comparison of the gap's MAGNITUDE is CONFOUNDED and I will label it NOT MEASURED
whichever way it comes out.** What the 51 targets *can* answer cleanly is P2 and P3 — whether the
qualitative conclusion and the matched primary hold on the half of the instrument where the
sampler was genuinely searching — and that is the question the task was set to answer.

**Also on the record before looking:** the coordinator flagged that folds 1 and 2 carried the
largest gaps and asked whether those folds skew long. That is instrument metadata rather than an
outcome, so checking it does not consume the prediction; I will check it and report it either way.

**Operator forks, enumerated before the run** (the rule adopted into `BRIEF` §7 this sprint):
*functional* — Bayes (primary, matches the pool arm) with squared reported beside it;
*basis* — the rebuild for both arms via `pool_argmin_rb`, with the window-basis row kept beside it;
*readout* — argmin on both sides, single structure, matched; *null* — `latent_mean`, the space's own
mean; *normalisation* — none, both arms are RMSD in Å. The alternative not taken on each is named.

**Budget.** 51 targets, 2.13 M configurations, ~30–60 min, one process.
**Promotion rule.** None. This removes a scope caveat; it does not create a method.

---

## D9 — ADDENDUM, registered 2026-09-07 before `s21/latentrank.py` produced any number

Appended after D8's data landed and **before** any output of the coordinator's percentile
experiment existed. Nothing above is edited.

**The question.** The best structure in the latent EXISTS (§1.5: the latent's ORACLE best is
1.18 Å better than the incumbent on 104/126) and the objective's exhaustive argmin does not return
it (+1.72 to +1.81 Å, **0 of 122 targets reversed**). *Where in the objective's own ranking does the
ORACLE-best sit?* If it is in the extreme tail a better **readout** can reach it; if it is in the
bulk, no monotone reranking of this objective ever can.

**My prediction, registered before looking.**

* **Q1.** The ORACLE-best's percentile in the objective's ranking will be **in the good tail but not
  at the argmin — my interval is the 0.1st to the 5th percentile**, median target. Basis: over the
  K = 500 pool the per-target `ρ(E_bayes, true)` is **+0.5678 mean / +0.6935 median** (§1.4), and a
  rank correlation of that size places the extreme of one variable well into the tail of the other
  but nowhere near its extreme. It also has to be consistent with the measured facts that the
  argmin misses by 1.72 Å while beating the space's mean by 1.37 Å — a bulk placement (≈50th) is
  inconsistent with the second, and a top-1 placement is inconsistent with the first.
* **Q2.** The **top-M ceiling curve will fall steeply and still not reach the ORACLE at M = 512.**
  At `N = 8192` the 1st percentile is rank 82, inside M = 512; the 5th is rank 410, just inside; the
  40th is rank 3277, far outside. So Q1 and Q2 are the same statement read two ways, and M = 512 is
  the right place to cut.
* **Q3 (the decision).** If Q1 lands inside my interval, a reranker route is **alive but expensive**:
  it needs M in the hundreds, not a top-1 readout. If it lands past the 20th percentile the route is
  **closed** and B's −0.697 Å readout lever is buying something other than reaching the ORACLE.

**Falsifier.** If the median percentile is **worse than the 20th**, Q1 is REFUTED and I was wrong
about the ordering being tail-informative. If it is **better than the 0.01st**, Q1 is also REFUTED
in the other direction and a top-1 readout was already nearly enough.

**A HAZARD IN THE NEW ARM, raised before its numbers exist, and it is the same class as the
min-of-N confound I raised on D8.** The **top-M ceiling is a minimum over M structures**. A curve
that falls from M = 1 to M = 512 therefore falls *partly because 512 draws beat 1 draw*, with no
ordering skill involved. **The matched control is a min-of-M over a RANDOM window of the same size
M**, plotted on the same axes. Without it, "a reranker could reach X Å at M = 512" is not separable
from "the minimum of 512 arbitrary latent configurations is X Å". With it, the vertical gap between
the two curves at each M **is** the objective's ordering skill in exactly the currency a reranker
spends. I predict that gap is **positive and substantial at every M** — the objective does order —
and that **it does not close the 1.72 Å discrimination gap at any M ≤ 512.**

**Promotion rule.** None; this is a bound and a control specification.

---

## D10 — ADDENDUM, registered 2026-09-07 before `s21/d_lrank.py` produced any number

Appended before this experiment ran. Nothing above is edited.

**The question the coordinator identified as the highest-value one left, and it is the right one.**
`latentrank.py` shows an ORACLE ceiling of **1.986 Å inside the objective's top-512** against the
deployed argmin's **3.281 Å** — 1.30 Å of headroom in a 512-candidate window. But every top-M number
there is an **ORACLE**: *"picks perfectly inside the top-M"* consumes the native. **The finding is
not the ceiling; it is whether any NATIVE-FREE readout can traverse that window.** That is measured
on the **pool** in this project (in-band ordering 0.600 across targets; 2.0 Å needs 0.638) and is
**not** measured on the **latent** — and my own D7 showed a pool-measured relationship failing to
transfer, so those bounds do not carry in either direction.

**The design, and it is deliberately the same operator decomposition I ran on the pool** (§1.2), so
the two are directly comparable rather than merely analogous. Within the objective's top-M window of
the enumerated latent, four readouts, each against its **matched-count control**:

    argmin        the deployed readout, M = 1                       the incumbent operator
    member        a uniformly random member of the top-M            worst case for a converged law
    medoid        the top-M's medoid                                single structure, set-informed
    avg           the top-M's coordinate average                    POINT CLOUD basis, labelled

    rand_member / rand_medoid / rand_avg     the SAME readout on a RANDOM M-window
    ORACLE_topM                              the ceiling, labelled, native-consuming

**Predictions, registered.**

* **T1.** `medoid` and `avg` beat `argmin` on the latent's top-M — the same direction as the pool,
  where `tail_medoid|disto` is −0.447 and `tail_avg|disto` −0.377 against matched random (§1.2d).
* **T2 (the one that matters, and the one I expect to be closest).** Against the **matched-count
  random window** — the only comparison that means anything — the gain will be **smaller on the
  latent than on the pool**, and I am not confident it clears the MDE. Mechanism, from R1's
  falsification (§4.1): what makes averaging work is **error incoherence**, not set spread, and the
  latent's top-M are all built from the *same per-residue basin mixtures*, so their errors should be
  far more **coherent** than those of independently retrieved pool windows. Coherent error survives
  averaging.
* **T3.** No native-free readout reaches the ORACLE top-M ceiling at any M ≤ 512. The residual is
  the part that needs a signal this project does not have.

**Falsifier.** If a native-free readout beats its matched-count random control by **more than that
comparison's own MDE with a CI excluding zero at every M**, then the latent's top-M *is* traversable
native-free and T2 is REFUTED — which would be a **positive** result and the most valuable outcome
available tonight. If instead every readout ties its matched control, the top-M window carries no
native-free structure beyond what the objective already used, and the reranker route is closed for
*this* objective's window even though its ORACLE ceiling is 1.30 Å below the argmin.

**Operator forks, declared per §7 rule 0.** *Functional* — Bayes (matches `pool_argmin_rb` and the
`latentrank` primary); squared not taken, and D7 licenses that as immaterial (median ρ 0.973).
*Basis* — built chain from `build_ca_exact` on **both** sides, so no window-vs-rebuild mismatch;
the `avg` arm is a POINT CLOUD and is labelled on its row. *Readout* — all four reported, none
privileged. *Null* — matched-count random window at the same M and the same readout, 8 draws.
*Normalisation* — none; both arms are RMSD in Å. *Tie handling* — stable sort throughout, and where
an argmin could land on tied scores the outcome is averaged over the tied set (project memory:
`np.argmin` on a tied signal once read the ORACLE sort order and invented a 1.386 Å winner).

**Budget.** `n ≤ 13` targets only (75 of 126, `2^n ≤ 8192`), M ∈ {1, 8, 64, 512}, one process, no
AMBER. **Declared limitation in advance:** restricting to `n ≤ 13` is the same length-defined
subset D8 was run to remove, and I am accepting it here because the experiment needs the top-M
*structures* rather than a summary and `n ≥ 14` costs 8× more. The result will be labelled as a
`n ≤ 13` result and **not** generalised to n = 126.

**Promotion rule.** None. A surviving native-free gain earns a replication on the pool's own
instrument and on `n ≥ 14`, not a claim.

---

## D11 — ADDENDUM, registered 2026-09-07 before the containment column was read

Appended before I looked at any AUC. The coordinator asked what I would conclude if the best
native-free containment predictor lands between 0.55 and 0.65 AUC. **The answer is arithmetic, and
it makes my own pre-registered 0.65 threshold TOO LENIENT.**

**The null for a best-of-K AUC, simulated at this experiment's own n and split** (n = 75, roughly
44 contained / 31 not, K = 10 candidate predictors, labels **pure noise**, folded |AUC| so
direction is free):

    a SINGLE predictor's |AUC|      median 0.546     95th percentile 0.633
    the BEST of K = 10              median 0.625     95th percentile 0.690

So **a best-of-ten AUC of 0.625 is the MEDIAN of pure noise** at this sample size, and the whole
0.55–0.65 band the question asks about sits **inside** the noise for a maximum over ten predictors.
My earlier "nothing below 0.65 is usable" was stated by assertion; it is wrong, and it is wrong in
the permissive direction.

**Registered conclusion rule, replacing the 0.65 threshold.**

1. Every predictor's AUC is reported with a bootstrap CI **and** beside a **permutation null on the
   label** — shuffle containment, recompute all K predictors, take the maximum, repeat. That is the
   correct null for *"the best of K"*, which is what a reader will look at.
2. **The bimodal route is DEMONSTRATED only if the best predictor's AUC exceeds the 95th percentile
   of that permuted best-of-K null** (≈ 0.690 here) **and** its own CI excludes 0.5.
3. **Anything in 0.55–0.65 is reported as NOT MEASURED and the route as NOT DEMONSTRATED — never
   as "promising", "suggestive", or "worth pursuing".** That range is what noise produces.
4. If a predictor clears the bar, it is still a single-instrument result on `n ≤ 13` and earns a
   replication on `n ≥ 14` and on the pool, not a claim.

This also retroactively raises the bar on anything else in this document that quoted a best-of-many
statistic without a maximum-order null. I have checked: **§1.2d, §4.1 and §5.2 each pre-specified
their comparison rather than selecting the best of several**, so none of them needs it. The one
place a maximum was taken over arms — §1.3e's four encoding cells — reports every cell rather than
the best, and its disposition rests on the *pooled* target-level number, not on a maximum.
