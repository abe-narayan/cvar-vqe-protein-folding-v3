# SPRINT 19 — SHARED AGENT CONTRACT
## Why are the predictor's errors structurally maladaptive, and can that be fixed?

Read in full before running anything. Binding.

---

## 1. The mission

Primary metric: **mean full-chain Cα-RMSD, frozen implementation, 126 cluster-disjoint targets.**

    incumbent                     3.204 A
    best structure built          3.048 A   (coordinate average of the shipped top-75)
    primary target              < 2.500 A
    stretch target              < 2.000 A

Minimum detectable effect at 80% power on this instrument: **0.084 A.** A null below that magnitude
is uninformative, not negative.

---

## 2. THE SPRINT-18 RESULT THIS SPRINT EXISTS TO ATTACK

> **The distogram's errors are systematically worse than random errors of the same magnitude.**

Permuting the residuals across pairs — keeping their sizes and the correct weights — takes the
refinement from **3.610 A to 2.560-2.635 A**. That survived separation-stratified permutation,
isotropic replacement, weight permutation, weight flattening, and an independent reimplementation in
another lane. It is the most robust result in the programme.

Everything downstream of the predictor is now measured to be near-optimal: the `1/sd^2` weighting is
**validated** (deleting it costs +0.153 [+0.045, +0.272], and the full sd^-p curve is monotone below
the deployed exponent), the functional form is sound (with true distances it reaches the **0.083 A**
parameterisation floor), the readout is optimal, and the optimiser is not the constraint (no
optimiser at any budget changes a number).

**So the question is no longer how to consume the distogram. It is why the predictor makes
structured errors, and how to change where they land.**

---

## 3. THE COORDINATOR'S OPENING MEASUREMENT — read this before designing anything

The predictor emits a **full distribution over 17 distance bins for every pair**
(`core/predict.py::Distogram.prob`). The *selection* path consumes the whole distribution through a
Bayes-risk lookup. **The objective path does not.** `s15/distcal.gather` collapses it:

```python
"dhat": np.asarray(dg["expected"], float),      # first moment
"sd":   np.maximum(np.asarray(dg["sd"], float), 1e-3),   # second moment
```

and **every objective in Sprints 17 and 18 minimised `((d - dhat)/sd)^2`.** Measured on 12 targets:

| quantity | value |
|---|---|
| probability mass the model puts **within 1 A of its own reported mean** | **0.641** |
| fraction of pairs with **two or more modes** (>2% mass) | **0.217** |
| mean predictive entropy | 1.43 nats (uniform = 2.83) |

> **CORRECTED 2026-09-06, BEFORE ANY LANE BUILT ON IT -- the mechanism stated here was WRONG, and
> it was the coordinator's error, not a workstream's.** The paragraph that stood here asserted that
> the collapse error acts *"coherently across correlated pairs, because a whole region flipping
> between two conformer families moves many pairs the same way at once."*
>
> **AGENT D refuted the coherence claim at n = 126** (artefacts complete):
>
> * multimodal pairs are **scattered, not clustered** -- adjacency 0.072 against a
>   separation-stratified label-permutation null of 0.070, z = +0.27 [+0.08, +0.46];
> * the implied correction has **zero spatial coherence** -- adjacent-pair sign agreement 0.540 vs
>   non-adjacent 0.537, diff +0.003 [-0.001, +0.008];
> * **the dissociation that kills it** -- the *real* residual sign(d_true - mean) **is** spatially
>   coherent (+0.0385 [+0.0320, +0.0457]), but **to the same degree on unimodal pairs (+0.0429
>   [+0.0336, +0.0537]) as on multimodal ones (+0.0489 [+0.0268, +0.0731])**. The coherence of the
>   error -- which is what the Sprint-18 result is about -- **is not carried by multimodality. It is
>   present in full where there is no multimodality at all.**
>
> **The "aims at unlikely values" harm does not survive its confound.** Unmatched, multimodal pairs
> look +0.410 A worse at the aim point, but carry mean sd 2.16 against 1.25. Matched within
> (separation x sd-decile) cells the sign **reverses** (-0.174 [-0.408, +0.072], NOT MEASURED); under
> regression adjustment on sep x log sd it reverses **significantly**, -0.422 [-0.595, -0.248].
> **Conditioned on spread and separation, the mean of a multimodal predictive distribution is at
> least as close to the truth as the mean of a unimodal one.**
>
> **The census is real but was overstated.** At n = 126 the n = 12 figures reproduce (0.241
> multimodal, mass-within-1A 0.634, entropy 1.44), and the bin-edge attack **failed** -- correcting
> mass to density for the non-uniform `BIN_EDGES` changes nothing (0.240). But **60% of flagged pairs
> are shoulders, not modes**: 0.097 under a prominence gate, 0.155/0.165 after 2x coarsening.
> **Quote ~0.10, not 0.22.**
>
> **What survives, and it is not nothing.** ORACLE nearest-mode beats the mean by **-1.068 A [-1.186,
> -0.953]** per pair -- but **61% is pure min-of-N** (matched-|offset| random-sign null: -0.654). The
> mode information proper is **-0.414 [-0.502, -0.324]**. **Both native-free mode rules are WORSE
> than the deployed mean**: highest-density +0.315, uniform-random +0.529.
>
> **A coordinator counter-argument is void, on an EXACT point.** This section previously argued
> "selection already consumes the full distribution, so why did it not show there?"
> `Distogram._risk[.,x] = sum_b p_b |x - c_b| w` is a non-negatively-weighted sum of absolute values,
> hence **convex in x with minimiser the weighted median of CENTRES** -- an identity, not a
> measurement. The selection path consumes the distribution in the letter but is **structurally
> incapable of expressing multimodality**; its silence says nothing either way.

**THE STANDING LEAD THIS LEAVES, and it is now the sprint's central target.** The moment collapse is
real and the modes carry real ORACLE information, but **the coherent error Sprint 18 found is not
multimodality** -- it lives on unimodal pairs in full.

> **What makes the residual spatially coherent on pairs where the model is confident and unimodal?**
> That is where the Sprint-18 mechanism actually is, and it is AGENT A's lane.

*The coordinator's `s19/distobj.py` arm still runs as the STRUCTURAL (RMSD) test of the same
question, pre-committed to close if `mode` does not beat `moment` -- which AGENT D's per-pair
measurement now predicts it will not.*

---

## 4. WHAT IS CLOSED — do not reopen without a materially new mechanism

Degree-1 / residue-additive / Walsh truncation · encoding-dependent objective tricks · objective
tempering (the deployed sd^-2 is at the optimum) · the torsion-prior rescue · Legacy as an objective
term (harmful with a sign, +0.722 [+0.541, +0.893]) · Legacy/AMBER as candidate rankers (**every
physics filter loses to a random gate of the same size**) · `leg_torsion` as a gate (recall converts
*negatively*) · harder search, classical or quantum · treating the optimiser as the bottleneck ·
"the objective's functional form is wrong" · Jacobian/coordinate steering · disagreement features ·
blindly widening K · generic rankers over existing signals.

**Standing roles**: AMBER is stereochemical repair only, with a measurable displacement tax
(+0.133 [+0.112, +0.165] at k=30) and an **irreducible +0.164 A** repair cost. The coordinate average
contracts the backbone 22.4%. The coordinate average is **not a privileged start** (a
zero-information helix start matches it at α=1).

---

## 5. THE DESIGN COMPASS — where error lands beats how much there is

At **identical residual RMS**, outcome spans **2.145 → 2.697 A**:

    fix_confident   2.145   -0.308 [-0.462,-0.150]      fix_long          2.618  +0.165 [+0.028,+0.303]
    fix_short       2.159   -0.294 [-0.424,-0.167]      fix_unconfident   2.697  +0.244 [+0.097,+0.386]
    uniform         2.453    0.000

**Improving confident, short-range pairs beats improving the worst, longest, least-confident ones —
and the latter is significantly WORSE than a uniform improvement.** Any new predictor, loss or
weighting is judged against this compass, not against MAE.

**Aggregate error metrics are not the objective.** A lower MAE/RMS that does not lower Cα-RMSD is a
negative result. Report the **spatial pattern** of the error, always.

---

## 6. FOUR THINGS THAT DISSOCIATE — report them separately, always

**objective quality** (does the number go down) · **structural quality** (does RMSD go down) ·
**alignment** (does lower objective mean lower RMSD — and *name the axis*: pairwise ordering
accuracy, Spearman, Pearson or copula ρ, which are not interchangeable) · **search quality** (can it
be found, at what budget).

Sprint 17 drove the objective down 71% and RMSD **up** 0.561 A. Never assume these move together.

Likewise separate, for every candidate-set intervention: **generation ceiling** (pool best) ·
**selection ceiling** (best a native-free selector can reach) · **repair ceiling** (what survives
valid terminal processing) · and the **realised** number. Never collapse them.

---

## 7. HARD RULES

**Sealed benchmark.** 60 targets. Do not read, probe, derive from, tune against, or use for any
architecture or hyperparameter decision. Clean through Sprint 18; keep it that way. **Do not open
`results/benchmark_manifest.json`.**

**Frozen metric.** Full-chain Cα-RMSD, existing implementation, all residues including termini, no
trimming, proper rotations only, model 1 of the native. No alternative metric reported as primary.
New diagnostics are welcome and must be labelled as diagnostics.

**Native information** is for evaluation, labelled **ORACLE** diagnostics and post-hoc interpretation
only. Native RMSD may be a training label **inside training folds only**; inference must be
native-free. Anything selected using native information is ORACLE and is not a predictive result.

**Seeding.** `s15/seed.py`'s `stable_rng` / `stable_seed`. Bare `hash()` is salted per process and
has destroyed results here.

**Legacy weights** are `DEFAULT_WEIGHTS`, never fitted.

**AMBER discipline.** Convergence gate declared before use and reported with its exclusion count;
the rotated-frame null reported with its **maximum**, not its mean; **a gated arm compared to its own
gated input.**

**The subset trap.** `s15/distcal.fit_correction` fits the leave-fold-out debias on whatever target
list it is handed — running on a subset silently changes `dhat`. Always gather and fit on the full
`I.targets()`, then iterate your subset.

**Artefacts.** Unique immutable result directory per experiment, deterministic seeds, completion
flags, persisted configs and hashes. Never treat a partial file as complete. Never overwrite a prior
sprint's files. Smokes named so they cannot be misread as results.

---

## 8. CONTROLS ARE MANDATORY, AND UNIFORM IS NOT A CONTROL

Every mechanism carries a **plausible zero-information null** and a **matched-random** operation of
the same magnitude/count.

**Uniform is not a zero-information control.** Uniform-on-the-torus places mass on impossible
backbone conformations, making it a *worse* measure rather than an uninformative one. Use
plausible-but-uninformative nulls: a constant ideal α-helix, matched empirical torsion marginals,
label permutations, shuffled assignments, matched-random candidates, isotropic perturbations,
classical thermostats.

**Four zero-information controls in Sprint 18 matched their informative arms.** A constant α-helix
start matched the coordinate average at α=1; a constant α-helix prior matched the retrieval prior at
every λ; and the MJ table's residue identities were worth nothing against a label permutation. **If
your arm does not beat a constant helix, it has not demonstrated conditioning.**

---

## 9. STATISTICS

**TARGET is the unit.** Paired target-level bootstrap CIs, fold-aware, medians and win/loss beside
every mean, per-target distributions, stratified by length and fold. Never a mean alone.

**A zero-spanning CI with insufficient power is NOT MEASURED — never "matched."** A match claim needs
a CI tight enough to exclude the effect size of interest, or n=126. A Sprint-18 "match" published at
n=14 reversed sign by n=30.

**A result below the 0.084 A MDE is not a validated improvement.** An argmin-on-the-mean over the
tuning instrument is a hyperparameter chosen there, and must be labelled as such.

No min-of-N ceiling without its min-of-N null, and the null must name its band. No directional
conclusions from smokes — record the smoke as a smoke.

---

## 10. LABEL EVERY STATEMENT

**EXACT** (a theorem or identity — not a discovery) · **ORACLE** (needs the native) · **ESTABLISHED**
· **SUPPORTED** · **PLAUSIBLE** · **OPEN** · **INCONCLUSIVE** · **REFUTED** · **RETRACTED**.

**Never present an algebraic identity as an empirical finding.** This programme has done it twice.
For every statistic ask first whether it is an identity, a theorem, an implementation consequence, an
observation, a causal hypothesis, or a learned relationship.

---

## 11. READ THE CODE BEFORE BELIEVING THE CLAIM

The source has repeatedly contradicted the table. Sprint 18's sharpest catch was an adversarial lane
reading `objceil.py` line 163 and finding a control that permuted the weights as well as the
residuals — it was not the functional it claimed to be, and a lead built on it had already been sent
to another lane.

For every decisive result inspect: argument order, weights, masks, seeding, frame conventions,
normalisation, whether the random control is genuinely matched, whether "same budget" is the same
budget, and whether evaluation is target-level.

---

## 12. WHEN TO KILL A BRANCH

Immediately downgrade when it loses to a plausible zero-information null · works only on the tuning
subset · has a CI spanning the required effect · dissolves under target-level analysis · is
reproduced by a simpler classical baseline · needs native information · rests on an arbitrary
encoding · improves an intermediate metric while worsening RMSD · needs post-hoc tuning · has a
mathematically invalid mechanism.

**A clean falsification is a successful result.** Interesting is not enough.

---

## 13. THE PILLARS

Genuine **VQE**, genuine **CVaR**, genuine **Legacy**, genuine **AMBER** — never replaced by
placeholders or classical approximations while still being described as quantum/physics components.
**Their architectural role is determined experimentally, not assumed**, and may be: accuracy
component, sampler, proposal mechanism, rejection gate, validity operator, diagnostic, or a
demonstrated non-contributor.

**Do not force a pillar into an accuracy role because the project needs it to exist.** A negative
result is acceptable. A false positive is not. Never call something a quantum improvement because
VQE was somewhere in the pipeline: it needs matched classical controls at matched budget, surviving
target-level statistics, not reproducible by a trivial thermostat, with a mechanism you can defend.

---

## 14. LEDGER

Every branch records: hypothesis · mechanism · prediction · primary metric · controls · budget ·
targets · result · CI · W/L · interpretation · status · why continue or stop. Status is one of
OPEN / SUPPORTED / INCONCLUSIVE / REFUTED / CLOSED / PROMOTED. Pre-register the primary outcome
before the run and **do not edit a pre-registration after seeing results** — if it was mis-specified,
say so, keep it, and record the untested regime as OPEN.

---

## 15. COMPUTE

8 cores (~6.43 core-equivalents), ~11 GB usable RAM. Target ~90-95% CPU; **never deliberately exceed
97%.** One heavy process per workstream; cap BLAS threads (`OMP_NUM_THREADS=1` etc.) on anything you
launch; check load before launching anything long; checkpoint per target. Cache deterministic
quantities; never cache anything that changes scientific randomness without explicit seed/version
control.

---

## 16. THE DEEPEST QUESTION

> The model's pairwise predictions contain real information, but the **distribution of its errors is
> structurally maladaptive.**

Attack that. Try to falsify it. Is it true for every target family? Is it the architecture, the loss,
insufficient joint modelling, multimodality, retrieval bias, a sequence-information limit, the
training distribution, the probability-to-distance conversion, downstream geometry, averaging
incompatible predictions — or a fundamental identifiability limit?

**Do not assume the answer. Find out.**
