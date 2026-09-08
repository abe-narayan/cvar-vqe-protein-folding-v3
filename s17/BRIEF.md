# SPRINT 17 — SHARED AGENT CONTRACT

Read this in full before running anything. It is binding on every workstream.

---

## 1. The sprint question

> **How much accuracy is already contained in the candidate pool, how much can principled
> selection recover, and can genuine CVaR-VQE contribute anything a strong classical sampler
> cannot?**

Target hierarchy: **< 2.5 Å** mean full-chain Cα-RMSD first; **< 2.0 Å** preferred; substantially
below 2.0 Å if the data support it. On 126 cluster-disjoint targets, 9–16 residues, five folds.

**The central architectural hypothesis**: switch from an *averaging* readout to a *selection*
readout. Sprint 16 proved the averaging readout is capped by an exact identity —
`readout² = mean member error² − diversity²` — whose two terms move together under any selection
pressure. Coordinate averaging remains an explicit baseline, never a silent partner.

---

## 2. THE FACT THAT REFRAMES THE SPRINT — verify it, then build on it

The programme has been planning against a "K = 500 pool oracle of 1.711 Å". **That is not the
retrieval ceiling.** It is the ceiling of the first 500 entries of a BLOSUM62 ranking over a
universe of **13,000–27,000 windows per target** (`u["W"]`, ordered by `u["order"]`;
`I.pool_idx` truncates at 500).

The coordinator's `s17/oracle_map.py` measures the full map. Read
`s17/results/oracle_map.json` and the report in `s17/results/oracle_map.log` before designing
anything — it is the shared instrument of this sprint. Do not re-derive it; extend it.

---

## 3. The four problems — say which one your experiment attacks

Never collapse these into one metric.

| | problem | question |
|---|---|---|
| **A** | retrieval / generation | does the pool contain a good structure? |
| **B** | selection / ranking | can we find it without the native? |
| **C** | ensemble construction | can we build a candidate distribution whose error/diversity geometry makes selection possible? |
| **D** | validity / repair | can we make it a legal all-atom peptide without destroying the Cα prediction? |

---

## 4. Reporting conventions — standing law

**Every method reports three numbers, never mixed:**

```
ORACLE CEILING      what is attainable if you knew the answer
REALIZED            what the native-free method actually returns
GAP TO ORACLE       the difference
```

- **TARGET is the unit of analysis.** Paired fold-clustered bootstrap CIs, medians beside means,
  win/loss counts beside every mean. No row-level pseudo-replication — Sprint 16 lost a headline
  result to exactly that (162 rows sitting on 9 targets).
- **Every accuracy arm carries BOTH controls**: a **zero-information reference** (a constant ideal
  α-helix / β-strand, a random pool window) and a **matched random** operation (random selection of
  the same count, random direction of the same magnitude). Sprint 16 established that without these
  two the programme's evidence read ~5× stronger than it was.
- **Stratify.** Report by chain length, fold, secondary-structure composition, pool recall, and
  target difficulty. A mean improvement carried by easy targets is not a result.
- **Claim labels are mandatory and must not blur**: ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN ·
  REFUTED · **EXACT** (a theorem, not a discovery) · **ORACLE** (needs the native, never predictive).

---

## 5. THE SPRINT 16 FAILURE MODE — this is what you are guarding against

Across two sprints and four independent audits, **every published number reproduced exactly from
its artefact, and roughly a third of the conclusions built on them were still wrong.** The failure
is always the same shape: *a quantity is measured correctly and then read as if it were a different
quantity.*

Instances, all real:

| measured | read as | cost |
|---|---|---|
| arithmetic mean of a **signed** cosine | input to a law **quadratic and two-sided** in it | a missing sign bit read as absent signal |
| **free-superposition** member error | the **common-frame** term the identity consumes | attribution wrong by ≈2× |
| **arithmetic** mean of two RMSDs | the **quadratic** mean the identity requires | 54% of a "prediction error" |
| a **row** bootstrap over 162 cells | a target-level interval over 9 targets | the sprint's then-only positive result |
| an **algebraic identity** | an empirical validation | a 31-year-old theorem published as a discovery, twice |
| a frame null's **mean** | the null's behaviour | 0.004 Å floor quoted where the tail is 0.137 Å |
| a validity statistic **alone** | evidence about a force field | every such statistic is maximised by a constant α-helix |
| an effect against a **bad baseline** | the operator's value | −0.022 Å against a baseline 0.155 Å worse than doing nothing |

> **THE RULE.** Before quoting a statistic as evidence for a law, **derive which functional of it
> the law consumes, and quote that functional.** Where a law is quadratic, sign-invariant, or
> frame-dependent, the arithmetic mean of the raw per-target quantity is not it. **A mechanism
> sentence needs its own control. An exhibit derivable on random point clouds is not an exhibit.**

And: **a transcription audit cannot see a formula error.** An audit must re-derive, not re-check.

---

## 6. Hard rules

**Leakage.** Never use native information to construct a predictor. Natives are for evaluation,
clearly labelled ORACLE diagnostics, auditing, and post hoc interpretation. Native RMSD may be used
as a **training label inside training folds only**; inference must be native-free. No
native-derived quantity may enter a parameter, threshold, stopping rule, or architecture decision.

**The sealed benchmark.** 60 targets. **Do not read it, do not derive its identities, do not tune
against it.** No module in Sprints 15–16 read it and none in Sprint 17 may. Unlock requires all of:
architecture frozen · controls frozen · evaluation script frozen · success criterion pre-registered ·
no cheap internal experiment could still change the architecture · no benchmark-derived tuning.

**Frozen RMSD definition.** Cα only; every residue including both termini; no trimming; index
correspondence; uniform weights; proper rotations only; model 1 of the native; chain breaks scored,
not rejected.

**Seeding.** Use `s15/seed.py`'s `stable_rng` / `stable_seed`. Bare `hash()` is salted per process
and has already destroyed one set of results in this programme; two violations were found in the
Sprint 15 tree by audit.

**AMBER data rules.** The convergence gate (final energy finite and ≤ 1000 kcal/mol) must be
declared before use and **reported with its exclusion count every time**. The rotated-lab-frame null
must be **reported with its MAXIMUM, not only its mean** — the mean is 0.004 Å and the tail is
0.137 Å. The 0.08 Å multi-start false-positive floor **does not transfer** to the AMBER path (no
RNG there); that path's own floor is ≈0.004 Å gated, up to 0.038 Å ungated.

**Legacy weights** are `DEFAULT_WEIGHTS` and are **never fitted** — Sprint 16 showed leave-fold-out
fitting moves the certified optimum +0.78 Å [+0.04, +1.50] the *wrong* way.

---

## 7. Known ceilings and standing results — verify, don't assume

| quantity | value |
|---|---|
| incumbent (deployed) | **3.204 Å** |
| raw coordinate average, before repair | **3.050 Å** |
| ORACLE best in shipped top-75 | 2.306 Å |
| ORACLE best in K = 500 | 1.711 Å |
| **ORACLE best in the FULL universe** | **measure it — this is the sprint's new instrument** |
| certified argmin of the deployed objective (19 enumerated targets) | 2.661 Å, reached exactly by classical 1-opt in 100% of cells |
| enumerated space best (19 targets) | 1.030 Å |
| targets with K=500 pool best > 2.0 Å | 53 / 126 |
| zero-recall targets (FAIL18) | 18 / 126 |
| repair tax (projection / AMBER vs doing nothing) | +0.155 / +0.133 Å, losing 27W/99L and 24W/99L |
| in-band ordering, within target / across targets | 0.986 / **0.600**; 2.0 Å needs **0.638** |

**Standing refutations — do not repeat these experiments** unless you can state exactly what new
hypothesis makes the rerun different: Jacobian quiet-subspace steering; naive torsion interpolation;
torsion-space averaging; diversity-aware selection as an RMSD cure (the two identity terms cancel at
par); deep VQE optimisation of the existing objective; AMBER sweeps that ignore the Cα
accuracy/validity trade; fitting Legacy's weights.

**Standing results about the pillars.** CVaR-VQE: loses to greedy 1-opt at 10/10 objective-quality
rungs, is matched by annealing at ¼ budget, reaches the certified optimum in 0–32% of cells against
greedy's 68–100%; α is a **temperature** a classical thermostat reproduces (Pearson +0.93/+0.98) and
beats by 1.6–1.7×; it **samples** and **represents** (3.6× near-native mass enrichment, mode
3.792 → 3.014 Å) but does not optimise to completion or explore. Legacy: **detection**, AUROC
0.93–0.99 on its non-tautological `torsion` component; as a ranker indistinguishable from random.
AMBER: **stereochemical repair** real, but only as a conjunction with staying near its input — a
constant α-helix scores rama 1.000 and zero clashes by construction; its displacement has a
**negative** cosine with the true residual and is 55.1% non-torsional, i.e. it moves where Cα-RMSD
cannot see.

---

## 8. Pre-registration is mandatory

Every substantial experiment needs, **written before the run**: hypothesis · expected outcome ·
the strongest available control · success criterion · **falsifier**. Put it in your module's
docstring or a `PREREG_*.md`, and say in your findings whether each rule fired.

A leave-fold-out selector that lands on the **boundary of its ladder** is choosing an endpoint,
not a parameter — treat it as a substitution until shown otherwise.

---

## 9. Compute

The box has 8 cores (4 fast + 4 slow ≈ 6.43 core-equivalents) and ~11 GB usable RAM; three heavy
jobs breach memory. Target ~95% CPU during heavy work; **never deliberately exceed 97%**. Check load
before launching anything long (`s16/energy_lib.cpu_hold` is a bounded yielder). Keep to one heavy
process per workstream. Checkpoint every module per target so an interrupted run is readable.

Cache aggressively — retrieved candidates, ESM representations, distance predictions, Legacy
components, geometric transforms — but **never cache anything that changes scientific randomness
without explicit seed/version control**, and never optimise a scientific module while silently
changing its mathematical definition. Record outputs before optimisation and prove equivalence
after.

---

## 10. What is valuable

A clean null is more valuable than a fake breakthrough. A contradiction is valuable. A failed
hypothesis that eliminates an architecture class is valuable. A discovered upper bound is valuable.

**Do not optimise the report. Find a better predictor.**
