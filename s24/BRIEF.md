# SPRINT 24 — LEARNED CANDIDATE GENERATION + CVaR-VQE. DRIVE THE DEV MEAN BELOW 3.0 Å.

**PRIMARY ENDPOINT: mean full-chain Cα-RMSD, 126 cluster-disjoint dev targets, point-cloud basis.
Incumbent 3.0483 Å. TARGET < 3.0 Å. The 60-target benchmark is SEALED — no inspection, no tuning,
no threshold or model selection against it, ever.**

---

## 0. THE THREE PILLARS — NON-NEGOTIABLE

1. **CVaR-VQE remains the selector.** Genuine VQE, genuine CVaR, genuine variational state
   preparation, genuine Hamiltonian energy evaluation, genuine CVaR-tail selection. Classical methods
   may generate, propose, refine, control and ablate. They may **not** be the selector and then be
   called the architecture.
2. **`H_Legacy` genuine and independently evaluable** — the 11-term potential at `DEFAULT_WEIGHTS`.
   Never fitted, never a proxy, never a regression target.
3. **`H_AMBER` genuine and independently evaluable** — ff14SB/GBn2 through the existing OpenMM
   pathway. Any surrogate is labelled a control or it does not exist.

The NN changes the **source of information**. It does not replace the machinery.

---

## 1. WHAT THE COORDINATOR ALREADY MEASURED — READ BEFORE PROPOSING ANYTHING

`s24/LEDGER.md` L1 and L2, both n=126, both complete, both pre-registered. Three results bind every
lane in this sprint.

### 1.1 The directive's mandated screen is inverted. Do not use it as a decision rule.

Directive §15 requires `fcommon` as a screen and calls a lower value promising. **Every alternative
source measured has a lower f than the incumbent and every one is worse**; a blind library draw scores
f = 0.4896 against the incumbent's 0.6758 and is **0.76 Å worse.** `f = |ebar|²/(|ebar|² + mean|d_k|²)`
falls whenever a source is merely more scattered, and `|ebar|² = n·RMSD²` exactly, so the numerator is
not independent evidence either. **Report f descriptively. Decide on the bias cosine (§1.2).**

### 1.2 The biases are ~31% independent, and the generator therefore has a hard quality bar.

    cos(incumbent bias, blind-library bias)        +0.6467
    cos(library bias, an independent library draw) +0.9330   <- the WITHIN-SOURCE control
    ratio                                           0.693

Confirmed geometry-free: fitting the quadratic error model to the matched-size mixture curve recovers
`c_eff = 0.655` against the directly measured `0.647`, and the curve bows **0.226 Å below** the
straight line between its endpoints. Parallel biases sit on that line. **The angular headroom is real.**

Inverting the same model, a second source pays at small mixing fraction exactly when its own error
ratio `q < 1/c_eff`, median **q\* = 1.274**:

> **THE SPEC. A learned source must reach ≤ ~3.9 Å standalone mean on the 126-target instrument with
> a bias cosine ≤ ~0.65 against the incumbent pool merely to BREAK EVEN, and must be materially
> better than that on at least one axis to move the mean. Diversity on its own is already available
> for free, already measured, and is worth +0.76 Å.**

### 1.3 Pool composition at fixed score and fixed readout is a FLAT lever.

Directive §17's central experiment, run with a zero-cost stand-in source: merge the shipped K=500 with
500 uniform library windows, score all 1000 with the **same** functional, take the **same** top-75.

    union vs incumbent  +0.0022  SE 0.0174  MDE 0.0487  fold[-0.027,+0.031]  61W/65L   NULL
    share of the merged top-75 taken from the NON-retrieved half:  0.355

**The scorer already takes a third of its set from candidates retrieval never proposed, and nothing
moves.** So: *the generator must beat the pool, not merely differ from it.* A model whose samples the
existing score likes about as much as what it already has will land on +0.002 Å.

---

## 2. DO NOT REDO — closed, with evidence

| closed | evidence |
|---|---|
| Any router over per-target quantities | 5 constructions, ~0% capture; finite-sample bound needs 280–1045 targets at n=126 |
| Tightening the retained set | clustering, medoid, sharp weights, small m — all at or worse; outliers are load-bearing (+0.142 Å to remove 4 deviant vs +0.005 at random) |
| Weighted averaging | leakage-oracle grid selects the UNIFORM point; no headroom exists |
| Global or per-target rescaling | nested CV +0.0002; ceiling −0.3403 is a function of the invisible common mode, unreachable in principle |
| Probability-weighted CVaR readout | 4 temperatures, never beats the size-matched uniform bar, significantly worse at 3 |
| Deeper ansatz / larger χ / other optimisers / encoding | layers 1→20, 5,760 encoding cells, none beats best-of-N |
| Cα-restrained or unrestrained AMBER repair of the average | 17 of 17 settings at or worse than no repair; the ladder bottoms out at the no-op |
| Naive torsion averaging | the prior test used a genuine circular mean; +1.024 Å is a result, not a bug |
| Diversity selection within the existing pool | predicted order held, all effects < 0.25× their own MDE |

**Anything on this list needs a fundamentally new mechanism, not a retry.**

---

## 3. LANES

Compute discipline first: **the box is 15.6 GB and 6.43 core-equivalents; it fits two heavy jobs.**
Never run more than two training/AMBER jobs at once. The AMBER/OpenMM lane is **serialised** — one
holder at a time, announced on open and on release.

**Lane A — data, corpus and leakage.** What corpus actually exists, what a generator may legitimately
train on, and proof that it is clean. Cluster separation stricter than random splitting; sequence
homology, PDB overlap, fragment overlap, source-organism overlap against all 126 dev targets **and**
against the sealed benchmark's identifiers without opening its results. Deliver the trainable corpus
as a materialised, versioned artefact with the exclusion audit attached. **Nothing else starts until
this lane's exclusion list exists.**

**Lane B — the residual generator (PRIORITY).** Directive §21/§39. `retrieved candidate + learned
stochastic torsional residual (Δφ, Δψ) → new candidate`. This is the lane with the strongest prior:
§1.2 says the corpus and retrieval already share 69% of their bias, so a model that RESAMPLES the
corpus is predicted to inherit it, while a model that learns a **correction** is not bound by that
argument. Sine/cosine residual parameterisation; conditioning on sequence embedding, local window,
distogram, and the retrieved structure itself. Must emit **many** samples per candidate, and must be
audited for mode collapse before its RMSD is ever quoted.

**Lane C — the from-scratch generator.** Directive §7. `p(φ,ψ | sequence, distogram)` → sample →
build. Compare representations empirically (torsions vs sin/cos vs latent), compare model classes on
feasibility and sample diversity rather than prestige, and start with the smallest model that shows
genuine multimodality. **Its first deliverable is not RMSD — it is the §1.2 spec: standalone mean and
bias cosine.** Report those two numbers before anything else.

**Lane D — quantum integration and Hamiltonian roles.** Whether a changed candidate manifold changes
what CVaR-VQE can contribute; the classical-ranking control on every promising pool; whether Legacy
and AMBER disagreement predicts useful complementarity (§28); AMBER strictly as feasibility filter,
energy measurement and constrained diagnostic — never as unconstrained refinement (17/17 closed).

**Lane E — statistics, replication and audit.** Independently re-derives every promoted number from
the artefacts, not from the claiming lane's code. Owns the fork enumeration for other lanes' primaries
(Rule 0 requires the enumerator to have no stake — the coordinator's own two Tier-1 files are marked
as violating this and want re-checking here).

---

## 4. HARD RULES

- **MDE = 2.8016 × SE, per comparison.** Report SE and the effect as a multiple of its own MDE.
  0.7–1.3× is the Type-M zone and is not a result.
- **Paired, fold-clustered CIs beside iid.** W/L on everything. Report worst-target degradation.
- **Basis discipline.** Point-cloud and built-chain RMSD are never compared. State the basis on both
  sides of every contrast. The gap is 0.156 Å of pure operator choice.
- **ORACLE vs ACHIEVABLE vs PRODUCTION**, labelled at every appearance. An oracle number is never a
  system result.
- **Best-of-K nulls are the distribution of the MAXIMUM.** A "−0.077 Å oracle" last sprint was 101%
  accounted for by its own best-of-K null.
- **A control must match the operator's space**, and a zero-information control must be **plausible**,
  not degenerate.
- **Nested CV for anything fitted.** A hyperparameter fitted on the folds it is scored on is leakage
  and must be reported separately as the optimism if shown at all.
- **Rule 0 — operator forks.** For any directional hypothesis, enumerate SIX axes in the module
  docstring and NAME THE ALTERNATIVE NOT TAKEN: functional, basis, readout, normalisation, null, and
  **THE LABEL**. Send the fork list to the coordinator BEFORE the run. Enumeration by someone with a
  stake is a weakening and must be declared as one.
- **Leakage.** The NN never trains on the 126 dev natives, never on benchmark structures, never on a
  loss that consumes native RMSD, native contacts, native compactness or native scale. "The sequence
  is not identical" is not evidence of no leakage — audit homology.
- **Mode-collapse audit precedes any RMSD claim** from a generator: unique-structure count, torsion
  entropy, pairwise RMSD, effective sample size, duplicate fraction.
- **Geometric validity is reported with every generated pool:** Ramachandran, ω, bond length, bond
  angle, clashes, chain breaks. RMSD bought with chemically absurd structures is not bought.
- **Completion flags demand the full key set**, not a row count. Atomic writes (tmp + `os.replace`),
  config-derived filenames, seeds from `s15.seed.stable_rng`.
- **Ablate anything that wins.** A win without an ablation is not promoted. Replication on independent
  seeds, full 126, pinned folds, and an implementation audit by Lane E.

## 5. LOGGING
experiment id · git commit · seed · target ids · config · corpus version · exclusion-list hash ·
model class · parameter count · representation · training objective · samples per target ·
Hamiltonian · ansatz · qubits · optimiser · α · shots/exact · ensemble size · post-processing ·
AMBER settings · mean/median/worst RMSD · SE · MDE · iid and fold CI · W/L · per-target delta.

## 6. THE STANDARD
If <3.0 Å cannot be reached honestly, **do not manufacture a win.** Identify the closest defensible
improvement and state exactly what remains limiting. A tiny numerical fluctuation is not a
breakthrough, and an elegant architecture that does not move the endpoint is not a result.
