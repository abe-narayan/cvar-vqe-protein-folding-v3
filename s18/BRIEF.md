# SPRINT 18 — SHARED AGENT CONTRACT

Read this in full before running anything. It is binding.

---

## 1. The sprint question

> **Does the degree-1 truncation of the distance objective actually solve the objective-alignment
> problem on the real 126-target continuous-torsion instrument — or was the 19-target result a
> small-instrument artefact?**

Targets: **< 2.5 Å** mean full-chain Cα-RMSD first, **< 2.0 Å** preferred. Incumbent **3.204 Å**;
the best structure the pipeline actually builds is the coordinate average at **3.048 Å**.

---

## 2. READ THIS BEFORE ANYTHING ELSE — Phase 0 is done, and it weakens the premise

`s18/PHASE0.md`. The 19-target numbers **reproduce exactly** (full 2.661, degree-1 2.411, space
best 1.030, weight-1 variance 0.613, cumulative ≤2 0.930). **But the inference drawn from them was
never statistically supported:**

    deg1 - full   mean -0.249   median +0.000   W/L 7/5   95% CI [-0.650, +0.093]

The interval includes zero, the **median is exactly zero** (7 of 19 targets share an identical
argmin), and the mean is carried by three targets. And there is a tension pointing the other way:
**ρ(objective, RMSD) is +0.264 for the full objective and +0.153 for degree-1** — the truncation is
a *worse* global correlate of RMSD even though its argmin is better on average.

**So degree-1 enters this sprint as a hypothesis at the noise floor of its own instrument, not as a
result to be translated.** Falsifier F4 is live before we start. Do not treat 2.411 as established.

---

## 3. The finding this sprint exists to act on

Sprint 17's closing measurement (`s17/refine.py`, n = 126): starting from the coordinate average
(3.048 Å) and refining toward the deployed distance objective,

| | |
|---|---|
| objective | **192.8 → 56.5 (−71%)** |
| RMSD | **3.048 → 3.610, +0.561 [+0.407, +0.712], 31W/95L** |
| window shape (w = 4/6/8) | within 0.035 Å of the full fit — the enumerated instrument's 0.7 Å local-refinement gain **does not transfer** |
| shuffled-distogram control | +2.346 |
| matched-magnitude random move | +1.499 |

> **The objective's optimum is 0.56 Å worse than the structure the pipeline already builds. The
> pipeline works because it does not optimise its own scoring function.** The objective is not
> weak — it is *mis-aimed*.

That is why an objective experiment is the only branch left: selection is closed at five
independent levels, retrieval is exhausted (+1.5 targets over a random 500), the readout is optimal
(averaging beats alternatives 39/39), and the repair tax (0.155 Å) is irreducible.

---

## 4. The mathematical object, and the trap in defining it

**Do not port a lattice Walsh truncation into continuous torsion space by analogy.** The rigorous
bridge is the **first-order functional ANOVA decomposition under a product reference measure**
μ = ⊗ᵢ μᵢ:

    E₀     = E_μ[E]
    fᵢ(θᵢ) = E_μ[E | θᵢ] − E₀
    E_{≤1}(θ) = E₀ + Σᵢ fᵢ(θᵢ)

> **CORRECTED 2026-09-06 — the sentence that stood here was FALSE, and it was the coordinator's
> error, not a workstream's.** It read: *“When μ is uniform on the enumerated lattice this **is**
> the Walsh weight-≤1 projection.”* It is not. Verified EXACT at 1e−12 by ADVERSARIAL and
> independently by MATH: the residue-additive ANOVA object equals the projection onto Walsh
> coefficients supported **inside one residue** — weight 0, weight 1, **and intra-residue weight 2**.
> The two objects differ by 1.4 objective sd. Point 1 below already anticipated the reason; the
> equivalence claim contradicted it and should never have been written. Anything derived from the
> equivalence is void; the residue-additive (RA) object is the gauge-invariant one and is the object
> the sprint should have been porting. See `s18/LEDGER.md` L5.

Three things to be careful about:

1. **Per-qubit ≠ per-residue.** The lattice encodes k = 4 states in 2 qubits, so a *per-residue*
   field needs weight-1 **and intra-residue weight-2** terms. Sprint 17 measured 95.7% of the
   weight-2 mass as intra-residue — meaning the strict weight-≤1 truncation **discards part of the
   per-residue field itself**. The proper "degree-1 in residues" object and the strict Walsh
   weight-≤1 object are *different*, and both should be built.
2. **The choice of μ is a modelling decision and must be native-free.** Candidates: uniform over
   the torsion library; the **retrieval pool's empirical per-residue marginal** (native-free and
   target-conditioned — the most informative legitimate choice); a generic Ramachandran prior.
   Report sensitivity to it.
3. **Computability.** E = Σ_p w_p (d_p(θ) − d̂_p)², and by the recorded locality theorem d_ij
   depends on exactly the j−i−1 residues between i and j. So E[E|θᵢ] decomposes over pairs and is a
   straightforward Monte Carlo over the other torsions.

---

## 5. Objective quality, structural quality, alignment, and search are FOUR different things

Sprint 17 proved they dissociate: a 71% objective reduction bought +0.561 Å of RMSD. For every
objective report all four separately:

- **objective quality** — does optimisation lower the mathematical objective?
- **structural quality** — is the resulting RMSD lower?
- **alignment** — does lower objective correspond to lower RMSD? (name the axis: pairwise ordering
  accuracy, Spearman, Pearson, or copula ρ — they are not interchangeable)
- **search quality** — can classical/quantum optimisation find the optimum, and at what budget?

---

## 6. Reporting conventions — standing law

- **TARGET is the unit.** Paired fold-clustered bootstrap CIs, medians and win/loss beside every
  mean, per-target distributions, stratified by length and fold. Never the mean alone.
- **Every arm carries both controls**: a **zero-information reference** and a **matched-random**
  operation of the same magnitude/count.
- **ORACLE CEILING / REALIZED / GAP** reported separately and never mixed.
- **Claim labels**, never blurred: ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN · REFUTED ·
  **EXACT** (a theorem, not a discovery) · **ORACLE** (needs the native).
- **No min-of-N ceiling without its min-of-N null, and the null must name its band.**
- **No directional conclusions from smoke runs.** Record the smoke as a smoke and wait.
- Instrument note: on the 126-target instrument the **minimum detectable effect at 80% power is
  0.084 Å**. A null below that magnitude is uninformative, not negative.

---

## 7. Hard rules

**Sealed benchmark.** 60 targets. Do not read, probe, derive, or tune against it. Verified clean in
Sprints 15–17; keep it that way.

**Frozen metric.** Full-chain Cα-RMSD, existing implementation, all residues including termini, no
trimming, proper rotations only, model 1 of the native. No alternative metric may be reported as
primary.

**Native information** is for evaluation, labelled ORACLE diagnostics, auditing and post-hoc
interpretation only. Native RMSD may be a **training label inside training folds only**; inference
must be native-free. Anything selected with native information is ORACLE and is not a predictive
result.

**Seeding.** `s15/seed.py`'s `stable_rng` / `stable_seed`. Bare `hash()` is salted per process and
has destroyed results here before.

**Legacy weights** are `DEFAULT_WEIGHTS` and are never fitted (leave-fold-out fitting moved the
certified optimum +0.78 Å [+0.04, +1.50] the wrong way).

**AMBER discipline.** Convergence gate (final energy finite, ≤ 1000 kcal/mol) declared before use
and reported with its exclusion count. The rotated-frame null reported with its **maximum**, not
only its mean. A gated arm must be compared to **its own gated input** — comparing it to the
ungated instrument-wide input is how Sprint 17 nearly published a 10× overstatement.

**Artefacts.** Unique immutable result directory per experiment, deterministic seeds, completion
flags, persisted objective coefficients and config hashes. Never treat a partial file as complete
(three workstreams once read a 60/126 artefact as final). Never overwrite a prior sprint's files.

---

## 8. Standing refutations — do not repeat

Jacobian/coordinate steering · torsion interpolation · disagreement features · target-level
calibration on these feature classes · AMBER as a ranker · Legacy argmin ranking · blindly widening
K · generic rankers over the existing signals · optimising the *current* distance objective harder.
Revisit only with an explicit statement of what new hypothesis makes the rerun different.

**Known roles** (hypotheses, not immunities): `leg_torsion` is the one gate that raises near-native
recall (+0.159 [+0.073, +0.239]) without damaging the set best; `leg_contact` (Miyazawa–Jennings)
carries **+0.080 [+0.032, +0.128]** of real in-band information yet **selects 0.136 Å worse than
random**; AMBER is stereochemical repair only (ranking ρ −0.027, AUROC 0.499, argmin worse than
random); the coordinate average contracts the backbone **22.4%** and that — not AMBER — causes the
cis-peptide defect (Spearman(spacing, cis) = −0.949).

---

## 9. Falsifiers for the whole sprint

**F1** degree-1 lands ≈3.6 Å on the real instrument · **F2** degree-1 is no better than the full
objective · **F3** the advantage dies under matched controls · **F4** the effect exists only on the
19-target instrument · **F5** the continuous mapping is mathematically invalid.

If any fires, **close the branch**. Do not rescue it with tuned weights, do not tune against RMSD,
do not expand the λ ladder beyond what was pre-registered.

**A clean falsification is a successful sprint.**

---

## 10. Compute

8 cores (~6.43 core-equivalents), ~11 GB usable RAM. Target ~90–95% CPU; **never deliberately
exceed 97%**. One heavy process per workstream; check load before launching anything long;
checkpoint per target. Cache deterministic quantities; never cache anything that changes scientific
randomness without explicit seed/version control.
