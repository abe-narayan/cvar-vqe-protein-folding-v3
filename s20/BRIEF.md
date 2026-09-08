# SPRINT 20 — SHARED AGENT CONTRACT
## Landscape geometry, decorrelated generation, and the Legacy–AMBER question

Read in full before running anything. Binding.

---

## 1. THREE IMMOVABLE PRIORITIES

    PRIORITY 1   minimise final full-chain Ca-RMSD
    PRIORITY 2   preserve GENUINE CVaR-VQE
    PRIORITY 3   rigorously understand and compare LEGACY vs AMBER

    best POINT CLOUD          3.048 A   coordinate average -- NOT a buildable structure
    best BUILT STRUCTURE      3.204 A   lambda=0 projection -- THE INCUMBENT
    deployed emission         3.236 A   after restrained AMBER
    primary target          < 2.500 A
    aggressive target       < 2.000 A

> **CORRECTED 2026-09-07, FIRST HOUR — the label "best built 3.048" was WRONG and it was the
> coordinator's, propagated from Sprint 19.** The coordinate average has a mean virtual Ca-Ca bond of
> **2.961 A against the physical 3.804 A -- a 22.2% contraction, with a global minimum bond of
> 0.649 A** across the 126 targets. **It is a shrinkage estimator of Ca positions, not a backbone**:
> no side chain can be placed on it and AMBER cannot score it. Verified independently by the
> coordinator against `bench_results/cache/1fc9f2dcf489e2fb`.
>
> **Projecting it onto the pipeline's own ideal-geometry manifold gives 3.204 A. So the incumbent
> already IS the best structure the pipeline builds**, and the 0.156 A gap is the price of
> stereochemical realisability (1.9x MDE). This is **not** a metric exploit other arms could copy:
> the best constant global rescale of `fit_ca` buys only 0.034 A, below MDE.
>
> **EVERY ARM MUST STATE ITS BASIS.** point cloud / built chain / repaired emission. Comparing a
> built chain against the 3.048 point cloud is a basis mismatch, and it is a **recurring** error in
> this project (`rmsd-reporting-basis-mismatch`).

**MDE at 80% power on the 126-target instrument: 0.084 A.** A null below that is uninformative, not
negative.

> **RMSD FIRST.** A method with beautiful Hessian spectra, lower energy, better convergence, higher
> validity or stronger quantum metrics **and worse RMSD is not an improvement.** Report it as such.

---

## 2. WHAT SPRINT 19 ESTABLISHED — start here, do not re-derive

**The mechanism.** The harm in the distogram's error is its **cross-pair sign coherence**. At
residual RMS matched to machine precision, destroying only the sign pattern is worth **−1.202 A
[−1.408, −1.004]**, 128% of the Sprint-18 gap. A residual borrowed from a *real, perfectly realisable
alternative structure* at matched magnitude is the **worst arm on the board**.

**Why.** The weighted fit projects onto the tangent space of the ideal-geometry manifold, rank
**2n − 5** (EXACT, verified: mean rank 20.92 vs 20.92; generic in-tangent fraction predicted 0.329,
measured 0.330). Real field **0.695** in-manifold vs **0.397** for a sign-randomised twin.
**The fit denoises the orthogonal component and has no power against the in-manifold one.**

**The source.** The harmful coherent component is **shared**, and the shares must be quoted
**null-subtracted** — the alignment statistic has a **shared-referent floor of 0.505** because both
arguments are scored against the same native (L5). Raw -> above-null share: same-seed ceiling 0.898 ->
1.00; same-arch 0.833 -> **0.83**; **retrieval pool 0.784 -> 0.71**; cross-architecture 0.731 ->
**0.57**; sequence-blind references 0.575-0.598 -> **0.18-0.24**. The *incoherent* component does not
share (0.00 against zero-information). **The ORDERING is what the mechanism needs and it survives
intact; "zero-information reproduces two thirds" does NOT — above the null the component is mostly
sequence-DEPENDENT.**

**The wall.** The harmful mode is the **per-target separation profile** (+52.5% of the gap). The best
native-free estimator of it *is* the retrieval pool — which shares 0.87 of the component being
estimated. **The estimator is made of the bias.**

**Selection.** Gate damage is **exactly** `readout² = ‖b_pool‖² + 2‖b_pool‖·ALIGN + ‖Δ‖²`. Score gates
raise ALIGN (+0.0209), score-free gates do not (−0.0048). **Any score preferring compact,
well-formed, pool-typical geometry selects toward the pool's own error** — which is why a
zero-information constant-helix gate is the *worst* arm (+0.1253, twice Legacy's). **But the
target-conditioned distogram filter beats matched-random by +0.401** — "any pool-typicality ordering
is harmful" is right; "any score ordering is harmful" is wrong.

---

## 3. WHAT IS CLOSED — do not reopen without a materially new mechanism

Predictor architecture and training loss (8 models, none past MDE, ORACLE MAE −11.7% for nothing) ·
joint consistency / metric realisability (more realisable ⇒ more coherently wrong) ·
distribution-shape consumption (mode information real, **zero transfers**, split-half inside a fold)
· physics ranking and steric rejection (every physics gate loses to a random gate) · Legacy as
objective term (+0.722, harmful with a sign) · `leg_torsion` as a gate (converts negatively) ·
softening the top-*m* cut (best T −0.028, a third of the MDE) · training "in the right basis" ·
degree-1/Walsh truncation · objective tempering · Jacobian steering · blindly widening K.

**Standing roles — AMBER's CORRECTED 2026-09-07, and the correction inverts it.** Turning the point
cloud into a valid backbone costs **+0.166 [+0.130, +0.202]**, 16W/110L, 5/5 folds; **that is the
PROJECTION's tax and it is irreducible** (s18 §6.2, whose title already said *"the AMBER pass is never
spent"*). **AMBER's own cost, on the deployed path** (`amber_k = 10.0`, `amber_steps = 0`, applied to
the **built chain**), is **+0.021 [+0.014, +0.028]** — *a quarter of the MDE*. And as a repair
operator on the point cloud, **AMBER at k=30 BEATS the projection by −0.022** [−0.036, −0.009]. All
four verified by the coordinator against the production cache at n=126.

> **The earlier clause hung the projection's number on AMBER**, handing every lane the prior that
> AMBER is expensive when the deployed operator costs a quarter of the MDE and is the *better* of the
> two repairs. This matters for Q1 specifically, which warns against assuming "AMBER is more
> physical" explains its behaviour — a lane primed with "AMBER costs 0.164" is primed wrongly.

Legacy is a
**diagnostic/labeller, never a selector**; m = 75 sits in a broad two-dimensional optimum; the
coordinate average beats alternatives 39/39 and is **not** a privileged start.

---

## 4. THE SPRINT'S TWO CENTRAL QUESTIONS

### Q1 — Is the AMBER landscape genuinely harder, and if so *what property* makes it harder?

**Do not assume "AMBER is more physical" is the explanation.** Separate experimentally:

    objective roughness | optimizer mismatch | CVaR concentration | ansatz limitation
    encoding limitation | higher-order coupling | poor initialisation | candidate-distribution bias

> **CORRECTED 2026-09-07 — Q1 CANNOT BE RUN AS FIRST BRIEFED, and the coordinator wrote it wrong.**
> "Change only `H_Legacy` <-> `H_AMBER`" presumes two energy functions on one landscape. **They are
> not.** Verified in source by the coordinator:
>
>     core/quantum.py  FoldingHamiltonian.energy(bitstring)  -> build coords, evaluate. NO relaxation.
>     amber_hamiltonian.py AmberHamiltonian.energy(bitstring) -> build coords, THEN
>         set k_rest = 100 kcal/mol/A^2 positional restraint on every heavy atom,
>         openmm.LocalEnergyMinimizer.minimize(...),
>         set k_rest = 0 and report the UNRESTRAINED energy at the MINIMISED coordinates.
>
> **So `H_AMBER = E o Relax`.** Swapping the Hamiltonians also swaps in a **relaxation operator**, and
> **every landscape metric Q1 names is confounded by it** — roughness, Hessian condition number,
> gradient noise, basin width. A relaxation *smooths* over most of the domain **and** injects
> *discontinuities* where the minimiser's basin assignment flips; both are properties of the operator,
> not of "AMBER physics". `_is_collapsed` additionally returns **+inf**, so `H_AMBER` has an
> infinite-valued region `H_Legacy` does not — `n_collapsed` is counted and **must be reported**
> (rule 3: report how many times the guard fired).
>
> **THE MISSING ARM is `E_amber o Relax_1`** — the relaxation reduced to its smallest honest setting
> (OpenMM's `maxIterations = 0` means *unbounded*, so **1**, not 0, is the floor). That is the only
> true "change only H" arm.
>
> **UPDATED 2026-09-07 after F-D2 fired on both clauses (L6): the cap binds on 192/192 calls, and
> AMBER's energy is NOT FINITE on 42% of the register without the relaxation.** The relaxation is
> **constitutive, not cosmetic** — without it there is no objective on most of the space, which kills
> the simple "add a `Relax_1` arm" fix. Q1 is therefore re-specified: **(b) PRIMARY — give `H_Legacy`
> the SAME `Relax_50` preprocessing, so the relaxation is a shared operator and cancels; this is the
> only genuine "change only H" and costs one wrapper. (a) DIAGNOSTIC — restrict both to the 58%
> subset where `E_amber o Relax_1` is finite and run three arms, so `Relax_50 - Relax_1` isolates the
> relaxation's own effect. (c) Otherwise label the conclusions RELAXATION-CONFOUNDED.**
>
> **And note before designing anything: `Spearman(E_legacy, E_amber o Relax_50) = +0.229`, of which
> about a third appears only after the relaxation** (`Relax_1`: +0.140). **The two Hamiltonians barely
> agree on ordering at all.**

> **Whoever holds Q1 must adopt (b), or (a), or take the label.** Pre-registered by WORKSTREAM D as
> Block E with falsifier F-D2: *if Spearman(E o Relax_50, E o Relax_1) >= 0.95 AND the iteration cap
> binds on < 10% of calls, the objection is withdrawn; otherwise Q1's conclusions are
> RELAXATION-CONFOUNDED and must be labelled so.*
>
> **AND A REPRESENTATION COLLISION**: both Hamiltonians take a **bitstring**. `AmberHamiltonian`
> refuses a lattice representation and requires `TorsionStateRepresentation`; production runs it at
> `torsion_window = 8` — **8 states, 3 bits per residue, not the k=4 lattice.** So Q1 as briefed is a
> **discrete-register** experiment and collides with section 5's continuous-torsion rule. **There is no
> continuous AMBER objective in the codebase**, and building one does not escape the relaxation:
> `energy_from_coords` goes through the same `_evaluate`.

Hold constant: sequence, torsion representation, candidate set, ansatz, qubit count, initialisation,
optimiser, evaluation budget, CVaR rule, measurement budget, convergence threshold — **and the
relaxation depth, which must be stated for every arm.**

**Every landscape metric must be connected to RMSD** (§42 of the directive). Hessian condition
number, negative curvature, gradient noise, basin width, CVaR concentration — for each, ask *does it
predict final Ca-RMSD?* **A landscape metric with no predictive relationship to structural outcome
must not be over-interpreted.**

### Q2 — Can any generator's error DECORRELATE from the retrieval pool's?

This is Priority 1's only live lever. Sprint 19 closed every route that *consumes* the pool better;
the pool is the bias, and it is simultaneously generator and training source.

**For every candidate source measure `rho(e_coherent, e_pool)` — not member RMSD.** A generator is
valuable if **quality is acceptable AND error correlation is low**, even when its mean member RMSD is
not the best available, because the incoherent component does not share and averaging suppresses it.

**Measure error diversity separately from geometric diversity.** Two populations can be
geometrically different and make the same structural mistake.

---

## 5. HARD RULES

**Sealed benchmark.** 60 targets. Do not read, probe, derive from, tune against, infer properties of,
or use for model selection. **Do not open `results/benchmark_manifest.json`.** Clean through Sprint
19; keep it that way.

**Frozen metric.** Full-chain Ca-RMSD, existing implementation, all residues including termini, no
trimming, proper rotations only. New diagnostics welcome and must be labelled as diagnostics.

**Continuous torsion representation.** Do **not** silently regress to the k=4 lattice or a binary
encoding because an experiment is easier there. If you use any discretisation, **prove
gauge-invariance of your conclusion or report it as encoding-dependent** — Sprint 18's headline was
an arbitrary 2-bit encoding artefact.

**Genuine pillars.** Genuine VQE, genuine CVaR, genuine Legacy (`DEFAULT_WEIGHTS`, never fitted),
genuine AMBER. Never a classical imitation described as quantum. **Their role is determined
experimentally**; a demonstrated non-contributor is an acceptable outcome, a false positive is not.

**Native information** is for evaluation, labelled **ORACLE** diagnostics and post-hoc interpretation
only. Native RMSD may be a training label **inside training folds only**.

**Seeding.** `s15/seed.py`'s `stable_rng`/`stable_seed`. Bare `hash()` is salted per process and has
destroyed results here.

**The subset trap.** `s15/distcal.fit_correction` fits the leave-fold-out debias on whatever target
list it is handed. Always gather and fit on the full `I.targets()`, then iterate your subset.

**Artefacts.** Unique immutable directory per experiment, deterministic seeds, completion flags,
persisted configs and hashes. Partials named `_PARTIAL_*` and quoted at their own n or not at all.
Superseded runs named `_SUPERSEDED_*` and never quoted. Never overwrite a prior sprint's files.

**Cached starts.** `s19/cache/start_<pdb>.npz` holds the deterministic projected start, identical for
every arm and lane. Read it; do not recompute `I.project` per target.

---

## 6. CONTROLS — and the five rules this programme learned the hard way

Every mechanism carries a **plausible zero-information null** and a **matched-random** operation of
the same magnitude/count.

1. **A control must be matched in the space the OPERATOR actually works in.** Five instances in two
   sprints, none self-caught: weights permuted with residuals; offsets orphaned from weights;
   isotropic in *raw* rather than *whitened* space; free-superposition mixed with common-frame; and a
   2.0 A **clip** that silently moved a "matched-magnitude" arm *after* the match. **The clip, the
   whitening and the projection are all part of the operator.**
2. **When an analytic null and a measured null disagree, the measurement is the null.** Two algebraic
   nulls were offered and wrong; a review of the derivation would have passed both.
3. **A gate can pass VACUOUSLY.** Report how many times the bound, guard, threshold or fallback
   actually **fired**; a pass with zero firings is not evidence. **And a bound must be verified on
   the calls it actually bound, not only on the ones it did not.**
4. **Uniform is not a zero-information control.** Uniform-on-the-torus places mass on impossible
   backbone conformations, making it a *worse* measure rather than an uninformative one. Use
   plausible-but-uninformative nulls (constant alpha-helix, matched empirical marginals, label
   permutations). **Note the role-dependence**: as a *reference* the helix is a null; as a *gate* it
   is the worst arm measured.
5. **Pre-register the FALSIFIER, not the hypothesis.** The control that refuted the coordinator's
   briefed mechanism in Sprint 19 existed *only* because a falsifier forced the lane to name in
   advance what would refute itself.

---

## 7. PRE-REGISTRATION FORMAT — mandatory before any large experiment

    Hypothesis:
    Expected mechanism:
    Primary endpoint:
    Falsifier:
    Null:
    Matched control:
    Budget:
    Promotion criterion:

**Do not move the goalpost after seeing the data.** If a pre-registration turns out mis-specified,
say so, keep it unedited, and record the untested regime as OPEN. Restoring a pre-registration is
more expensive than explaining a deviation and is the correct choice.

---

## 8. STATISTICS

**TARGET is the unit.** Paired target-level bootstrap CIs, fold-aware, medians and win/loss beside
every mean, stratified by length and fold. Never a mean alone.

**A zero-spanning CI with insufficient power is NOT MEASURED — never "matched."** A Sprint-19 "match"
published at n=14 reversed sign by n=30. **A result below the 0.084 A MDE is not a validated
improvement**, and an argmin over a grid on the tuning instrument is a hyperparameter chosen there.

No min-of-N ceiling without its min-of-N null, naming its band. No directional conclusions from
smokes.

---

## 9. LABELS

**EXACT** (theorem/identity, not a discovery) · **ORACLE** (needs the native) · **ESTABLISHED** ·
**SUPPORTED** · **PLAUSIBLE** · **OPEN** · **INCONCLUSIVE** · **NOT MEASURED** · **REFUTED** ·
**RETRACTED**.

**Never present an algebraic identity as an empirical discovery.** For every statistic ask first
whether it is an identity, a theorem, an implementation consequence, an observation, a causal
hypothesis, or a learned relationship. **Derive the operator before interpreting its statistic.**

---

## 10. READ THE CODE BEFORE BELIEVING THE CLAIM

The source has repeatedly contradicted the table, and this applies to **performance claims exactly as
to scientific ones** — Sprint 19 recorded a "pathological minimiser" that was CPU starvation from
concurrent lanes, amplified into three documents without an independent check.

For every decisive result inspect: argument order, weights, masks, seeding, frame conventions,
normalisation, whether the random control is genuinely matched, whether "same budget" is the same
budget, and whether evaluation is target-level.

---

## 11. COMPUTE

8 cores (~6.43 core-equivalents), ~11 GB usable RAM, and **the box is shared**. Target ~90-95% CPU;
**never deliberately exceed 97%.** One heavy process per workstream; cap BLAS threads
(`OMP_NUM_THREADS=1` etc.); check load before launching anything long; checkpoint per target.

**Timing measured under contention is not a property of the code** — rule 10 above exists because
that mistake was made. Measure performance on a quiet box or not at all.

---

## 12. WHEN TO KILL A BRANCH

Immediately downgrade when it loses to a plausible zero-information null · works only on the tuning
subset · has a CI spanning the required effect · dissolves under target-level analysis · is
reproduced by a simpler classical baseline · needs native information · rests on an arbitrary
encoding · improves an intermediate metric while worsening RMSD · needs post-hoc tuning · has a
mathematically invalid mechanism.

**A clean falsification is a successful result. Interesting is not enough.**
