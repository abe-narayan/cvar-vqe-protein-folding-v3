# SPRINT 21 — SHARED AGENT CONTRACT
## CVaR-VQE as the selector · Legacy and AMBER as competing Hamiltonians

Read in full before running anything. Binding.

---

## 1. THE MISSION AND THE LADDER

    PRIORITY   minimise final mean full-chain Ca-RMSD
    ARCHITECTURE (non-negotiable)   genuine CVaR-VQE as the SELECTOR
    HAMILTONIANS (non-negotiable)   genuine H_Legacy and genuine H_AMBER, always separably evaluable

    best POINT CLOUD        3.048 A   coordinate average -- NOT a buildable structure
    best BUILT STRUCTURE    3.204 A   lambda=0 projection -- THE INCUMBENT
    deployed emission       3.236 A   after restrained AMBER
    primary target        < 2.500 A
    stretch target        < 2.000 A

> ### MDE IS AN OPERATOR, NOT A CONSTANT — CORRECTED 2026-09-07 (L8), and this corrects four sprints
>
>     MDE = (z_.975 + z_.80) * SE = 2.8016 * sd(paired differences) / sqrt(n)
>
> **`sd(paired differences)` is a property of THE COMPARISON, not of the instrument.** Sprints 18, 19,
> 20 and this brief all quoted **"MDE = 0.084 A"** as a universal bar. Measured across **26 real
> comparisons** from this sprint, Sprint 20 and the production cache, **the per-comparison MDE ranges
> from 0.11x to 9.54x that constant — a factor of 84 — and only 12% fall within 1.5x of it.**
>
> **It errs in BOTH directions and both errors are already in this programme's record:**
>
> * **Too LARGE for low-variance comparisons.** The deployed AMBER repair (n=126) has paired
>   sd 0.0385, **SE 0.0034, MDE 0.0096 A — one ELEVENTH of the quoted bar.** Its effect of +0.0207 is
>   **6.09 SE, power 1.000.** §3 of this brief calls that *"a quarter of the MDE"*, **which invites
>   reading a six-sigma effect as negligible.** It is one of the most precisely measured quantities in
>   the programme.
> * **Too SMALL for high-variance comparisons.** Sprint 20's encoding cell AMBc/spsa (n=10, 2 seeds)
>   has **SE 0.2860, MDE 0.8014 A — 9.54x the bar** — while reporting an effect of −0.649 as
>   significant. **The study's own 80%-power MDE is LARGER than the effect it detected**: |d|/SE =
>   2.27, post-hoc power **0.62**, Type-M magnitude exaggeration **1.27x**. The LEG/spsa cell is worse:
>   power **0.18**, Type-M **2.39x**.
>
> **THE RULE, binding on every lane: report SE beside every mean, and quote the MDE that comparison's
> own paired sd implies.** Never label an effect against the 0.084 constant. *(The programme derived
> this correctly once — `s16/review_FINDINGS.md` Q4.1 — and it was the CONSTANT that generalised into
> later briefs rather than the METHOD.)*

**EVERY ARM STATES ITS BASIS** — point cloud / built chain / repaired emission. Comparing a built
chain to the 3.048 point cloud is a basis mismatch and it has happened four times in this project.

---

## 2. THE DISTINCTION THIS SPRINT TURNS ON

Sprint 20 measured **argmin-by-energy over the K=500 pool**:

    argmin RMSD   AMBER 4.990 | Legacy 5.487 | distogram 3.676     pool MEAN 4.739
    [n = 20, subset(20), REBUILD basis, squared-distance functional]

**Picking the lowest-energy structure is worse than picking at random, for both physics energies.**
That conclusion is UNAFFECTED by the provenance note below — its four numbers share one n, one
basis and one pool, so they are internally comparable.

> **PROVENANCE CORRECTED 2026-09-07 by WORKSTREAM D — the coordinator quoted this row without its
> operator, which is the FIFTH basis/provenance hazard this project has met and the fourth of mine.**
> These are **n = 20 targets on subset(20)**, on the ideal-geometry **rebuild** of each window's
> torsions, and the "distogram" column is `s20.qb2_lib.Ham("DIST")`, a **squared-distance**
> functional — **not** the deployed selector, which is `s12.instrument.shipped_score`, a Bayes-risk
> score off the full posterior. **Three differences at once, none of them stated where the row was
> quoted.**
>
> **DO NOT BENCHMARK AN n=126 ARM AGAINST 3.676 OR 4.739.** subset(20) is simply **harder** than the
> instrument: pool mean **4.732 vs 4.453**, pool best **1.888 vs 1.711**. The n=126, window-basis,
> deployed-selector row is:
>
>     argmin distogram 3.454 (= s12/instrument's pinned `shipped argmin 3.4540`)   pool MEAN 4.453
>
> **The 3.676-vs-3.454 gap is the SUBSET, not the functional.** D pre-registered the opposite
> hypothesis — that the squared-vs-Bayes functional explained it — and **refuted its own hypothesis**
> at n=126: per-target `rho(E_bayes, E_sq)` median **0.9729**, and `argmin_sq - argmin_bayes =
> +0.050 [-0.019,+0.156]` against its own MDE of **0.124**, a null.
>
> **Consequence adopted for §4:** the matrix's **Distance row is robust to which distance functional
> a lane builds it on** — the two are interchangeable at the argmin. Lanes need not re-run on the
> other one. But a lane mixing them **within a single paired comparison** must still say so, since
> the point estimate moves ~0.05 A in a known direction.

> **CORRECTED 2026-09-07, within the hour, by WORKSTREAM D — the coordinator's framing named an
> operator the pillar does not have.** The paragraph here claimed a CVaR-VQE emits a **tail mean**.
> **It does not.** Verified in source (`core/quantum.py`): the run returns `vqe_bitstring` (argmin
> over the final distribution's samples), `vqe_modal_bitstring` (the mode) and `best_seen_bitstring`
> (argmin over everything seen). **None is an average of a tail. CVaR is the TRAINING objective; the
> READOUT is an argmin.**
>
> **And the consequence is sharper than a correction.** Q6 is EXACT: CVaR's minimiser is a **face
> supported on the alpha-tail**, and *the alpha-tail of H contains H's global pool minimum.* So for a
> **pool-restricted** selector whose **training objective and readout energy are the same H**,
> CVaR-VQE and argmin have the **same optimal answer** — CVaR changes the *sampling distribution*,
> not the *selected point*.
>
> **Scope conditions, which are where the live questions now are**: it needs same-H readout and pool
> restriction, and it **breaks** if the readout is an **average**, or if the **readout H differs from
> the training H**. WORKSTREAM D is formalising it.
>
> **So the sprint's opening question is answerable ON PAPER for the argmin readout**, and what remains
> live is (a) the **averaged** readout — which is the Sprint-19 averaging lever wearing a CVaR label,
> and must be labelled as such — and (b) a **readout H different from the training H**, which is a
> genuinely new and untested design.

**Do not assume the Sprint-20 argmin result transfers to any other readout. Do not assume it does
not. Measure it — and state which readout every arm uses.**

---

## 3. WHAT SPRINT 19-20 ESTABLISHED — evidence, not commandments

**On search.** Six independent instruments: search is saturated. At 8192 evaluations **no sampler,
quantum or classical, beat the zero-evaluation retrieval pool**; across 36 matched-budget cells
**0 improved objective and RMSD together**, rho(optimisation gained, structure gained) =
**+0.003 [-0.068, +0.072]**.

**On the quantum pillar.** Genuine and correct. **It genuinely trains** on the continuous encoding
(-0.210 [-0.339, -0.082], 5/5 folds vs best-of-N from the *untrained* circuit). **Two EXACT limits**:
CVaR's minimiser is a face; and the deployed ansatz is **bond-dimension 4 — a <=16-state HMM,
classically samplable by construction.** Entanglement: -0.013 [-0.095, +0.077], NOT MEASURED.
**Ansatz-seed sensitivity is 0.200 A within-target sd, 2.4x MDE — any variational arm with fewer than
4 seeds is NOT MEASURED.**

**On the two Hamiltonians.** `Spearman(E_L, E_A) = -0.0886 [-0.1239, -0.0514]`, Pearson **+0.0312**
(opposite sign), `cos(grad E_L, grad E_A) = -0.3512 [-0.4719, -0.2221]`, 27/30 negative. **They agree
about catastrophes and disagree about ordering.** What separates them is **compactness**:
Legacy-preferred candidates are **0.45 A more compact** (-0.4477, 124W/2L, 5/5 folds). **Legacy is a
compactness/typicality model wearing a physics vocabulary.** Both ORACLE axes fail to separate —
**complementary about geometry, jointly blind to accuracy.**

**On AMBER's landscape.** Not more negative curvature — **concentrated** curvature: participation
ratio **0.0745 vs Legacy 0.4221 (30W/0L)**, anisotropy 18.8 vs 5.8, condition number three orders
larger, 44% near-zero modes. **~7% of modes carry all of it.** That is a **steric singularity**, and
it is the same object as the collapse census: on 9,450 real rebuilds unrelaxed AMBER is *finite
everywhere* with median **+16,062**, p99 3.1e13, max **5.5e23** kcal/mol, **53.5% above 1e4**.
**Finite-but-meaningless is more dangerous than +inf — nothing throws.**

**And `H_AMBER = E o Relax` in the deployed object.** The 50-iteration cap binds **192/192**;
`Spearman(Relax_50, Relax_1)` = 0.358-0.882; **AMBER is not finite on 42% of the register without the
relaxation.** The relaxation is **constitutive, not cosmetic.** A bare single point is available via
`core.amber._run(steps<0)` — verified bit-exact over 150 comparisons and it never calls the collapse
guard.

**On where the error lives.** Changing the **corpus** does not decorrelate (rho 0.968/0.923 vs a 0.85
bar); changing the **selector** does (rho 0.794). Corpus swap costs +0.109, selector swap **+0.387**.
**The carrier is the distogram-as-selector, not the pool.**

**On the readout.** The projection tax **+0.166 [+0.130, +0.202]** is the price of the **manifold**,
not of un-contracting: a manifold-constrained Frechet mean recovers **+0.0048 [-0.0117, +0.0218]**,
and five different manifold-constrained operators land within **0.016 A** of each other.

**On AMBER's real cost.** Deployed (`k=10`, `steps=0`, on the built chain) it is **+0.021 [+0.014,
+0.028]** — **SE 0.0034, so 6.09 SE at power 1.000, one of the most precisely measured quantities in
the programme** *(the earlier "a quarter of the MDE" gloss was wrong; see §1)* — and at k=30 it
**beats the projection by -0.022**. It is the
*cheaper* of the two repairs. Do not carry "AMBER is expensive".

**On validity.** `caonly_k300` is 0.110 A more accurate with an **indistinguishable clash count** and
**cis_frac +0.4242**. **A scalar validity score would have promoted a broken structure. Validity is a
VECTOR** — clash, cis, Ramachandran, omega, bond/angle.

---

## 4. THE MANDATORY MATRIX

For the same candidate problem, same representation, same budget, same seeds, with **genuine
CVaR-VQE**:

    Legacy | AMBER | Legacy+AMBER | Distance | Distance+Legacy | Distance+AMBER | Distance+Legacy+AMBER

Report for each: **final Ca-RMSD (mean AND median), target-level paired difference, CI, W/L, fold
behaviour** — in that order, before anything else — then objective convergence, gradient and Hessian
geometry, CVaR behaviour, candidate diversity, and the **validity vector**.

**Legacy and AMBER must remain separably evaluable at all times.** The system must be able to report
Legacy-only, AMBER-only, hybrid, and a structural-only baseline. **Do not bury either inside an
opaque learned score.**

> ### EVERY MATRIX ROW MUST STATE ITS READOUT. "Which Hamiltonian is best" is NOT WELL-DEFINED without one.
>
> Measured at n = 126 (L3), against matched-count random controls at alpha = 0.15 — same pattern at
> 0.05 and 0.30, all CIs excluding zero:
>
>     readout        disto                    legacy
>     tail_member   -0.90 [-1.08,-0.73]      -0.41 [-0.56,-0.26]   both BEAT their control
>     tail_medoid   -0.45 [-0.67,-0.24]      +0.30 [+0.20,+0.40]   legacy WORSE, 35W/91L
>     tail_avg      -0.38 [-0.55,-0.22]      +0.33 [+0.21,+0.45]   legacy WORSE, 41W/85L
>
> **Legacy's tail BEATS a random subset if you draw ONE member from it, and LOSES to a random subset
> if you take its medoid or its average.**
>
> **Mechanism — CORRECTED 2026-09-07 (L7), and the first version was wrong.** It is **not** the
> tail's *spread*; it is whether the tail's errors are **COHERENT**. Averaging cancels i.i.d. error
> and **preserves systematic** error (`error-coherence-decides-correctors`: at identical 0.688 sign
> accuracy, coherent mistakes emit +0.31 A and i.i.d. mistakes -0.14 A). Legacy's tail shares a
> **coherent compactness bias** (L2c: 0.45 A more compact, 124W/2L), so **averaging preserves the
> bias** — which is why Legacy costs +0.33 on `tail_avg` while *helping* `tail_member` by -0.41.
> **Diversity is NOT the lever**: maximising set spread is refuted on all three bands, and
> *minimising* it HELPS (-0.214 [-0.397, -0.040], 66W/60L).
>
> **A matrix run at `tail_avg` reports Legacy as harmful; the same matrix at `tail_member` reports it
> as helpful. Both are true.**

> ### AND THE `disto` ROW AT alpha = 0.15 IS NOT A NEW ARM — IT IS THE PRODUCTION PIPELINE.
>
> `tail_avg(0.15)` under the distogram Hamiltonian = **3.0483**, matching the shipped
> `rmsd_avg` at **max |diff| = 0.000000 on 126/126 targets, corr 1.000000**. It must: the shipped
> pipeline *is* "rank the K=500 pool by the Bayes-risk distogram score, keep the top 75 (= 0.15 x 500),
> coordinate-average". And `argmin|disto` = 3.454 reproduces the instrument's pinned
> `shipped argmin 3.4540`. **Both of the instrument's own constants return exactly** — the strongest
> soundness gate available here.
>
> **So a CVaR-VQE that selects an alpha-tail of this pool and averages it is RE-DERIVING THE
> INCUMBENT. Do not report 3.048 as a CVaR result.**

---

## 5. HIGH-VALUE HYPOTHESES, in the directive's own order

1. **Legacy -> AMBER continuation** (`H(lambda) = (1-lambda)H_L + lambda H_A`), flagged as potentially
   the highest-value idea. Measure gradient variance, Hessian spectrum, basin movement, CVaR
   concentration and **RMSD** as functions of lambda. **Do not assume smooth interpolation.**
2. **AMBER preconditioning** — pre-relaxation, bounded/smooth barrier transforms, local surrogates,
   staged or coarse-to-fine AMBER, trust regions, regularised sterics, torsion-space surfaces.
   **If a surrogate is used it must be validated against true AMBER and NEVER called "AMBER".**
3. **Encoding — DOWNGRADED 2026-09-07, do not build on it until the control reports.** Sprint 20's
   `theta` vs `(cos theta, sin theta)` result is **confounded by step count in every cell**: the
   embedded arm takes **0.49-0.77x the iterations** and reaches a **worse objective in 9/12 cells** —
   the signature of an arm that optimised *less*, on an instrument whose most reproduced law is that
   optimising harder makes structure worse. **SPSA's deficit is a `track_quality` DIAGNOSTIC**
   consuming **288/512 budget units embedded against 144 in theta** (d doubles). And **two of the
   twelve cells are the same measurement** — `AMB|nelder` and `AMBc|nelder` are bit-identical because
   `AMBc` is a strictly monotone transform and Nelder-Mead is comparison-only, so the count is
   **11 distinct cells, p = 0.0117 not 0.0063**. *An identity counted as an empirical trial.*
   **Exclude every comparison-only arm from any AMB-vs-AMBc table — those rows are zero by
   construction and are not nulls.**

   **What survives**: at target level over the 11 distinct cells, **-0.2167 [-0.3239, -0.1095],
   W/L 9/1, p = 0.0215**. *The measurement is real; the interpretation is not established.*

   **A mechanism for the confound**: in the embedding, `||u||` is a **pure gauge direction** —
   `arctan2` is scale-invariant, so the true gradient has **exactly zero radial component and half the
   embedded parameters are unphysical**. Per-coordinate normalisation gives a pure-noise direction a
   full-size step, and if `||u||` drifts upward the effective angular step shrinks as `1/||u||` — an
   **implicit step-size annealing schedule**.

   **The deciding control is matched EFFECTIVE STEP COUNT, in both directions** (rule 1: matched in
   the space the operator works in), plus angular displacement and gauge-radius drift as instruments.
4. **Population-aware Hamiltonians** `H(X_i | X_1..X_N)` — consensus, outlier detection, cluster
   membership, disagreement.
5. **Target-specific weights** `lambda(S)` from native-free inputs only, against fixed-global controls.

---

## 6. HARD RULES

**Sealed benchmark.** 60 targets. Do not read, probe, derive from, tune against, or use for any
decision. **Do not open `results/benchmark_manifest.json`.** Prove the seal with a hash, not a read.

**Frozen metric.** Full-chain Ca-RMSD, existing implementation. New diagnostics welcome and labelled.

**Continuous torsional space is the architecture.** Do not regress to a discrete lattice as the main
representation. Any discretisation must **prove gauge-invariance** or be reported as
encoding-dependent — Sprint 18's headline was an arbitrary 2-bit encoding artefact, and Sprint 20's
gauge check (`Z2^n`, percentile 0.500) is the standard to meet.

**Genuine pillars.** Genuine VQE. Genuine CVaR with real tail selection. Genuine Legacy at
`DEFAULT_WEIGHTS`, never fitted. Genuine AMBER ff14SB/GBn2. **The final candidate decision must come
from the VQE/CVaR.** A classical layer may generate, preprocess, diagnose or control — **it must not
secretly replace the selector.**

**Normalisation.** Legacy and AMBER live on wildly different scales. **Measure mean, variance, tails,
gradient scale, Hessian scale, outlier behaviour and numerical range before combining**, and build
physically defensible dimensionless forms. **Never pick a normalisation because it improved RMSD on
the same data without pre-registration.**

**Native information** is for evaluation and labelled ORACLE diagnostics only. Separate **pool
ceiling / selector ceiling / geometry ceiling / repair ceiling** and never present an ORACLE number as
achieved.

**Seeding.** `s15/seed.py`'s `stable_rng`. **4 seeds minimum on any variational arm.**

**Artefacts.** Completion flags; `_PARTIAL_` and `_SUPERSEDED_` prefixes; never overwrite a prior
sprint. **A completion flag must require the FULL configuration, not the subset it was called with.**

---

## 7. CONTROLS — the six rules this programme paid for

0. **ENUMERATE THE OPERATOR FORKS BEFORE THE RUN — and name the alternative you did NOT take.**
   *(Added 2026-09-07; proposed by WORKSTREAM D, adopted verbatim, binding.)* For any comparison
   with a **directional hypothesis**, list the five operator axes in the module docstring —
   **functional, basis, readout, normalisation, null** — and beside each declared choice write the
   alternative you rejected.

   **The evidence for making this binding rather than advisory is the coordinator's own arm tonight**
   (L14). Its primary carried **three** unstated operator differences — argmin vs average, squared
   vs Bayes functional, built chain vs retrieved window — worth **+0.215 Å between them against a
   true effect of +0.025 Å**, i.e. **89% of the reported effect was operator, and all three pointed
   toward the conclusion the author expected.**

   **D's diagnosis is the reason the rule takes this form:** the tell is *not* the direction of any
   single choice — each fork had a defensible answer on both sides — it is that **all three went
   unstated. A fork you have registered is a fork you cannot take silently.** "Be more careful" is
   not a mechanism; enumeration is. `s21/tailprice.py` already implements exactly this pattern on
   the **normalisation** axis alone (declared rank-normal, computed raw-sum and z-score alongside,
   with the rule *"if a non-declared arm wins, the declared choice is recorded as WRONG, not
   swapped"*) — and it caught nothing only because it was applied to **one of the five axes**.
   Applied to all five it would have caught all three of tonight's.

   **Corollary, from the same incident: never correct an ARGMIN comparison with a mean-shift
   constant measured per member.** Rebuilding or reweighting changes *which* item wins, an effect
   with no reason to equal the average displacement — measured here at **7× the per-member figure**.
   Same class of error as quoting an MDE as a constant (§1).

1. **A control must be matched in the space the OPERATOR works in** — the clip, the whitening, the
   projection and the relaxation are all part of the operator.
2. **When an analytic null and a measured null disagree, the measurement is the null.**
3. **A gate can pass VACUOUSLY** — report how many times the bound/guard/threshold actually FIRED,
   and **verify a bound on the calls it actually bound.**
4. **Uniform is not a zero-information control** — use plausible-but-uninformative (constant helix,
   matched empirical marginals, label permutations). Note the role-dependence: as a *reference* the
   helix is a null; as a *gate* it is the worst arm measured.
5. **THE SHARED REFERENT FLOOR** — two quantities measured as deviations from a common reference
   correlate by construction. Measured floor here: **0.505**. **Measure the floor before interpreting
   any alignment.**

**Matched classical controls on every quantum claim**: greedy, classical analogue of the optimiser,
simulated annealing, classical gradient, entropy-matched sampler, and **best-of-N from the UNTRAINED
circuit** — never an initialisation mean. **The goal is not to make classical lose; it is to identify
exactly what the quantum method contributes.**

---

## 8. PRE-REGISTRATION — mandatory before any large experiment

    Hypothesis / Prediction / Primary endpoint / Falsifier / Null / Matched control / Budget / Promotion rule

**Pre-register the FALSIFIER, not the hypothesis.** Do not move the goalpost after seeing data; if a
pre-registration was mis-specified, say so, keep it unedited, and record the untested regime as OPEN.

---

## 9. STATISTICS AND LABELS

**TARGET is the unit.** Paired target-level bootstrap CIs, medians and W/L beside every mean, fold
behaviour. **Quote the CI construction with the CI** — `s12.instrument.paired` is **i.i.d. over
targets, not fold-clustered**, across 293 call sites; use `s18.phys_lib.paired`'s `ci_fold` where the
disposition could depend on it.

**A zero-spanning CI with insufficient power is NOT MEASURED, never "matched".** A result below the
0.084 A MDE is not a validated improvement. An argmin over a grid on the tuning instrument is a
hyperparameter chosen there.

Labels: **EXACT · ORACLE · ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN · INCONCLUSIVE · NOT MEASURED ·
NOT SUPPORTED · REFUTED · RETRACTED.** **Never present an identity as an empirical discovery; derive
the operator before interpreting its statistic.**

---

## 10. READ THE CODE BEFORE BELIEVING THE CLAIM

The source has contradicted the table in each of the last four sprints, and **this applies to
performance claims exactly as to scientific ones** — a "pathological minimiser" turned out to be CPU
starvation from concurrent lanes. **Measure timings on a quiet box or not at all.**

---

## 11. COMPUTE

8 cores (~6.43 core-equivalents), ~11 GB RAM, shared box. Target ~90-95% CPU, **never deliberately
exceed 97%.** One heavy process per workstream; cap BLAS threads; checkpoint per target; read
`s19/cache/start_<pdb>.npz` rather than recomputing the deterministic start.

> **OpenMM/AMBER CONTEXTS ARE SERIALISED — ONE LANE AT A TIME.** `core.amber.memory_guard` refuses
> above 92% and the box has been sitting at 95-96%; **two lanes' runs have already been killed by
> it.** Priority tonight: **A's mandatory matrix** > **C's continuation** > **D's encoding control
> (AMBc half)**. If you need an AMBER context and another lane holds one, **do the Legacy or
> distogram half of your work first and queue the AMBER half** — every lane has one. Announce
> acquisition and release through the coordinator.

---

## 12. WHEN TO KILL A BRANCH

Loses to a plausible zero-information null · works only on the tuning subset · CI spans the required
effect · dissolves under target-level analysis · reproduced by a simpler classical baseline · needs
native information · rests on an arbitrary encoding · improves an intermediate metric while worsening
RMSD · needs post-hoc tuning · mathematically invalid mechanism.

**A clean falsification is a successful result. Interesting is not enough.**
