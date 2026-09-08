# Sprint 22 — WORKSTREAM B PRE-REGISTRATION

Written before any Sprint-22 number exists. Dated addenda only from here on; nothing above this
line is edited after data.

Date: 2026-09-07.

---

## 0. SCOPE FOR THIS DOCUMENT

BRIEF.md hands Workstream B four numbered priorities. This registration covers **Priority 1, the
controlled synthetic perturbation probe**, named the highest-priority, mechanism-defining
experiment of the workstream. Priorities 2-4 (surrogate search, normalised continuation,
Hamiltonian-vs-repair) are addressed in `agentB_FINDINGS.md` primarily by **citing and auditing
what Sprint 21 already established** (C5, C6/C6′/C7′, C8c — all in `s21/CLAIMS.md` section C),
because BRIEF section 1 says plainly **"WHAT SPRINT 21 CLOSED — DO NOT REDO"**, and continuation
(item 3) and preconditioning (the mechanistic half of item 4's AMBER-as-Hamiltonian role) are on
that list. Re-litigating a closed result without new information would burn compute BRIEF section
8 asks me to conserve. If time and the AMBER queue permit after Priority 1 is complete, a small
labelled addendum will be appended below the line — never a rewrite of this section.

---

## 1. HYPOTHESIS

**H1.** Legacy and AMBER respond to *different* physical variables under matched-magnitude
torsion-space perturbation, and the difference is not incidental: it follows from an already-measured
structural fact about the two Hessians (`s21` C5, reproduced and cited, not re-run):

    participation ratio   Legacy 0.497   AMBER 0.067   (curvature spread vs concentrated)
    condition number       Legacy 6.26    AMBER 7,756   (AMBER's spectrum is far more anisotropic)

A Hessian whose curvature is concentrated in ~7% of directions (AMBER) should show a **disproportionate**
energy response when a fixed torsion-space perturbation budget is *concentrated* onto few coordinates
(a single residue, or a targeted residue pair) relative to when the identical budget is *spread*
across the whole chain — more than a Hessian whose curvature is spread across ~50% of directions
(Legacy) would show under the same contrast. Separately, established mechanism (`s21` C7″, L33: "Legacy
prefers 0.45 Å more compact structures... compaction closes the contacts that ARE the [AMBER steric]
singularity") predicts Legacy's response is **larger, relative to AMBER's, for a compactness-directed
move** than for a concentration-directed one, and the reverse for AMBER.

## 2. PRIMARY ENDPOINT

**NOT RMSD** (Hard Rule 1 exception, stated explicitly per that rule's own requirement: "or is
labelled as not doing so"). The endpoint here is an **energy response**, because the question is
mechanistic (what does each Hamiltonian see), not evaluative (does it help RMSD) — RMSD-facing
consequences of this mechanism are for a downstream lane, not this probe.

For each of five perturbation **axes** (A: diffuse/random torsional deviation, B: compactness-directed,
C: contact-density-directed local window, D: steric-pair-directed, E: single-residue concentrated
kick), and for each of the two potentials (Legacy, AMBER — **AMBER here means the BARE SINGLE POINT,
`ConstrainedBox.energy_point`, exactly the object `s21/c_norm.py` and `s21/c_cont.py` declare and
use for the identical reason: a landscape needs a function of theta, and the deployed `H_AMBER =
E ∘ Relax_50` is not one, because the relaxation leaves the torsion manifold** — stated once here and
then in every table, per BRIEF's instruction to state which AMBER is in play every time):

    NORMALISED SENSITIVITY(axis, potential) =
        robust standardised slope of  z(dE_potential)  on  z(target-proxy_axis)
        pooled over all perturbation trials for that axis,
        z(dE) := (dE - median(pool dE for that potential, that target)) / (1.4826 * MAD(...))
        the "pool" is that target's own shipped top-75 real rebuilds, reusing the median/MAD
        already computed and persisted by `s21/results/c_norm.json` (bare single-point AMBER,
        genuine Legacy total) -- not recomputed, so the normalisation is identical to the one
        already audited by C5/C6.

This is dimensionless in both the energy axis and (via Spearman rho, reported alongside) the proxy
axis, so a Legacy number and an AMBER number are directly comparable — the thing BRIEF's compactness
clause under Hard Rule 1 (scale-free wherever possible) requires.

Per-axis target proxy:

    A  ||delta_theta||_active,rms      ("torsional deviation", the design/dose variable itself)
    B  delta Rg  (radius of gyration, core.geometry.radius_of_gyration on the rebuilt CA trace)
    C  delta contact_frac  (fraction of CB-CB pairs, sep>=3, within the 8.5 A contact-term cutoff
       -- purely geometric, NOT MJ-weighted, so it is an energy-model-independent proxy)
    D  delta d_ij  (CB-CB distance of the target-engineered pair, and delta steric (raw unweighted
       `core.energy.steric_term`) reported alongside)
    E  the single-residue kick magnitude (== axis A's dose m, but concentrated on one residue's
       (phi,psi) instead of spread over the whole active vector)

## 3. EXPECTED MECHANISM

Stated in section 1. The falsifiable, single-number form:

    CONCENTRATION RATIO  R(potential) =
        median( |dE| / ||delta_theta||_rms )  over axes {D, E}     [concentrated]
      / median( |dE| / ||delta_theta||_rms )  over axis  {A}        [diffuse]

Prediction: **R(AMBER) > R(Legacy)**, bootstrap CI on `R(AMBER)/R(Legacy)` excludes 1 on the side
predicted.

## 4. FALSIFIER

Registered before any trial is run, never moved afterward:

**F-B1** (concentration). If the bootstrap CI (percentile, over targets as the resampling unit,
2000 draws) on `R(AMBER)/R(Legacy)` **includes 1, or excludes 1 on the wrong side** (i.e. Legacy
turns out more concentration-sensitive), H1's concentration clause is **REFUTED**.

**F-B2** (compactness/contact specialisation). If the normalised sensitivity to axis B (compactness)
is **not** larger in magnitude for Legacy than for AMBER (bootstrap CI on the Legacy-minus-AMBER
difference of |normalised sensitivity| does not exclude zero on the Legacy-larger side), H1's
compactness clause is **REFUTED**.

Both clauses are reported regardless of outcome; a partial refutation is written up as a partial
refutation, not silently dropped.

## 5. NULL

**Within-target, within-axis permutation of the (proxy, dE) pairing** across the perturbation trials
generated for that (target, axis) cell: shuffle which trial's realised proxy value is paired with
which trial's dE, holding target identity and axis fixed, 2000 shuffles, giving a null distribution
for the standardised slope. This is deliberately **not** a cross-target shuffle: BRIEF section 2
records that a cross-target permutation null on this exact instrument family is mis-specified,
because targets are correlated through a shared difficulty factor and shuffling across targets lets
the null draw from that shared factor (`the routing-ceiling null I already got wrong`). Shuffling
within one target's own trials removes only the perturbation-response correspondence, not the
between-target structure, so it targets the right null hypothesis: "this axis's proxy predicts this
axis's dE beyond what the trial-to-trial noise of that target alone would produce."

## 6. MATCHED CONTROL

Every perturbation across all five axes is matched on **realised RMS torsion-space displacement**
over the full active coordinate set (`||delta_theta||_active,rms`), at three shared magnitudes
(0.02, 0.05, 0.10 rad), so magnitude is never a confound between axes — the design principle BRIEF
section 7 states as "a control must be matched in the space the operator works in." The random axes
(A, E) additionally control for direction idiosyncrasy with independent repeated draws; the directed
axes (B, C, D) are run at both signs (toward/away from the target reference) so a one-sided reading
of a symmetric potential cannot masquerade as a directional finding.

## 7. BUDGET

Reuses `s21/c_norm.SUBSET` — the same 30-target subset already carrying the C5 curvature numbers
this hypothesis is built from, so the mechanism claim and the curvature claim are about the *same*
targets rather than requiring a second, unstated generalisation. Per target: 3 starts (pool medoid +
2 deterministic random pool draws, `SD.stable_rng` seeded identically to `s20/c_land.starts_for`'s
convention) x ~43 perturbation trials (axis A: 3 magnitudes x 4 repeats = 12; axis E: 12; axis B:
3 magnitudes x 2 signs = 6; axis C: 6; axis D: 6; +1 baseline) = ~130 AMBER single points and ~130
Legacy evaluations per target, ~3,900 AMBER evaluations total. At the measured ~7-10 ms per bare
AMBER single point (module docstring, `core/amber.py`) plus per-target box setup, budget estimate is
under 10 minutes of wall clock, run as a single serialised block.

**AMBER/OpenMM SERIALISATION (BRIEF section 8).** This run is the only process in this session
touching an OpenMM context; announced here before the run starts, and each target's context is
closed (`Pot.close()`) immediately after that target's trials complete, before the next target opens
one — so at no point does this lane hold more than one live context, and the box never sees more
than this one heavy job from this lane during the run.

## 8. PROMOTION CRITERION

This is a mechanism probe, not a selector. Nothing here is promoted into the pipeline. Promotion
criterion is definitional: **F-B1 and F-B2 either survive intact (report as ESTABLISHED with
mechanism), survive with one clause refuted (report the surviving half only), or both refute (report
"Legacy and AMBER are not measurably specialised along these five axes, and the earlier curvature
asymmetry does not translate into differential axis-sensitivity" — a genuine negative, written up
with the same weight as a positive).**

## 9. RULE 0 — SIX OPERATOR FORKS, ALTERNATIVE NOT TAKEN NAMED

**FUNCTIONAL.** TAKEN: bare AMBER single point (`ConstrainedBox.energy_point` via `s20.c_land.Pot`).
NOT TAKEN: the deployed `E ∘ Relax_50` (not a function of theta — relaxation leaves the torsion
manifold, exactly `s21/c_cont.py`'s stated reason) and `E ∘ Relax_1` (+inf on 42% of the register per
Sprint 20 L6). A landscape mechanism probe needs a function of theta; this fork is forced, matches
`s21`'s convention, and is why no number in this file is a statement about the deployed operator —
see Priority 4 in `agentB_FINDINGS.md` for that distinction made explicit.

**BASIS.** TAKEN: real starting structures from the shipped top-75 real-rebuild pool (pool medoid +
2 deterministic pool draws), i.e. genuine candidate geometries, not synthetic decoys. NOT TAKEN:
starting from an idealised straight or helical chain, or from lattice states — a mechanism claim
about "what the deployed pipeline's energies respond to" should be measured where the pipeline
actually operates.

**READOUT.** TAKEN: signed raw dE (Legacy total, AMBER bare single point), standardised per target
by that target's own pool median/MAD. NOT TAKEN: an RMSD-based readout (explicitly declared as the
Hard-Rule-1 exception in section 2), and NOT an argmin-over-perturbations readout (there is no
selection step here; every trial is reported, not just the best one — Rule 0 clause 3's best-of-K
null hazard does not apply because nothing here takes a max or min over trials).

**NORMALISATION.** TAKEN: per-target robust z (median, 1.4826×MAD) from the already-persisted,
already-audited `s21/results/c_norm.json`. NOT TAKEN: a global cross-target normalisation (would
mix target difficulty ranges into the sensitivity estimate) and NOT the `Nt` asinh transform `s21`
declared for the continuation work — reasoned here rather than assumed: `Nt` exists to keep a
*mixture* `H(lambda)` well-behaved under summation; this probe never sums or mixes the two energies,
it only correlates each one separately against a proxy, so the additional nonlinearity would only
complicate the slope's interpretation for no protective benefit. If a reviewer wants the `Nt` version
computed alongside, it is cheap to add as an audit column; declared as NOT the primary here, before
any number exists.

**NULL.** TAKEN: within-target, within-axis shuffle of the (proxy, dE) pairing (section 5). NOT
TAKEN: a cross-target shuffle (mis-specified for this instrument family, per BRIEF section 2's
worked example) and NOT an isotropic-random-direction null in place of the constructed axes (BRIEF
section 7: a zero-information control must be plausible, and axis A already *is* the isotropic
random control for axes B-E — a second one would be redundant, not more rigorous).

**THE LABEL.** The one binary-flavoured quantity in this design is the falsifier's direction test
(F-B1, F-B2). Both are stated as **continuous ratios with a bootstrap CI test against a fixed
reference value (1 or 0)**, not a threshold whose cut point depends on any covariate of the data
(peptide length, difficulty, or otherwise) — the exact hazard Rule 0's sixth fork names. No length-
or difficulty-dependent binarisation appears anywhere in this design.

---

*(Addenda, if any, appended below this line with today's date, never editing the above.)*

---

## ADDENDUM 1 — 2026-09-07, same day, after F-B1/F-B2 were read on the main run

F-B1 (section 4) came back **REFUTED, and on the wrong side**: `R(AMBER)/R(Legacy) = 0.19, CI
[0.041, 0.514]` — the *opposite* of the predicted direction. Before writing that up as a clean
refutation, I re-examined the operationalisation and found it was not testing what section 1's
mechanism argument actually claims, for a reason specific to how `core.geometry.build_backbone`
works: it is a **sequential NeRF chain**, so perturbing one residue's (phi, psi) moves every atom
from that residue to the C-terminus, not "one residue's worth of geometry." **Torsion-INDEX
concentration (few nonzero coordinates) is not CARTESIAN concentration, and C5's participation
ratio describes concentration in the Hessian's own EIGENBASIS, which axes D and E never targeted.**
This is caught and reported here, before any further data is read on the corrected version, rather
than being folded into the original F-B1 with the mismatch left implicit.

**F-B1′ (replaces F-B1's operationalisation; does not touch F-B1's own recorded result, which
stands as measured, under the operator it actually used).** Perturb along each potential's OWN
top-|eigenvalue| Hessian eigenvector (computed once per target at the medoid start, identical
`hess_fd` machinery C5 used) instead of along a coordinate-index-concentrated direction. Falsifier:
bootstrap CI on `OWN-RATIO(AMBER) / OWN-RATIO(Legacy)` (defined in `s22/probe_eigen.py`'s docstring)
excludes 1 on the wrong side or includes 1 (REFUTED / NOT MEASURED); supported only if the CI
excludes 1 on the AMBER-larger side. Same matched RMS magnitudes (0.02, 0.05, 0.10 rad), both signs,
same 30-target subset, same medoid starts already used for the main run's start 0. Null and control
conventions are unchanged from sections 5-6. This is a **correction of the operator, registered
before its own data is read** — exactly the discipline BRIEF section 6/Rule 0 requires of a reviewer
who finds a fork after the fact, applied here by the same lane to its own design.

**Result, recorded here rather than only in FINDINGS**: F-B1′ **SUPPORTED**. Reproduced median
participation ratio on this 30-target subset (Legacy 0.4965, AMBER 0.0668) matches C5's 0.497/0.067
almost exactly, certifying the instrument before the new number is read. `OWN-RATIO(Legacy)=2.559`,
`OWN-RATIO(AMBER)=74.43`, ratio `29.09`, bootstrap 95% CI `[4.49, 542.68]` (2000 draws, targets as
the resampling unit) — excludes 1 on the predicted side by a wide margin. **A bootstrap coding
defect was found and fixed before this number was quoted**: the first version re-seeded
`np.random.default_rng(0)` inside the per-draw resampling closure, so all 2000 "draws" resampled
identically and returned a degenerate CI (`[12.72, 12.72]`) — caught because a CI with equal bounds
is itself a red flag, fixed by hoisting the RNG outside the loop, and reported here rather than
silently discarded.

## ADDENDUM 2 — 2026-09-07, same day, Priority 2 (AMBER surrogate)

**Hypothesis.** `s21.c_cont.MixPot.amber_bonded` — the bonded (bond+angle+torsion) force-group
subset of the SAME genuine ff14SB System, nonbonded/solvation off, already built and named (never
"amber") in `s21/c_cont.py` — tracks the genuine bare-single-point `H_AMBER`'s RESPONSE to a torsion
perturbation well enough to be a defensible cheap RANKING proxy, even though it structurally cannot
match its magnitude (nonbonded/solvation is ~98% of one AMBER evaluation's cost per `core/amber.py`'s
own docstring).

**Falsifier.** Pooled Spearman rho between `dE_amber` (genuine) and `dE_bonded` (candidate
surrogate) over the same axis-A (diffuse) and axis-B (compactness) trials, medoid start, full
30-target subset: rho > 0.7 supports a defensible ranking surrogate; rho <= 0.7 refutes it. Threshold
fixed here, before the run, as a round, pre-committed number rather than picked after seeing the
result.

**Basis/readout/normalisation/null forks**: BASIS is the same real top-75 pool medoid starts as
Priority 1 (not idealised geometry). READOUT is the raw signed dE of each object, unstandardised —
a ranking (Spearman) question does not need a robust z. NORMALISATION: none — Spearman is already
rank-based and scale-free by construction, so no transform is needed or applied. NULL: none required
— this is a direct correlation test with a fixed pre-registered threshold, not a comparison against
a zero-information control. THE LABEL: the 0.7 threshold is a fixed constant, not a function of any
covariate of the data (length, difficulty, or otherwise).

