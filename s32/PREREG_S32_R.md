# PREREG — S32 LANE R (RECONSTRUCTION)

Committed **before the first number exists**. Any arm added after this file is committed is
labelled EXPLORATORY in `s32/MULTIPLICITY.md`.

Basis for every RMSD in this lane unless stated otherwise: **built-chain Cα RMSD,
`tuning126`, n = 126**, chain = `s12.instrument.project(C, seq, fold)["ca"]` (ramah penalty,
λ = 0.3, multi-start, maxiter = 300, grad = "exact") scored against `u["nat_ca"]` by
`s12.instrument.ca_rmsd`. **Production = 3.2105 Å. Cloud = 3.0483 Å is a DIFFERENT OBJECT and
is never differenced against a chain number except as the explicitly named `price`.**

Input cloud: `s29/results/s29_O_structs/<pdb>.npz["prod"]`, the production 75-member
coordinate average. Rule 3 check: this array is asserted bit-identical to
`bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json["avg_ca"]` for all 126 targets, and the
check is printed and persisted. Both arms of every chain contrast are projected **in the same
job**.

---

## R1 — Is the +0.1622 Å projection price real, and is any of it recoverable?

**HYPOTHESIS.** The price `RMSD(chain) − RMSD(cloud) = +0.1622` is a real property of the
production average, and it decomposes into (a) a *magnitude* term forced by the cloud being
off the valid-chain manifold by `d = RMSD(chain, cloud)`, and (b) a *direction* term: whether
the projection displacement points toward or away from the native.

**MECHANISM.** Kabsch RMSD is a metric on shape space. If the projection displacement of
magnitude `d` were statistically orthogonal to the cloud's native error `e`, the chain error
would be `sqrt(e² + d²)`. This is the ISOTROPIC NULL, matched to the operator's own
displacement magnitude (contract rule 6).

**PREDICTION.**
- P1.1 the reproduced production price is +0.1622 ± the chain floor (mean floor 0.0134).
- P1.2 the measured price is **not below** the isotropic null; i.e. the projection does not
  systematically move toward the native.
- P1.3 the cloud is contracted relative to an ideal chain (virtual Cα–Cα < the builder's
  ideal), and a native-free isotropic rescale to the ideal bond reduces the price.

**FALSIFIER.** P1.2 is falsified if `mean(price) < mean(price_null) − MDE`. P1.3 is falsified
if the rescale arm's **chain** mean is not lower than production's chain mean by ≥ 1× its own
MDE, in which case scale repair is declared not the mechanism.

**CONTROL.** the isotropic null `sqrt(e²+d²) − e`, computed per target from that target's own
`d`. ORACLE rescale (best `s` on a grid, per target) prices the ceiling of the scale axis and
is labelled ORACLE / NOT DEPLOYABLE.

**TEST.** paired `s24.stats_lib.compare(arm, prod_chain, folds=pinned_folds())`.

**STATISTICAL RULE.** contract rule 1. < 0.7× MDE = NOT A RESULT; 0.7–1.0× = NOT MEASURED;
≥ 1.0× with fold CI excluding zero and ≥ 4/5 folds agreeing = RESULT.

**DEPLOYMENT CONDITION.** the rescale factor is computed from the cloud's own geometry only
(no native, no label), so a RESULT here is deployable as stated.

---

## R2 — How much RMSD is the multi-start branch choice leaving on the table? (ORACLE)

**HYPOTHESIS.** The λ-path argmin over branches is decided by an objective that is degenerate
to ~1e-7 among branches that are ~1e-1 apart in RMSD (S31-L6). Therefore the ORACLE-best chain
over a larger restart set sits materially below 3.2105, and the gap prices a **discrete,
target-dependent decision problem whose candidates the existing code already generates**.

**MECHANISM.** `core.project.fit_multi` starts L-BFGS-B from 4 fixed constant-torsion points
and keeps the strict objective argmin. Warm starts cannot cross between the near-degenerate
ideal-geometry torsion solutions (S9-2), so the emitted chain is determined by which basin each
start falls into, not by the objective.

**START SETS** (per target, all native-free to construct):
- `GEN4` — the 4 production generic starts. **This is production.**
- `MEM75` — the φ/ψ of the production top-75 filtered pool members (`sub` in the production
  cache record). Real torsions, target-specific, already computed by the pipeline.
- `RAND75` — 75 windows drawn uniformly at random (seed pinned) from the same target's window
  universe, **not** score-selected. The matched control: same object type (real torsions), same
  count, no selection. Its draw distribution is reported per contract rule 10 (≥ 3 independent
  draws, mean and draw-to-draw sd).

Each start `s` is run as production runs a branch: `fit_prior(C, s, lam=0)` then
`fit_prior(C, φ₀, ψ₀, pen, lam=0.3)`. Recorded per branch: λ=0 and λ=0.3 objective values,
`d_to_C`, the emitted Cα trace, and its RMSD to native (ORACLE **label only**).

**PREDICTION.**
- P2.1 branches within a target span ≥ 0.3 Å in RMSD on a majority of targets while their
  λ=0.3 objectives span < 1e-3.
- P2.2 `ORACLE-best over GEN4 ∪ MEM75` is lower than 3.2105 by ≥ 0.3 Å.
- P2.3 the ORACLE-best gain is **not** explained by best-of-K alone: `stats_lib.best_of_k_null`
  / `best_of_k_accounted` against the branch distribution, and split-half transfer where a
  criterion is involved.

**FALSIFIER.** P2.2 is falsified if ORACLE-best-over-branches is within 1× MDE of production.
P2.1 is falsified if the within-target RMSD spread is < 0.05 Å (then the branches are the same
structure and there is no decision to make).

**CONTROL.** (i) `RAND75` at matched N; (ii) the best-of-K order statistic priced explicitly;
(iii) the *mean* over branches, which is what a criterion with zero skill delivers.

**STATISTICAL RULE.** ORACLE arms are reported as ceilings with SE, never as effects that
clear MDE for deployment. The *deployable* arms (objective-argmin, native-free-criterion
argmin) are tested by `compare(..., folds=...)` against production.

**DEPLOYMENT CONDITION.** ORACLE-best is **NOT DEPLOYABLE** and is labelled so on every
occurrence. Only R3 can make any of it deployable.

---

## R2a — the branch census (registered, added after the coordinator's 2026-09-21 steer)

Before any ORACLE ceiling, count the object. Per target, over the full start set: the number of
**distinct converged branches** (clustered at Cα-RMSD < 1e-3, the same threshold
`core.project.degeneracy` uses for `n_same_point`), the λ=0 and λ=0.3 objective gap between the
best and runner-up **distinct** branch, and the **RMSD spread** across distinct branches.
`verify/project_degeneracy.json` and `verify/project_stability_partial77.json` are read first
as the existing artefacts of exactly this question, and my census must agree with
`project_degeneracy`'s `runner_up_dist` distribution on the 4 production starts or I say why.

**FALSIFIER for the whole lane.** If the median target has 1 distinct branch, R2/R3 are empty
and the lane reports that.

---

## R3 — Is there a NATIVE-FREE criterion that picks a good branch?

**HYPOTHESIS (revised, primary).** The branch degeneracy is a **torsion-branch** degeneracy with
a known physical asymmetry: `core.project`'s own docstring states "a CA trace admits two
ideal-geometry torsion solutions at near-equal objective distance, **one Ramachandran-plausible
and one not**". If that is true in practice, a Ramachandran log-likelihood — native-free,
target-specific through `res_classes(seq)` and the fold, and **chiral** — has in-band skill at
picking the branch.

**MECHANISM, and why it is not another scalar.** Every native-free ranker this project has
tested in band is a function of the distance map and is therefore **achiral**: a structure and
its mirror image score identically. A torsion-branch choice is exactly the thing an achiral
observable cannot see. The Ramachandran density over (φ, ψ) for L-amino acids is strongly
chiral (the α_L region is rare), so it is the first native-free in-band signal in this project
with a mechanistic reason to work rather than a fitted one.

**WHY THE FINAL RUNG DOES NOT ALREADY DO THIS.** The branch is fixed at the **λ = 0** rung of
`lam_path`, where the penalty has weight zero, so the branch is chosen by a coordinate distance
whose best-vs-runner-up gap is ~1e-7 on the majority of targets. At λ = 0.3 the penalty that
enters is `RamaHingePenalty`, which has **sparse support**: a residue already as plausible as
95% of real residues of its class contributes exactly zero. Among branches that are all above
the hinge threshold it is identically blind. The criterion proposed here is the **un-hinged,
dense** `RamaPenalty` log-density, which is a different function of the same table.

**HOW THIS DIFFERS FROM S31's BRANCH-CARRY NULL** (−0.0055 Å at 0.44× MDE): S31 *carried* a
branch down the λ ladder — it changed which trajectory was continued. This lane *scores the
final converged branches* with a denser, un-hinged, chiral likelihood and re-selects among
them. Different operator, different object, different selection rule. S31's null therefore
does not close this, and this sentence is the named mechanism required by contract rule 19.

**CRITERIA SCORED PER BRANCH** (all native-free, fixed here before any number exists):
1. `rama_nlp` — `RamaPenalty(seq, fold, kind="rama")(φ, ψ)`, the dense 4-class negative mean
   log-density. **PRIMARY.**
2. `rama20_nlp` — the same with `kind="rama20"` (per-amino-acid table). Secondary.
3. `ramah` — `RamaHingePenalty`, the production penalty. **Expected near-blind** — the
   sparse-support prediction, and a falsifier for the mechanism story if it is *not* blind.
4. `obj_lam03`, `obj_lam0`, `d_to_C` — the projection's own objective. **The NULL**: this is
   what production already uses and it should have in-band ρ ≈ 0 by construction.
5. `posphi_frac` — fraction of non-glycine residues with φ > 0. The crude chirality readout and
   the cheapest possible version of the mechanism.
6. `disto_risk` — `s12.instrument.shipped_score(distogram, pair_dists(chain))`; `disto_mae` —
   mean |d_ij(chain) − E_disto[d_ij]|. Achiral by construction; included precisely to test the
   achirality argument.
7. `legacy` — `lf.BatchLegacy` `lg_all` with `FITTED_WEIGHTS`, from the branch's own φ/ψ.
8. `typicality` — mean Cα-RMSD of the branch to the production top-75 members.
9. AMBER ff14SB/GBn2 single-point — **one process at a time**, on distinct branches only, with
   its own coverage reported. Lowest priority; skipped with a stated reason if the branch
   census makes it uneconomic.

**PREDICTION.**
- P3.1 `obj_lam0` / `obj_lam03` in-band Spearman ρ with true RMSD is indistinguishable from 0.
- P3.2 `rama_nlp` in-band ρ is **positive** (lower neg-log-density ↔ lower RMSD) with a
  fold-clustered CI excluding zero.
- P3.3 the achiral channels (`disto_risk`, `disto_mae`) have strictly smaller |ρ| than
  `rama_nlp`.
- P3.4 **directional** (contract rule 18): `rama_nlp`-argmin over branches beats the branch
  mean AND beats production by ≥ 1× its own MDE, n = 126, built chain.

**FALSIFIER.** P3.2 is falsified if the fold-clustered CI on in-band ρ covers zero — which
would say the degeneracy is **not Ramachandran-separable in practice**, contradicting
`core/project.py`'s own docstring, and that is itself a reportable result. P3.4 is falsified if
the argmin arm is within its own MDE of production; then the finding is "correlated, not
directional" (the S31 rule-18 shape).

**CONTROL.** (i) a random branch (its own draw distribution, ≥ 3 draws); (ii) the branch mean;
(iii) ties: `stats_lib.argmin_tied` everywhere, never `np.argmin` on a possibly-tied
criterion — averaging the outcome over the tied argmin set (standing lesson, S12).

**STATISTICAL RULE.** in-band ρ is aggregated over targets with a fold-clustered CI.
Directional arms use `compare(..., folds=...)`. A criterion selected by looking at the ORACLE
ranking is ORACLE; criteria are fixed in this file before any number exists.

**DEPLOYMENT CONDITION.** a criterion is deployable only if it is computed from the sequence,
the leave-fold-out model and the branch geometry — no native, no label, no tuning on RMSD.

---

## R4 — Can a readout land in the cheap projection regime without losing the averaging benefit?

**S31 already refuted the per-target `m` axis (ORACLE global m −0.0044 at 0.19× MDE;
leave-fold-out m +0.0075, 63W/63L). THIS LANE DOES NOT REDO THE `m` AXIS.**

**HYPOTHESIS.** The +0.1622 price is paid because the average is snapped back onto the
valid-chain manifold from four *generic* constant-torsion starts, which is a badly conditioned
route. Snapping back from a start that is already a valid chain near the average — the medoid's
own torsions, or any top-75 member's torsions — reaches the manifold without the branch lottery.

**MECHANISM.** the coordinate average is built by superposing the top-75 on the MEDOID and
taking the mean, so the medoid's torsions are, by construction, the nearest valid-chain
parameterisation of the object being projected. `core.pipeline.project()` calls
`lam_path(..., extra=None)`: the medoid start is **never offered to the production
projection**, although `lam_path` accepts exactly that argument and `project_single` already
computes it as a reporting arm.

**ARMS.**
- `PROD` — production, GEN4.
- `MEDOID_EXTRA` — production's own rule (strict objective argmin) over `GEN4 ∪ {medoid
  torsions}`. Fully deployable, one extra start.
- `GEN4∪MEM75_OBJ` — production's rule over all 79 starts. Fully deployable.
- `MEDOID_ONLY` — the medoid-start branch alone, no objective selection. Deployable.
- `SCALE+PROD` — R1's native-free rescale then production's rule.
- Each arm reports, from its first row (contract rule 15): virtual-bond mean and sd of the
  emitted chain, and its displacement (Cα-RMSD) from the production chain and from the cloud.

**PREDICTION.** P4.1 `MEDOID_ONLY`'s price (chain − cloud) is smaller than production's
+0.1622. P4.2 at least one of the deployable arms beats production's **chain** mean by ≥ 1×
MDE.

**FALSIFIER.** P4.2 is falsified if every deployable arm is within its own MDE of 3.2105. In
that case the conclusion is that the objective argmin over a larger, better start set is a
random redraw (P2.1's degeneracy), which is a *mechanism*, not "didn't work".

**CONTROL.** `RAND75_OBJ` — the same rule over the random-window starts, at matched N, with its
own draw distribution. This separates "better starts" from "more starts".

**DEPLOYMENT CONDITION.** no native anywhere in the start set, the selection rule or the
scale factor.

---

## Shared rules for this lane

- Every number is persisted to `s32/results/*.json` (stable keys) before it is written in prose.
- Every emitted comparison is appended to `s32/MULTIPLICITY.md` as it is emitted, with
  registered-vs-exploratory marked.
- Per-target chain floor 0.0134 mean / 0.0329 p90 / 0.2285 max. No per-target claim smaller
  than the p90 floor is made; no mean claim smaller than 0.0134 is made.
- `stats_lib.compare` is LOWER IS BETTER and is always passed `folds`.
- ORACLE arms carry "ORACLE / NOT DEPLOYABLE" in the artefact key AND in every table row.
