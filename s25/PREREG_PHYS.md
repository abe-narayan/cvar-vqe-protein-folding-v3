# SPRINT 25 — LEGACY / AMBER PHYSICS LANE. PRE-REGISTRATION.

Written and committed to disk **before any configuration was run**. Rule 0 six forks below,
each naming the alternative NOT taken. Falsifiers registered in §4.

---

## 1. DELIVERABLE 1 — THE 7-CONFIGURATION COMPARISON SUITE

### 1.1 The seven configurations

Three channels, all NATIVE-FREE, all "lower is better", all defined on the **same K = 500
candidates** of the **same target**:

| channel | symbol | provenance |
|---|---|---|
| Legacy | `x_LEG` | `s16.energy_lib.legacy_total_from(legacy_components_of_windows(seq, PHI, PSI))`, 11 terms at `core.energy.DEFAULT_WEIGHTS`, **never fitted** |
| AMBER | `x_AMB` | genuine ff14SB/GBn2 **single point** (`s20.qb2_lib.AmberSP`, `k_restraint=0`, `steps=-1`), read from `s24/cache_amber/` after a bit-exactness re-check |
| Distogram | `x_DIS` | the shipped Bayes-risk score, `s24.d_harness.score_shipped` = `I.shipped_score(distogram, pair_dists(W))` |

    C1 = {LEG}   C2 = {AMB}   C3 = {DIS}
    C4 = {LEG,DIS}   C5 = {AMB,DIS}   C6 = {LEG,AMB}   C7 = {LEG,AMB,DIS}

### 1.2 NORMALISATION, STATED MATHEMATICALLY (this is what makes "equal weight" mean something)

For a score vector `x ∈ R^K` on one target's pool, with `r = rankdata(x)` (average ranks):

    zrank(x)  =  ( r − mean(r) ) / sd(r)                                          (1)

With no ties `r` is a permutation of `1..K`, so `mean(r) = (K+1)/2` and
`sd(r) = sqrt((K²−1)/12)` **exactly** — constants that do not depend on `x`, on the target or
on the functional. Hence `zrank` maps **every** channel onto the *identical* empirical
distribution: the standardised discrete uniform, mean 0, sd 1, support `±sqrt(3(K−1)/(K+1))`
≈ ±1.7315 at K=500. With ties, (1) still returns mean 0 and sd 1 by construction.

**The combined energy of configuration S is**

    E_S  =  zrank( Σ_{c ∈ S} zrank(x_c) )                                         (2)

Two properties, both required and both proved rather than assumed:

* **Equal weight is exact.** Inside the sum every channel enters with unit coefficient onto a
  marginal with identical mean, variance, support and quantile spacing. "Equally weighted" is
  therefore a statement about matched marginals, not a hope about matched scales. *This is the
  step the brief demands and it is why a bare `E_LEG + E_AMB` would be meaningless: the raw
  Legacy total spans ~1.6 decades and raw AMBER spans 16.39.*
* **Every configuration hands the VQE an identically-scaled Hamiltonian.** Without the OUTER
  `zrank` in (2) the sum of `|S|` channels with mean pairwise correlation `ρ̄` has
  `sd = sqrt(|S| + |S|(|S|−1)ρ̄)`, which at ρ̄ ≈ 0 is 1.00 / 1.41 / 1.73 for |S| = 1 / 2 / 3.
  The CVaR objective trades energy against `T·H`, so an un-restandardised combination would
  silently change the temperature between configurations and confound "which functional" with
  "how hard the entropy term pulls". The outer `zrank` removes that exactly.
* `zrank` is **idempotent** (`zrank(zrank(x)) = zrank(x)`), so for |S| = 1 equation (2)
  reduces to `zrank(x_c)` and C1/C2/C3 are unmodified single channels.
* `zrank` is **strictly monotone**, so it changes no ordering, no argmin, no level set, and no
  CVaR tail membership. It is the currency `core.pipeline._zrank` already deploys.

**DECLARED SECONDARY ARM — raw moment.** `E_S^raw = Σ_{c∈S} (x_c − mean)/sd`, run to full
completion on all seven so the normalisation choice is visible rather than asserted. Registered
prediction in §4 F4.

### 1.3 What is held identical across all seven (the pairing discipline)

| | fixed value |
|---|---|
| targets | all **126** dev targets, `sorted(t["pdb"] for t in s12.instrument.targets())`, identical order |
| folds | the **5 pinned folds** read from the instrument (`ST.pinned_folds`), never recomputed |
| pool | the shipped **K = 500** BLOSUM retrieval pool, `d_harness.Candidates.from_universe(pdb, k=500)`, bit-identical to `s24/cache_amber`'s `universe_idx` (asserted per target) |
| basis | **point cloud** — `I.ca_rmsd(I.coordinate_average(W[S]), nat_ca)`. Never compared to built-chain. |
| metric | full-chain Cα-RMSD, Å |
| selector | **genuine CVaR-VQE**, `d_harness.arm_vqe` → `core.quantum.run_cvar_vqe`, exact `StatevectorCircuit`, 9 qubits (512 states, 500 real + 12 padding at `max + 10 sd`), 3 layers, 80 Adam iters on the **exact parameter-shift gradient**, seed 0, **α = 0.18, T = 0.5** |
| readout | **uniform coordinate average over the realised CVaR tail**, in the retained set's own medoid frame |
| pairing | every contrast is **paired per target**; `d_a − d_b` per target, never a difference of means |

**α = 0.18 is pinned on a NATIVE-FREE criterion and no other.** A 12-target probe on the
distogram channel alone (no RMSD read, no outcome inspected) gave median realised tail size
**75** at α = 0.18 — exactly the production rung `M_PROD = 75` — against 64 at α = 0.15 and 84
at α = 0.20. Pinning α there makes the VQE readout and the classical top-75 control consume the
same number of candidates, which is the only way the brief's "same readout" and "CVaR-VQE for
all seven" can both hold. The probe also recorded 8.82 bits of state entropy out of a maximum
of 9 and ESS ≈ 52, with realised m = 64 against a *nominal* α·2ⁿ = 77 at α = 0.15 — i.e. the
trained state puts ~18% more mass on the low-energy prefix than uniform, so **the tail
materially participates** (Pillar 1) rather than collapsing to an argmin.

### 1.4 The m-confound, and how it is removed rather than hidden

By the s24 set-equality theorem the realised tail is the classical top-m of the configuration's
own energy order, and `m` is data-dependent. Different configurations may therefore realise
different `m`, which would confound *which energy* with *how many candidates were averaged*
(`operator-consumes-set-mean`). Every cell therefore reports the exact decomposition

    RMSD_VQE(c)  =  RMSD_top75(c)  +  [ RMSD_top-m_c(c) − RMSD_top75(c) ]  +  ε_c
                    └ energy only ┘   └── the m-ladder movement ──────┘      └ ≈0 by theorem ┘

so any configuration-to-configuration difference can be split into an energy part and an m part.
`ε_c` is asserted ≈ 0; a non-zero value is a **bug**, not a quantum effect.

### 1.5 Statistics, per configuration and per pairwise contrast

`s24/stats_lib.py` verbatim (`ST.compare`, `ST.fmt`) — mean, median, SD, quartiles, best/worst
target, paired Δ, SE, **MDE = 2.8016 × SE computed per comparison**, effect/MDE, W/L/T, iid CI
**beside** fold-clustered CI, per-fold deltas and folds-same-sign, worst-target degradation, p90
degradation, retrodesign power/Type-M/Type-S, and the concentration test against a
uniform-effect null. All 21 pairwise contrasts plus all 7 against the incumbent. Fold-wise and
length-stratified breakdowns.

**"Best of the seven" is an ORDER STATISTIC** and is scored against `ST.best_of_k_null` with
`share_accounted` and `k_eff`, never quoted as a result.

---

## 2. RULE 0 — THE SIX OPERATOR FORKS, WITH THE ALTERNATIVE NOT TAKEN

**1. FUNCTIONAL.** *Taken:* equal-weight rank-standardised sums of three unmodified channels —
genuine 11-term Legacy at `DEFAULT_WEIGHTS`, genuine ff14SB/GBn2 single-point AMBER, the shipped
distogram Bayes risk. *Not taken:* (a) a **fitted** blend weight `w` — s24 L16 closed it, an
oracle `w` with full leakage on the very targets it is scored on is worth 0.0148 Å; (b)
**interaction-only AMBER** (`nonbonded + solvation`), which project memory records as the right
way to score AMBER but which requires `k_restraint=10, steps=0` minimisation and is therefore
forbidden by s23 L11 and by this lane's brief; (c) any Legacy term subset or re-weighting; (d)
any learned surrogate for either energy (none is defensible — the bonded subset tracks at
ρ = 0.211 against a 0.7 bar).

**2. BASIS.** *Taken:* point-cloud Cα-RMSD of a coordinate average against `nat_ca`. *Not
taken:* built-chain / ideal-geometry-projected RMSD (`d_harness.readout_projected`). The gap is
0.156 Å of pure operator choice; the two are never placed in the same table.

**3. READOUT.** *Taken:* uniform coordinate average over the CVaR tail, α pinned so the median
tail is the production rung of 75. *Not taken:* (a) the `p_theta`-weighted average — s23 L8
measured it worse at 3 of 4 temperatures; (b) the consensus medoid, which is the operator
`core/pipeline.quantum_stage` deploys but which would not anchor on the 3.0483 Å endpoint this
sprint is measured against; (c) an argmin / small-m readout, which `structural-objective-beats-
the-energies` shows would change *which* w is optimal and is a different question.

**4. NORMALISATION.** *Taken:* rank standardisation, equation (2), with an outer re-standardisation
of every combination. *Not taken:* (a) **raw moment z** — run to completion as a declared
secondary rather than silently dropped, because 57.5% of candidates exceed 1e4 kcal/mol, the
worst single point is 7.1e18, and 99.4% of a pool lands inside |z_raw| < 0.1, so moment-z
degenerates into "which candidate has the worst clash" (`pauli-spectrum-delta-spike-artefact`,
s24 L10); (b) **99th-percentile winsorisation**, which memory records as *not enough* — a
winsorised AMBER table still carried 0.986 of its variance in ten configurations; (c) signed-log
conditioning, monotone and therefore ordering-equivalent to rank on the classical arm but not on
the VQE's energy gaps, so it would be a third arm without a distinct question; (d) no
normalisation, which would make α and T incomparable across configurations.

**5. NULL.** *Taken, four of them, every one matched in the operator's space:* (a) classical
top-75 on the **identical** energy, run for **every** configuration as a full parallel table —
never as a substitute for a VQE row; (b) classical top-m at each configuration's **own realised
m** (the size-matched bar); (c) random m-subsets through the **identical** coordinate average,
16 draws, mean as the null and max/min kept so any best-of-K claim is priced against the
distribution of the extremum; (d) **rank-permuted energy** per configuration — AMBER's/Legacy's
marginal preserved exactly, its correspondence to candidates destroyed — which is the control
that caught s24 L16's real arm as indistinguishable from noise. *Not taken:* a uniform-on-the-
torus or synthetic-structure null (`zero-information-control-must-be-plausible`: uniform is a
*worse* measure, not an uninformative one); a comparison against an initialisation mean instead
of best-of-N (`concentration-is-wrong-when-discrimination-binds`).

**6. THE LABEL.** *Taken:* **mean full-chain point-cloud Cα-RMSD over the 126 dev targets, on
the VQE arm, paired per target.** Every RMSD is an ORACLE evaluation of an ACHIEVABLE
(native-free) selection; nothing is PRODUCTION-labelled unless it is the shipped configuration.
*Not taken:* (a) a correlation label — ρ(E, RMSD), in-band ρ, native percentile — which are
landscape diagnostics in Deliverable 2 and are explicitly *not* the label, because
`structural-objective-beats-the-energies` measured global ordering and argmin quality moving in
**opposite** directions along the same axis; (b) a "best of the seven" minimum, which is an
order statistic; (c) the bias cosine, which s24 D1 used because n = 30 could not clear an RMSD
MDE — at n = 126 that reason does not apply and using it would be avoiding the endpoint.

---

## 3. WHAT IS NOT DONE, AND WHY

* **No minimisation, no relaxation, no repair, anywhere in the production path.** s23 L11
  measured 17 of 17 restrained and unrestrained repair settings at or worse than doing nothing,
  with the ladder bottoming out at the no-op. AMBER is an **energy measurement** and a
  feasibility diagnostic in this lane and nothing else.
* **No re-derivation of AMBER energies at OpenMM cost.** `s24/cache_amber/` holds 126 × 500 =
  63,000 genuine single points with `amber_verify_max_rel = 0.0` on every target. Integrity is
  re-checked (§4 F0) and then the cache is used.
* **No candidate generation, no pool change, no prior change.** All three are on the
  DO-NOT-REDO list and none is this lane's question.

---

## 4. FALSIFIERS, REGISTERED BEFORE THE RUN

**F0 — CACHE INTEGRITY (gate, not a result).** Five fold-stratified targets × 4 candidates are
re-scored with a **fresh** `AmberSP` single point under `LOCK_AMBER` and compared to the cached
`e_amber`. Required: max relative difference **0.0**. Alongside, on all 126 targets: the pool's
`universe_idx` must reproduce bit-for-bit and `score_shipped` must reproduce the cached
`score_dist` bit-for-bit. If any of these fails, nothing downstream is read.

**F1 — THE ANCHOR.** Configuration 3 (Distogram) under the **classical top-75** control must
return **3.0483 Å** to four decimals. `zrank` is monotone, so w=0-style exactness is available
and is asserted, not approximated. If it does not, the harness is mis-wired and the suite stops.

**F2 — THE THEOREM.** `gate_set_equality` (subset-hood) must pass on all 7 × 126 = 882 cells.
Equality is *reported*, never asserted. Any subset-hood failure means a mis-wired label, a
padding leak, or the two arms seeing different energies — a bug, and the suite stops.

**F3 — MY REGISTERED PRIOR ON THE OUTCOME. I EXPECT EVERY AMBER CONFIGURATION TO FAIL.**
I predict, before looking:
  (i) C3 (Distogram alone) is the best single configuration or statistically tied for best;
  (ii) **every configuration containing AMBER is null-to-worse than its AMBER-free counterpart**
       (C5 ≥ C3, C6 ≥ C1, C7 ≥ C4), because ρ(AMBER, distogram) = −0.0270, because nested CV in
       s24 L16 chose a **negative** AMBER weight on 4 of 5 folds, and because the AGREE_BAD
       partition (both energies reject the candidate) came back *better* than AMBER_PREFERS;
  (iii) C4 (Legacy+Distogram) is worse than C3, because ρ(Legacy, distogram) = +0.3875 — Legacy
       is a compactness/typicality prior partly aligned with the thing that already selects, so
       it adds correlated noise rather than a new direction;
  (iv) C1, C2, C6 (no distogram) are 0.5–1.5 Å worse than C3.
**If any AMBER-containing configuration beats its AMBER-free counterpart past its own MDE with
the fold CI excluding zero, my registered prior is refuted and that is the sprint's positive
result.** It would also contradict s24 L16's oracle ceiling and would have to be reconciled with
it before being believed.

**F4 — THE NORMALISATION ARM.** Under the raw-moment secondary I predict configurations 2, 5, 6
and 7 move materially (top-75 set overlap with the rank arm well below 1.0) and get **worse**,
reproducing s24 §4 on a new instrument. If raw-moment is BETTER, the rank choice is wrong and
must be revisited before anything else in this lane is believed.

---

## 5. DELIVERABLE 2 — THE LANDSCAPE COMPARISON

Reuse-and-cite where the pool has not changed (s20's Hessian, gradient, anisotropy, ruggedness
and participation-ratio measurements on the continuous torsion manifold; s24's conditioning
census). **Re-measure only what the pool or the score can plausibly have changed**: energy
scale and candidate-energy distribution on this K=500 pool, ρ(E, distogram), ρ(E, RMSD) global
and in-band, native percentile, compactness response with an Rg control, CVaR concentration and
realised tail size, and structural usefulness (the endpoint itself). Each row carries a one- or
two-sentence physical interpretation. Every standing fact the brief lists is **verified on this
instrument**, not assumed.

---

## 6. PROVENANCE

Every artefact written with `ST.save_atomic(..., complete_keys=..., rows=..., n_expected=126,
module_file=__file__)`. `LOCK_AMBER` announced on open and on release, held for the F0 spot check
only. ORACLE / ACHIEVABLE / PRODUCTION labelled at every appearance.

---

# ADDENDUM — FORM 5. PRE-REGISTERED IN FULL BEFORE THE RUN. Authorised by the coordinator.

## A1. The question

The functional lever has failed as a **filter** (s24 D1-C), a **partition** (s24 D1), a **fitted
score** (s24 L16, oracle ceiling 0.0148 A) and an **equal-weight score** (s25 suite SS1.3-1.5).
Form 5 is the one remaining class: the physics never chooses a candidate, never splits the pool
and never enters the ranking — it supplies a **per-target SIGN at inference**, which is the only
leverage `in-band-ordering-is-per-target` leaves open.

## A2. THE RULE — PRIMARY, exactly as written in agentPHYS_FINDINGS SS5 and unchanged

For target `t`, let `T_t` be the DISTOGRAM's own top-75 set (configuration C3's retained set):

    s_t  =  mean_{i in T_t} [ zrank(E_LEG)_i  -  zrank(E_AMB)_i ]          NATIVE-FREE

Both channels are rank-standardised over the full K=500 pool before the mean is taken, so `s_t`
is dimensionless, scale-free and comparable across targets. The arm switches between **C3
(Distogram)** and **C5 (AMBER+Distogram)** by thresholding `s_t`.

**FITTING: nested leave-one-fold-out, threshold AND direction fitted on the OTHER FOUR FOLDS
ONLY.** No per-target fitting. No target is ever scored under a threshold or a direction its own
fold helped choose. The candidate thresholds are the training folds' own `s` quantiles at
5%..95% in 5% steps, plus the two degenerate ends (always-C3, always-C5) so the fit can decline
to switch. Both directions (`s > tau -> C5` and `s < tau -> C5`) are on the grid.

## A3. THE FALSIFIER — unchanged, and it is the whole test

> The switched arm must beat C3 by **MORE THAN ITS OWN MDE (2.8016 x SE)** with the
> **fold-clustered CI excluding zero**, under nested LOFO CV. If it does not, **Form 5 is closed
> and the functional lever is closed in all FIVE of its forms.**

## A4. NULLS AND CEILINGS, all registered

* **N1 — permuted signal.** `s_t` shuffled across targets with a stable RNG: marginal preserved
  exactly, correspondence to targets destroyed, the identical nested CV re-run. This is the
  control that decides whether any movement is about the signal or about the fitting procedure.
* **N2 — degenerate arms.** always-C3 and always-C5 reported, so a "gain" that is really just
  the grid drifting toward one endpoint is visible.
* **C1 — ORACLE ceiling.** per-target `min(C3, C5)`, full leakage. Reported as ORACLE, never as
  a rule, and priced with `k_eff` because a best-of-K null is INAPPLICABLE when the K arms are
  correlated (the suite's own seven arms have `k_eff = 1.13` against a nominal 7).
* **M1 — mechanism check, ORACLE DIAGNOSTIC.** the correlation between `s_t` and the per-target
  advantage `C5_t - C3_t`. This reads the native and is therefore diagnostic only; it can
  explain a result but can never produce one.

## A5. DECLARED SECONDARY, with its multiplicity stated

A second native-free signal from the same mechanism — the scale-free compactness contrast

    r_t  =  mean_{i in T_t} Rg_i  /  mean_{i in pool} Rg_i  -  1

— is run through the IDENTICAL nested CV. **It is NOT the registered falsifier.** Two signals
are tested, so a secondary that clears the plain bar must clear a Bonferroni-corrected one
(alpha/2) before it is believed, and that is stated here rather than decided afterwards.

## A6. MY PRE-STATED EXPECTATION: IT FAILS. Three quantitative reasons, all on the record.

1. The ORACLE ceiling on ANY per-target switch among the suite's seven arms, with full leakage,
   is **-0.3349 A**, and `k_eff = 1.13` says the seven are effectively one arm — so most of that
   is noise on nearly-identical draws.
2. `decorrelated-errors-exist-but-are-unusable` prices a fusion gain as the **SQUARE** of the
   weaker channel's skill, and AMBER's in-band rho is **+0.0913**.
3. Four prior forms failed, and the equal-weight form failed with **no free parameter at all**.

**Registering a loud expectation of failure is the point, not a reason to skip:** a closure claim
in five forms is worth materially more than one in four with the fifth untried.

## A7. CONSTRAINT

The architecture is **FROZEN**. Form 5 cannot change what ships unless it clears A3 decisively,
in which case the result goes to the coordinator and the freeze is theirs to reconsider. Nothing
here assumes either way.
