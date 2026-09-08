# SPRINT 15 — PAPER-DRIVEN QUANTUM PEPTIDE FOLDING PROGRAM. Shared agent brief.

Read this completely before writing code. It is the contract.

---

## 0. WHAT THIS IS

Not a tuning sprint. A research program whose output should support a **peer-reviewed paper**.

The final work must answer:

> **What is newly learned about quantum variational conformational optimisation that another
> researcher did not know before this study?**

Structural goals: **< 2.5 A** is a major success, **< 2.0 A** exceptional, on the 126-target
development instrument. 3.0-3.2 A is **not** a success unless it carries a major mechanistic
discovery. A sub-2 A result resting on an oracle, leakage, target removal or a physically
invalid structure is worth nothing. **A rigorous negative that closes an important open
question is a legitimate success tier.**

---

## 1. NON-NEGOTIABLES

**The final system must contain all four as real scientific variables, not decoration:**
genuine **VQE** (real parameterised state, real variational optimisation, real sampling),
genuine **CVaR** (correct quantiles, ties, finite shots, gradients, bias), a genuine
**Legacy** non-all-atom energy, and genuine **AMBER** ff14SB all-atom.

**Benchmark protection is absolute.** `results/benchmark_manifest.json`, `peptide_folds.json`,
`peptide_clusters.json` are pinned. The 60-target benchmark is **inaccessible until the
architecture and protocol are frozen**. No peeking, no benchmark-adjacent tuning, no "just
checking one target".

**NO NATIVE INFORMATION AT INFERENCE.** Native coordinates, torsions, distances, RMSD,
energies, contact maps and secondary structure may be used **only** for training under
leave-fold-out discipline, for experiments explicitly labelled **ORACLE DIAGNOSTIC**, and for
post-hoc evaluation. Any quantity derived from native data inherits `ORACLE = TRUE` and may
never enter a predictive headline table. Use `s14.ladder.audit_emitter`, which poisons the
native trace and demands bit-identical output — it is mechanical and it works.

**No fake components.** A classical search is not VQE. Enumeration is not VQE. A surrogate is
not AMBER. If you run a classical control, say so in the arm name.

**Preserve history.** Never delete a superseded claim. Write the correction beside it. The
corrections ledger is the most valuable thing this project produces.

---

## 2. WHAT SPRINT 14 ESTABLISHED. These are binding constraints, not doctrine.

Re-testing one requires a specific new hypothesis and a stated scope violation, not a rerun.

### The arithmetic

| quantity | value |
|---|---|
| incumbent retrieval pipeline | **3.204 A** (`synthesis_fit`) |
| in-band ordering accuracy required for 2.0 A through a top-100 operator | **0.638** |
| measured: Legacy / learned objective / 1-local prior / AMBER | 0.539 / 0.523 / 0.501 / **0.463** |
| per-torsion sigma required for 2.0 A (i.i.d. error) | **15.1 deg** |
| best torsion channel available (retrieval top-75 circular mean) | **phi 33.6, psi 59.2 deg** |
| sigma-equivalent of the incumbent | 27.1 deg |
| torsion-space content at k=4, ORACLE descent from random start | 1.982 A |
| retrieval pool best member | 1.711 A |

**Both the generation side and the selection side are roughly a factor of two short.**

### The five results that constrain everything

**1. Running VQE is WORSE than not running it.** Against the only fair control — best of the
same number of shots from the **untrained** circuit — optimisation loses by +0.65 to +1.32 A
at **0/12**, including from ORACLE warm starts, and it **erases its own initialisation**
(after-VQE lands 2.42-2.76 A regardless of a start spanning 2.97-4.13 A). On the objective
axis VQE reaches the certified optimum in 0-12% of cells against greedy 1-opt's 71%.

**2. The mechanism.** Optimisation concentrates the distribution onto the objective's
low-energy tail, and **every objective is at chance or worse inside its own tail** (0.512 for
the prior; the best *bulk* ranker anti-ranks at **0.390** in its own lowest 0.1%). The argmin
sits at the **41st-48th percentile of its own tail** — a random draw.
**`tail mean - tail best` = 1.93 A**, against a certified selection gap of 1.876 and a sampled
2.040. Three measurements, one quantity. **Concentration is the wrong move when
discrimination binds.**

**3. The deepest limit.** In-band ordering is learnable to **0.986 WITHIN a target** with an
overfitting gap of **-0.0005 [-0.0011,+0.0000]**, and collapses to **0.600 ACROSS targets**
(transfer gap **-0.386 [-0.437,-0.332]**, 770x larger). A **linear** pair potential already
saturates the within-target problem, so no capacity, nonlinearity or equivariance helps.
Per-target skill correlates **+0.909** with the native radius of gyration. **The ordering axis
is real, learnable, and per-target; its sign is a property of the target that inference cannot
observe.** The learning curve **declines** — one or two training targets beat twelve.

**4. Selection decomposes exactly.** `sel = pool + FILTERING + ORDERING`. Legacy: filtering
-0.124, ordering **+0.263**, net **+0.139** — exactly the Sprint 13 certified-optimum-minus-random
figure. **The ORDERING term's sign IS "optimise harder, get worse", measured directly.**
Legacy is a clash gate (98.8% one term, constant across its own decile). AMBER is the best
filtering term measured (-0.141) and an inert refiner (+0.014 A [-0.001,+0.028]).

**5. Aggregation beats the objective 3.4 to 1.** Coordinate averaging beats torsion averaging
of identical windows by **+1.024 A**; the whole contribution of any objective is 0.171 A.
Averaging **contracts the backbone by 25.8%** (rho(contraction, projection cost) = -0.552),
  **[CORRECTED 2026-09-05: measured directly on the production consensus, the contraction is
  3.5% against the TRUE distances and 9.5% against the predicted ones — not 25.8%. The direction
  and the rho(contraction, projection cost) = -0.552 relationship stand; the magnitude does not.
  See `coord_FINDINGS.md` K11. Do not build on the 25.8% figure.]**
and restoring valid bonds costs 0.157 A by any route. A *perfect* ranker emits 2.474 A through
a decile and **1.219 A through a top-100** — the terminal operator is worth 1.26 A on identical
rankings.

### Closed, with reasons

- **Chemical shifts as a generation channel.** 54/126 targets have usable heteronuclear
  shifts; **ORACLE-perfect torsions on all 54 still leave the instrument at 2.021 A** and 55
  are needed. Closed by arithmetic, independent of predictor quality. Torsions from
  *predicted* shifts are 1.0 A **worse** than using the sequence directly.
- **A learned many-body objective.** -0.166 A, CI crosses zero, reverses on drop-top-3, both
  nulls worse than random, leaked-label control loud at -1.405 A on 19/19.
- **Sequence conditioning of a learned objective is HARMFUL** — the sequence-blind twin wins.
- **Fusing Legacy and AMBER.** Truth-partialled error correlation +0.096 (genuinely
  decorrelated) but fusion is +0.032 rho under a per-target ORACLE weight and **loses 0.240 to
  leave-one-target-out**. Pareto discards 72% of the truly-best 1%.
- **The published log-scale continuous-torsion encoding** (arXiv:2609.02113, QTF): its energy
  is never a qubit observable and `P ~ 2N+6n` — 132 classical parameters for 43 torsions.
  A reparameterisation, not a compression. Adopting it would breach rule 1.
- **In-band ranking over retrieval pools** (Sprint 12): flat learning curve against a leaked
  label loud at n=8. Signal-limited, not sample-limited.

### Traps this project has actually fallen into. Do not repeat them.

- **A weak control.** A zero-information constant alpha-helix (phi=-63, psi=-42) emits 4.065 A
  and **beats uniform random by 0.457 A**. "Beats random" proves nothing. Clear the helix.
- **The proposal confound.** An in-decile rank correlation of +0.370 collapsed to **+0.046**
  when the proposal stopped being drawn from the prior being scored. Always state the proposal.
- **In-decile rho is NOT monotone in objective quality** — it peaks at moderate signal and
  falls where the objective is much better, because a better decile is narrower. **Report
  global rho, in-decile rho, argmin and decile mean together, always.**
- **Bulk pairwise accuracy does not predict argmin quality.** The best bulk ranker has the
  worst optimum. Only **tail-restricted** accuracy tracks it — read against a **random-tail
  null at 0.524-0.527**, never 0.500 (a gap-matched/binned variant returns to a true 0.500;
  do not conflate the two).
- **A raw drop-top threshold is NOT a valid concentration test.** When mean/sd is small,
  discarding the ten best removes much of the total **even if every target carries an identical
  effect**. Compare every concentration statistic to a **simulated uniform-effect null**, and
  print **mean/sd** so a reader can see when the test has no power. This rule caught one real
  artefact and misfired on one real result in the same sprint.
- **Reading past a check.** `I.paired` returns `drop_top10_mean_diff` and `top10_share` in the
  same dict as the W/L. A claim was published on the reassuring field. **Emit concentration
  statistics as one PASS/FAIL verdict**, never as fields a reader can select from.
- **A near-even W/L with a CI excluding zero** is suggestive of a *concentrated* effect and a
  prompt to run the null-calibrated check. The **median-vs-mean ratio** flags it for free.
- **n <= 4 samples from the enumerated nine have reversed a conclusion four times.** 2MK7 is a
  recurring outlier.
- **The cached AMBER subset in `s13/results/qarch_enum_*.npz` is 21.7% ORACLE-CONDITIONED**
  (0.401 A better than its space). **Only `amber_kind == 0` rows are an unbiased sample.**

### Corrections to older records, already applied

- **FOUR torsions per chain are inert** for the CA trace (`phi[0]`, `psi[0]`, `phi[n-1]`,
  `psi[n-1]`) — both terminal residues invisible, dead block `2*log2(k)`, **14 of 18 live** at
  n=9 k=4, **21.9 mean live qubits at k=4, not 25.9**. Follows from the locality theorem.
- **The locality theorem**: `d_ij` depends on exactly the `j-i-1` residues strictly between i
  and j. A separation-8 pair is a **14-qubit interaction at k=4**; **no 2-local Ising form of a
  distance-based objective exists in this encoding**. All-atom: CA `i<m<j`, N `i<=m<j`,
  C/O `i<m<=j`, CB `i<=m<=j`.
- **QNG is refuted ONLY at depth 1**, for every entangler. At depth 2 the condition number
  reaches 155.2 for a torsion-aware block entangler and 478.9 at depth 3. **QNG is open again
  at depth >= 2** — whether it *helps* was never measured and is a HYPOTHESIS.
- **CVaR**: value estimator correct to 1.7e-15. Three defects: the known `baseline="tail"`
  gradient bias now has a **closed form** `-c * grad P(E<q)` (the recorded "cosine -0.023" is
  one draw from a +0.06 to +0.96 spread — quote the closed form); a **sampled-CVaR upward bias
  at non-integer `alpha*N`** (+0.134 sd at N=13, alpha=0.1, decays as 1/N); and
  **`dCVaR/dp` identically zero iff `p(argmin E) >= alpha`** (0/3000 counterexamples) — at
  alpha=0.01 one initialisation in four to twenty **starts dead**.
- **`alpha` is not a learning rate**, and at small alpha CVaR **collapses to an argmin-finder**.
  The recorded "small alpha collapses two AMBER variants" has **the sign backwards**: small
  alpha discriminates better; **concentration** collapses them.
- **Binary encoding**, for four independent reasons: it attains the information-theoretic
  qubit bound at power-of-two k; needs **510 Pauli terms where one-hot needs 454,463**; is
  surjective so needs no penalty; and has the highest gradient variance tested. Gray buys
  nothing. **A non-surjective encoding's Pauli spectrum measures its CONSTRAINT, not its
  objective.**
- **The budget trap is a property of BAD OBJECTIVES.** Searching harder makes structures worse
  on a bad objective, is neutral on a mediocre one, and **helps monotonically on a good one**
  (it inverts by signal 0.50). Always state the objective-quality condition.

### The one open thread Sprint 15 inherits

**Anything that supplies the PER-TARGET SIGN of the in-band ordering axis at inference.** The
axis is compactness. Native-free proxies for the native radius of gyration reach
**0.24-0.37** (length-residualised, all CIs excluding zero) against the ORACLE's 0.909 — real,
free to compute, weak. **But the end-to-end chain does not reproduce**: the ORACLE ceiling
measures +0.416 (p=0.177) on one skill definition against +0.909 on another. **Resolving that
discrepancy on a common definition is the first job**, before building on it.

---

### PHASE 0 AUDIT RESULTS — BINDING ON EVERY AGENT, added 2026-09-05

The Phase 0 provenance and instrument audit is complete (`s15/audit_FINDINGS.md`, 25 JSON
artefacts, nine `s15/audit_*.py` modules). **Phase 0 gate: PASS**, with one retraction, one
correction to this brief, one new defect, and four disclosures. Everything below overrides any
earlier statement in this file or in a Sprint 14 record.

**What reproduces bit-identically from primary inputs.** `peptide_db.npz` rebuilt from 1,524 PDB
files (787/787 records, max diff exactly 0.0); all 126 window universes (9 arrays × 126 targets,
2.35 M windows); `pool_best = 1.7108244199364904` from the PDB files; `top75_best` and
`n_zero_recall` without the pipeline cache; `shipped = 3.4540004952559396` at torch 1, 2 and 8
threads with zero argmin flips; `synthesis_fit = 3.2040761603809194` end-to-end including
`amber_ca`; the AMBER golden −489.9138948277905 with ff14SB/GBn2 parameters verified against an
independently built ForceField at max diff 0.0; Legacy's total as the weighted sum of its 11
components (worst 2.7e-5 over 1.28e7 structures); and CVaR against Rockafellar–Uryasev (max error
4.2e-14 over 400 cases).

**`python -m s12.instrument` is a CACHE READ, not a reproduction.** In its 2.2 s it recomputes
only the Kabsch step, and it asserts only three of its five constants — `shipped` and
`synthesis_fit` are printed unchecked. Do not describe running it as verifying the pipeline.

**The distogram is NOT bit-stable under `torch.set_num_threads`.** 1/2/8 threads give three byte
patterns (`prob` up to 2.0e-6, `risk` up to 1.3e-4). The cache was written at torch's default 8;
this brief mandates 2, so *following the brief changes the bits*. It moves no headline constant,
but any claim of bit-identical reproduction must state the thread count.

**RETRACTION — every AMBER row in `s14/results/obj_floor.json`.** It computes on the full
labelled subsample with **no `amber_kind` mask** (0.42 Å oracle-conditioned), carries no
stratification string, and is flagged `"complete": true`. Corrected values are in
`s14/obj_FINDINGS.md` §1.3; the artefact was never regenerated. Do not read that file. Every
other consumer (s13 `qarch_robust`/`sepfit`/`validity`/`locality`, s14 `ener_*`, `vqe_collapse`)
is correctly masked.

**CORRECTION TO THIS BRIEF — the oracle-conditioning share is TRANSPOSED above.** Measured on all
19 enumerated files: oracle-conditioned (`kind == 2`) is **21.7%**; **40% is the *unbiased*
share**; 39% is prior-conditioned. The −0.401 Å effect size is right (−0.4009 on the 9, −0.4243
on 19), and "use `amber_kind == 0`" remains correct. The Sprint 14 modules and findings label the
strata correctly — only the one-liner in §2 of this brief was wrong.

**NEW DEFECT, in no prior record — `amber_kind == 0` is unbiased on the mean but not in the
tails.** `a_idx` force-includes `snap_idx = ORACLE_snap()`, which lands in the "unbiased" stratum
on 16/19 files. There it is the **minimum-RMSD member on 6/19**, in the **RMSD-lowest 1% on
10/16**, and **inside the <1.5 Å in-band set on 14/19** — touching 40% / 25% / 17% / 15% / 11% of
all in-band pairs on 2MJQ / 7VI4 / 2MK7 / 8HVS / 5V5B, whose bands hold only 5, 8, 12, 13 and 19
members.

> **BINDING RULE: use `amber_kind == 0 AND amber_idx != snap_index` for every tail, in-band,
> argmin or top-k statistic, and print the per-target n beside it.**

The direction is safe — the headline 0.463 for AMBER can only have been flattered, so Sprint 14's
negative stands and is if anything understated — but the number is not clean.

**`BAND = 1.5 Å` has no derivation anywhere, and `I.FAIL18` is a threshold artefact.**
`n_zero_recall` runs 45 → 2 across BAND 0.5 → 3.0. On Sprint 12's own 99-combination sweep,
**exactly one of the 18 FAIL18 targets (9KAR) appears in every version of the set**; mean Jaccard
0.498. FAIL18 has been treated as an object across three sprints. **Keep reporting the FAIL18
column — it is comparable to every prior sprint — but never argue from set membership, and never
describe FAIL18 as "the hard targets" without this caveat attached.**

**Other load-bearing constants.** `K = 500` is load-bearing for `pool_best` (1.970 → 1.504 across
K = 100 → 2000) but *not* for `shipped` (3.425 → 3.520, flat). `m = 75` moves `top75_best` 2.609 →
2.106 across m = 25 → 150. Both were fitted on the same 126 targets they are reported on — state
this. `IDENTITY_THRESHOLD = 0.6` is reused without justification, and the repo already records 573
library members that contain a target at ≥ 0.6 while passing the filter. Undocumented and
load-bearing: `min_sep = 2` (no sweep exists), `torsion_window = 8`, the AMBER tolerance of 1.0,
and the `bond + angle > 1000` strain gate. **Legacy's 11 weights were never fitted** — declared "a
variance-balanced starting point, not a fit", which is a strong claim to defend in a paper.
Pre-registration mismatch: `s8/project_devarm.json` pins `phip@0.03` while production ships
`ramah@0.3`; worth 0.004 Å on RMSD, but every backbone-plausibility claim rests on the post-hoc
arm.

**RMSD is frozen, and independently reimplemented.** Horn's quaternion method (no SVD, proper
rotation by construction) agrees with `kabsch_rmsd_batch` to **2.04e-13 Å over 63,000 real pool
structures**. All 15 edge cases pass, including the mirror test (both give 2.9602 where the same
code without the determinant fix gives 0.0). **Frozen definition: CA only; every residue including
both termini; no trimming; correspondence by index; uniform weights; proper rotations only;
ångström; model-1 native; chain breaks scored, not rejected.** Full-chain versus `[1:-1]` differs
by up to **2.4 Å** on one structure, so the choice must be stated explicitly in the paper.

**AMBER costs less than this brief says.** Memoisation defeated and hit-counted: **8.3 ms at 127
atoms rising to 23.3 ms at 237 atoms**, so **8–14 ms on the n = 9–10 enumerated targets**, not
28 ms. Budget-matching that charged 28 ms **over-charged AMBER by 2–3×**. A repeated structure
costs 0.33 ms (a 33× memo speedup). `s13/qarch_lib.amber_energies` still says "~6 ms" in its
docstring; both figures are wrong.

**Artefact hygiene, for the reproducibility document.** 325 of 363 result JSONs carry no
completeness flag; of the 38 that do, **all 12 "incomplete" flags are false alarms** — `I.write`
takes `len()` of a dict of *arms*. A flag with a 12/12 false-alarm rate is worse than none; do not
cite it as evidence a run finished. Five Sprint 14 results exist twice (`s14/results/X.json` and
`s12/results/s14_X.json`) with no canonical copy. `verify/amber_audit.json`'s
`5_translation_invariant: 0.5446` measures the **minimiser**, not the physics (single point
1.3e-6; minimised 0.148 kcal/mol) and has no PASS/FAIL. **37.6% of cached AMBER single points
exceed 1e6 kcal/mol**, 3.2% exceed 1e12, max 2.6e20 — the 99th percentile is itself 1e12–1e15, so
winsorising there is not enough. Sprint 14's flagship ablation ladder (4.072 → 3.485 → 3.314 →
3.204) traces cleanly but to four different files, one misleadingly named
`vqe_classical_limit.json`, with no artefact holding the ladder itself.

**Leakage disclosures, not retractions.** 16/126 universes contain a window that is verbatim an
own-fold database member's full sequence, and on **4 targets the top BLOSUM hit is the target's
own sequence**. Measured impact on `pool_best`: **exactly zero on all four** — they are different
structural determinations, 0.59–4.13 Å from native. Disclose in the paper; do not remove targets.

### CORRECTION TO A COORDINATOR MODULE, 2026-09-05

**`s15/distml.py`'s `LogPTable` had a uniform-grid bin lookup on a NON-UNIFORM grid.** Bin index
was `floor((d − c[0]) / (c[1] − c[0]))`, but the distogram's centres run 0.5 Å apart at short
range and up to 4.0 Å apart at long range, so every query above ~5 Å read the wrong bin (8.50 Å
read the 7.25 Å bin; 15.00 Å read the 17.50 Å bin; everything above 15.25 Å collapsed into one).
**Fixed**: `searchsorted` on the real centres with per-interval widths, derivative re-verified at
cosine 1.000000000000000 against central differences. **Void**: every `ml_*` arm in
`s15/results/distml.json`, the `pool_hist` / `pool_sim_weighted` / `combined` arms of
`s15/pooldist.py`, and the first `s15/cascade.py` run (killed at 40/126 and relaunched).
**Unaffected**: `s15/distgeo.py` and `s15/distcal.py`, which are weighted least squares and never
touch the class. See `s15/coord_FINDINGS.md` K4-RETRACTED. **If your module imports `LogPTable`,
`ShiftedLogP` or `SumLogP`, re-run anything you computed before this correction.**

---

## 3. THE MACHINE

15.6 GB RAM, often only **2-3 GB free** (shared with a browser and editor). 8 logical cores =
**6.43 core-equivalents** (4 fast + 4 slow). Sprint 14 saturated it: at the worst point agent
processes were measured at **3-5% of one core**, and one workstream reported 0.07 cores.

- Cap threads at import: `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2`.
- Gate large allocations with `s13.qarch_lib.wait_for_memory` or `core.quantum.mem_gate`.
  **Wait for headroom — never catch a MemoryError and return NaN.** A NaN column silently
  becomes "no signal", which is a fabricated negative.
- **Checkpoint incrementally.** Sprint 14 lost an entire ansatz study because it wrote its
  JSON only on completion. Write partial results every N targets.
- Do NOT load `esm_cache.npz` (1.5 GB); use `s12/esm_bank.py` (100 MB).
- Costs: AMBER single point **8-14 ms distinct on n=9-10 targets** (0.33 ms memoised), `refine_coords` **9.1 s**, Legacy 0.3 ms,
  a projection ~5-20 s.
- Two heavy jobs fit. Three do not. Say so if your job is heavy.

---

## 4. MACHINERY THAT EXISTS. Use it; do not rebuild it.

**`s12/instrument.py`** — the shared instrument. `python -m s12.instrument` must reproduce
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`. **Run it at start and end.**
`I.targets()`, `I.load_univ(pdb)` (W/PHI/PSI/S/org/sim/order/rr/nat_ca), `I.build_ca`,
`I.ca_rmsd`, `I.kabsch_rmsd_batch`, `I.project`, `I.coordinate_average`, `I.distogram`,
`I.shipped_score`, `I.pair_index`, `I.pair_dists`, `I.paired`, `I.summary`, `I.write`,
`I.FAIL18` (report the column, but see the BAND caveat in the Phase 0 audit results below).

**`s13/qarch_lib.py`** — `Space(pdb, k)` (sequence-conditioned, target held out), `PHI`/`PSI`
(n,k) tables, `empirical_prior`, `legacy_components`, `amber_energies`, `wait_for_memory`.

**`s13/results/qarch_enum_<PDB>.npz`** — 9 targets fully enumerated at k=4 (262,144 configs),
plus **10 more n=10 targets from Sprint 14** in `s14/cache/` — **19 targets, 1.28e7 exactly
labelled structures**. Keys: `rmsd`, `legacy`, `prior`, eleven `leg_*`, and `amber_total` +
five `amb_*` on a subset that is **21.7% oracle-conditioned** (use `amber_kind==0` AND
`amber_idx != snap_index`; see the Phase 0 audit results below).

**`core/quantum.py`** — `StatevectorCircuit`, `cvar_exact`, `grad_cvar_paramshift/fd/score`,
`run_cvar_vqe`, `MPSAnsatz`, `Adam`, `mem_gate`. **`qansatz.cvar_gradient` with
`baseline="tail"` is the known-biased one; `baseline="const"` is correct.**

**Sprint 14 modules worth reusing:** `s14/ladder.py` (emitters + the leakage guard),
`s14/retprior.py` (retrieval-conditioned torsion prior, the best channel measured),
`s14/hamil.py` (native-free structural objective, `Terms` with `e_prior`/`e_disto`),
`s14/ener_lib.py` (`decile_rho`, `pair_accuracy`, `argmin_rmsd`, `tail_accuracy`,
`selection_decomposition`, tie-averaged), `s14/ener_decoy.py` (cached decoys + matched nulls),
`s14/obj_ceiling.py`, `s14/vqe_tailacc.py`, `s14/signpred.py`.
Normalisation constants: `s14/cache/ener_norm.json` — **use `grad_sd`, never raw AMBER's sd**.

---

## 5. STATISTICS AND REPORTING

Every comparison: **paired mean difference, bootstrap 95% CI, median, mean/sd, W/L, per-fold
values, and a NULL-CALIBRATED concentration check** emitted as one PASS/FAIL verdict. Use
`I.paired`. Every attractive correlation gets a null — permutation, label shuffle, matched
random arm, or sequence-blind twin.

**Report the cascade for every architecture, never a single final number:**

    G (generation) -> S (selection) -> A (aggregation) -> F (final)

and the gaps `G-S`, `S-A`, `A-F`. **The largest gap is the next research target.**

**Budget matching: report BOTH** (A) equal objective evaluations, and (B) equal computational
cost. State which is matched. **If the conclusion changes between them, report the conflict.**
Never pick the favourable convention after seeing the result.

**Tier every claim**: DEMONSTRATED / ORACLE DIAGNOSTIC / LITERATURE-SUPPORTED / HYPOTHESIS /
REFUTED. Distinguish **exploratory** from **confirmatory**. Do not present the best of dozens
as a pre-specified primary result.

**Independent verification** for every headline claim: an independent implementation or a
mathematically independent check for numbers, an independent falsification experiment for
mechanisms. **Two copies of the same code are not replication.**

---

## 6. HOW TO WORK

Write reusable modules under `s15/`, prefixed with your workstream tag. Write findings
**continuously** to `s15/<TAG>_FINDINGS.md` — the coordinator reads them to schedule.

Use the Write tool for file creation; heredocs break on apostrophes and on `\n` inside Python
strings, and it has cost this project hours. Set `PYTHONIOENCODING=utf-8` and keep console
output ASCII.

**Do not send interim reports to the user.** Report to the coordinator when you have a result
that changes a decision, or when you finish.

Your final message: what you measured with numbers and CIs, what you refuted **including your
own hypotheses**, what remains open, and exact reproduction commands.
