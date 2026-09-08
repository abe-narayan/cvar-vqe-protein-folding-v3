# AGENT D — THEORY / LITERATURE / ADVERSARIAL REVIEW
## Sprint 20 findings

Pre-registered in `s20/PREREG_D.md`, written before the corresponding results were inspected and
not edited after. Every claim carries a `s20/BRIEF.md` §9 label. **TARGET is the unit** throughout.

**Every number below states its basis** — *point cloud* / *built chain* / *repaired emission* —
which is §1 of the brief as this lane corrected it in its first hour.

The sealed 60-target benchmark was not read, probed, derived from, or inferred about, and
`results/benchmark_manifest.json` was not opened.

Artefacts, each with a `COMPLETE` flag:

    s20/results/D_QB_CLOSE/        Sprint 19 section 5 closed from Agent B's own artefact  n=126
    s20/results/D_C_ALIGNNULL/     the null the shared-bias statistic never had            n=126
    s20/results/D_E_AMBEROP/       H_AMBER as an operator: cap, ordering, definedness      3 targets x 64 bitstrings x 3 caps

Code: `s20/d_qb_close.py`, `s20/d_alignnull.py`, `s20/d_amberop.py`.

---

# 0. EXECUTIVE SUMMARY

> ## Two mission numbers were mislabelled, one of them inverting a standing role. Sprint 19's undelivered §5 is a clean falsification on all three of its own pre-registered endpoints. And the sprint's Q1 cannot be run as briefed, because the two Hamiltonians are not two energies on one landscape.

| # | statement | label |
|---|---|---|
| D1 | **"Best built 3.048 Å" is not a structure.** The coordinate average has a mean virtual Cα–Cα bond of **2.961 Å against the physical 3.804 Å** — 22.2% contraction, global minimum bond **0.649 Å**. Projecting it onto the pipeline's own manifold gives **3.204 Å**: the incumbent already IS the best structure the pipeline builds. Not a metric exploit — the best constant global rescale of `fit_ca` buys 0.034 Å, below MDE. | **ESTABLISHED** — coordinator error #7, corrected in five places |
| D2 | **The "irreducible +0.164 Å AMBER tax" is the PROJECTION's tax**, measured against the point cloud. AMBER's own cost on the deployed path is **+0.0207 [+0.0143, +0.0276]**. And the sign inverts: as a repair operator on the point cloud, AMBER at k=30 **beats** the projection by −0.022. | **ESTABLISHED** — a standing role was inverted in the live brief |
| D3 | **Sprint 19 §5 is a clean falsification.** `PREREG_B` P1 **REFUTED** (+0.079 vs best classical), P2 **REFUTED** (every trained quantum point dominated, CI excluding zero), P3 **REFUTED** (+0.067 generation). **Four of nine kill rules fire.** | **REFUTED** (n=126, pre-registered) |
| D4 | **The VQE genuinely trains, for the first time in this programme** — on the continuous basin-latent encoding, against the mandatory untrained control: −0.210 [−0.339, −0.082] realised, 5/5 folds; −0.299 on generation. This **reverses a project memory** that was measured on the k=4 lattice. | **ESTABLISHED** (ORACLE scoring, native-free decisions) |
| D5 | **Entanglement contributes nothing here.** Deleting the CNOTs changes realised RMSD by **−0.013 [−0.095, +0.077]**. Latent nearest-neighbour mutual information 0.045 bits. | **NOT MEASURED** — the honest label, not REFUTED |
| D6 | **No sampler, quantum or classical, at 8192 objective evaluations beats the zero-evaluation retrieval pool.** And it is not min-of-N: pool500 selects 75 from 500 while every sampler selects 75 from 8192. `c_lbfgs`, which optimises the deployed objective hardest, is second-worst on the board. | **ESTABLISHED** — `search saturates, discrimination binds` on a fifth instrument |
| D7 | **Every circuit-side landscape metric collapses when target difficulty is partialled out** — entropy +0.144→+0.065, max_prob −0.241→−0.099, ESS +0.204→+0.082. Only M survives, and M is nearly the outcome. | **ESTABLISHED** — the brief's §4 warning, measured |
| D8 | **The CVaR tail parameter is worth nothing and every setting of it is on the wrong side.** α=0.05/0.25/anneal are +0.012/+0.013/+0.025 against α=1.00 (the plain mean, no truncation), all CIs spanning zero. The machinery is genuine and correct; its measured contribution is zero. | **NOT MEASURED**, three ways, all adverse |
| D9 | **Q1 is not executable as briefed.** `H_AMBER = E ∘ Relax₅₀`: a 50-iteration restrained OpenMM minimisation lives *inside* `AmberHamiltonian.energy`. Swapping Hamiltonians swaps in a relaxation operator, plus a `+inf` region `H_Legacy` does not have. | **ESTABLISHED** (from source) — see §5 for the measurement |
| D10 | **The relaxation is what makes AMBER's objective DEFINED.** Finite on **58%** of the register at `Relax₁` and **97%** at `Relax₅₀`; the 50-iteration cap binds on **192/192** calls; `Spearman(Relax₅₀, Relax₁)` = 0.358–0.882. **F-D2 fires on both clauses.** | **ESTABLISHED** (operator property, 192 calls) / **INCONCLUSIVE** (the Spearman magnitude at n=3) |
| D11 | **The shared-bias alignment statistic had an unsubtracted null of 0.505** [0.469, 0.540], agreed by three independent constructions. **"Two thirds of the harmful component from zero information" becomes about a fifth** (0.64→0.18, 0.67→0.24). S1's ordering, S3, S7 and W5 are untouched; **S5 is strengthened three- to four-fold**. | **ESTABLISHED** — F-D1 fires at its weak clause |
| D12 | **A new methodological class — the SHARED REFERENT FLOOR.** Two quantities measured as deviations from a common reference have a correlation floor set by that reference, and it must be measured before the correlation is interpreted. Beside Z1 and Z2; costs one arm; unrun for four sprints. | **ESTABLISHED** — methodological |

---

# 1. PRIORITY 1 — THE STARTING AUDIT

## 1.1 The two mission numbers, reproduced and one of them relabelled

Reproduced end to end from `bench_results/cache/1fc9f2dcf489e2fb`, the production run pinned as
`s12/instrument.PROD_KEY`, 126/126 records present:

| record field | mean | median | what it is | **basis** |
|---|---|---|---|---|
| `rmsd_avg` | **3.0483** | 2.8373 | coordinate average of the shipped top-75 | **point cloud** |
| `rmsd_fit` | **3.2041** | 2.9661 | λ=0 arm of the projection — *the incumbent* | **built chain** |
| `rmsd_arm` | 3.2148 | 2.9661 | λ=0.3 ramah arm | built chain |
| `rmsd_full` | 3.2355 | 2.9757 | after restrained AMBER — *what the pipeline emits* | **repaired emission** |

Both quoted mission numbers reproduce exactly. `3.2041` is `s12/instrument`'s own asserted
`synthesis 3.2041` and `core/pipeline.Config`'s documented multi-start `lam_path[0.0]`; the
codebase itself flags that a *second* number, 3.2005/3.2007, is the single-start projection from
the medoid's torsions and is a different construction. **No discrepancy.**

**The discrepancy is in the label.** Virtual Cα–Cα bond lengths, n=126:

    avg_ca      2.9614 A      min bond over all targets  0.6492 A      22.2% contraction
    fit_ca      3.8040 A
    ca (arm)    3.8040 A
    amber_ca    3.8666 A

`avg_ca` is a **shrinkage estimator of Cα positions, not a backbone**. No side chain can be built
on it and AMBER cannot score it. Restoring physical bond length by uniform rescale gives
**4.091 Å**; projecting it onto the pipeline's own ideal-geometry manifold gives **3.204 Å**.

> **The incumbent already IS the best structure the pipeline builds. 3.048 Å is the best point
> cloud it emits, and the 0.156 Å between them is the price of stereochemical realisability —
> 1.9× the MDE.**

**And it is not a metric exploit other arms could copy.** The best *constant* global rescale of
`fit_ca` buys 0.034 Å (3.2041 → 3.1700 at s = 0.94), below the MDE. The averaging operator's
contraction is *per-target adaptive shrinkage*, which is real. Per-target ORACLE-optimal isotropic
rescaling buys 0.375 Å on `fit_ca` and 0.340 Å on `avg_ca`, but that is an oracle and not
available.

**Where the label came from, and where it was used as a bar.** Sprint 12 had it right —
`s12/dir_FINDINGS.md:241` reads *"the shipped Bayes incumbent (`bayes` = 3.048 avg / 3.204
projected)"*. It drifted by Sprint 19 and was then used as the baseline for arms measured as built
chains:

- `s19/FINAL_REPORT.md` §4 and `CLAIMS.md` O6 — *"a perfect ORACLE mode-picker lands at 3.472 Å
  against the 3.048 Å the pipeline already builds, +0.424 Å worse than doing nothing."* I verified
  in `s19/results/D_P5_headroom/headroom.json` that every arm there is a **built chain** via
  `I.build_ca`, while `avg` = 3.048 is the point cloud and `proj` = 3.213 is the built start.
  Like-for-like the number is **3.472 vs 3.213 = +0.259**. Direction unchanged; **magnitude
  overstated by 64%**. That is this lane's own predecessor's arm.
- `s19/a_sepfix.py:164` prints *"coordinate average (the incumbent)"* at 3.048 while its arms are
  `I.build_ca` chains — `CLAIMS.md` W4 carries the same mismatch.
- `s18/exp_report.py:352` prints *"pool → coordinate average (INCUMBENT)"* and computes all seven
  `[H5]` rows against it, including `pool_best`/`pool_mean`/`sel_full`/`sel_le1`, which are **real
  fragments with valid geometry**.

**The sweep of Sprints 16–19 for further instances is reported as NEGATIVE where it is negative.**
I recomputed the whole `[H5]` table on the built baseline from `s18/results/exp_main_pool_512.json`
(n=126, `proj` = 3.2126):

| route | RMSD (built) | vs `avg` (published) | **vs `proj` (correct)** | W/L |
|---|---|---|---|---|
| pool ORACLE best | 2.306 | −0.742 | **−0.906 [−1.041, −0.780]** | 124/2 |
| pool mean | 3.551 | +0.502 | **+0.338 [+0.267, +0.413]** | 14/112 |
| select by `E_full` | 3.531 | +0.483 | **+0.318 [+0.160, +0.480]** | 43/83 |
| select by `E_le1` | 3.600 | +0.551 | **+0.387 [+0.234, +0.549]** | 40/86 |
| avg → refine on `E_le1` | 4.335 | +1.287 | **+1.123 [+0.891, +1.366]** | 22/104 |
| avg → refine on `E_full` | 3.610 | +0.561 | **+0.397 [+0.240, +0.563]** | 41/85 |
| avg → `E_le1` exact argmin | 4.346 | +1.297 | **+1.133 [+0.885, +1.365]** | 26/100 |

**No Sprint-18 conclusion changes.** The offset is the constant +0.164 and every arm loses to both
baselines by far more than it. The drift is confined to the label — with one exception, which is
§1.2.

`s15/expand.py:147` already sweeps an isotropic scale on the average. The exploit was noticed in
Sprint 15 and never connected to the reporting basis.

## 1.2 The "+0.164 Å AMBER tax" is the projection's, and the sign inverts

`s20/BRIEF.md` §3 as briefed: *"AMBER is stereochemical repair with an irreducible +0.164 Å tax."*
Measured on the production cache, n=126, target-level paired bootstrap:

| operator | cost | CI | W/L | folds |
|---|---|---|---|---|
| projection λ=0.3, **vs the point cloud** | **+0.1664** | [+0.1307, +0.2019] | 16/110 | 5/5 |
| projection λ=0, vs the point cloud | +0.1557 | [+0.1222, +0.1889] | 18/108 | 5/5 |
| AMBER emission, vs the point cloud | +0.1871 | [+0.1490, +0.2253] | 20/106 | 5/5 |
| **AMBER on the BUILT chain — the deployed operator** | **+0.0207** | **[+0.0143, +0.0276]** | 40/86 | 5/5 |

**+0.1664 is the projection's tax.** Sprint 18 said so: `s18/FINAL_REPORT.md` §6.2 is titled *"The
repair tax is irreducible **and the AMBER pass is never spent**"*, and `s18/PHYS_STATUS.md:232`
calls it *"the incumbent **PROJECTION's** +0.163 [+0.132, +0.205] tax (L30's +0.164)"*.
`s19/BRIEF.md:127` still keeps them apart. `s20/BRIEF.md:83` collapsed them.

**The deployed AMBER is also not the arm any of those figures describes.** `core/pipeline.Config`:
`amber_k = 10.0`, `amber_steps = 0` (minimise to convergence). Every "AMBER k=30" figure in
Sprints 16 and 19 is a different restraint constant applied to a different input.

**And the sign inverts.** `s16/CLAIMS.md` J1 measured both repair operators against the same point
cloud: projection **3.2052 (+0.1554)**, AMBER k=30 **3.1831 (+0.1333)**. Against the only operator
it competes with, **AMBER is the cheaper of the two repairs by −0.022** — which is what project
memory records, and the opposite of what the brief's standing role implies.

> Corrected standing role: **turning the point cloud into a valid backbone costs +0.166
> [+0.131, +0.202]; that is the PROJECTION's tax and it is irreducible (s18 §6.2). AMBER's own
> cost, on the deployed path, is +0.021 [+0.014, +0.028] — a quarter of the MDE.**

This matters this sprint rather than as bookkeeping: Priority 3 is Legacy-vs-AMBER and Q1 warns
lanes not to assume "AMBER is more physical". A lane starting from "AMBER costs 0.164 Å" has been
handed a prior that AMBER is expensive, when the deployed operator costs 0.021 and is the better
of the two repairs.

## 1.3 The pipeline, pinned from source

| stage | implementation | pinned parameters |
|---|---|---|
| **retrieval** | BLOSUM62 sum over every length-n window of the leakage-safe library; `order` = stable argsort of −sim | `k = 500`, `tie_break = "stable"` (PINNED — changes the pool) |
| **filter** | shipped Bayes-risk score `I.shipped_score`, top-m | `m = 75`, `min_sep = 2` |
| **synthesis** | superpose all 75 on the pairwise-RMSD **medoid**, then arithmetic mean | → `avg_ca`, **point cloud** |
| **projection (3b)** | `core.project.lam_path(C, ramah, (0.0, λ))`, multi-start, exact gradient | `λ = 0.3`, `maxiter = 300`, `project_grad = "exact"`, `multi_start = True` |
| **repair (4)** | `core.amber.refine_coords` on the full backbone built from the λ=0.3 torsions | `amber_k = 10.0`, `amber_steps = 0` |

`reference_precision = True` round-trips the pool through float32 to reproduce `s9/final.py`
bit-for-bit. `project_grad = "exact"` vs `"analytic"` is documented as a *different pipeline* —
126/126 targets move, median 0.042 Å, worst 1.90 Å — and is hashed as one. Both facts are
correctly documented in the source and both are the kind of thing that would silently invalidate a
cross-sprint comparison; neither is violated in anything I checked.

## 1.4 The metric

`I.ca_rmsd` → `I.kabsch_rmsd_batch`: centre both, SVD of the cross-covariance, sign-correct the
smallest singular value by `det`, no scaling, all residues, no trimming. **Proper rotations only,
verified**: `rmsd(A, mirror(A)) = 1.28 ≠ 0` on a random 12-mer. It agrees with `s7.audit`'s
implementation (which is what `core/pipeline` calls) at **max |Δ| = 0** over 200 random pairs.

**One hazard, not a live defect.** `core.geometry.rmsd` is in the geometry backend's contract and
does **not** superpose — it is documented as "RMSD of two already-superposed point sets" and
differs from the frozen metric by up to 5.5 Å on unsuperposed input. Every live call site
(`core/data.py:424,759`, `peptide_db.py:85`, `s12/coord_benchN.py:151`, `s8/consensus2.py:189`,
`s8/inband.py:433`) superposes first. `fragment_db.py:98` and `verify/projection_divergence.py:59`
do not, and both are same-frame diagnostics. **Clean, but one careless import away from a silent
metric change.**

## 1.5 Legacy

`energy_terms.DEFAULT_WEIGHTS` is 11 terms, and `core.energy.FITTED_WEIGHTS` is **the identical
dict** — steric 4.0, contact 1.0, hbond_local 1.0, hbond_longrange 3.0, coop_helix 2.0,
coop_sheet 2.0, solvation 0.5, electrostatic 1.0, aromatic 0.8, torsion 0.15, compactness 0.4.
**Legacy is genuine and never fitted.** The *name* `FITTED_WEIGHTS` is a hazard and should be
aliased or renamed; nothing in the numbers is wrong.

## 1.6 AMBER, and the thing that blocks Q1

`amber_hamiltonian.AmberHamiltonian` is **not an energy function**. `energy(bitstring)`:

1. builds the full heavy-atom structure from the bitstring (hydrogens placed by **frozen local
   frames**, calibrated once — not relaxed per configuration);
2. runs `openmm.LocalEnergyMinimizer.minimize(ctx, tolerance = 2.0, maxIterations = 50)` under a
   **k = 100 kcal/mol/Å² positional restraint on every heavy atom**;
3. sets `k_rest = 0.0` and reports the **unrestrained** potential at the **minimised** coordinates;
4. returns `+inf` if `_is_collapsed` fires (`collapse_mode = "geometry"`, min heavy–heavy contact
   < 1.05 Å), counted in `n_collapsed`.

`core.quantum.FoldingHamiltonian.energy(bitstring)` builds coordinates and evaluates. No
relaxation, no infinite region.

    H_AMBER  =  E_ff14SB/GBn2  o  Relax_50        H_Legacy  =  E_knowledge

**§4 Q1's instruction — "change only `H_Legacy` ↔ `H_AMBER`" — is therefore not executable.** Every
landscape metric Q1 names is confounded by a relaxation operator that a *smoothing* argument and a
*basin-boundary-discontinuity* argument both apply to, and neither is "AMBER physics". The measured
half of this is §5 and falsifier **F-D2**.

Two further notes for Q1. First, both Hamiltonians take a **bitstring**; `AmberHamiltonian`
explicitly refuses a lattice representation and requires `TorsionStateRepresentation`, which
production instantiates at `torsion_window = 8` — 8 states and **3 bits per residue, not the k=4
lattice**. Q1 as briefed is a discrete-register experiment and collides with §5's continuous-torsion
rule; there is no continuous AMBER objective in the codebase, and `energy_from_coords` does not
escape the relaxation because it goes through the same `_evaluate`. Second, the deployed *pipeline*
AMBER (`core.amber.refine_coords`, k=10, steps=0) is a **third** object, different again from both
`AmberHamiltonian` and the k=30 arms.

## 1.7 CVaR and SPSA

**CVaR is genuine and the implementation is better than the reports imply.** `tail_indices`
replaces `argsort` with `argpartition` and proves the mask is *identical*, ties included, by
constructing the stable-sort tie rule explicitly. `cvar_from_samples` sorts before the mean
specifically so the floating-point sum is bit-identical to the shipped `np.sort(e)[:k].mean()`.
`cvar_gradient`'s default `"const"` baseline subtracts the mean weight over **all** samples — a
genuine control variate, measured at cos **+1.000000** against `cvar_gradient_exact`. The recorded
defect (`"tail"`: centre within the tail, accumulate over the tail only, biased at cos +0.655634)
is **reproduced verbatim as a named arm** rather than quietly deleted, so it stays measurable.
That is the correct disposition of a known defect and it is the project memory
`cvar-gradient-baseline-defect` handled properly.

**SPSA.** `_spsa` is Spall's, with three documented departures: a **progress-driven** gain
schedule (because the runner sets `maxiter` far above what the budget affords, which had made the
first step five times too small), **common random numbers** across the ± pair, and two evaluations
per iteration rather than three. All three are improvements and all three are documented.

Two observations the reports do not make.

- **`best_f` is an argmin over noisy estimates.** `_spsa` returns the parameter vector whose
  *midpoint estimate* `(fp+fm)/2` was lowest. That is a winner's-curse selection over a noisy
  sequence and biases the reported objective downward. It does not affect the returned *structure*
  (which is re-scored), so it is a reporting hazard, not a result hazard. **Diagnostic only.**
- **Common random numbers interact with CVaR's tail discontinuity in a way nobody has checked.**
  arXiv:2608.09810's SPSA-under-barren-plateaus analysis assumes the ± evaluations use
  *independent* measurement outcomes; this implementation shares the stream. CRN is the right
  choice for variance reduction, but with a discrete sample space and a tail selector, a small
  `ck` can leave both the sampled multiset *and* the tail membership unchanged, giving
  `fp − fm` **exactly zero** and a zero gradient estimate. `MIN_USEFUL_SPSA_ITERS = 100`'s own
  docstring records the shipped configuration returning "a bit-identical answer across 3 seeds and
  5 hyperparameter settings — the signature of an optimiser that never moved", attributed to
  budget. **The exactly-zero-gradient rate is a one-line instrument** (count `fp == fm`) and has
  never been measured. Recorded as **OPEN**, cheap, and for whoever owns the global VQE path.

## 1.8 Seeds and the target split

`s15/seed.stable_seed` is `blake2b(salt|repr(parts)…, digest_size=4)` — stable across processes,
platforms and Python versions. **No `hash()` anywhere in the live path.** Sprint 19 verified 630
arrays across a fresh process with 0 mismatches and distinct streams per target; I did not
re-verify that and cite it.

**The 126 tuning targets are cluster-disjoint from both the dev set and the benchmark by
construction**, not by audit: `s7/debias.tuning_targets` computes
`bad = {clusters of dev_set(24)} ∪ {clusters of benchmark()}` and skips any peptide whose identity
cluster is in it, taking one representative per remaining cluster with 9 ≤ n ≤ 16.
`tests/test_data.py` carries the disjointness assertions and is the authorised route to
re-verifying them; I did not run it, because it reads the manifest.

**Benchmark integrity certificate, established without opening the file.** A digest of the bytes
carries no information about which targets are in it, so this is not a probe:

    results/benchmark_manifest.json
      sha256  a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d
      size    8002 bytes        mtime  2026-09-01 18:52:50

Unmodified since before Sprint 12, and `git log` shows one commit touching it. **Future sprints can
check byte-identity against this line without ever reading a target name**, which is a cheaper and
stronger discipline than each sprint asserting "we did not open it".


---

# 2. BLOCK B — SPRINT 19 §5 IS DELIVERED, AND IT IS A CLEAN FALSIFICATION

`s19/FINAL_REPORT.md` records §5 as **NOT DELIVERED**: Agent B's main sweep complete on disk,
analysis and gauge check unreported at freeze, no number quoted. The artefact
`s19/results/qb_main.json` is **126/126 targets × 17 arms with zero missing cells**. I executed
`s19/PREREG_B.md` §§6–9 exactly as its author wrote them. **No structure was rebuilt; every number
is Agent B's own.**

**Soundness gate first, per that pre-registration:** the (M, D) identity `readout² = M² − D²`
asserts at **max |residual| = 5.3e-14 Å²** across all 126 × 17 cells. PASS.

**Basis is clean throughout** — `realised_ORACLE` is the projected **built chain** for every arm
including the incumbent, `avg_ORACLE` is the point cloud, both recorded per arm. Agent B's readout
got this right where three other lanes did not.

## 2.1 The table

| arm | realised (built) | sel_best | gen_best | gen@500 | M | D | cov<2.5 | evals |
|---|---|---|---|---|---|---|---|---|
| **pool500** — *the shipped generator* | **3.230** | 2.285 | 1.721 | 1.721 | 3.713 | 1.901 | 73.5 | **0** |
| c_marg — *zero-info, matched marginals* | 3.326 | 2.491 | 1.662 | 2.160 | 3.918 | 2.236 | 291.9 | 8192 |
| c_metroH — *best classical sampler* | 3.407 | 2.682 | 1.670 | 2.030 | 3.824 | 1.975 | 647.9 | 8192 |
| c_cem0.05 | 3.419 | 2.621 | **1.605** | 2.084 | 3.849 | 1.957 | 555.7 | 8192 |
| c_helix — *zero-info constant helix* | 3.426 | 2.623 | 1.924 | 2.372 | 4.068 | 2.419 | 164.9 | 8192 |
| c_chain0.25 — *classical learned proposal* | 3.427 | 2.622 | 1.636 | 2.139 | 3.912 | 2.102 | 244.1 | 8192 |
| **q_a1.00** — *best quantum* | **3.486** | 2.650 | 1.745 | 2.209 | 3.968 | 2.102 | 325.8 | 8192 |
| q_a0.05 / q_a0.25 / q_anneal | 3.497 / 3.499 / 3.511 | | 1.672 / 1.690 / 1.698 | | | | | 8192 |
| c_prod0.25 — *CNOTs deleted* | 3.498 | 2.684 | 1.739 | 2.113 | 3.863 | 1.876 | 486.7 | 8192 |
| c_lbfgs — *the incumbent optimiser* | 3.658 | 3.311 | 1.903 | 2.314 | 3.796 | 1.231 | 1132.8 | 8193 |
| q_untrained | 3.696 | 2.755 | 1.971 | 2.493 | 4.315 | 2.581 | 69.4 | 8192 |

## 2.2 All three pre-registered endpoints fail

**P1 (PRIMARY — realised).** `q_a1.00 − c_metroH` = **+0.079 [+0.001, +0.163]**, fold-clustered
[−0.039, +0.203], 60W/66L, 4/5 folds. The quantum arm is *worse*, and nowhere near the required
≤ −0.084 with a CI excluding zero. **REFUTED.**

**P3 (generation).** `q_a0.05 − c_cem0.05` = **+0.067 [+0.013, +0.121]**; at matched count 500,
+0.041 [−0.029, +0.111], **NOT MEASURED**. Against `pool500` at matched count: **+0.404
[+0.270, +0.546]**, 5/5 folds. The one axis Agent B pre-registered as having a chance does not.
**REFUTED.**

**P2 (ε-dominance).** Every *trained* quantum point is dominated by the classical point set with a
CI excluding zero: q_a1.00 −0.066 [−0.096, −0.037], q_a0.05 −0.097, q_a0.25 −0.111, q_anneal
−0.078. The leave-one-out control does exactly the job Sprint 17's recorded trap demanded: the only
non-dominated points are `pool500` (+0.204), `c_helix` (+0.157) and `q_untrained` (+0.150) — and
`q_untrained` is non-dominated **only because it is the most diverse** (D = 2.581), which C3 already
established is an accounting identity and not a value. **REFUTED.**

**Four of the nine `PREREG_B` §8 kill rules fire:**

    loses to c_marg (zero-information)   +0.160 [+0.084,+0.240]  5/5 folds
    loses to pool500 (the incumbent)     +0.256 [+0.131,+0.390]  5/5 folds
    matched by c_metroL (a thermostat)   +0.001 [-0.089,+0.094]
    matched by c_anneal                  -0.001 [-0.097,+0.100]

## 2.3 Two results worth keeping out of the wreckage

**(a) The VQE genuinely trains — first time in this programme, and it corrects a memory.**

    q_a1.00 - q_untrained   realised    -0.210 [-0.339,-0.082]  69W/57L  5/5 folds
    q_a0.05 - q_untrained   generation  -0.299 [-0.407,-0.195]  82W/44L  5/5 folds

Project memory `concentration-is-wrong-when-discrimination-binds` records *"running the VQE is
WORSE than not running it, 0/12 cells, +0.65 to +1.32 Å"* against best-of-N from the untrained
circuit. **That was measured on the k=4 lattice register.** On the continuous basin-latent encoding
at matched budget the sign **reverses**, decisively, on the mandatory control. The memory's
methodological half — *always control against best-of-N, never against an initialisation mean* — is
untouched and is exactly what made this measurable. **Its empirical half is encoding-scoped and
should be annotated, not deleted.**

**(b) Entanglement contributes nothing.** `q_a1.00 − c_prod0.25` — the identical pipeline with the
**CNOTs deleted**, same CVaR estimator, same seed — is **−0.013 [−0.095, +0.077]**. Latent
nearest-neighbour mutual information is 0.045 bits (mean) against the product latent's exact 0.000.
The correlated basin latent — the specific mechanism the lane was built to test, *"a whole region
flipping between two conformer families"* — is worth 0.013 Å with a CI spanning zero. The honest
label is **NOT MEASURED**, not REFUTED: the instrument has 0.084 Å of resolution and the effect, if
any, is inside it.

## 2.4 The finding that matters most for Priority 1

> **No sampler on the board — quantum or classical, learned or thermostatted — at 8192 objective
> evaluations, beats the zero-evaluation shipped retrieval pool on realised RMSD.**

And it is **not** a min-of-N artefact in the pool's favour: `pool500` selects its 75 from **500**
candidates while every sampler selects its 75 from **8192**. The sampler with the deepest search
(`c_lbfgs`, multi-start L-BFGS on the deployed objective) is the **second-worst arm on the board**
and has the lowest diversity (D = 1.231).

This is `search saturates, discrimination binds` reproduced on a **fifth** instrument, and it
prices Q2 before Q2 starts: **fifteen independent generators, every one of them worse than
retrieval after the terminal operator.** Any Q2 arm must clear `pool500` at 3.230 Å on a built
chain, and thirteen 8192-evaluation samplers did not.

## 2.5 The landscape metrics do not predict the outcome

The brief's §4 demands each landscape metric be connected to RMSD. Rank-based partial Spearman,
control = `pool500`'s realised RMSD (an ORACLE difficulty proxy), n=126:

| arm | metric | raw ρ | **partial ρ** |
|---|---|---|---|
| q_a0.05 | circuit entropy | +0.144 | **+0.065** |
| q_a0.05 | max probability | −0.241 | **−0.099** |
| q_a0.05 | effective sample size | +0.204 | **+0.082** |
| q_a0.05 | latent nn mutual information | +0.035 | **+0.167** |
| q_a0.05 | objective best | +0.335 | **+0.234** |
| q_a0.05 | D (diversity) | +0.499 | **+0.356** |
| q_a0.05 | M (member error) | +0.982 | **+0.930** |
| q_a1.00 | entropy / max_p / ESS | +0.222 / −0.268 / +0.248 | **+0.051 / −0.066 / +0.051** |
| c_prod0.25 | entropy / max_p / ESS | +0.156 / −0.136 / +0.146 | **+0.055 / −0.052 / +0.050** |

**Every distributional metric is a proxy for "this target is hard".** Only M survives, and M is
member error — nearly the outcome itself, an identity rather than a predictor. The pattern is
identical on the CNOT-free classical control, so it is not a property of the quantum circuit.
**The Q1 lane must partial out difficulty on every landscape metric it reports; a raw ρ on this
instrument is a difficulty measurement.**

## 2.6 The CVaR tail parameter is worth nothing, in the direction theory predicts

Paired against `q_a1.00`, which is **α = 1.00 — the plain sample mean, no tail truncation at all**:

    q_a0.05  - q_a1.00   +0.012 [-0.084,+0.105]  56W/70L  folds 3/5
    q_a0.25  - q_a1.00   +0.013 [-0.066,+0.089]  56W/70L  folds 2/5
    q_anneal - q_a1.00   +0.025 [-0.057,+0.108]  52W/74L  folds 3/5

All three **NOT MEASURED**, and all three point estimates on the wrong side. **Priority 2 is
preserved and its contribution is measured at zero**: the CVaR machinery is genuine and correct
(§1.7), and on this instrument the tail parameter buys nothing.

There is 2026 theory for exactly this regime, and it is new since Sprint 19's literature pass —
**arXiv:2605.02850**, *Quantum Tilted Loss in Variational Optimization*. It places CVaR in a
tilted-loss family and shows (i) it **cannot** eliminate barren plateaus, only reshape local
geometry; (ii) CVaR's hard α-quantile cutoff *"creates optimization discontinuities as parameters
shift samples in/out of the tail set"* where the smooth tilted loss does not; (iii) the sample
complexity of resolving the sharpened gradient scales as `Õ((e^{|γ|Δ} − 1)² / γ²ε²)` —
**exponential in the tilt**. Landscape sharpening is paid for in shots, exponentially. At 512 shots
an α = 0.05 tail is 26 samples. **That α = 1 wins is the predicted regime.**

> **Corollary the Q1 lane needs.** "CVaR concentration" is not a clean landscape property of the
> Hamiltonian. The tail-membership discontinuity belongs to the CVaR **rule**, identically under
> `H_Legacy` and `H_AMBER`. Any difference in measured CVaR concentration between them is a
> difference in energy **spectra** passed through a shared discontinuous operator. **Derive the
> operator before interpreting its statistic** (§9).

---

# 3. BRANCH ENUMERATION, SPRINTS 16–19

Branch = a research direction, not a claim. Status is mine, argued from the artefacts; where I
disagree with the sprint that produced it, the disagreement is stated.

| # | branch | status | basis |
|---|---|---|---|
| 1 | Objective functional redesign (weights, robust losses, Bayes risk) | **CLOSED** | s18 §3.1 ceiling 0.083 Å; s19 O9 `riskw` is the only arm ever to beat the deployed functional (−0.134) and lands at 3.476, worse than not refining |
| 2 | Degree-1 / Walsh truncation | **REFUTED** | s18, six independent closures; headline was a 2-bit gauge artefact at percentile 0.00 |
| 3 | Objective tempering (`sd^−p`) | **REFUTED** | s18 G6c |
| 4 | Legacy as an objective term (`leg_contact`) | **REFUTED** | s18 §5, +0.722, harmful *with a sign* |
| 5 | `leg_torsion` as a gate | **REFUTED** | s18 §5.5 — converts negatively |
| 6 | Legacy as selector / ranker | **REFUTED** | s16/s17 |
| 7 | Legacy as a steric rejection gate | **REFUTED** | s19 C7 — worse than random, than diversity-preserving rejection, and than *not rejecting*; monotone from r=2 |
| 8 | Physics filters/gates in general | **REFUTED** | s18 §6.1 — every physics score loses to a random gate of the same size |
| 9 | AMBER as a ranker | **REFUTED** | memory `physics-ranks-real-geometry-not-lattice`: on matched pools with an rg control AMBER puts the native at the 51st percentile |
| 10 | AMBER as stereochemical repair | **SUPPORTED — standing role**, but **relabelled this sprint**: its cost is +0.021 on the built chain, not +0.164 (§1.2) | production cache, n=126 |
| 11 | AMBER Pareto (k-ladder, blends, Ca-only) | **OPEN / PARTIAL** — s19 §6 is PROVISIONAL from a superseded capped run; the clean rerun is at 41/126 | see §3.1 |
| 12 | AMBER restraint constant `k` as a tunable | **INCONCLUSIVE** | memory: chosen on dev-set RMSD with **no native-free rule that selects it** |
| 13 | Quantum VQE as an optimiser (k=4 lattice) | **CLOSED** | s18 C1–C5; greedy certifies the optimum at 36 evaluations in 100% of cells |
| 14 | **Quantum CVaR-VQE as a structured sampler (continuous torsion)** | **REFUTED** — all three pre-registered endpoints, 4/9 kill rules | **this sprint, §2** |
| 15 | Quantum in-band selection via `log q_θ` | **NOT MEASURED** — s17 carries it as OPEN at +0.083 [+0.001, +0.179]. That is *at* the 0.084 MDE with a CI lower bound of +0.001. By §8 it is not a validated effect and **must not be carried forward as support** | s17 F10 — **my downgrade** |
| 16 | Predictor architecture / training loss | **CLOSED** | s19 P1 — 8 leave-fold-out predictors, none past MDE; ORACLE MAE −11.7% for nothing |
| 17 | Joint consistency / triangle / metric realisability | **REFUTED** | s19 P2 — more realisable ⇒ more coherently wrong |
| 18 | Distribution-shape / multimodality consumption | **REFUTED** | s19 O5/O7 — mode information real, zero transfers, split-half inside a fold |
| 19 | Training loss "in the right basis" (short-favouring) | **OPEN — BLOCKED, not declined** | s19 P5 — ESM cache/memory; pre-registration standing unedited |
| 20 | Post-hoc distogram correction / calibration | **CLOSED** | memory `error-shape-not-mae-decides-ranking` |
| 21 | Native-free estimation of the separation profile | **DOWNGRADED** | s19 W4 — the ladder is interpolation toward `poolfull`, and it fails its own pre-registered trap |
| 22 | Error-subset detection and repair | **REFUTED** | s19 W6/W7 — detecting is partly possible, repairing what you detect is worth nothing |
| 23 | Selection: score gates / in-band ranking | **REFUTED** (as a lever) / **ESTABLISHED** (as a mechanism) | s19 C2/C4 |
| 24 | Diversity-preserving gate design | **ESTABLISHED and worth zero** | s19 C5 — recovers the whole Sprint-18 deficit, then ties matched-random |
| 25 | Consensus / medoid selection | **ESTABLISHED** — the only in-band discriminator, −0.172 | memory |
| 26 | Retrieval: widening K | **REFUTED** | s17 — worse realised answer; composition confound never controlled |
| 27 | Shortlist construction for coverage (k-means + score) | **CLOSED BY THEORY, never run** — see §3.2 | s17 §11 named it "the cheapest outstanding experiment"; three sprints later it has not been run at n=126 |
| 28 | Softening the top-*m* cut | **CLOSED** | best T −0.028, a third of the MDE |
| 29 | Averaging space (coordinate vs torsion) | **ESTABLISHED** | memory — coordinate beats torsion by 1.024 Å |
| 30 | Jacobian steering | **CLOSED** | s16/s18 |
| 31 | Learned representations / ESM features | **CLOSED as features** (s17 §7); ESM worth +0.288 Å on *selection* | memory |
| 32 | **Local window refinement of a candidate** | **DOWNGRADED — not comparable to the mission number**; see §3.3 | s17 F12 |
| 33 | Chemical-shift / torsion restraints | **CLOSED** | memory `torsion-restraints-reach-the-target` |
| 34 | Why post-fit κ predicts per-target outcome | **OPEN**, no explanation, and I decline to offer one (§6) | s19 M7 |
| 35 | "Do not consume the argmin at all" | **OPEN — nobody is working on it** | s19 §9 item 2 |
| 36 | **The repair operator itself** — projection vs restrained relaxation vs de-shrunk projection | **OPEN, and newly so** | §1.1/§1.2 and §7 |

## 3.1 Unresolved and partial artefacts, by inspection rather than by report

A mechanical scan of `s1{6,7,8,9}/results/**.json` for `complete: false` under a clean name:

| file | rows | expected | hazard |
|---|---|---|---|
| `s19/results/agentC_pareto.json` | **41 → 48 while this was written** | 126 | the clean AMBER-Pareto rerun, **still writing**; the headline name holds a growing partial |
| `s19/results/agentC_pareto_report.txt` | — | — | **was generated from the SUPERSEDED capped-2000 run** and sat under a clean name; its header read "n = 126" (the pre-registration's n) while its rows read "n=62" (post-gate). **FIXED during this sprint** — now `_SUPERSEDED_agentC_pareto_report.txt`. |
| `s18/results/exp_budget.json` | 85 | — | `complete: false`, clean name, no `n_expected` |
| `s18/results/exp_leverage.json` | 88 | — | as above |
| `s18/results/exp_main_uniform_512_sub30.json` | 17 | — | as above |

I have not found any of the three Sprint-18 files quoted at the wrong n. They are readable traps,
which is the hazard `s12/instrument.write` was written to document, and the fix is one rename each.

Claims still carrying an unresolved label: s17 F10 (§15 above), s18 B10 (a smoke *stopped* at n=29
with both point estimates going the wrong way), s18 D12 (withdrawn, not measured), s19 M7, s19 P5,
s19 Block F (F-C2 UNTESTED at freeze).

## 3.2 The "cheapest outstanding experiment" is closed by the operator, not by neglect

`s17/FINAL_REPORT.md` §11 names the shortlist as **the** gap to attack (0.791 Å of ceiling destroyed
by the top-75 cut) and records that a 6-target smoke of a coverage construction reached
−0.275 [−0.473, −0.091], *"the full comparison was stopped for compute and is the cheapest
outstanding experiment."* Three sprints later it has not been run.

**It should be marked CLOSED, and the reason is the operator, not the effort.** Project memory
`operator-consumes-set-mean` fits `d_out = 1.16·d_set_mean + 0.04·d_set_best` at R² 0.89. A
shortlist construction that raises the **ceiling** (the set best) without also improving the set
**mean** converts at the 0.04 coefficient: 0.791 × 0.04 = **0.032 Å**, well under the MDE. Sprint 17
in fact wrote the correct instruction — *"build the shortlist for the readout you intend to use"* —
and then filed it as an experiment about the ceiling. The deployed readout consumes the mean.

That leaves a genuinely open sub-branch, which is *not* the one that was filed: a shortlist
construction that improves the set **mean**. `legacy_clust` (s19 C5) is the nearest thing tried and
it ties matched-random.

## 3.3 s17's "one positive lever" is not comparable to the mission number

`s17/FINAL_REPORT.md` §11 item 4: *"the one positive lever measured anywhere in the sprint is local
refinement — 3.636 → 2.95 Å over 6-residue windows."* Reading `s17/quantum_FINDINGS.md` §5, the
comparison is:

- **n = 9 targets** × 3 seeds × 4 windows, not 126;
- the readout is the **argmin**, not the m=75 coordinate average the pipeline uses;
- the baseline 3.636 Å is *"the retrieval torsion prior's mode"*, **not** the 3.204 Å incumbent.

So it is a 0.7 Å improvement over a baseline 0.43 Å worse than the incumbent, through a readout the
pipeline does not use. Memory prices the difference between those readouts directly: a perfect
rank-1 is worth **−1.74 Å through argmin and −0.03 Å through the m=75 average.** The lever is real
on its own instrument and there is no evidence it survives to the deployed one. **This is the same
class of error as §1.1 — a number quoted against a baseline and a readout other than the mission's
— and it is the reason nobody pursued the "one positive lever" for three sprints without anyone
writing down why.**

## 3.4 A cheap unmeasured instrument on the SPSA path

`core.quantum._spsa` uses **common random numbers** across the ± pair while the CVaR objective
selects a **tail**. With a discrete sample space, a small `c_k` can leave both the sampled multiset
and the tail membership unchanged, making `fp − fm` **exactly zero**. `MIN_USEFUL_SPSA_ITERS`'s own
docstring records the shipped configuration returning *"a bit-identical answer across 3 seeds and 5
hyperparameter settings — the signature of an optimiser that never moved"*, attributed to budget.
**Counting `fp == fm` is one line and has never been run.** OPEN.

---

# 4. BLOCK C — THE SHARED-BIAS ALIGNMENT STATISTIC HAD AN UNSUBTRACTED NULL OF 0.505

Pre-registered in `s20/PREREG_D.md` Block C, falsifier **F-D1** fixed before the run. Artefact
`s20/results/D_C_ALIGNNULL/` (COMPLETE, n=126); code `s20/d_alignnull.py`.

## 4.1 The operator, derived before its statistic is interpreted

`s19/a_source.py` reports cross-family alignment of the coherent component as

    corrC(A, B)  =  corr( d(X_A) - d_nat ,  d(X_B) - d_nat )

where `X_F` is family F's own terminal fit. Both arguments are deviations from **the same
native**. Write, for any reference field `d_typ` that both families share,

    d(X) - d_nat  =  (d(X) - d_typ)  -  (d_nat - d_typ)

The second term is **identical for every family by construction**. So `corrC` carries a floor set
by the ratio of the native's own deviation from typicality to the families' deviations. **Any two
realisable length-n peptide structures scored against the same native must correlate.** Sprint 19
never measured how much.

## 4.2 The null, three ways, and they agree

Eight independent draws per target, averaged; N1 and N2 use **the same draws**, so their difference
isolates the fit operator and nothing else.

| null | construction | value | 95% CI |
|---|---|---|---|
| **N1** pool × pool | two independent windows from this target's own K=500 pool, no fit | 0.507 | [0.472, 0.542] |
| **N2** fit × fit — **MATCHED** | the same `align_lib.fit`, the same cached start, the same `1/sd²` weights, driven to those same two windows' distance fields | **0.505** | [0.469, 0.540] |
| **N3** cross-target | windows of the same length drawn from *different* targets | 0.529 | [0.489, 0.568] |

**The three agreeing to within 0.024 is the internal check.** The null is a property of "two
realisable peptides scored against one native" — not of the fit operator, and not of this target's
own pool. §6.2's rule applies in the lane's favour here: the measurement settled it and no
derivation was needed.

## 4.3 F-D1 fires at the weak clause, not the strong one

I pre-registered that if N2 ≥ 0.575 with a CI containing it, the zero-information references are
*at* the null. **They are not.** N2 = 0.505 and its CI excludes 0.575. Both zero-information
references carry genuine excess, 5/5 folds, CIs excluding zero. **S2 is not refuted.**

## 4.4 What breaks is the "two thirds" framing

Excess over the matched null, and the share of the ceiling recomputed as
`(corrC − null) / (ceiling − null)`:

| pair | raw | **excess [95% CI]** | W/L | folds | share **as published** | share **null-subtracted** |
|---|---|---|---|---|---|---|
| same arch, different seed — **the CEILING** | 0.898 | +0.394 [+0.357, +0.432] | 123/3 | 5/5 | 1.00 | 1.00 |
| same arch, different regularisation | 0.833 | +0.329 [+0.290, +0.366] | 122/4 | 5/5 | 0.93 | **0.83** |
| **retrieval pool mean** | 0.784 | +0.279 [+0.244, +0.316] | 117/9 | 5/5 | 0.87 | **0.71** |
| different architecture, matched reg (`d20~pairnet`) | 0.754 | +0.250 [+0.207, +0.292] | 105/21 | 5/5 | 0.84 | **0.63** |
| different architecture (`deployed~pairnet`) | 0.731 | +0.226 [+0.183, +0.270] | 101/25 | 5/5 | 0.81 | **0.57** |
| ZERO-INFO separation prior | 0.575 | **+0.070 [+0.028, +0.112]** | 83/43 | 5/5 | 0.64 | **0.18** |
| ZERO-INFO constant α-helix | 0.598 | **+0.093 [+0.037, +0.146]** | 84/42 | 5/5 | 0.67 | **0.24** |
| *(zero-info vs the joint architecture: `pairnet~sepprior`)* | 0.564 | +0.060 [+0.015, +0.104] | 78/48 | 5/5 | 0.63 | **0.15** |

*Both cross-architecture pairs are shown because `s19/FINAL_REPORT.md` quotes **0.754**, which is
`d20~pairnet` (matched regularisation), while `s19/CLAIMS.md` S1's "cross-architecture 0.81" is the
share for `deployed~pairnet` at 0.731. They are different pairs and I report both rather than
introduce a third mismatch of my own.*

> **The zero-information references reproduce about a FIFTH of the same-architecture ceiling, not
> two thirds.**

Three consequences:

- `s20/BRIEF.md` §2 "The source": *"0.64–0.67 for sequence-blind references"* should read
  **"0.64–0.67 raw, 0.18–0.24 above a 0.505 shared-target null"**.
- `s19/FINAL_REPORT.md` §0 and §2: *"Two zero-information references each reproduce ~two thirds of
  the same-architecture ceiling"* and the bolded *"**And the harmful component is largely
  sequence-independent**"* **do not survive**. Above the null the component is mostly
  sequence-*dependent*.
- **S5 is strengthened, not weakened.** Sprint 19 stated it defensively — *"sequence conditioning is
  not worthless; 0.833 and 0.784 sit well above 0.598"*. After null subtraction the gap between
  conditioned (0.83 / 0.71 / 0.57) and zero-information (0.24 / 0.18) is three to four times larger
  than it looked. **The lane that wrote S5 as a caveat was writing the main result.**

## 4.5 What is undamaged, said as loudly as the correction

**S1's ordering survives intact, and the ordering is what the mechanism needs**: ceiling >
same-architecture > retrieval pool > cross-architecture ≫ zero-information, before and after.
**S3** (the incoherent component does not share) is untouched — it is a separate statistic.
**S7 and W5 — the pool shares the harmful component and is simultaneously the estimator of it — are
untouched**: the pool's null-subtracted 0.71 is the second-highest number in the table and still
above cross-architecture. **THE WALL IS WHERE AGENT A PUT IT.**

I recorded in `PREREG_D` before the run that I expected the ordering to survive and the "two
thirds" framing to break. That is what happened, and it is recorded here so the outcome cannot be
re-narrated either way.

## 4.6 Limitation, volunteered

The null is built from pool windows at ≈4.45 Å mean member RMSD while the predictors' fitted
structures sit at ≈3.5–3.7 Å. The null trends **mildly upward** with member distance (N1 0.507 on
this target's own pool; N3 0.529 on further-away windows), so a null matched to the predictors' own
distance would be slightly **lower** and every excess slightly **larger**. That moves all six rows
the same way and cannot change the ordering, but **0.18 / 0.24 is a lower bound on the
zero-information share** and I do not claim it to two decimals. A quality-matched null is the clean
follow-up and is cheap.

A second choice worth naming: `(corrC − null)/(ceiling − null)` is one normalisation; a Fisher-z
excess is another and would give different second digits. The **direction and ordering** are
invariant to that choice; the exact fractions are not.

## 4.7 A new methodological class — THE SHARED REFERENT FLOOR

This is §9 in a form the register does not yet have. Not an identity published as a discovery, but
**a statistic published without deriving its operator**: `corrC` is a correlation between two
vectors that share a large additive term by construction.

> **Whenever two quantities are both measured as deviations from a common reference, their
> correlation has a floor set by that reference, and the floor must be measured before the
> correlation is interpreted.**

It sits beside **Z1** (a control matched in the wrong space) and **Z2** (a measured null beats an
analytic one) as a third cheap, mechanical check. It costs one extra arm. Nobody in four sprints
ran it, and it moved a headline from "two thirds" to "a fifth".

---

# 5. BLOCK E — `H_AMBER` IS NOT AN ENERGY, AND THE RELAXATION IS WHAT MAKES IT DEFINED

Pre-registered in `s20/PREREG_D.md` Block E, falsifier **F-D2** fixed before the run. Artefact
`s20/results/D_E_AMBEROP/` (COMPLETE); code `s20/d_amberop.py`.

Three targets (1CS9 n=9, 1A1P n=13, 1A13 n=14), 64 random bitstrings each on the **production
k=8 torsion register** (`torsion_window = 8`, 3 bits/residue), **192 calls per arm**. Identical
`AmberHamiltonian` object, identical bitstrings; **only the iteration cap differs**. OpenMM's
`maxIterations = 0` means *unbounded*, so 1 is the smallest honest floor and the two arms then
differ in nothing else — the §6.1 rule applied to an operator rather than to an arm.

| pdb | n | bits | finite @50 | `E₅₀ − E₁` mean | ρ(50, 1) | ρ(100, 50) | **cap bound** | collapsed @50 |
|---|---|---|---|---|---|---|---|---|
| 1CS9 | 9 | 27 | 48/64 | −3.2e+05 | 0.358 | 0.835 | **1.00** | 1 |
| 1A1P | 13 | 39 | 29/64 | −3.5e+08 | 0.827 | 0.797 | **1.00** | 4 |
| 1A13 | 14 | 42 | 35/64 | −3.3e+07 | 0.882 | 0.903 | **1.00** | 1 |

## 5.1 The cap binds on 100% of calls

Doubling the cap to 100 iterations still moves the energy by more than 1e-3 kcal/mol on **every one
of the 192 calls**. By **Z6a**, `H_AMBER` is not "AMBER energy" — it is *"50 iterations of
restrained relaxation, then evaluate"*, a **different operator**, and it is a different operator on
every call rather than on an unlucky one. **Z6b's operational half applies in full**: a bound must be
verified on the calls it actually *bound*, and it bound all of them.

## 5.2 The ordering changes, and the ordering is what CVaR consumes

`Spearman(E ∘ Relax₅₀, E ∘ Relax₁)` = **0.358 / 0.827 / 0.882**, mean 0.689, against the
pre-registered 0.95 threshold. And it has not converged *at* 50 either:
`Spearman(Relax₁₀₀, Relax₅₀)` = 0.797–0.903. **The ordering is still moving when the cap stops it.**
A CVaR tail selected under `Relax₅₀` is a different subset from the one selected under `Relax₁` or
`Relax₁₀₀`.

## 5.3 The result I did not predict, and it is the decisive one

    finite at Relax_1     112/192   (58%)
    finite at Relax_50    186/192   (97%)
    finite at Relax_100   190/192   (99%)

> **The relaxation is what makes the objective DEFINED.** On **42% of the register** AMBER's energy
> at the as-built ideal-geometry coordinates is not a finite number, and the 50 minimisation steps
> are what turn it into one.

Mean shifts of 3e5 to 3.5e8 kcal/mol say the same thing: the unrelaxed structures are in hard steric
overlap. So `Relax₅₀` is **not a smoothing convenience bolted onto an energy — without it there is
no objective on most of the space.**

**This half-kills my own pre-registered fix**, and I say so rather than quietly narrowing it. "Just
add the `E_amber ∘ Relax₁` arm" yields a function that is `+inf` on 42% of the register and cannot
be CVaR-optimised as it stands. The pre-registration is left unedited and the defect recorded here.

The `+inf` collapse guard fired **1 / 4 / 1** times per target at `Relax₅₀`, reported per **Z6**
because a guard's firing count is the evidence that it is not vacuous.

## 5.4 The relaxation manufactures part of the Legacy–AMBER agreement

    Spearman(E_legacy, E_amber o Relax_50) = +0.229
    Spearman(E_legacy, E_amber o Relax_1)  = +0.140

**The two Hamiltonians barely agree on ordering at all** — +0.229 — and about a third of even that
agreement appears only *after* the relaxation. **Any Q1 statement of the form "Legacy and AMBER rank
configurations differently" is measuring `E_legacy` against `E_amber ∘ Relax₅₀`, not against
`E_amber`.**

## 5.5 What Q1 should do instead

The comparison the brief wants — one register, one ansatz, one optimiser, two energies — is
available only if both energies are defined on the same domain. Three options, in my order of
preference:

**(b) HOLD THE RELAXATION FIXED ON BOTH SIDES — recommended as primary.** Give `H_Legacy` the same
`Relax₅₀` preprocessing: evaluate the knowledge-based energy at the AMBER-relaxed coordinates. The
relaxation becomes a *shared* operator and cancels. This is the cleanest "change only H" available
and it costs one wrapper.

**(a) RESTRICT THE DOMAIN — recommended as the diagnostic beside it.** Run both Hamiltonians on the
subset of the register where `E_amber ∘ Relax₁` is finite (58%), report that restriction as part of
the operator, and compare `E_legacy`, `E_amber ∘ Relax₁` and `E_amber ∘ Relax₅₀` as **three** arms.
The `Relax₅₀ − Relax₁` contrast is then the relaxation's own effect, cleanly separated from the
physics, on a common domain.

**(c)** Failing both, report every landscape conclusion as **RELAXATION-CONFOUNDED** in its label.

## 5.6 Label and limits, stated plainly

**n = 3 targets, 192 calls per arm. This is not an RMSD claim and I make none.** The cap-binding
rate (**192/192**) and the definedness rates (**58% → 97%**) are properties of the operator counted
over calls and are not marginal statistics. The mean Spearman **0.689** is a three-target mean
spanning **0.358–0.882** and must be quoted as that range, not as a point.

> **Label: ESTABLISHED as a property of the operator; the Spearman magnitude INCONCLUSIVE at this n.**

**A self-correction, volunteered.** The smoke run of this module wrote a `COMPLETE` flag for a
1-target, 8-draw configuration — **precisely the partial-under-a-clean-name hazard this lane flagged
in `agentC_pareto.json` and three Sprint-18 files in §3.1**. I caught it on the next read, removed
the stale flag, and added a guard so `d_amberop.py` writes `COMPLETE` only for the full
pre-registered configuration. Recorded because a lane that flags a hazard and then commits it should
say so.

---

# 6. PRIORITY 3 — LITERATURE, and what is genuinely usable

Sprint 19's pass (`s19/agentD_FINDINGS.md` §4) was thorough and I do not repeat it. Everything
below is either **new since that pass** or **used differently**.

## 6.1 Quantum — four items that bear on this sprint's questions

**arXiv:2605.02850, *Quantum Tilted Loss in Variational Optimization* (2026) — the single most
usable quantum item.** It places CVaR inside a tilted-loss family (Theorem 7:
`CVaR_α ≥ L_γ + (1/γ)log(1/α)` for γ<0) and establishes three things this programme needs:

1. **CVaR cannot eliminate barren plateaus.** It reshapes local geometry when structure and tilt
   schedule align; it does not remove the flatness.
2. **CVaR's hard α-quantile cutoff "creates optimization discontinuities as parameters shift
   samples in/out of the tail set"**, where a smooth tilt does not. That discontinuity is a property
   of the **rule**, present identically under `H_Legacy` and `H_AMBER`.
3. **Sharpening is paid for in shots, exponentially**: sample complexity
   `Õ((e^{|γ|Δ} − 1)² / γ²ε²)`. At 512 shots an α=0.05 tail is 26 samples.

This is the theory for §2.6: **α = 1.00 winning is the predicted regime**, not an anomaly. It also
kills "CVaR concentration" as a clean landscape *property of a Hamiltonian* — the operator must be
derived before its statistic is read (§9).

**arXiv:2608.09810, *From Barren Plateaus to SPSA Optimization in VQE* (2026).** Iteration
complexity `T = Ω(κ²/ε²·(E‖∇f‖²)²)` and measurement complexity with exponent 5/2 under BP; explicit
gain schedules `μ_t = μ₀T^{−1/2}`, `c_t = c₀T^{−1/8}`, `M_t = M₀T^{1/4}`. **Its Lemma 1 assumes the
± evaluations use INDEPENDENT measurement outcomes.** `core.quantum._spsa` uses **common random
numbers**. CRN is the right engineering choice, but it puts the implementation outside this
analysis, and combined with CVaR's tail discontinuity it admits an exactly-zero gradient estimate.
See §3.4 — one line to measure, never measured.

**Geo-ADAPT-VQE and the ADAPT family (QFI-driven operator selection).** **"Landscape-aware ansatz
selection" already exists**, and is done with the quantum Fisher information rather than a Hessian
spectrum. Relevant to the novelty analysis below, and relevant to Q1: if a landscape-aware choice is
to be made, the published machinery is QFI-based.

**arXiv:2606.18580 and the expressibility literature.** *Expressibility alone does not determine
trainability*; highly expressive circuits are more prone to flat landscapes but poorly-scaled
gradients and unfavourable parameter geometry are a separate axis. **Directly relevant to Q1's
"ansatz limitation" candidate**: a Q1 arm that measures expressivity and infers trainability is
making an inference the literature does not support.

**arXiv:2609.00235, *Optimization Landscape Geometry in VQE for Frustrated Quantum Spin Models*
(2026)** — the closest published methodology for the landscape-geometry half of Q1.

## 6.2 Structure — the item that externally corroborates §1.1

**PMC2662860, *Improving consensus structure by eliminating averaging artifacts* (COMBO/MCORE).**
This is the external version of §1.1 and it should have been found four sprints ago:

- coordinate-averaged consensus structures have **unrealistic local geometry** — **63% of atoms in
  clashes** in the averaged model against **1% after refinement**;
- **the averaged model has BETTER RMSD than its refined version**: COMBO **3.28 Å** vs MCORE
  **3.36 Å** — a **0.08 Å** penalty for restoring geometry;
- the fix is **Monte Carlo refinement driving a structure *toward* the average under harmonic
  bond/angle potentials** — i.e. restrained relaxation, **not** a least-squares projection.

**Two things follow.** First, this programme's penalty is **+0.156 to +0.166 Å, roughly twice
theirs** — consistent with peptide length: fewer residues, fewer constraints, proportionally more
contraction. Second, **this programme independently reproduced MCORE's central design finding and
never connected it**: `s16/CLAIMS.md` J1 measures the projection at **+0.1554** and AMBER k=30 —
restrained relaxation toward the average — at **+0.1333** against the same point cloud. Restrained
relaxation beats least-squares projection, in this codebase, by the amount the published method
claims. The deployed pipeline uses the projection.

**Krogh–Vedelsby ambiguity decomposition** (`E_ens = Ē_member − Ā`). Sprint 19's C1 —
`readout² = M² − D²`, labelled **EXACT** — is that decomposition specialised to squared error. The
label is right and the programme is not overclaiming, but **the citation is owed**, and it means the
ensemble-diversity literature's conclusions transfer: ensemble gain goes as the *ambiguity*, and
reducing member correlation is the classical lever. That is Q2, in the classical form, with fifty
years of results behind it.

**PMC13370982, *Benchmarking AI Protein Structure Predictors Reveals a Persistent Bias in
Multi-State Proteins* (2026).** Predictors default to the dominant PDB state, and **MSA-based
approaches recapitulate the same bias** — so it is not architectural. That is the external analogue
of S4 (*"every predictor emits a typical peptide of that length"*) on a different problem, and it
supports S4 while being independent of it.

**Conformer generation.** RINGER and torsional diffusion remain the mature comparison; new since
Sprint 19 is **PPI-Diff (2026)**, internal-coordinate manifold diffusion on a Riemannian manifold of
dihedrals, explicitly balancing local stereochemical validity against flexible conformational
capture. That is the generator R1 below would eventually reach for — *after* R1's measurement says
whether the lever exists.

**The Moreau envelope** is the right theoretical frame for §1.6/§5: partial minimisation of an
energy is a smoothing operator, and minimising the envelope has the same minimisers as minimising
the original. **But `E ∘ Relax₅₀` is not a Moreau envelope**: it is an *inexact, iteration-bounded,
start-dependent* partial minimisation, which is neither smooth nor a proximal map, and can be
discontinuous where the minimiser's basin assignment flips. The literature gives the intuition
(smoothing) and denies the guarantee.

## 6.3 Novelty analysis (directive §40) — what is genuinely new, and what is not

The goal is not to force a novelty claim. Item by item, from the methods rather than the abstracts:

| combination | verdict |
|---|---|
| quantum conformational proposal generation | **NOT NOVEL.** arXiv:2609.02113 (Cleveland Clinic, Sept 2026), QuPepFold (PLOS One, Feb 2026), arXiv:2608.05491. |
| CVaR-driven conformer ensembles | **NOT NOVEL as a concept** — QuPepFold explicitly uses CVaR-tuned VQE to generate IDR ensembles. |
| landscape-aware ansatz selection | **NOT NOVEL** — Geo-ADAPT-VQE and the ADAPT family select operators from quantum Fisher information. |
| **AMBER vs a coarse knowledge-based potential inside one VQE, everything else held** | **NOT FOUND anywhere.** Genuinely new — and this sprint has just found it is **not executable as briefed** without an `E ∘ Relax₁` arm (§5). |
| **error-decorrelated structure generation, measured as ρ(e_coherent, e_pool)** | **NOT FOUND** for structure generators. The classical ambiguity theory is old; the specific measurement is not in this literature. |
| the 2n−5 tangent rank | standard distance-geometry rigidity (n−2 virtual angles + n−3 virtual dihedrals). Correctly labelled **EXACT** in Sprint 19. |
| the "shared referent floor" (§4.7) | elementary statistics, but **I could not find it named** as a check in the structure-prediction or ensemble literature. |

**What is genuinely new here, ranked by how much I would defend it:**

1. **`search saturates, discrimination binds`, with matched controls** — externally reproduced twice
   (arXiv:2609.02113 at a 2.18–3.26 Å selection gap; arXiv:2606.21241 on cost-Hamiltonian ranking),
   and now on a **fifth internal instrument** (§2.4: fifteen generators at 8192 evaluations, none
   beating zero-evaluation retrieval).
2. **Cross-pair sign coherence as the harmful property of a distance-prediction error**, with a
   control matched to machine precision on residual RMS. I know of nothing equivalent.
3. **The control discipline itself.** The matched-random gate; the zero-information gate; best-of-N
   from the *untrained* circuit; the min-of-N null with its band named; the gauge test over
   encodings; the leave-one-out ε-dominance control. **Nothing in the quantum-folding literature
   carries any of these**, and three of the papers above would not survive one of them. Sprint 19
   said this and I confirm it after reading the methods: the closure of the quantum branch in this
   programme is the most controlled negative result in the area. That is a publishable contribution
   in its own right and it is a *methods* contribution, not a folding one.

**What the field is doing that this programme is not, and should note:** every 2026 quantum-folding
paper I read reports **energies**, **convergence** or **fidelity** and not structural accuracy
against a held-out set. `arXiv:2607.02749` (Graph-VQE) reports *"lower final energies"* as its
result. That is the RMSD-FIRST rule of §1, violated, in four independent groups.

---

# 7. PRIORITY 4 — TWO RECOMMENDATIONS

Judged by expected information value per unit compute, not by novelty. Each carries a mechanism, a
first experiment, a matched control and a pre-registered falsifier. **Priors against each are stated
before the case for it.**

## R1 — Measure Q2's own statistic on the SIXTEEN generators that already exist, before building a seventeenth

**The defect this exposes.** Q2 asks for `ρ(e_coherent, e_pool)` on every candidate source and calls
it "Priority 1's only live lever". **Sprint 19's fifteen samplers cannot answer it**, because they
are all conditioned on the pool: `s19/qb_lib.fit_basins` fits the per-residue von Mises mixtures to
`RP.windows(pdb, "pool")` — **the retrieval pool's own torsions**. Every quantum arm, every
classical learned proposal, every thermostat draws from a law estimated on the pool. **The study
never contained a pool-independent generator**, so its uniform failure to decorrelate is a property
of its construction, not evidence about Q2.

The two exceptions are the zero-information arms, and they are informative in the wrong direction:
`c_helix` carries no pool information at all and realises **3.426 Å** against `pool500`'s **3.230**.
So at the one point on the board where decorrelation is guaranteed, quality collapses.

**First experiment.** One readout pass over the sixteen generators already on disk, storing the
emitted distance fields, and computing per generator: (i) member quality, (ii)
`ρ(e_coh_g, e_coh_pool)` using the **same** coherent/incoherent decomposition `s19/a_source.py`
uses — each family split by **its own** terminal fit — and (iii) realised RMSD. Then plot the
(quality, ρ) frontier. **No model is trained and no structure is invented.**

**Matched control.** `pool500` is the maximally-correlated end and `c_helix`/`c_marg` the
zero-information end; the frontier is read *between* them, and the null for ρ is §4's shared-referent
floor, which must be subtracted here exactly as it had to be in §4.

**Pre-registered falsifier.** *If no generator occupies the region {member quality within 0.3 Å of
the pool's} × {ρ below 0.5 after null subtraction}, Q2's premise is empty on everything measurable,
and Q2 should be closed by measurement rather than left open by effort.* If one does occupy it, it
is the first Priority-1 candidate in four sprints.

**Why this is the highest information value per unit compute available.** It is hours of readout
against a sprint of generator-building, it uses artefacts that already exist, and it can *close* a
question the brief calls the only live lever. A negative here is worth more than a new generator.

## R2 — Attack the repair operator, the last large measured loss nobody has attacked in its actual form

**Priors against, stated first.** (a) `s18` F3 already tested a minimum-displacement spacing
restoration and it cost **+0.023 more** than the projection. (b) `s16` J1 shows AMBER k=30 already
recovers **−0.022** of the tax, so ~13% is known-recoverable by a different repair — and that
restraint constant was **picked on dev-set RMSD with no native-free rule selecting it**. (c)
`s12/forensics_rggate` shows the deployable size channel has **zero skill** against a shuffled null
on the targets where size matters (FAIL18 **+0.058**), so any variant that estimates the right
*scale* natively starts from a measured null. **Expected recovery is modest — plausibly 0.02–0.05 Å
— and I will not sell it as more.**

**The case anyway.** Turning the point cloud into a structure costs **+0.166 [+0.131, +0.202]**,
2× MDE, 16W/110L, 5/5 folds — one of the most reliable effects in the project — and the *operator*
has never been varied, only its output patched. Sprint 18's F3 changes the **spacing**; it does not
change **which manifold point is chosen**. The projection minimises `‖X − avg‖²`, the distance to a
**shrunk, biased estimator**, and inherits that bias in every direction the manifold can represent.

**The untried operator: the manifold-constrained Fréchet mean.**

    argmin over (phi, psi) of   sum_k  min_R  || X(phi,psi) - R . W_k ||^2

with each member's rotation **re-optimised jointly** with the torsions. The pipeline currently fixes
every rotation to the medoid superposition, averages, and *then* projects. For a curved constraint
set those are different objects, and the Fréchet mean **cannot contract**, because the manifold
forbids it. It is native-free, it is one L-BFGS per target over 2n parameters with 75 inner
Procrustes solves, and it has never been run. External support: MCORE (§6.2) reaches the same
conclusion from the other direction — refine *toward* the average under geometric potentials rather
than project onto them.

**Matched control**, per §6.1: the projection of the **same** 75 members' fixed-rotation average —
the incumbent — plus an arm in which the rotations are re-optimised but the manifold constraint is
dropped, which isolates "re-optimising rotations" from "constraining to the manifold". Without that
second arm the comparison cannot attribute any effect.

**Pre-registered falsifier.** *If the Fréchet mean's built-chain RMSD is within ±0.03 Å of
`project(avg)` at n=126, R2 is dead* — the two operators reach the same object and the +0.166 is
irreducible exactly as `s18` §6.2 says.

**Cost.** ~1–2 min/target on a quiet box; ~3 h at n=126, subsampleable to a pre-declared 40-target
prefix. Cheap enough that its negative is affordable.

## What I am NOT recommending, and why

- **Building any new generator before R1's measurement.** Fifteen already lost to zero-evaluation
  retrieval (§2.4).
- **Any further work on the CVaR α parameter, on the ansatz, or on entanglement.** §2.6 and §2.3(b)
  price all three at zero with the right controls, and §6.1 explains why α should be expected to be
  worth nothing at this shot budget.
- **Q1 as briefed, until the `E_amber ∘ Relax₁` arm exists** (§5).
- **Anything that consumes the argmin.** A perfect rank-1 is worth −1.74 Å through argmin and −0.03 Å
  through the m=75 average; every argmin arm in four sprints has lost.
- **Re-opening the predictor.** §4 *reduces* the zero-information share of the harmful component, so
  more of it is sequence-dependent than Sprint 19 thought — but P1 already measured eight
  leave-fold-out predictors and none moved past the MDE. §4 changes the *interpretation* of the
  wall, not its height.

---

# 8. WHAT I LEAVE OPEN, AND WHAT WOULD CLOSE IT

| # | open question | what would close it |
|---|---|---|
| 1 | **M7 — why post-fit κ predicts the per-target gap when start-point geometry does not.** I pre-committed in `PREREG_D` Block F **not** to offer a story for it, and I am keeping that. | A **native-free predictor of κ**. Until one exists, κ and the gap share two fits and a common cause is not excluded. Anything else is narrative. |
| 2 | A **quality-matched** version of §4's null (pool windows restricted to the predictors' own distance-from-native band) | One extra arm; the direction of the bias is known (§4.6) and the ordering cannot change, so this refines a number rather than a conclusion. |
| 3 | The **exactly-zero SPSA gradient rate** under CRN + a tail selector (§3.4) | Count `fp == fm` in `core.quantum._spsa`. One line. |
| 4 | Whether **`E_amber ∘ Relax₁`** — AMBER's energy with no relaxation — is a usable objective at all, given §5's definedness result | It is the missing Q1 arm; measuring it *is* the answer. |
| 5 | Q2's premise, per **R1** | R1's falsifier, pre-registered above. |
| 6 | The **repair operator**, per **R2** | R2's falsifier, pre-registered above. |
| 7 | s19 **P5** (short-favouring training loss) — **BLOCKED, not declined** | ESM cache headroom; the pre-registration stands unedited. |
| 8 | The clean **AMBER Pareto** rerun at 41/126 (§3.1) | Time; it is running. Until it finishes, `s19` §6 is PROVISIONAL and `agentC_pareto_report.txt` must be renamed. |
