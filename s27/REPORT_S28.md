# S28 REPORT: THE CVaR-VQE SPINE UNDER THREE NEW FORMULATIONS, AND WHERE THE ERROR LIVES

Sprint 28 of the peptide structure-prediction programme. Repository branch `s26`, work under
`s27/`. Run 2026-09-14 19:00 to 23:40 Pacific (paused on the user's usage limit) and 2026-09-19
20:08 to 22:00 (resumed from `s27/RESUME_S28.md` with zero recomputation). Coordinator plus four
lanes (A, B, C, D) throughout, two second-wave assignments per lane. The ledger is
`s27/LEDGER.md` S28-L0 to S28-L50; every number below carries the entry that holds its artefact
path, and no number is quoted from prose. Style follows `s26/REPORT.md` Parts VII and VIII:
what was tested and found, what moved, what closed, what is newly open, the retractions, and the
one-paragraph answer on the quantum component.

## 1. The question and the answer, stated first

The brief asked for real movement on the built-chain RMSD, or a real quantum result, or both,
with the CVaR-VQE kept as the spine and everything around it (what it sees, how it is encoded,
what it optimises, how its output is consumed) open. The sprint built the three formulations
the record had never tried, measured each on the built chain against production with the full
control set, and had an Adversary lane attack every entry within the hour.

The answer. Nothing moved the built chain. Every deployable arm of every formulation is either
inside its own MDE of production or worse than production beyond it, and every one of those
verdicts survived the Adversary. What the sprint bought instead is the mechanism, measured on
both sides:

- The circuit family the deployed CVaR-VQE uses, read as SIGNED affine weights over the pool,
  contains a structure within 0.25 A of the native on every one of the 126 targets (ORACLE,
  S28-L26b). Expressivity is not the barrier.
- The optimiser reaches the objective's optimum (lam = 0 reproduces `run_cvar_vqe` bit for bit;
  the objective falls on 126/126; S28-L18b, S28-L26b). Optimisation is not the barrier.
- The objective is the barrier, three ways. Its minimiser is away from the native (every
  signed arm worse, S28-L26b). Its gradient at the pipeline's own output carries no directional
  information about the native (cosine -0.03, inside the random-direction null, S28-L23b), and
  a step along it is worth what a random step is worth (S28-L39). And it is not alone: 20 of
  31 native-free scorers in the whole S27 library prefer the projected production average to a
  0.25 A structure (18 once the 5 tie-dominated scorers are set aside, S28-L49), and the one
  two-clause pass is at the mean of the multiplicity null (S28-L48, L49).
- On the quantum side, the first non-diagonal Hamiltonian in the project's history was built,
  differentiated exactly and run end to end. Its hopping term is invisible to the circuit's
  gradient at the deployed width because the pool's similarity graph is a near-rank-one
  typicality projector (S28-L8b, S28-L11). A spread-spectrum graph halves that decay and the
  circuit then finds the coherent basin on 25 to 39% of (target, seed) cells, and the emitted
  structure does not move on either class (S28-L25, S28-L29, S28-L46, S28-L47).

## 2. How S28 was run

- Contract `s27/S28_CONTRACT.md` (12 non-negotiables, the S26 discipline) with addendum 1 at
  19:10 on the user's steer (built chain is the only verdict basis; no cosmetic variants; every
  positive answers "does it survive the averaging bottleneck", "what is new vs S26/S27", and the
  three-way quantum split; lane D gates every positive before it is built on) and addendum 2 at
  the resume (reuse everything, recompute nothing).
- Briefs `s27/briefs/S28A.md` to `S28D.md`; second wave `S28A2.md` (local behaviour of the
  objective), `S28B2.md` (spread-spectrum graph, conditional), `S28C2.md` (recognition audit).
- Pre-registrations `s27/PREREG_S28_A.md` (base + 4 addenda), `PREREG_S28_B.md` (base + 3),
  `PREREG_S28_C.md` (base + 3), `PREREG_S28_C2.md` (base + 1); every falsifier written before
  its number and checked by lane D before the run (S28-L1, L2, L3, L23, L34).
- Every job through `s26/jobrun.py` under the S26 governor (`s26/governor.py` v2.4); 60 job
  records in `s26/jobs_done/s28*.json`. The user's own programs held the box at 84 to 98% RAM
  for most of both sessions; the governor suspended and killed 11 lane jobs as collateral
  (S28-L5, the STATUS log), every one of which resumed from a per-target checkpoint. No result
  was lost or recomputed. The pause and resume are `s27/RESUME_S28.md`.
- Statistics as the contract: `s24.stats_lib.compare` (paired, SE, MDE = 2.8016 SE, iid and
  fold-clustered CIs, W/L never used as evidence), `best_of_k_within` on every grid, matched
  controls in the operator's own space, random-subset nulls for every stratum claim, NaN-poison
  for every deployable operator (run by lane D, not the lane).
- Suite at the close (S28-L45): every test file in the tree ran in the quiet window, on the same
  flag basis as the S26 close: 440 passed, 2 skipped (the same two absent artefacts as S26), 0
  failed, 0 errors on 442 unique tests; all 11 opt-in OpenMM/pipeline checks run and pass.
  S27's artefacts reproduce to the last digit on six random rows across the three bases (seeds
  101 to 106; S28-L14, L38). `s26/examine.py` at the close: 836 modules, 101 pinned hashes with
  no drift, 21 of 21 recorded claims check (S28-L45).

## 3. What was tested and found (built chain unless stated; production 3.2126 A)

### 3.1 Lane A: the amplitude readout (a CVaR-VQE whose output is a signed combination)

The formulation. The deployed circuit's real amplitudes psi_i are signed; C(theta) = sum psi_i
W_i / sum psi_i is an affine combination of the pool windows. Objective F = CVaR_alpha(E; p) - T
H(p) + lam S(C), with S the shipped distogram Bayes-risk of the emitted structure. The first
non-convex consumer of the pool on record (S10-5's affine bound was ORACLE only).

| Measurement | Result | Entry |
|---|---|---|
| ORACLE expressivity of the 27-parameter family (theta chosen against the native) | 0.288 A point cloud best-of-5 (single start 0.36 to 0.42; the best-of-5 is an order statistic, k_eff 3.4); 0.252 A EMITTED on the built chain, 126/126 under 2 A; against 0.608 for a random 27-dim affine subspace | S28-L1b, L13, L26b |
| Recognition, point cloud: circuit at lam 0.3 / 1 / 3 vs production 3.0483 | +0.335 / +0.337 / +0.390, all WORSE (2.0x to 2.4x MDE, 5/5 folds) | S28-L18b, L20 |
| Recognition, BUILT CHAIN: lam 0.3 / 1 / 3 vs production | +0.236 (1.48x) / +0.226 (1.29x, Type-M zone, same sign) / +0.264 (1.59x), fold CI above zero; F1 silent; lam = 0 read as signed weights 6.14 A | S28-L26b, L27b |
| Classical controls at matched budget (unconstrained signed weights on 75 and 500; random 27-dim subspace; softmax simplex; untrained circuit best-of-16) | every signed family worse in proportion to how well it minimises S; the convex simplex lands at production (3.056) | S28-L18b |
| What the objective does with sign freedom | S at the circuit optimum 1.34 < production 1.67 < the ORACLE structure 2.09: it uses the freedom to move AWAY from the native; the native sits at the objective's 36.9th percentile | S28-L18b item 2, L26b |

Verdict: REFUTED for accuracy, as registered. Expressivity is not the barrier; the objective is.

### 3.2 Lane A2: is the objective locally informative at the production point?

| Measurement | Result | Entry |
|---|---|---|
| ORACLE cosine between the objective's steepest-descent direction at the production cloud and the direction to the native (rigid body removed) | DIS -0.034 (SE 0.021), 56/126 positive; FAIL18 -0.143 (SE 0.066); the other 108 -0.016; random-direction reference 0.140; Spearman with production RMSD -0.37 | S28-L23b, L24 |
| The same cosine for every S27 CA channel with a gradient | RG_LAW +0.034, DISTPOT -0.002, CONTACT +0.030, ENV -0.006, CAGEO +0.010: none above the random reference | S28-L23b |
| Deployable step ladder e in {0.1, 0.3, 1.0} A along the descent direction, BUILT CHAIN, vs production | -0.004 (0.21x) / +0.018 (0.74x, sign only) / +0.105 (2.10x, WORSE) | S28-L39, L40 |
| The same step vs the mean of the SAME two projected random directions | -0.005 / -0.000 / -0.014 (0.22x / 0.00x / 0.26x): the descent direction is worth a random direction | S28-L39 |
| One Adam step within the circuit family from the production-reproducing theta | +0.000 / +0.008 / +0.122 (e = 1 WORSE, 1.34x); the family's FAIL18 excess is a real set property and it is harm (random-18 null p 0.001) | S28-L39, L40, L42 |

Verdict: REFUTED, as registered. The objective is blind at production; a trust region around
production has nothing to trust.

### 3.3 Lane B: the first non-diagonal Hamiltonian, H = diag(E) - J A

The formulation. A_ij = exp(-d_ij^2 / 2 sigma^2) on the pairwise CA-RMSD of the 500 windows
(native-free), sigma the pool median; F = CVaR - T H - J <psi|A|psi> with the hopping term's
exact parameter-shift gradient added to the deployed gradient; readouts R1 (the deployed CVaR
tail), R2 (p-weighted average), R3 (p-top-75); controls PERM (graph destroyed, spectrum kept),
RAND (degree sequence kept), the exact ground state by eigensolver; J in {0, 0.1, 0.3, 1, 3}.

| Measurement | Result | Entry |
|---|---|---|
| Gradient variance of the hopping term alone, n = 4..9 | decays -1.7 to -1.8 per qubit (4.06e-3 at n = 4 to 4.16e-6 at n = 9); at n = 9 it is 7,300x below the CVaR-plus-entropy terms at J = 1 (820x at J = 3); the full objective's variance is flat in n at every J | S28-L8b, L11 |
| Mechanism | the rank-one part of A carries 97% of the variance at n = 9 and reproduces the slope; the Perron vector is 95 to 99% the uniform state; a diagonal observable with the same spectrum decays at the same rate: the decay is the near-rank-one SPECTRUM's, not off-diagonality's | S28-L8b, L11 |
| Set-equality departure | `gate_set_equality` passes 4,914/4,914 cells; the deployed tail equals the classical top-m of E at its realised m to 1.1e-13 A on every cell at every J: any R1 movement at J > 0 is the m-ladder term | S28-L21, `s28_B_mladder.json` |
| Point cloud, 39 arms x 126 | R1 within 0.013 A of production at every J; R2 / R3 +0.24 to +0.33 A worse; the eigensolver at J = 1 re-selects the DIS top-75 (Jaccard 0.960) | S28-L21, L22 |
| BUILT CHAIN, 18 arms x 126 (comparator re-projected in-lane, equal to S27's rows at 0.0 on 126/126) | R1 at every J, graph, seed inside 0.46x MDE of production and 0.22x of J = 0; R2 / R3 +0.25 to +0.37 A WORSE (1.28x to 1.58x, fold CI above zero); F1 does not fire (the one BETTER cell is a 0.28 A-worse readout climbing 0.04 A, unreplicated on seed 0, matched by PERM); F3 / F4 silent; the exact ground state ties production (-0.003 / -0.006) | S28-L41, L43, L44 |
| Three-way split, middle leg | the J = 3 ground state is not representable by the circuit (max overlap 0.80 to 0.88 over 16 starts, saturated at 80 iterations), but a representable state has F below the ground state's on 10/12 targets and below every reached F on 24/24 cells: the gap is OPTIMISATION by basin selection; the trained state is BIMODAL on all 126 (43 / 32 cells by seed fully coherent, the rest at the J = 0 state, none between) | S28-L41, L43, `s28_B_represent.json` |

Verdict: REFUTED as registered (prior WORSE or null held). The coupling is real, exact and
correctly differentiated (16 tests); it changes the emitted structure by nothing measurable.

### 3.4 Lane B2: a spread-spectrum graph (the one follow-up the mechanism licenses)

| Measurement | Result | Entry |
|---|---|---|
| kNN graph (k 5, 10; degree-normalised, unit spectral norm), hop-only gradient variance | decays -1.00 / -1.03 per qubit (Gaussian -1.84), 30 to 60x the variance at n = 9; the pre-registered -1.0 line is straddled, not cleared; the remainder (not the rank-one part) now carries 96 to 99% of the variance | S28-L25, L26, L29 |
| Coherence classes at J = 1 / J = 3 (class counts, never a mean; retraction R4) | at J = 1 the circuit is at the J = 0 state on nearly every cell; at J = 3, k = 10 it is fully coherent on 49 / 33 of 126 cells by seed (Gaussian 43 / 32) | S28-L29, L46 |
| Endpoint at the one licensed cell (k = 10, J = 3; prereg addenda 2 and 3), point cloud | R1 inside 0.42x MDE of J = 0 and 0.37x of production; R3 +0.28 A worse; REAL vs PERM silent; the kNN eigensolver (PR 53 cluster) ties production (0.54x); no measurable endpoint difference between the coherent and incoherent classes | S28-L46, L47 |
| The two R3 cells at 0.7x MDE, BUILT CHAIN | +0.31 A worse than production (1.45x, 5/5), within 0.35x of their J = 0 twin | S28-L46 |

Verdict: REFUTED as registered; B2 closed at its scope.

### 3.5 Lane C, Part 1: a FAIL18 detector on a never-used feature class

The 18 targets on which production fails (5.83 A against 2.58 A on the other 108, S27) are
where every alternative Hamiltonian helps; a native-free switch would be worth 0.03 to 0.06 A.
Seven router constructions had failed (S22 L7, S23 L7, S26 L110/L115). New inputs: the pool's
own statistical-potential distribution and its agreement with the distogram (21 features).

| Measurement | Result | Entry |
|---|---|---|
| Nested leave-fold-out ridge logistic, held-out AUROC vs a 500-draw label-permutation null | SP block 0.468 (null p95 0.629); SP + controls 0.608 (p 0.10); the S22 set through the same harness 0.591: no block clears its null | S28-L6, L8 |
| Best single feature priced as an order statistic over 26 | DISTPOT_spread 0.706 vs max-over-26 null p95 0.720 (p 0.076); it is a length proxy (corr -0.43 with n; 0.593 after residualising) | S28-L6 |
| ORACLE ceiling: switch exactly the 18 to DIS+DISTPOT / DIS+ENV, BUILT CHAIN | -0.020 (0.28x MDE) / -0.047 (0.62x): even a perfect detector's prize is under detectability at n = 126 | S28-L6 |

Verdict: CLOSED on both grounds; the switched arm was never run.

### 3.6 Lane C, Part 2: readouts that consume ranking information without averaging it away

Medoid-plus-k-neighbours, a ranker-based outlier trim, and a diversity-preserving weighted
average, each with a permuted-ranker control in the operator's own space and grid pricing.

| Measurement | Result | Entry |
|---|---|---|
| Every cell, BUILT CHAIN, vs production | on the harmful side, +0.015 to +0.149 A (0.29x to 0.99x MDE); none clears its MDE either way | S28-L15, L16, L32, L33 |
| Ranker trims vs random trims of the same size | worse, fold CI above zero (sign measured, size not) | S28-L32 |

Verdict: CLOSED. All three are convex combinations of the same 75 members; none escapes the
averaging bottleneck of S27 section 6, and consuming the ranking makes them worse.

### 3.7 Lane C2: which native-free scorer, if any, recognises the structures the family can express?

ORACLE diagnostic with a closure falsifier. Per target: production, the ORACLE circuit
structure (best-of-5 and single start), the 0.61 A random-subspace structure, and two matched
controls (a random signed combination at the same distance; a Gaussian perturbation at the same
RMSD). Statistic: the fraction of targets on which a scorer prefers the ORACLE structure to
production, ties at 0.5, fold-clustered CI, one max-over-scorers sign-flip null, the pool-member
control (does the scorer prefer any real trace to the contracted average?).

| Measurement | Result | Entry |
|---|---|---|
| CA level, 15 scorers, two control seeds | 9 of 15 prefer the contracted production average to a 0.29 A ORACLE structure with the fold CI below 0.5 (DIS 0.206, DIS_MEAN 0.151, SS_MATCH 0.393, CONS_POOL 0.024); the learned linear combination prefers random signed combinations | S28-L35 |
| The one CA-level candidate, CAGEO (0.611) | VETOED by the pool-member control: it prefers any real protein trace to the 22%-contracted average exactly as often (-0.007, 0.18x) | S28-L36, L37, R3 |
| BUILT CHAIN, 31 scorers (16 backbone + 15 CA on the projected chains), one max-over-31 null | 20 of 31 prefer the projected production average to a 0.25 A ORACLE structure with the fold CI below 0.5 (18 informative once the 5 tie-dominated are set aside, S28-L49; DIS on 93% of targets, LEG on 79%, RAMA on 59%); 6 coin tosses; CAGEO collapses to anti-recognition on ideal geometry (0.611 to 0.421) as registered; the one two-clause pass, CONTACT@chain 0.583, sits at the max-over-31 null's mean (0.579, p95 0.627, p_max 0.39) | S28-L48, L49 |

Verdict: the closure claim STANDS over the whole library. No native-free scorer in the project
recognises a near-native structure when it is offered one.

### 3.8 The instrument

| Measurement | Result | Entry |
|---|---|---|
| The built chain's reproducibility under a floating-point input difference | lane A's re-projection of the production cloud (input equal to 5.7e-14) differs from S27's rows by 0.006 A mean with a 0.5 A tail on one target (a branch flip of the multi-start projection); lane B's re-projection equals S27's rows at 0.0 on 126/126: the floor belongs to the input difference, not the projection; every S28 MDE is above it | S28-L18, L27b, L43; R1, R2 |
| The B chain's J = 0 anchor | equals S27's `vqe_rows.jsonl :: DIS` at 0.0 on 126/126, both seeds | S28-L9, L38 |
| Lane A's own NaN-poison test | was a determinism test; lane D's `tests/test_s28_D.py` runs the real one (natives NaN'd through a monkeypatched loader; every deployable field bit-identical) and it passes | S28-L10 |

## 4. What moved

Nothing on the built chain. The largest deployable movement in the sprint is lane B's deployed
readout at J = 0.1, seed 0, -0.0024 A against its J = 0 twin (0.22x MDE, S28-L41), and the
largest against production is A2's e = 0.1 step, -0.004 A (0.21x, S28-L39). Both are inside
their own noise and neither is called anything.

What did move is the account of why. Before S28 the record said "ranking information is
anti-useful on an averaging readout" (S27 section 6) and "the objective does not rank the
native" (S8-9). After S28 it says: the deployed circuit can express a near-native structure
on every target and the optimiser can reach the objective's optimum, and the objective (and
every native-free scorer the project owns) points away from that structure, globally (its
minimiser), locally (its gradient at production) and comparatively (its preference between
the two). The convexity of the readout was not the binding constraint; sign freedom makes the
output worse, not better, in proportion to how well the objective is minimised.

## 5. What closed

1. The signed-amplitude readout as an accuracy mechanism (S28-L26b, L27b).
2. The shipped objective as a locally informative function at production (S28-L23b, L24, L39,
   L40); with it, any trust-region or gradient-step refinement from production.
3. The Gaussian-graph hopping Hamiltonian as a selector or readout mechanism, at every J, on
   every readout, on the built chain (S28-L41, L43).
4. The spread-spectrum follow-up at the one cell the mechanism licenses (S28-L46, L47).
5. The FAIL18 detector on the statistical-potential feature class, and with it the eighth
   router construction; the ORACLE prize of a perfect detector is itself under MDE on the
   built chain (S28-L6, L8).
6. Ranking-consuming readouts (S28-L32, L33): the "consumer of ranking information" question
   S27 opened.
7. Recognition of the ORACLE structures by any native-free scorer in the S27 library, at the
   CA level and on the built chain, with the multiplicity priced (S28-L35, L36, L48, L49).
8. Two instrument questions: the built chain's floor (input difference, not projection; S28-L43)
   and the set-equality theorem's behaviour under a non-diagonal H through the tail readout
   (identity to 1e-13, S28-L21).

## 6. What is newly open

- The objective. Every route in this sprint died at the same place: the native-free objective
  and the library of native-free scorers built around the same distogram all prefer the
  contracted average to a near-native structure (S28-L48). What would recognise one is not in
  the project's channel library. The prior's derivative was already the only steep lever
  (S23); S28 says the same from the readout's side: the family that can express the answer
  exists, and nothing native-free can pick it out. A scorer trained to prefer the ORACLE
  structure to production would be trained on the native and is not deployable; the open
  question is whether any native-free signal outside the distogram's own error structure
  (S23 L9's 68% common mode) exists at peptide length. S28 measured the library and found none.
- The coherent basin. At J = 3 the circuit either finds the sign-coherent ground-state basin or
  stays at the J = 0 state, with nothing between; the coherent share is 25 to 39% of cells by
  seed and graph (S28-L41, L46). Which
  cells is not a function of production's RMSD or of FAIL18. It is a property of the
  optimisation landscape of a non-diagonal H under the deployed ansatz, measured for the first
  time; it has no accuracy consequence here because both basins emit the same structure, and
  it is the one quantum-side object the sprint leaves for a study of its own.
- The 0.5 A branch-flip tail of the projection (S28-L27b): a 1e-14 input difference moves one
  target in 126 by half an angstrom. Inside every MDE quoted, but a reproducibility property
  worth a pinned tolerance in the projection's own tests.

## 7. Did the quantum component earn anything this sprint?

No accuracy, and it was given every chance: its state was read as signed weights (the one
readout that could have escaped the convex bottleneck), it was handed the first non-diagonal
Hamiltonian in the project's history with an exact hopping gradient, and that Hamiltonian was
then re-built on a spread-spectrum graph at the one cell the mechanism licensed. On the built
chain every deployable arm is inside its MDE of production or worse, every verdict survived the
Adversary, and the classical spectral counterpart (the exact ground state) ties production
wherever the circuit does. What the quantum component earned is three measured properties
that stand on their own and were not known before: the pool's similarity graph is a near-rank-
one typicality projector whose hopping term decays at -1.7 to -1.8 per qubit and is 7,300x
below the diagonal terms at the deployed width (S28-L8b, L11); the CVaR tail read from a
non-diagonal H is still exactly the classical top-m prefix to 1e-13 at every J (S28-L21), which
is the set-equality theorem surviving its first genuine test; and under a spread spectrum the
variational optimiser is bimodal, finding the coherent basin or not on a per-cell basis with
nothing between (S28-L41, L46). The first two say why the circuit cannot matter through this
readout; the third is a landscape fact about the ansatz that a future non-diagonal formulation
would have to design around. The honest summary is that the quantum component's ceiling in
this pipeline was measured to be the objective's ceiling, and the objective's ceiling is zero.

## 8. Retractions (`s27/RETRACTIONS_S28.md`)

| | Retracted | By | Kind |
|---|---|---|---|
| R1 | S28-L7: the built chain's floor is 1e-5 A (one target) | S28-L18 | scope |
| R2 | S28-L18: the floor is 0.02 A per target, 0.0025 mean (22 targets) | S28-L27b; then S28-L43 (the floor is lane A's input difference) | scope |
| R3 | S28-L35: the closure claim is falsified at the registered bar by CAGEO | S28-L36, L37 (pool-member control) | veto of a candidate |
| R4 | S28-L11 / L8b / L21 / L25: "the circuit collects a third of its hopping bound" (a mixture of coherent and J = 0 cells read as partial alignment) | S28-L43, L44 | reading; no number changed |

Also struck: lane B's "slightly worse in the coherent basin" (post-hoc subgroup, 0.26x, S28-L43
/ L44); lane A2's FAIL18 "regime" for the raw step (18-target noise by the random-18 null,
scoped to the circuit arm, S28-L40 / L42).

## 9. What damaged the pre-registered expectations

- The ORACLE ceiling of a 27-parameter curved family (0.25 A emitted, 126/126) was expected
  far below 2.954 and came in far below even that; the affine-hull rows are exactly zero by
  dimension counting (3n <= 48 coordinates, 75 members), which is why S10-5's affine bound was
  0.035 (S28-L1b).
- The Gaussian graph was expected to give the circuit "something to do"; it gave it a
  projector on the uniform state (S28-L8b).
- Lane D's own rule "overlap near 1 means optimisation shortfall" was wrong in form: the ground
  state is only 85% representable AND the gap is optimisation, because a representable state
  with lower F than every reached state exists (S28-L43, R-D2 in `s27/s28_D_FINDINGS.md`).
- CAGEO looked like the first native-free scorer to prefer a near-native structure; it prefers
  any real trace to a contracted one (S28-L36).

## 10. What was not done, and why

- No AMBER scorer in the C2 audit (one AMBER process at a time; the CA and backbone channels
  closed the question without it; S28-L48 "what I did not do").
- No lam sub-grid, no finer J or sigma, no k beyond {5, 10}: contract addendum 1 forbade
  cosmetic variants, and the mechanism licensed exactly one B2 cell.
- No second-seed re-run of any null (none was owed: nothing positive).
- The switched FAIL18 arm (F1 did not fire; addendum 2 of `PREREG_S28_C.md`).
- No AMBER relaxation of any S28 output: nothing reached the built chain as a candidate.

## 11. Operating record

- Two sessions; one pause on the user's usage limit with a full resume note; zero recomputation
  on resume (checkpoints B 103, A2 94, C2 80 of 126 verified before relaunch; every job resumed
  at its first undone target; `s27/RESUME_S28.md`, `s27/STATUS.md`).
- 11 governor kills of lane jobs under the user's own memory load (chrome 3 to 6.4 GB, Zoom, a
  Windows servicing worker); all resumed from checkpoints; one coordinator kill of lane D's
  forking test job to protect the lane jobs (S28-L5).
- Ledger numbering collided five times (two lanes appending in the same minute); resolved with
  b-suffixes and correction entries (S28-L12, L19, L27, L28). Lanes were told to number from
  the tail in one process.
- Three coordinator timestamp errors in STATUS, corrected in place; the rule "run the clock in
  the same command as the stamp" is in memory.
- A handed-back lane stays idle until messaged (lane C lost 25 minutes at the resume); the
  coordinator owns every wake-up.

## Appendix. The ledger (S28-L0 to S28-L50), one line each

| Entry | Lane | Content |
|---|---|---|
| L0 | coord | the brief, the lanes, where the coordinator placed the bets, the registered priors |
| L1, L2, L3 | D | prereg checks A, B, C: STANDS WITH CAVEAT |
| L4, L10, L38, L45 | D | suite status at start, 19:45, resume, close (440 / 2 / 0) |
| L5 | coord | the test-job kill; heavy files deferred |
| L6, L8 | C, D | FAIL18 detector null; ORACLE switch ceiling under MDE; STANDS |
| L7, L7b | C | S27 chain rows reproduced bit-exactly on 6 targets |
| L8b, L11 | B, D | F5 hopping trainability; rank-one mechanism; STANDS WITH CAVEAT |
| L9 | D | B's J = 0 anchor bit-exact |
| L1b, L12, L13 | A, D | ORACLE expressivity 0.288 / 0.252 emitted; STANDS WITH CAVEAT |
| L14 | D | S27 reproduction seed 102 |
| L15, L16 | C, D | Part 2 point cloud (intermediate); STANDS |
| L17 | A | partial note, superseded |
| L18, L27b | D | the built chain's floor, scoped (R1, R2) |
| L18b, L19, L20 | A, D | recognition point cloud: every signed arm worse; STANDS |
| L21, L22 | B, D | hopping point cloud; m-ladder identity; STANDS |
| L23 | D | second-wave prereg checks (A2, B2) |
| L23b, L24 | A, D | A2.1 ORACLE cosine -0.03; STANDS |
| L25, L26 | B, D | B2 kNN trainability; middle leg; STANDS WITH CAVEAT |
| L26b, L27, L28 | A, D | THE A VERDICT: worse at every lam on the built chain; STANDS |
| L29 | B | F5-B2 second clause fails; kNN decomposition |
| L30, L31 | A, D | A2.2 ladder on the point cloud; STANDS |
| L32, L33 | C, D | THE C PART 2 VERDICT: closed on the built chain; STANDS |
| L34 | D | prereg check C2 |
| L35, L36, L37 | C, D, C | C2 CA level; CAGEO vetoed (R3); accepted |
| L39, L40, L42 | A, D, A | THE A2.2 VERDICT: refuted on the built chain; STANDS; wording adopted |
| L41, L43, L44 | B, D, B | THE B VERDICT: refuted on the built chain; STANDS; R4; wording adopted |
| L46, L47 | B, D | THE B2 ENDPOINT: refuted at the licensed cell; STANDS |
| L48, L49, L50 | C, D, C | THE C2 BUILT-CHAIN AUDIT: closure stands over 31 scorers; STANDS WITH CAVEAT; caveats adopted |

Findings files: `s27/s28_A_FINDINGS.md`, `s27/s28_B_FINDINGS.md`, `s27/s28_C_FINDINGS.md`,
`s27/s28_D_FINDINGS.md` (with the qualifier table for every number above). Results:
`s27/results/s28_*.json`, `s28_*.jsonl`; job records `s26/jobs_done/s28*.json`; logs
`s26/logs/s28*.log`. Second-wave briefs `s27/briefs/S28A2.md`, `S28B2.md`, `S28C2.md`. No
further waves are planned; this report closes the sprint.

## 12. What to work on next: making cost reduction mean RMSD reduction

The full programme, with the meter's four numbers, falsifiers, priors and costs, is in the
published page `s27/REPORT_S28.html` (https://claude.ai/artifact/5iE59KicGARNF5CBjofZFw),
section 8. In one paragraph: the shipped cost's Spearman with RMSD along the ladder from
production to the native is -0.40 (`s27/results/s28_C2_chain_summary.json :: ladder_rho ::
DIS@chain`), its gradient at production is blind (cosine -0.03, S28-L23b) and the native is at
its 37th percentile (S28-L18b). The programme, in dependency order: (1) package the meter
(`s27/cost_audit.py`: ladder rho, cosine at production, native percentile, preference vs
production with the pool-member control; all four pieces exist in `s27/`); (2) calibrate the
2x over-confident posterior leave-fold-out and re-read the meter (hours; the identified
mechanism, S25 L2 and S28-L35); (3) score the de-contracted structure (the average is 22%
contracted and two scorers' preference was the contraction, S28-L36); (4) the pool-residual
prior: train, leave-fold-out on `s8/generate_univ`, a head that predicts the pool's 68%
common-mode error (S23 L9) from ESM features, with the S19 error-coherence check on its
mistakes; (5) an RMSD-supervised cost on the ORACLE ladder structures with candidate-only
features, prior null. Not to be run: any further readout, graph, J, lam, k or sigma variant;
any global rescaling; a ninth router; harder optimisation of the current cost.
