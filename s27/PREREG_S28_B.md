# S28 LANE B PRE-REGISTRATION -- THE FIRST NON-DIAGONAL HAMILTONIAN: ENERGY PLUS HOPPING ON THE POOL GRAPH

Written 2026-09-14 19:10 before any endpoint number (any RMSD to a native) was read. Brief:
`s27/briefs/S28B.md`. Contract: `s27/S28_CONTRACT.md`. Code: `s27/s28_B_hop.py`, tests
`tests/test_s28_B.py`, results `s27/results/s28_B_*.json`.

## 0. What was read before this file, and the one native-free probe

The graph was probed on three targets (1A13, 2N9M, 9WXK) with NO native quantity read
(`s26/logs/s28B_probe_graph.log`, job `s28B_probe_graph`, peak RSS 0.329 GB, 10 s), to put
the J grid on a measured scale rather than a guessed one. Findings that shape the design:

- pairwise CA-RMSD over the 500 windows: median 3.97 to 5.00 A, p10 1.4 to 3.0, p90 5.8 to 6.9;
- with sigma = median, A_ij = exp(-d^2 / 2 sigma^2) has mean off-diagonal 0.61: the graph is
  DENSE and SMOOTH, lambda_max 311 to 321 against lambda_2 37: A is close to rank one;
- the weighted degree (row sum of A) is strongly ANTI-correlated with E = zrank(DIS score):
  corr(degree, E) = -0.47, -0.89, -0.90. The pool's most typical members are its best-scored
  members. The Perron vector of A is therefore a typicality direction, the object S27 L3
  measured as harmful when it enters selection (CONS +0.29 A);
- with A scaled to unit spectral norm, the exact ground state of diag(E) - J A stays
  localised (participation ratio 1.0 to 1.5) at J <= 0.3, delocalises to PR 52 to 80 at J = 1,
  and to PR 277 to 335 at J = 3. Its 75 most probable states have Jaccard 0.92 to 0.97 with the
  DIS top-75 at J = 1 and 0.67 to 0.97 at J = 3.

## 1. Held fixed (identical to S27 / S25)

The 126 dev targets sorted, 5 pinned folds, the shipped K = 500 pool
(`s24.d_harness.Candidates.from_universe`), E = zrank(shipped distogram Bayes-risk score) from
`s27/cache/<pdb>.npz :: DIS` through `s27.run_pool.zr`, the S25 selector settings (9 qubits over
500 + 12 padding states, padding energy max + 10 sd, depth 3, 80 Adam iterations at lr 0.15 on
the exact parameter-shift gradient, alpha 0.18, T 0.5, seeds 0 and 1), the deployed readout
(uniform coordinate average in the retained set's medoid frame, `s24.d_harness.readout_uniform`),
exact ties in any top-m broken by the S27 stable per-target key (`s27.run_pool.topm`), statistics
`s24.stats_lib.compare` (paired, MDE = 2.8016 x SE, iid CI beside the fold-clustered CI; the
decision is on the fold CI). Point cloud is the screening basis; the BUILT CHAIN
(`s24.d_harness.readout_projected` = `s12.instrument.project`) is the reporting basis.

Anchors that must reproduce before any other number is read: the J = 0 VQE arm at seed 0 must
equal S27's DIS VQE arm bit-for-bit per target (`s27/results/vqe_rows.jsonl`, config DIS, seed 0:
mean 3.058016; seed 1: 3.056181); the DIS top-75 point cloud 3.048338.

## 2. The Hamiltonian and the objective

- Graph: d_ij = CA-RMSD after optimal superposition of windows i and j
  (`core.backend("numerics").kabsch_rmsd_batch`, one row per call), symmetrised, zero diagonal.
  sigma = median of the 124,750 off-diagonal d_ij (native-free). A_raw_ij = exp(-d_ij^2 / 2
  sigma^2), A_raw_ii = 0. NORMALISATION: A = A_raw / lambda_max(A_raw), so that <psi|A|psi> <= 1
  for every state and the hopping term -J <psi|A|psi> is bounded by J in the units of E, whose
  sd is exactly 1 (zrank) and whose consecutive-rank gap is 2 sqrt(3)/499 = 0.00694 and range
  2 sqrt(3) = 3.464. Padding states (12 of 512) have zero rows and columns.
- H = diag(E) - J A on the 512-state register. Real symmetric, off-diagonal for J > 0.
- Objective, keeping the spine: F(theta) = CVaR_alpha(E; p_theta) - T H(p_theta)
  - J <psi_theta|A|psi_theta>. The first two terms and their gradient are
  `core.quantum.free_energy` unchanged; the hopping term is evaluated on the real statevector
  (`StatevectorCircuit.state`) and its gradient by the two-term parameter-shift rule
  (pi/2 shift, factor 1/2, exact for the RY generator) on the batched shifted states, added to
  the deployed gradient. At J = 0 the loop is `core.quantum.run_cvar_vqe` line for line
  (same RNG, same Adam), so J = 0 is bit-identical to `s24.d_harness.arm_vqe`.
- A statement that matters for the quantum-side reading: for J > 0 the objective is NO LONGER
  A FUNCTION OF p. Two states with identical p differ in F through the signs of their
  amplitudes (A >= 0 entrywise, so aligned signs maximise <psi|A|psi>). Every earlier objective
  in this project was a function of p alone.
- Grid: J in {0, 0.1, 0.3, 1.0} (the brief's) plus J = 3.0 (the delocalised limit the probe
  showed; added here, before any endpoint number). Any minimum over the grid is an order
  statistic and is priced with `ST.best_of_k_within` and the split-half transfer.
- Graph scale: sigma = median is PRIMARY. sigma = median/2 (sharper) is a SECONDARY graph run
  on the same grid and reported separately with its own pricing; it does not change the
  primary decision.

## 3. Readouts (each a separate arm; every arm consumes the same p or state)

Three consumption operators applied to a probability vector p over the 512 states (padding
mass discarded, p renormalised over the 500 real candidates where a renormalisation is needed):

- R1 TAIL (the deployed readout): the alpha-tail support of p along the E order
  (`s22.qcand_lib.exact_face`), real candidates, uniform average. m is the realised tail size.
- R2 WEIGHTED: the p-weighted coordinate average over all 500 real candidates in the
  p-weighted consensus-medoid frame (`core.pipeline.consensus_medoid` with weights, the S23 L8
  convention, applied to the full p rather than the tail because at J > 0 the mass outside the
  E-prefix is the object of interest). At J = 0 this is a near-uniform average over the whole
  pool (the trained state carries 8.8 of 9 bits) and is expected to be far worse than R1; it is
  kept because it is the only readout that consumes p without a prefix truncation.
- R3 PTOP75: the 75 most probable real candidates (ties by the stable key), uniform average.
  This holds the production rung m = 75 fixed and lets ONLY THE SET change, which is where an
  off-diagonal term can act; R1 cannot see it (section 5).

Two sources of p:
- VQE: p_theta at the Adam optimum of F, seeds 0 and 1.
- GS: the exact ground state of H by `numpy.linalg.eigh`, p = v^2 (the classical spectral
  counterpart; at J = 0 it is the one-hot argmin, so R1 gives m = 1 and R3 is degenerate: 74
  exact-zero ties broken by the key; both are reported and labelled degenerate, never
  compared as a result).

Arms: {VQE(s0), VQE(s1), GS} x {R1, R2, R3} x J in {0, 0.1, 0.3, 1.0, 3.0} x graph in {REAL,
PERM, RAND} x sigma in {median, median/2}. All on the point cloud. Built chain (section 7) for
J = 0 R1, the best J R1 (by the point cloud, priced), and any VQE arm at 0.7x MDE or better on
the point cloud against its J = 0 counterpart.

## 4. Controls

- J = 0: reproduces S25/S27's VQE arm (3.0580 seed 0) bit-for-bit; the comparator for every
  hopping arm of the same readout, seed and graph.
- PERM: A_perm = P A P^T with P a stable random permutation of the 500 real candidates
  (`s27.run_pool.rng_for(pdb, "s28B_perm")`). Same spectrum, same degree multiset, same
  normalisation constant; the correspondence between graph position and E destroyed. The
  analogue of S27's rank-permuted channel control.
- RAND: a random geometric graph with the same degree SEQUENCE: 500 points iid N(0, I_3), the
  same Gaussian kernel with its own median-distance sigma, then symmetric Sinkhorn scaling
  D R D so that row i's sum equals the real graph's row-i sum to 1e-8 (candidate i keeps its
  own degree; only who-is-similar-to-whom is randomised). Normalised to its own unit spectral
  norm. This separates "degree is anti-correlated with E" (kept) from "the specific
  similarity structure" (destroyed).
- The random-75 null: `s27/results/pool_rows.jsonl :: rand_mean` (16 draws per target, mean
  3.4209), not recomputed.

## 5. What the set-equality theorem does and does not say at J > 0 (registered reading)

The S24/S25 theorem is a property of the tail-reading operator `cvar_from_probs(E, p, alpha)`:
for ANY p the alpha-mass is allocated along the E order, so the tail's support is a subset of an
E-prefix. That holds at J > 0 by construction, so `gate_set_equality` passing 126/126 at every J
is NOT a finding and will not be reported as one. What CAN depart at J > 0, and is measured:

- the realised rung m (the prefix the alpha mass reaches) against its J = 0 value, per target;
- holes (exact zeros inside the prefix; a generic RY state has none);
- the p-ranked set: Jaccard(R3 set, E top-75) and the p-mass on the E top-75 (real states);
- total-variation distance 0.5 sum |p_J - p_0| for the same seed and graph;
- participation ratio 1/sum p^2; the hopping value <psi|A|psi>, its same-sign upper bound
  sum |psi_i| A_ij |psi_j|, and the sign coherence (sum psi)^2 / (sum |psi|)^2.

The "first measured departure from the set-equality regime" is therefore the R3/TV/Jaccard
column, not the gate.

## 6. Falsifiers (registered before any endpoint number)

F1 "HOPPING HELPS": some VQE arm (readout R, J > 0, graph REAL, sigma median) beats the
  SAME readout at J = 0 on the BUILT CHAIN beyond its own MDE, with the fold-clustered CI
  excluding zero on 5/5 folds, on BOTH seeds, AND beats its PERM control at the same J beyond
  MDE. Screening on the point cloud at 0.7x MDE or better decides which arms are projected.
  Registered prior: WORSE or NULL at every J. Reason: the probe shows the hopping direction is
  the typicality direction (corr(degree, E) -0.5 to -0.9; the Perron vector of A), which S27 L3
  measured as anti-useful on an averaging readout (+0.29 A for DIS+CONS, worse than its own
  permuted control at 1.7x MDE); at J <= 0.3 the ground state is localised (PR <= 1.5) and the
  trained p at T = 0.5 is near-uniform, so R1 is expected to sit inside its MDE of J = 0 (a
  null the instrument cannot distinguish from "identical"); at J >= 1 the retained set moves
  toward the pool's mode and the arms are expected WORSE with the fold CI above zero.
  If R3 at J = 1 lands inside its MDE of J = 0 R1, that is the Jaccard 0.92 to 0.97 of the
  probe reading back (the graph re-selects the DIS top-75), not a result.

F2 "THE CIRCUIT MATTERS": at the same J and graph, VQE-R beats GS-R beyond MDE with the fold
  CI excluding zero, for R in {R2, R3} (R1 at GS is m = 1 and is not compared). Registered
  prior: the eigensolver wins or ties on R3 at J >= 1 (its R3 set is nearly the DIS top-75);
  on R2 the VQE's T = 0.5 entropy makes its p near-uniform, which is a regularisation effect
  and would be reported as such, never as a quantum effect.

F3 "THE PERM CONTROL": at each J, REAL vs PERM. If REAL is not better than PERM beyond MDE,
  the specific similarity structure contributed nothing measurable and any REAL-vs-J=0 effect
  is a property of the amplitude of the coupling, not the graph.

F4 "THE RAND CONTROL": at each J, REAL vs RAND. If REAL ties RAND, the degree sequence (which
  is anti-correlated with E) carries everything the graph does.

F5 TRAINABILITY VS J (property measurement, no RMSD): gradient variance of F over
  theta ~ N(0, 0.6^2), depth 3, exact parameter shift, at n = 4..9, with E the standardised rank
  ladder over the DIS top-2^n candidates and A their real graph (renormalised), J in {0, 0.1,
  0.3, 1.0, 3.0}, 120 draws per cell, on S27's 12 trainability targets (`P.targets()[::11][:12]`),
  median over targets. Beside it: the HOPPING-ONLY linear cost <psi|A|psi> (J = 1, no CVaR, no
  entropy: the textbook linear cost with a NON-diagonal observable) and the S25 linear diagonal
  control (alpha = 1, T = 0, J = 0). Reported: Var[dF/dtheta_0] and E|g|^2/P per (n, J), the
  fitted log2 slope per qubit, and the ratio Var(J)/Var(0) at each n. Registered falsifier for
  "the off-diagonal term changes trainability": the slope per qubit at J = 3 differs from the
  slope at J = 0 by more than 0.3 (the S27 T11 falsifier was a 2x variance ratio at fixed n;
  here the slope is the object). No prior on the sign. No slope will be called a barren
  plateau or its absence.

## 7. What is projected to the built chain, and when

After the point-cloud screen: J = 0 R1 (both seeds; seed 0 reproduces S27's chain DIS row
3.2126 only up to the VQE's realised m, so it is projected afresh), the best J for R1 (by the
point-cloud mean; priced as an order statistic over the 4 non-zero rungs), and every VQE arm
whose point-cloud effect against its J = 0 counterpart is at or beyond 0.7x its MDE in either
direction (a WORSE arm is projected too, so the chain verdict is not one-sided). GS-R3 at the
same J as the best VQE R3 is projected as the F2 comparator. At most eight arms x 126
projections (about 3.4 s each, `s26/jobs_done/s27_chain.json`).

## 8. Diagnostics (ORACLE where they read a native, labelled at every appearance)

ORACLE: the per-target RMSD of every arm; FAIL18 (`s23/d_scale_globalcapture.py :: FAIL18`)
vs the 108 others, strata reported as in S27 `strata.json`. Native-free: participation ratio,
overlap with the DIS top-75, mean pairwise RMSD of the retained set (the cluster it selects),
the realised m, entropy, ESS, the hopping value and sign coherence, cost per target.

## 9. Integrity

- NaN-poison: `select_all(W, E, key, ...)` (the native-free half) is run with nat_ca and
  oracle_rr replaced by NaN and its outputs (sets, p, coordinates) asserted bit-identical
  (`tests/test_s28_B.py`). Natives are read only in `oracle_rmsd`.
- No native quantity chooses J, sigma, the readout, a subset or a stopping rule. The "best J"
  is chosen on the point-cloud endpoint over the grid and PRICED, not believed.
- Memory: the probe peaked at 0.329 GB (graph 500 x 500 float64 = 2 MB, statevector 512, the
  instrument's universe cache is the bulk). Estimate 0.4 GB per job; one job at a time under
  `s26/jobrun.py --agent S28B`. Checkpoint: one JSONL row per (pdb, arm) appended atomically;
  resumable.
- Cost (from S27's 0.32 s per VQE run and the probe): about 5 J x 3 graphs x 2 sigma x 2 seeds
  = 60 VQE runs per target plus 30 eigh calls and one 500 x 500 graph: about 40 s per target,
  90 min for 126. Chain: 8 arms x 126 x 3.4 s = 57 min.

## 10. Rule 0 forks, the alternative not taken

- Sigma: the median pairwise RMSD (brief) rather than a fixed Angstrom value; median/2 is the
  registered secondary. A learned or oracle sigma is not run.
- Normalisation of A: unit spectral norm rather than unit mean degree. Under mean-degree
  normalisation the same J values are 1.03 to 1.05x larger for this dense graph (lambda_max /
  mean degree = 1.02 to 1.05 on the probe), so the grid would be indistinguishable.
- The CVaR term is over the diagonal E (the on-site energies) and the hopping enters as a
  linear expectation, per the brief. The alternative (CVaR over the spectrum of H in its
  eigenbasis) has no readout in the candidate basis and is not run.
- R2 is over the full p, not the tail (S23 L8's object). The tail-restricted weighted average
  is a fourth readout I did not add, because at J > 0 the mass outside the prefix is the point.
- No AMBER. No projection beyond the arms in section 7.

## ADDENDUM 0 (2026-09-14 19:15, coordinator steer received; still before any endpoint number)

1. The verdict basis is the BUILT CHAIN. Point-cloud numbers screen; every J arm at 0.7x MDE or
   better on the point cloud goes to the chain before it is called anything; J = 0 and the best
   J go to the chain regardless (section 7 stands).
2. No cosmetic variants. The SECONDARY graph scale sigma = median/2 (section 2, last bullet) is
   WITHDRAWN from the endpoint runs; sigma = median only. J = 3.0 is KEPT: it is not a finer
   grid but the one rung at which the exact ground state is delocalised over more than 75
   states (PR 277 to 335 on the probe), and the coordinator's question (i) is about the hopping
   GROUND STATE consumed by the deployed readout. R3 is kept as the readout form of the
   departure measurement (section 5), not as a readout tweak: R1 cannot see an off-diagonal
   term at fixed m by construction.
3. Order of work: code and tests, then F5 (trainability vs J, no RMSD) and the departure
   diagnostics, posted to the ledger first; then the endpoint run.
4. Any positive is answered in its own ledger entry with (a) the built-chain survival and the
   FAIL18 / 108-other split, (b) what is new against S26/S27 with the ledger line cited,
   (c) the three-way split: exact ground state (Hamiltonian quality) vs VQE state (optimisation
   quality) vs emitted structure. A gain that survives only with the eigensolver is not a
   quantum result and will be written as such. Nothing is built on a positive before lane D
   posts STANDS.

## ADDENDUM 1 -- B2, A HOPPING GRAPH WITH A SPREAD SPECTRUM (2026-09-14 21:40, written after S28-L21 (point cloud, intermediate) and lane D's S28-L11, before any B2 number and before the S28B built-chain verdict)

Brief: `s27/briefs/S28B2.md`. Conditions: (b) holds (S28-L11 confirms the near-rank-one
mechanism: lambda_2 / lambda_1 = 0.11 to 0.14, the rank-one part of A carries 97% of the
hop-only gradient variance, `s27/results/s28_B_rank1.json`); (a) is read on the point cloud
(S28-L21: the deployed readout within 0.013 A of J = 0 at every J, the other readouts worse
than production) and is confirmed or not by the built-chain verdict entry, which gates the
B2 ENDPOINT arms. The B2 TRAINABILITY measurement (no RMSD) runs before that gate.

The new angle, stated up front (contract rule 10): S28-L8b/S28-L11 found the hopping term's
gradient variance decays at -1.7 to -1.8 per qubit BECAUSE the Gaussian similarity graph is a
typicality projector (one eigenvalue 1.0, the next 0.11 to 0.14), and that a diagonal
observable with the same spectrum decays at the same rate. That is a statement about the
SPECTRUM. A graph with a spread spectrum (many eigenvalues of comparable size) is a different
Hamiltonian in the one respect the mechanism names. Nothing else (finer J, other sigma) is run.

Design, held fixed from the base prereg unless stated:
- Graph: the symmetric k-nearest-neighbour graph on the same pairwise CA-RMSD matrix
  (i ~ j if j is among i's k nearest or i among j's; binary weights), k in {5, 10};
  degree-normalised, D^-1/2 A_knn D^-1/2, then scaled to unit spectral norm (the normalised
  adjacency of a connected graph already has lambda_1 = 1; the rescale is the identity then and
  is kept so the code path is one). Zero diagonal, zero padding rows. k = 10 is the ARM,
  k = 5 its replication; no third k.
- Reported per target beside the Gaussian graph's: lambda_2 / lambda_1, the top eigenvector's
  participation ratio / dim, its overlap with the uniform state, the number of connected
  components, the degree range after symmetrisation (>= k), and corr(degree, E).
- Objective, circuit, settings, readouts (R1, R2, R3), seeds (0, 1), PERM and RAND controls
  (RAND: a random symmetric k-regular-in-expectation graph is NOT used; RAND is the same
  Sinkhorn degree-matched construction as the base prereg applied to the binary kNN degrees,
  which for a near-regular graph is a random graph of the same density), as in the base.
- J in {0.3, 1, 3} (the 0.1 rung dropped: invisible at every readout in S28-L21).
- F5-B2, the trainability measurement FIRST (it does not depend on RMSD): hop-only gradient
  variance at n = 4..9 (E the rank ladder over the DIS top-2^n, A_knn built on that sub-pool
  with the same k) and the full objective at J in {0, 0.3, 1, 3}, 120 draws, S27's 12
  trainability targets, median; the sign coherence / share of the same-sign bound the circuit
  collects at each J on the endpoint rows when they exist. FALSIFIER for "the spectrum was the
  mechanism": the kNN hop-only slope is SHALLOWER than -1.0 per qubit (Gaussian: -1.7 to -1.8)
  AND the circuit collects more than half of its same-sign bound at J = 1 (Gaussian: 1% at
  J = 1, 25 to 34% at J = 3). If the slope stays at or below -1.7 the mechanism reading was
  wrong or incomplete and is said so. Registered prior: the slope is shallower (a sparse
  graph's normalised adjacency has a spread spectrum, so the hopping term is no longer one
  squared overlap) but the endpoint is null or worse (see next).
- Endpoint (ONLY after the S28B built-chain verdict is posted and lane D has checked it): the
  built chain for J = 0 R1, each J's R1, and any arm at 0.7x MDE on the point cloud, paired
  against production (S28-L2(a)) and against J = 0, PERM and RAND, both seeds, FAIL18 / 108,
  the best J priced as an order statistic. Registered prior: WORSE or null. Reason: a kNN
  hopping ground state delocalises over the argmin's structural neighbours, which is the
  consistency mechanism S27 section 6 measured as harmful (+0.29 A for DIS+CONS) and which the
  Gaussian eigensolver at J = 1 turned into a re-selection of the DIS top-75 (S28-L21); a
  sharper graph makes the retained cluster tighter, not more native.
- Files: `s27/s28_B2_knn.py` (a graph builder reusing `s28_B_hop` for everything else),
  `tests/test_s28_B2.py`, results `s27/results/s28_B2_*.json|jsonl`, ledger entries `(date,
  B2)`, findings under a B2 heading in `s27/s28_B_FINDINGS.md`.
- Memory and cost: as the base (0.35 GB; the kNN graph is sparser but stored dense).

## ADDENDUM 2 -- B2 ENDPOINT SCOPE (2026-09-14 23:10, written after S28-L25 and S28-L29 (property entries), before any B2 endpoint number and before the S28B built-chain verdict)

F5-B2 as registered was not passed (first clause straddled, second clause failed, S28-L29).
The one cell the mechanism licenses for an endpoint look is the one where the circuit
collects its hopping: k = 10 (connected on every target measured; k = 5 is disconnected on
some targets and is NOT run), J = 3, both seeds, on the point cloud: VQE R1 and R3, GS R3
(labelled DEGENERATE on any target whose k = 10 graph has more than one component; the count
is written per target), the PERM control at the same cell, and J = 0 (already in
`s28_B_rows.jsonl`, bit-identical). No other k or J. Comparators: the same readout at J = 0
(F1), PRODUCTION (DIS top-75 uniform; S28-L2(a)), PERM (F3); both seeds; FAIL18 / 108. The
built chain only if a cell reaches 0.7x MDE on the point cloud against production or J = 0,
and then paired against production on the chain. Registered prior: WORSE or null (the kNN
ground state is a tight cluster of PR 53 around the argmin, the consistency mechanism). The
job runs only after the S28B built-chain verdict entry is posted and lane D has checked it.
