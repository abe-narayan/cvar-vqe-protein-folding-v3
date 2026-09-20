# S29 STATE (coordinator; the running written state the charter requires; updated as results arrive)

Last update: 2026-09-20 00:02 Pacific.

## Leading hypothesis (H1: the typicality axis)
The missing information is not in any scorer; it is in the SIGN of the pool's systematic error.
Every predictor in the record emits "a typical peptide of that length" (S18/S19); the coherent
error is the gap between typical and this native; the pool's error is 68% common-mode (S23 L9);
sequence conditioning moves the answer partway from typical toward the native (blind pipeline
3.989 vs shipped 3.213, S19). If the conditioned answer lies BETWEEN the typical answer and the
native, then the native-free direction (conditioned minus blind) points toward the native, and
extrapolating along it with one leave-fold-out scalar is a deployable operator that uses the
common-mode error instead of averaging it in. A CVaR-VQE whose Hamiltonian's off-diagonal term
is the CENTERED (zero-Perron-mode) agreement of candidates' deviation-from-typical selects a
coherent set along that axis; its readout is the tail average extrapolated along the axis.
Attacks findings 3, 4, 5, 8, 11 (the anti-correlation is systematic, hence informative; the
common-mode error is the signal, not the noise). Accepts 1, 2, 6, 7, 9 as binding.

## Strongest objection
S24 L2/L3 (`s24/LEDGER.md` lines 115 to 160): the unselected blind-library source has a bias
cosine of 0.647 with the incumbent (31% angularly independent) but q = 1.231 (worse), so its
projection on the shared axis is about 0.8 of the incumbent's: the blind answer is NEARER the
native along the shared axis and much worse orthogonally. Then u = (conditioned minus blind)
has a +0.2 component ALONG the shared error and a large component that is minus the blind's
orthogonal error, and extrapolating along u increases the shared error. Mixtures under the same
score sit on a straight line (cos 0.943 for a score-selected retrieval-free pool; the union is
worth +0.002). H1's sign is therefore expected WRONG by the record's own geometry; the probe is
kept because it is 20 minutes, decisive, and the extrapolation (negative weight) was never run.
If it dies, the leading hypothesis becomes H0 below.
S16 (native-free error-direction steering: "every arm chose do nothing") and S24 L13 (the
prior's slope holds only along the native's own direction; a real operator at cos 0.5 buys
+0.024 A). If the ORACLE cosine between (conditioned minus blind) and (native minus conditioned)
is near zero, H1 is dead before any build. Also: the blind pool's bias may be the same
direction as the conditioned pool's, in which case the contrast is zero.

## What falsifies it (pre-registered, lane O, first measurement of the sprint)
ORACLE cosine between u = (shipped average minus blind-pool average) and v = (native minus
shipped average), rigid body removed, per target: falsified if mean cos is inside the random
reference (|cos| ~ 0.14 at 3n-6 dof) with the fold CI including it, or if the leave-fold-out
best global step along u does not clear 0.7x MDE on the point cloud. Prior: cos 0.2 to 0.4 on
the 108, higher on FAIL18 (where conditioning was harmful the sign may reverse: that is the
regime question and the random-18 null decides it).

## H0 (the fallback that is probably the truth, and the sprint's ceiling measurement)
The ceiling is imposed by the absence, anywhere in the system, of per-target information
orthogonal to the typicality axis (the "typical peptide of this length and sequence"). S29's
job under H0 is (a) to measure it decisively across every recognisable operator class (the
ORACLE ladder, lane O), (b) to search the literature and the data for any native-free source of
such information (lanes L, T, M), and (c) to build the CVaR-VQE formulation in which the quantum
stage consumes the pool JOINTLY (a correlated state over mutually compatible hypotheses, an
interacting Hamiltonian with a spread spectrum by construction) so that, if such information is
found, the quantum stage is where it acts and where its removal degrades the result. Under H0 a
clean "the quantum stage contributes a specific classically irreproducible quantity" at
unchanged RMSD is the reachable major result (charter section 17).

## Divergent direction (H2, lane X): the state space
The candidate pool may be the wrong state space. A quantum state over per-residue structural
hypotheses (fragment / torsion-bin assignment) with an interacting Hamiltonian (posterior
pairwise log-likelihood as 2-body terms, torsion prior as 1-body, a transverse mixer as the
off-diagonal), CVaR over the sampled configurations' posterior NLL, readout = coordinate average
over the CVaR tail's configurations (a coherent basin, not a pool average). New angle vs
S13/S15/S21: the readout is the tail ENSEMBLE, not the argmin, and the posterior is consumed
jointly, not as marginals. Prior: ties the pool (S21: search saturates).

## Standing measurements before any build
- The ORACLE ceiling ladder by recognisable-operator class (lane O): best basin average for
  k clusters; best sparse convex combination with s members; best m-subset; best single; the
  typicality-axis cosine and its leave-fold-out step. This says what class of operator could
  reach 2.5 A at all.
- The data-path map and the convenience-choice list (lane M); the evaluation-harness audit
  (lane D).
- The cost-RMSD meter (lane D) for every objective proposed.
- Theory (lane T): why the L1 Bayes risk of a 2x over-confident posterior contracts; the class
  of objective that is locally informative given only marginals; the spectral condition for a
  non-degenerate compatibility Hamiltonian (centering); the reachable set of the deployed ansatz.
- Literature (lane L, permanent): the two areas of charter section 6, with "what information
  does it contain that we do not" per paper.

## Agents (wave 1, 2026-09-19 23:05)
| lane | role | brief |
|---|---|---|
| L | literature (permanent) | s29/briefs/S29L.md |
| M | data-path map, convenience choices, harness audit | s29/briefs/S29M.md |
| T | theory and derivations | s29/briefs/S29T.md |
| O | ORACLE ceiling ladder and the typicality-axis probe | s29/briefs/S29O.md |
| D | adversary (permanent, rotating), cost meter, suite | s29/briefs/S29D.md |
| X | divergent: configuration-space CVaR-VQE | s29/briefs/S29X.md |


## Integration note 1 (2026-09-19 23:55, after S29-L1)
Lane L's topic 1 closes the "import a QA method" route from outside: no published native-free
quality method has been trained or benchmarked below 40 to 50 residues (ProQ3 filters under 50,
VoroMQA trains above 99, DeepAccNet 50 to 300), the four signal classes they reduce to are each
absent at this length, already ours, or measured dead here, and AlphaFold2's own pLDDT has no
within-target ranking skill on 588 peptides of 10 to 40 aa (McDonald 2023: rank-1 costs 0.2 to
1.1 A against best-of-5). Charter finding 8 is the field's position at this length, not a defect
of the S27 library. Consequence for the sprint: H0 rises; the two published native-free
selectors that work at 9 to 25 aa (PEP-FOLD 2.6 A, APPTEST 1.96 A) both rank an ensemble their
OWN energy generated, which is a convergence diagnostic over a self-consistent set, not a
ranker over a foreign pool. That is a structural hint for the architecture: a selector is only
known to work when it scores what it generated. Lane X's configuration space is the only S29
direction with that property (the Hamiltonian that scores the configurations is the Hamiltonian
that generates them), which raises X's priority from "divergent" to "the second real candidate".

## Wave-2 candidate probes (coordinator's list; not yet assigned)
P1 THE PROJECTION PRICE. Production pays +0.159 A to go from cloud (3.0483) to chain (3.2126),
and S28-L39 measured that a de-contracted cloud pays only +0.026. The cloud is contracted 22%
(mean virtual bond 2.96 vs 3.77). S23 L6 closed the RMSD-optimal global scale as unreachable in
principle, but that is a different quantity from a GEOMETRIC consistency rescale (make the
cloud's mean CA-CA bond equal the ideal 3.8 A), which is fully determined and native-free. The
probe: rescale the production cloud to ideal bond length, project, compare on the built chain;
then the same with the per-target scale that equalises the bond (not the RMSD). Cheap (126 x 3 s),
decisive, and it attacks the one stage of the pipeline no S28 lane touched. Falsifier and the
S23 L6 distinction go in its prereg.
P2 THE TWO-SOURCE CONTRAST at the level lane O is already measuring (rung 6).
P3 The centered compatibility Hamiltonian's spectrum (lane T section 3c) if T finds the
gradient-variance decay flattens.


## Integration note 2 (2026-09-20 00:00): H1 IS DEAD; THE S28 HAMILTONIAN CLOSURE IS SCOPED; EIGHT LANES
1. H1 (the typicality axis) is FALSIFIED on its own ORACLE ceiling, before any deployable arm:
   in `s29/results/s29_O_lfo.json` the ORACLE global step t along u = (shipped average minus the
   blind average) is EXACTLY 0.0 on every fold and the ORACLE global mean equals production to
   the last digit (3.048338). Even with the native, one global step along that axis buys nothing;
   the per-target step (2.742) is an order statistic over 31 choices. This reproduces S16's
   "every arm chose do nothing" and confirms the S24 L3 geometry the objection cited. The
   leading hypothesis is now H0. (Lane O posts the entry with the controls lane D required.)
2. CHARTER FINDING 7 IS SCOPED, AND IT IS THE SPRINT'S FIRST REAL OPENING. Lane T measured, on
   12 real pools (`s29/results/s29_T_spectra_rows.jsonl`, 72 cells): lambda_2/lambda_1 at n = 9
   is 0.138 for the raw Gaussian similarity A (S28's matrix, near rank one, the mechanism that
   killed the hopping term), 0.465 for the double-centered A_c, and 0.634 for the signed
   agreement matrix G = D D^T of deviations from the pool mean. The S28 closure was about ONE
   DEGENERATE SIMILARITY MEASURE, not about off-diagonal Hamiltonians. G is also the charter's
   "Hamiltonian encoding the disagreement between the prior and the pool" and it removes the
   68% common mode by construction instead of averaging it in. Lane B is spawned to build it.
3. Eight lanes now (the charter's maximum): L literature, M map and audit, T theory, O ladder,
   D adversary and meter, X configuration space, P projection price, B compatibility Hamiltonian.
   Roles kept distinct: B and X are different state spaces, not variants; D attacks both; L
   keeps reading; O finishes the ceiling ladder that tells us what any of them could reach.


## Integration note 3 (2026-09-20 00:02, after S29-L6): THE FUNCTIONAL, NOT THE INFORMATION
Lane D metered lane X's pair log-score and found the first cost in the record that is not
ANTI-informative on the near-native ladder: ladder rho +0.018 [-0.069, +0.123] against the
shipped cost's -0.182 [-0.308, -0.053], a paired difference of +0.200 at 1.45x MDE with 5/5
folds and power 0.98 (`s29/results/s29_D_cost_audit_X_cost_nll_ca.json`). The mechanism lane D
states is the sprint's first structural insight: the shipped cost is a BAYES RISK (expected L1
distance error under a posterior that is about 2x over-confident, S25 L2), whose minimiser is a
CONTRACTED structure; the pair log-score is a PROPER SCORING RULE of the SAME posterior, whose
minimiser is not driven to contract. Same information, different functional, +0.200 of ladder
rho. Lane D also scoped it correctly: the native still sits at the 37.8th percentile of its own
pool (36.8th under the shipped cost, difference NOT MEASURED), so finding 8 survives and this
cost can only stop making things worse, not recognise; and the memory `better-matrix-worse-
ranking` is on the record for exactly this shape.
CONSEQUENCE. Two independent stages now attack ONE mechanism, contraction under an
over-confident posterior:
  P1 (lane P, running): the contraction is removed downstream, at the projection.
  F1 (queued, lane M on release): the contraction is never created, by swapping ONLY the
     functional at the selection stage of the SHIPPED pipeline (L1 Bayes risk -> the log score
     of the same 17-bin posterior, same pool, same K, same readout, same projection), as a
     deployable leave-fold-out arm on the built chain with the matched control (a monotone
     re-ranking of the shipped score, which changes the functional's SHAPE but not its
     information) and the S24-L3 parallel-bias check. If both fail, contraction is closed as a
     lever and the ceiling is information, not functional form. If either moves the chain, it is
     the sprint's first movement and lane D attacks it the same hour.

## Closed in S29
- H1, the typicality axis (S29-L<O's entry>): the ORACLE global step is exactly zero; the axis
  carries no deployable signal, and the per-target step is an order statistic.
- Importing a native-free QA method from the literature (S29-L1, lane L): no method exists at
  this length; the four signal classes are absent, ours, or measured dead. The route is closed
  from the outside as well as from the inside (S28-L48).

## Comparison count (multiplicity)
0 endpoint comparisons run. (Every lane reports its count per entry; the coordinator sums here.)

## Budget plan
Reading, theory, literature: ~35% of the sprint. Probes (12 targets, pre-registered): ~35%.
Full-instrument runs: ~30%. Revised consciously at each STATE update.
