# S29 STATE (coordinator; the running written state the charter requires; updated as results arrive)

Last update: 2026-09-19 23:05 Pacific (sprint start).

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

## Closed in S29
(none yet)

## Comparison count (multiplicity)
0 endpoint comparisons run. (Every lane reports its count per entry; the coordinator sums here.)

## Budget plan
Reading, theory, literature: ~35% of the sprint. Probes (12 targets, pre-registered): ~35%.
Full-instrument runs: ~30%. Revised consciously at each STATE update.
