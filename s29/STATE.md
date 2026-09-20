# S29 STATE (coordinator; the running written state the charter requires; updated as results arrive)

Last update: 2026-09-20 00:24 Pacific.

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


## Integration note 4 (2026-09-20 00:06, after S29-L7 and S29-L9): TWO GUARDS AND A CONSTRAINT
1. THE THEOREM (S29-L7, lane T) is the sprint's spine result so far: an objective is locally
   informative iff its per-pair force coefficients are negatively correlated with production's
   own signed error against the native; under the record's error model the expectation carries
   NO term in the native's deviation from typical, is second order, and vanishes when the
   posterior's median map and production deviate from typical alike. Corollaries: 45% of pair
   space is in ker(Jc^T) and invisible to any marginal objective's gradient; non-separability
   buys nothing; the sign is -sign(beta - 1) with beta the over-confidence, which predicts the
   measured -0.034 and the -0.143 on FAIL18. This is H0 DERIVED, not merely observed.
2. GUARD 1, now contract rule 20: the cosine is gameable by shrinking the target map toward
   typicality (positive cosine, zero information, worse structure). Every cosine gain must ship
   with its implied shrink, the native percentile, and the emitted bond and Rg. Lane D is
   adding the shrink to the meter and running the theorem's three predictions as a measured
   guard.
3. GUARD 2, from lane M's audit (S29-L9): the distogram MEMORISES its training peptides by 8x
   (in-fold vs out-of-fold NLL delta +2.075, SE 0.152, 4.88x MDE, correct sign on 5/5 folds;
   372,881 parameters against 787 peptides + 6,003 fragments with dropout deliberately 0).
   Consequence for every S29 lane: any diagnostic computed on the corpus (calibration,
   sharpness, MAE, a fitted residual) describes MEMORISATION and means nothing about the
   deployed model until it is recomputed out of fold. This bounds the one route lane T's
   theorem leaves open (a learned residual that sees the native's deviation from typical): it
   must be trained and read strictly leave-fold-out, and its in-sample fit is worthless.
   The harness itself is sound (nine checks, all pass; one declared non-bit-exactness in the
   s12 score cache that affects only operators reading the order below the top-75 cut).
4. The sprint's remaining live routes, in order of what would move the endpoint:
   B (a non-degenerate off-diagonal Hamiltonian: does the quantum state select a better set),
   P (contraction removed at the projection), F1 (contraction never created, at the selection
   functional), X (the configuration state space where the scorer scores what it generated).
   All four are gated, pre-registered, and attack mechanisms the theorem names.


## Integration note 5 (2026-09-20 00:08, after S29-L11): I WAS WRONG ABOUT WHY LANE B WAS WORTH SPAWNING
Lane T's theory section 3 corrects my integration note 2, and the correction is mine to own.
DERIVED LAW: for any unit-spectral-norm observable, Var[dF/dtheta] = r_stable/D^2 and nothing
else (r_stable = ||A||_F^2 / ||A||_2^2), which predicts S28's Gaussian graph to 7%, S28-B2's kNN
to 46x and lane D's J* = 85.7 as 88, with NO free parameter. Consequences:
(a) Note 2's trainability half is WRONG. Centering removes the lambda_2/lambda_1 degeneracy
    (0.138 -> 0.465 -> 0.634) but makes the gradient decay WORSE (-2.305 for A_c and -1.900 for
    G against A's -1.830), because r_stable rises only 1.04 to 1.6. "The spectrum is no longer
    degenerate" is the wrong justification for a build, and I gave it.
(b) DESIGN RULE, general: an off-diagonal term is gradient-visible at the deployed width iff its
    STABLE RANK grows with the register (parity at n = 9 would need r_stable ~ 8000). No dense
    kernel can, centered or not; a k-regular graph reaches M/k; ANY Gram matrix of structural
    deviations is capped by rank at r_stable <= 3N_res - 6 <= 42 here. This closes the whole
    "make the coupling matrix better conditioned" family by derivation, not by one more run.
(c) WHAT SURVIVES is the MEANING of the ground state, not its trainability: <v|G|v> =
    |sum_i v_i delta_i|^2 / N_res, so the top eigenvector is the signed combination whose
    deviations from the pool mean add to the largest displacement -- the pool's principal
    contrast, positive on one pole and negative on the other. Every deployed readout is a
    function of p = psi^2 and is blind to that sign, so the contrast SELF-CANCELS and the
    emitted structure is the pool mean again. A signed readout (S28 lane A's, refuted under the
    shipped objective but not as a readout) is required, and under it the family collapses to
    ONE parameter: production +- eta PC1(pool). Everything then rests on the SIGN, which theory
    section 2 says the marginals cannot supply. Lane B is redirected accordingly; lane O is
    measuring that family's ORACLE ceiling as rung 8 (T predicts under 0.15 A better than
    production; above 0.30 A would mean the build is worth much more than the theory says).
This is the charter's method working: theory before the build killed a justification I had
already acted on, and replaced it with a sharper, cheaper, falsifiable claim.


## Integration note 6 (2026-09-20 00:09, after lane D's first turn): TWO CORRECTIONS THAT BIND EVERY LATER NUMBER
1. THE LADDER'S SIGN IS LADDER-DEPENDENT, AND THE CHARTER'S RUNGS ARE THE WEAKER TEST. For the
   shipped cost, ladder rho is -0.182 (CA) / -0.402 (chain) on S28's five rungs, but +0.260 (CA)
   / -0.092 (chain) on the charter's six, and +0.118 / -0.236 on all nine. The cost orders the
   BULK correctly and anti-orders the NEAR-NATIVE half; a ladder weighted toward the bulk
   therefore reports a positive rho for a cost that is adversarial exactly where it matters
   (the memory `decoy-bank-not-a-pool-proxy` is this shape). STANDING RULE for the rest of the
   sprint and for the report: a candidate cost must beat the shipped cost ON S28's RUNGS; a
   positive charter rho beside a negative S28 rho is the expected failure mode and is quoted as
   such. Every "-0.40" in any S29 text names its ladder.
2. MY CONTRACTION STORY IS HALF WRONG, AND THE HALF THAT IS WRONG IS THE GRADIENT. Lane D's
   shrink signature (contract rule 20, now mechanical in the meter) measures what a 0.3 A
   descent along -grad f does to the emitted geometry: for the SHIPPED cost the descent EXPANDS
   (bond x1.0438, Rg x1.0248; only 18/126 contract). So the shipped cost's local blindness is
   NOT a contraction artefact at the gradient level, though contraction remains true of its
   RANKING (production sits at the 12.6th percentile of its own pool, S28-L36). Integration
   note 3's framing survives for the ranking and for F1, and dies for any gradient-level story.
   P1 (the projection rescale) and F1 (the functional swap) are unaffected as experiments, but
   their entries must not claim a gradient mechanism.
3. H1's GRAVE, ONE DETAIL WORTH KEEPING: lane D reports the per-target ORACLE argmin sits at
   the grid's LEFT edge (t = -1, i.e. AT the blind average) on 29/126 targets. Where the axis
   has any ORACLE content at all, it points TOWARD the sequence-blind answer on a quarter of
   targets, the direction opposite to H1's hypothesis. Lane O states it in the rung-6 entry.
4. Lane D's own two defects, recorded because the sprint's standard applies to the Adversary:
   an all-NaN cosine axis crashed the renderer after a completed 126-target run, and a
   partially defined cosine was printed as a measurement. Both fixed with regression tests.


## Integration note 7 (2026-09-20 00:14, after S29-L13): THE CONDITION FOR NON-CLASSICALITY, AND A RECORD CORRECTION
1. RECORD CORRECTION (emphasis, not fact). The project's "set-equality theorem" (the CVaR tail's
   support is a prefix of the energy order; S25, S28-L21) is equation (12) of Barkoutsos et al.,
   Quantum 4:256 (2020): CVaR is DEFINED on sorted samples and the paper scopes itself to
   diagonal Hamiltonians. It is a definition restated, not a discovery. Every future citation of
   it in this project, including the S28 report and its published page, cites Barkoutsos eq (12)
   beside it. The S28 measurement (the identity holds to 1e-13 at every J on a NON-diagonal H,
   where the definition alone does not guarantee it) stands as a measurement.
2. THE CONDITION FOR A GENUINELY NON-CLASSICAL FORMULATION, which is the charter's "real quantum
   result" target stated precisely for the first time in this project: BOTH
   (C1) the Hamiltonian's terms do not commute, so the eigenbasis is not the computational
        basis (S28's hopping H satisfied this, with a degenerate off-diagonal), AND
   (C2) the prepared object is NOT an eigenvector, so an eigensolver is not the classical
        counterpart either -- satisfied by a thermal/Gibbs state, by a state whose role is to be
        a sampling distribution, or by a free-energy objective.
   THE PROJECT HAS NEVER SATISFIED BOTH. Everything diagonal is a sort; S28's non-diagonal arm
   targeted a ground state, so its counterpart was an eigensolver, and it tied one. Caveat from
   lane L: for a DIAGONAL H the Gibbs state is a classical Boltzmann distribution over
   candidates and S21 enumerated the latent exhaustively on 75/126, so (C2) alone buys nothing.
   The untested cell is (C1) AND (C2) together: a free-energy / Gibbs objective over a
   NON-COMMUTING Hamiltonian. Lane T is asked whether that cell can contain anything measurable
   here before any lane builds it.
3. A CANDIDATE MECHANISM FOR THE SPRINT'S CENTRAL PUZZLE, from the same paper and new to this
   record: for any theta* whose state has overlap rho with the best candidate, theta* is a
   GLOBAL minimum of CVaR_alpha for alpha <= rho. The global-minimiser set is therefore
   {theta : overlap with the best candidate >= alpha}: large and flat, and the objective is
   indifferent to exactly the freedom an averaging readout consumes (which OTHER candidates
   populate the tail). That is a candidate explanation for "the optimiser reaches the optimum
   on 126/126 and the emitted structure does not move" (S28-L18b, S28-L26b). Lane T checks it;
   lane D adds the clause that any accuracy change attributed to CVaR optimisation must be shown
   not to be a tie-break inside that flat set.
4. THE WARNING WE INHERIT: Cerezo et al., Nat Commun 16:7907 (2025), argue that provable absence
   of barren plateaus often implies classical simulability. The project's one genuine quantum
   positive (the optimiser trains, no plateau at any measured width) sits in that regime. It is
   not retracted; it is scoped, and rule 9 already forbids the reading that would be wrong.


## Integration note 8 (2026-09-20 00:17, after S29-L14): THE COMPARISON SPLITS, AND IT VALIDATES THE SPRINT'S TARGET
Lane L's topic 5: there is NO published ceiling for native-free peptide prediction at 9 to 16
residues; the field reports method scores on small curated sets (PEP-FOLD 2.6 A on 25 NMR
peptides, APPTEST 1.96 A on 42, AF2 best-of-5 by class 2.2 to 4.5 A, MD folds 10 to 20-mers at
1e5 to 1e6 CPU-hours). None is like-for-like with this instrument on three counts: composition
(curated NMR peptides with regular secondary structure, versus 126 identity-clustered PDB targets
whose hard stratum is 56% steric-zipper amyloid and lasso peptides that no linear-window
retrieval can represent), reporting (best-of-N there, one deployable answer here), and regime
(every published method selects inside an ensemble ITS OWN energy generated, which S29-L1 showed
is the only regime where native-free selection works at this length; we rank 500 real windows
from other proteins with an independently constructed objective).
THE USEFUL OUTPUT, and it is the report's framing: the comparison SPLITS.
  GENERATION: our ORACLE ceilings, 2.31 A (top-75) and 1.71 A (pool best), sit INSIDE or below
  the published band. Generation is not this project's problem, and the charter's 2.5 A target
  is inside the pool.
  SELECTION: the 0.9 to 1.5 A between 3.2126 and those ceilings is the ENTIRE gap, and the
  field's own best selector has no in-band skill at this length either.
Also recorded: reporting 3.21 against 1.96 without those three caveats would be misleading in
the project's own disfavour, and lane L says so explicitly. And if S29 produces a MEASURED bound
on native-free selection for 9 to 16-mers with controls, that is a contribution to the field and
not only to the project -- which is the shape the sprint's H0 outcome would take.

## Resource note (2026-09-20 00:17)
CPU 94%, the band the user asked for, with 8 governed jobs on 8 cores (lane O's chain ladder in
4 shards plus PC1, lane P's probe, lane X's probe, lane D's checks). RAM 71.6%: this work is
compute-bound, not memory-bound, and manufacturing memory pressure to reach 94% would risk the
real jobs for a number. If a memory-heavy step becomes scientifically justified (an ESM re-embed,
a full 500x500x126 tensor), RAM rises then.


## Integration note 9 (2026-09-20 00:24, after S29-L15): THE SPRINT HAS CONVERGED, AND THE MECHANISM IS DERIVED
THE MECHANISM FOR THE CENTRAL PUZZLE, complete and derived (lane T, Q1). By the envelope theorem
the deployed CVaR is EXACTLY constant along 437 of 511 simplex directions at the realised tail
(85.5%); the entropy term resolves those by flattening, so the exact optimum is uniform ABOVE
the VaR and exponentially enhanced below it, with p* a function of (alpha, T) ALONE because
E = zrank is the same rank ladder on every target to 1.18% of range. The readout consumes only
the tail SET; the set is a prefix fixed by the ordering; the one scalar left is where the prefix
cuts. THE DEPLOYED CVaR-VQE IS EQUIVALENT AT THE ENDPOINT TO CHOOSING ONE NUMBER m -- and the
m-ladder has been priced three times. That is "the optimiser reduces the objective on 126/126
and the structure does not move", derived rather than observed, and it supersedes the
flat-minimiser-set story from lane L (Barkoutsos's overlap condition has probability e^-46 here;
the envelope argument needs no overlap assumption and holds at every point).
THE CONTROL IT IMPLIES (now lane D's M6, and the sprint's sharpest): replace the entire quantum
stage by the target-independent rank-weight profile p*(alpha, T), no circuit, no optimiser, no
per-target computation, and run it to the built chain. T predicts agreement with the deployed
arm on >= 120/126. This is the charter's "removing the quantum stage must degrade the result"
at its sharpest and the report carries it whatever else happens.
THE ARCHITECTURE THE THEORY ENDORSES, and three lanes now agree on it. Lane L: non-classicality
needs non-commuting terms AND a non-eigenvector target. Lane T Q2: DERIVED NO for that cell in
the candidate-index encoding (the Duhamel first-order term vanishes for zero-diagonal couplings,
so the thermal state's leading quantum content is a classical reweighting by squared-similarity
degree -- exactly what S28's degree-matched control held fixed), and the only operator class that
escapes the trainability obstruction is a sum of LOCAL PAULI terms: for a transverse field
r_stable = D/n, the slope is exactly -1 per qubit and J* is 11.8 instead of 90. A local mixer is
meaningful only where basis states have local structure, i.e. in a CONFIGURATION-SPACE encoding.
SMALLEST QUALIFYING FORMULATION: basis state = per-residue configuration assignment; H = 1- and
2-body posterior terms + Gamma sum_q X_q; a thermal / free-energy target, never an eigenvector;
readout = the CVaR tail's coordinate average. That is lane X's encoding plus a mixer, and lane X
is redirected to build exactly it, running T's falsifier ladder cheapest first: TV > 0.45 on the
sampled distribution (below it the readout provably cannot resolve the difference, S25 L15),
then the correct classical counterpart (a thermal sampler or SA over the same space at matched
evaluations, NOT an eigensolver), then the endpoint.
HONEST CEILING ON IT, stated now so no entry overclaims: by theorem 2 this cell creates no
information about the native's deviation from typical, so its upside is the charter's
"classically irreproducible contribution at unchanged RMSD" unless the configuration space's own
posterior carries more than the pool's marginals -- lane X's premise, to be measured by the
meter, not asserted.

## Closed in S29
- The non-commuting free-energy cell IN THE CANDIDATE-INDEX ENCODING (S29-L15 Q2, derived):
  the thermal state's leading quantum content is a classical degree reweighting.
- The "better-conditioned coupling matrix" family (S29-L11, derived): gradient visibility needs
  stable rank growing with the register; no dense kernel and no Gram matrix of deviations can.
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
