# S29 STATE (coordinator; the running written state the charter requires; updated as results arrive)

Last update: 2026-09-20 00:38 Pacific.

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


## Integration note 10 (2026-09-20 00:25, after S29-L16): THE IN-BAND EXPERIMENT IS NOVEL, AND ITS NULL IS WEAKER THAN IT LOOKS
Lane L's topic 6 supplies the design, the bound, the prior and one scoping caution I have adopted
as the sprint's position.
DESIGN, published and human-validated: the PIRM 2018 challenge built this construction in
transpose (fix distortion, rank by a NO-REFERENCE realism index validated by 35 raters who never
saw the ground truth) and found the index well correlated with human judgement ACROSS bands
(Spearman 0.83) and unreliable WITHIN them. That is this project's own +0.653 global / +0.091
in-band split reproduced independently in another field on this exact design. Expect a small
in-band effect and power for it; and because the tradeoff is steepest at the low-distortion end
where we operate, the band must be narrow, which is the design's central tension.
THE BOUND, exact: for jointly Gaussian (scorer, accuracy, realism) the in-band correlation is the
PARTIAL correlation rho_SY.R = (rho_SY - rho_SR rho_RY)/sqrt((1-rho_SR^2)(1-rho_RY^2)), constant
across the band, and its zero is exactly "the scorer's whole association with accuracy is
MEDIATED by realism". So the answer is predictable in closed form before any run, and computing
the three correlations IS the pre-registration. Report skill as a CURVE over band widths (global
correlation at wide, partial at thin): one pre-registered object instead of k bands.
THE PHYSICS, which makes it interpretable either way: banding on a statistic destroys that
statistic's own discriminating power by construction, so an in-band experiment measures exactly
what is ORTHOGONAL to realism -- a positive is direct evidence of the quantity H0 says the system
lacks; a null says the library carries nothing beyond realism.
THE CAUTION, ADOPTED AS THE SPRINT'S POSITION: every scorer in the S27 library was FITTED for the
between-band task (ANDIS says it in our own field: native recognition and decoy discrimination
"cannot be optimized simultaneously with the same parameter sets"). Therefore a null from the
current library does NOT close the in-band question; it says the library is the wrong instrument,
which is a weaker and different claim. Lane D states which claim it is making before running, and
lane L is reading topic 7 (training a within-group ranker; the free-energy class) so the
follow-up is designed before the result lands.
NOVELTY, stated because it bears on the report: no published QA evaluation conditions on a
native-free realism statistic before correlating with accuracy. The field's most careful
methodological paper (Hamelryck et al.) documents the confound across 139 of 149 decoy sets and
explicitly declines the matching fix. D's measurement is novel rather than derivative.
HONEST PRIOR on its value: the record's squared-skill law means a small partial correlation buys
nearly nothing in Angstroms, so a positive is a MECHANISM result first and a candidate operator
second.


## Integration note 11 (2026-09-20 00:27, after S29-L17): THE SPRINT'S THESIS, AND THE BUILD THAT TESTS IT
1. CONTRACTION IS AN EXACT VARIANCE IDENTITY, not Jensen loosely and not a posterior bias:
   d_ij(C)^2 = mean_k d_ij(W_k)^2 - s_ij^2 and Rg(C)^2 = mean_k Rg(W_k)^2 - Delta^2, where
   Delta^2 is EXACTLY S23 L9's idiosyncratic term. The operator cannot take the 68% common mode
   and cannot avoid paying the 32% in contraction. The fractional contraction falls from 22% at
   the bond to 6% at the envelope because the spread does not grow with separation while the
   distance does, which derives S23 L1's "averaging smooths".
   CORRECTION TO MY OWN PUBLISHED PROGRAMME (`s27/REPORT_S28.md` section 12 item 2 and the
   published page): "calibrate the posterior and re-read the meter" is DEAD as written. A width
   error does not move the median map (so not the minimiser, deriving S25 L2's null), a uniform
   over-confidence multiplies the metric by a constant and moves no minimiser at all, and
   calibration cannot touch the contraction because that is the pool's dispersion. The only
   live version is a SEPARATION-BAND RE-WEIGHTING with one or two parameters: the shipped
   objective over-weights mid-range pairs by about 2.7x against a calibrated one (z_sd 1.23 /
   2.05 / 1.87 / 1.28 across separations 2-2 / 4-5 / 6-8 / 9-15). Lane M or D runs T's
   12-target prediction on it; the artifact update at the sprint's end must carry this.
2. THE SPRINT'S THESIS, NOW BUILDABLE. T's section 4 gives the lift that fixes the flatness:
   CVaR over a STRUCTURAL observable, "tail then aggregate", F = CVaR - T H + lam f(R_alpha(p))
   with R_alpha the TAIL'S OWN coordinate average (the deployed readout as a function of p).
   The objective then SEES which candidates populate the tail -- precisely the freedom it is
   currently blind to. Derived and ready: the envelope gradient dR/dp_y = (W_y - W_x_q)/alpha,
   the same 2P parameter-shift cost as the deployed objective, convex cells indexed by the tail
   SET, non-smooth only where the scalar CVaR already is, and it is the MORE device-realisable
   lift (a quantile plus a mean structure, not the full 2^n distribution). It is NOT S28 lane A,
   which put the same term on the signed-amplitude readout.
   MY REGISTERED PREDICTION, on the record before lane B runs it: the flatness FALLS (the
   mechanism works) and the endpoint gets WORSE, because f is a marginal-class objective bound
   by theorem 2 and the shipped score's near-native ladder correlation is -0.182. If that is
   what happens, the sprint has demonstrated rather than asserted its thesis: THE FLATNESS WAS
   FIXABLE AND THE INFORMATION WAS THE BARRIER. If the endpoint improves, lane D attacks it the
   same hour and the sprint has its first movement.
3. The two questions now hang together: lane B's lift makes the objective see the set, and lane
   D's band experiment asks whether ANY scorer can order structures of equal realism. If D finds
   one, it is the f that lane B's lift should carry.


## Integration note 12 (2026-09-20 00:30, after lane M's three deliverables): THE ANCHOR IS CLASSICAL
1. THE HEADLINE, AND IT REFRAMES THE CHARTER'S CONSTRAINT. The 3.2126 A production anchor NEVER
   PASSES THROUGH THE QUANTUM STAGE. `core/pipeline.py:179` has `quantum: bool = False` as the
   production default (I verified the line independently), `s27/run_vqe_chain.py --chain` is the
   classical tie-safe top-75, `arm_vqe` runs only under `--vqe`, and the production cache record
   carries `quantum: null, n_top: 75`. So every S28 and S29 endpoint contrast "against
   production" is against a CLASSICAL pipeline, and the CVaR-VQE arm is a parallel arm whose
   S25 point-cloud number is 3.0580 against the classical 3.0483. This is consistent with, and
   explains, S28-L21: the CVaR tail equals the classical top-m prefix to 1e-13, so running the
   quantum stage reproduces the classical set. It is not a defect and nothing is retracted, but
   the final report must say it in the first paragraph: the charter's "CVaR-VQE remains the
   central component" is a requirement about the FUTURE architecture, because in the deployed
   present the quantum stage is switched off and would change nothing if switched on.
2. LEAVE-FOLD-OUT CONSTRAINS ONLY 9.5% OF THE DISTOGRAM'S TRAINING DATA: `fold_fragments`
   removes 0 to 15 of 6,003 fragments per fold, so 90.5% of every fold model's chains are shared
   across all five folds. With the 8x memorisation (check 9) this bounds how independent the
   five fold models are, and therefore how much a fold-clustered CI can protect against a
   corpus-level artefact. Lanes quoting fold CIs keep this in view; it does not invalidate them.
3. T = 0.5 HAS NO RECORDED CRITERION ANYWHERE in the repository. Neither does the 17-bin edge
   set (the outer centres 4.0 and 25.0 are invented), the soft-bin sigma 0.6, or the medoid
   frame (the one readout choice with no measurement, and the frame in which S23 L9's
   decomposition is defined). 17 of lane M's 31 convenience choices are untested.
4. ASSIGNED, one line of code and it bounds every production quantum arm: the ORACLE ceiling of
   the TOP-128 set (the prefix the quantum stage actually sees), which is absent from the record.
   Lane O takes it as a ladder rung.
5. F1's 12-target probe (NOT evidence for the instrument): PROD 3.3816, LOG 3.6161 (+0.2345,
   0.35x MDE), L2RISK 3.5632. The mechanism corroboration is the useful part and it is the THIRD
   independent one: LOG contracts LESS (bond 3.106 vs 2.987) and is worse; L2RISK contracts MORE
   (2.868) and is worse. CONTRACTION AND RMSD DO NOT TRACK ACROSS FUNCTIONALS. Lane M also
   declined my literal matched control (a monotone re-ranking) because through a top-m readout
   it is the IDENTITY by algebra, and replaced it with a zero-information member swap. That was
   the right call and I record it as my error.


## Integration note 13 (2026-09-20 00:32, after S29-L20 and S29-L21): TWO FAMILIES CLOSED AT THEIR ORACLE CEILINGS
1. RUNG 6 (H1, the typicality axis): the falsifier fires on both clauses and lane O reports the
   kill is stronger than the one I registered. H1 is closed.
2. RUNG 8 (lane T's one-parameter family, the signed-readout collapse of any centered or
   agreement-matrix Hamiltonian): T predicted the ORACLE ceiling of production + eta*PC1 would be
   under 0.15 A better than production. IT IS EXACTLY ZERO: the best GLOBAL eta over a 61-point
   grid is eta = 0.0 and the mean is 3.0483, production to the last digit. The deployable
   leave-fold-out arm is +0.0071 (worse). And the family is not degenerate: the pool's spread
   along PC1 is 1.068 A and PC1 carries 36.3% of the members' deviation variance, so this is an
   informative-LOOKING direction that is empty.
   THE MEASUREMENT THAT MATTERS IS THE SIGN. Best per-target eta with a free sign buys -0.4543
   (an ORDER STATISTIC over 61 values); forcing the sign positive buys -0.2142 and negative
   -0.2455. So most of the apparent per-target gain IS the freedom to choose the sign per
   target, which is exactly what lane T's theorem 2 says the marginals cannot supply. Theory
   predicted the family, predicted its collapse to one parameter, predicted the sign would be
   the binding quantity, and predicted the ceiling's magnitude; measurement confirmed all four
   at the floor.
   CONSEQUENCE: lane B's REDIRECTED build (the signed readout over a centered/agreement
   Hamiltonian) is dead at its ceiling and must not be run. Lane B's NEW build (the
   tail-then-aggregate lift, note 11) is a different object -- it is about the objective's
   flatness, not about PC1 -- and stays live.
3. Standing tally of what the sprint has closed by MEASUREMENT AT AN ORACLE CEILING rather than
   by a null endpoint run: the typicality axis (rung 6), the PC1 family (rung 8), the
   better-conditioned coupling matrix (derived, S29-L11), the non-commuting free-energy cell in
   the candidate-index encoding (derived, S29-L15), the QA import route (literature, S29-L1),
   a scorer that prefers the near-native answer (perception-distortion, S29-L12), and
   "average more or differently" (bounded at 0.008 A, S29-L8). Seven routes, none of which cost
   a 126-target endpoint run to close.


## Integration note 14 (2026-09-20 00:32, after S29-L19 and S29-L21): THE SPRINT'S UNIFIED FINDING IS THE SIGN
Three independent lines converged on the same quantity within one hour, and it is now the
sprint's thesis rather than a conjecture.
(a) RUNG 8 (measurement, lane O): along the pool's first shape mode the ORACLE per-target gain
    is -0.4543 with a FREE sign and only -0.2142 / -0.2455 with the sign forced positive /
    negative. Most of the apparent per-target signal IS the freedom to choose the sign.
(b) THEOREM 2 (derivation, lane T): no objective built from the posterior's marginals carries a
    term in the native's deviation from typical, so the marginals cannot supply that sign.
(c) S14 VIA LANE L (the record, topic 7): a linear 4,125-parameter pair potential already
    saturates the within-target ordering problem at 0.986 with an overfitting gap of -0.0005,
    and cross-target transfer is 0.600 against the 0.638 needed -- a transfer gap 770x the
    overfitting gap. No loss, architecture, capacity or equivariance touches transfer. S14's own
    conclusion, which the sprint has now re-derived from two other directions: the only open
    direction is "a conditioning signal, not a better objective".
THE UNIFIED STATEMENT: what the system lacks is a PER-TARGET SIGN (equivalently, a conditioning
signal supplied at inference), not a better objective, not a better ranker, not a better
Hamiltonian, not a better readout. Every S29 route that died this sprint died for that one
reason, and the routes still running are each a test of whether some quantity can supply it.
ALSO, AND IT ARRIVED BEFORE THE EXPERIMENT RATHER THAN AFTER: lane L found that the in-band
ordering axis IS compactness (S14: per-target in-band skill correlates +0.909 with the native's
z-scored Rg and +0.951 with rho(contacts, RMSD)) and that every realism statistic in the S27
library is compactness-like. So banding on a compactness-loaded statistic would remove the very
axis the experiment is looking for and make a null self-fulfilling. Lane D must report
rho(R, Rg) for its band statistic beside the result, and re-run S28-L48's 0.960 pairwise learner
INSIDE the band as a validity check (it should collapse toward chance if the band is real).
This is the second time this sprint that the literature lane has saved an experiment from being
uninterpretable before it ran.


## Integration note 15 (2026-09-20 00:38, after S29-L22): THE DISTORTION IS NOT A CONTRACTION
Lane P measured the averaging distortion as a function of sequence separation on all 126 clouds
and it is NOT a contraction. The ratio of the cloud's distances to the native's is monotone in
|i-j|: 0.773 at the virtual bond, CROSSING 1.00 near |i-j| = 8, and 1.10 by ratio of means
(1.05 by median ratio) at |i-j| = 13. The NATIVE-FREE profile against the posterior's own median
map has the same shape and the same crossing, so this is not an ORACLE-only statement.
CONSEQUENCES, and they tidy up three standing items.
(a) S23 L1's two numbers (bond 22% short, envelope 6% short) are the TWO ENDS OF ONE CURVE and
    the middle is where the sign changes. The project has been describing a shape distortion as
    a contraction for several sprints, including in my own integration notes 3 and 11 and in the
    published S28 page; the correct statement is "averaging shortens local geometry and
    lengthens long-range geometry", which is lane L's Jensen mechanism and lane T's variance
    identity seen per separation: averaging shrinks a distance by more when the pool disagrees
    more about it RELATIVE to that distance's size, and the pool disagrees most about local
    geometry.
(b) THE BOND-RESCALE ARM IS PREDICTED HARMFUL BY MEASUREMENT RATHER THAN BY ARGUMENT: setting
    the bond right requires g = 1.29, which inflates every separation beyond 8 by about 29% on
    top of distances that are ALREADY too long. The best single scalar must sit near the middle
    of the curve, about 1.05 to 1.08 (the lane's SPAN and ISO factors are 1.102 and 1.077).
(c) CORRECTING THE PROFILE IS REFUTED, AND AN ORACLE-FITTED PROFILE IS NO BETTER -- another
    closure at an ORACLE ceiling rather than by a null endpoint run. This bears on lane T's
    "separation-band re-weighting" as the only live form of the calibration item: P corrected
    the CLOUD's separation profile and it did not help, which is a different operation from
    re-weighting the OBJECTIVE's pairs but close enough that T should state whether its
    suggestion survives P's measurement before anyone builds it.
Also from the same entry: the projection price is NOT a flat toll (mean +0.1643, sd 0.2014,
positive on 111/126) and correlates +0.604 with the native-free contraction factor, +0.481 after
partialling out the cloud's own RMSD and the length. And S23 L2/L3 reproduce exactly on these
clouds (42.9% of targets want expansion against S23's 42.1%), with every native-free factor
inside the +-0.11 band S23 measured and the only one outside it carrying the WRONG sign.
Lane P also withdrew its own first statistic (a mean of per-pair ratios, 1.37, dominated by small
denominators) before using it. The shape conclusion is identical under all three statistics.

## Closed in S29
- The PC1 / signed-readout family (S29-L21): ORACLE best global eta is exactly 0.000; the
  per-target gain is the SIGN, which the marginals cannot supply.
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
