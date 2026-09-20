# S29 STATE (coordinator; the running written state the charter requires; updated as results arrive)

Last update: 2026-09-20 01:32 Pacific.

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



## THE SPRINT'S HEADLINE RESULT (2026-09-20 00:44, S29-L23; provisional until lane D's attack)
EVERY NATIVE-FREE OPERATOR IS A DISPLACEMENT AND ITS VALUE IS ONE COSINE.
    RMSD_achievable = RMSD_prod * sqrt(1 - rho_max^2)
with rho the cosine between the displacement and the direction to the native. Selection,
re-weighting, a gradient step, a basin average, a signed readout, a projection and a quantum
tail are all displacements, so every operator class -- and every rung of lane O's ladder --
collapses onto one number. Inverted from 3.2126: 3.00 A needs rho 0.358; 2.50 A needs 0.628;
2.31 needs 0.695; 1.71 needs 0.847. The random-shape-field reference is 0.140, so the charter's
2.5 A is 4.5 random directions' worth of alignment.
MEASURED rho FOR EVERY FIELD THE PROJECT HAS BUILT (all ORACLE): the shipped objective's descent
-0.034; every other S27 channel with a gradient +0.034 to -0.006; the typicality axis -0.058;
PC1 per-target |rho| 0.37 to 0.39 but with the sign correct on 52% of targets, hence a SIGNED rho
of 0.015 and a leave-fold-out arm +0.0071 A worse than production.
THE BOUND: every native-free field ever built here (|rho| <= 0.04) gives >= 3.210 A, a 0.003 A
gain. A field as strong as a RANDOM one with the sign right everywhere (0.140) gives >= 3.181 A.
The best structured field known, PC1, with a PERFECT ORACLE per-target sign (0.37) gives
>= 2.98 A. CENTRAL STATEMENT: no native-free operator over the present information reaches below
about 3.18 A on the built chain; the honest central estimate is 3.21 A, production itself; and
the nearest margin, 2.98 A, requires a per-target sign measured at chance.
THE GAP TO 2.5 A IS NOT a search gap, NOT expressivity (a 27-parameter family holds a 0.25 A
structure on every target), NOT aggregation (the hull of the shipped top-75 contains a 2.00 A
point and the pool's hull a 1.12 A point) and NOT optimisation. IT IS ONE NUMBER.
THE BOUND'S END-TO-END CHECK, WITH ITS PROVENANCE CORRECTED (S29-L28, lane D, certified from
git): T's sign formula gain = RMSD(1 - sqrt(1 - rho^2 (2q-1)^2)) reproduces lane O's PC1 triple
-- per-target magnitude, 52% sign accuracy, +0.0071 A leave-fold-out -- exactly, BUT IT WAS NOT
DERIVED BEFORE THOSE NUMBERS WERE READ. The term `2q-1` first enters THEORY.md at commit
53c42bd4, 00:38:49, which is 7.5 minutes AFTER lane O posted S29-L21 at 00:31:19, and
`git log --all -S "2q-1"` finds nothing earlier. The DERIVATION is unaffected and the agreement
is real; the agreement is POST HOC and must be labelled so wherever it is quoted, including in
the final report. I repeated the "derived before" claim in this file and to the user before it
was checked, and both are corrected here. The provenance rule (contract 27) caught its author.
LOAD-BEARING ASSUMPTION (B2): rho_max <= 0.14 for every field constructible from the present
information. Derived for the marginal class by theorem 2; measured outside it across 31 scorers,
38 signals, eight routers, the typicality axis and PC1. Lane D's attack is to find a field the
survey missed (the AMBER step, the Legacy gradient, the projection's own displacement, top-25 vs
top-150 differences, medoid minus mean, the ESM contact-weighted displacement).
THE THREE EXITS, each of which must break B2 and only B2: a better distance prior (moves e
itself, slope -2.15 A per unit); a learned residual with errors decorrelated from the
predictor's (blocked by S19's error coherence); and a physics term on the emitted structure in
its FREE-ENERGY form, which sits outside the marginal class so theorem 2 does not bind it. The
third is the queued S8 free-energy item and is now the only structurally live exit in the sprint.


## Integration note 16 (2026-09-20 00:45, after S29-L24): A WITHDRAWAL, THE REPORT'S SENTENCE, AND THE ONE LIVE EXPERIMENT
1. LANE T WITHDREW ITS OWN SECTION 1.4 RECOMMENDATION. The separation-band re-weighting it
   offered as the only live form of the calibration item is NOT the same operation as lane P's
   cloud rescale (re-weighting the objective changes WHICH 75 real windows are retained, and a
   difference of real windows is realisable by construction, where a per-separation correction
   of a cloud is not). But it is closed anyway by a measurement that predates the sprint:
   `s12/obj_FINDINGS.md` section 4, 126 targets -- leave-fold-out per-shell weighted L1 gives
   3.058 against the shipped 3.048 (worse), the deployable shell-profile-only arm gives 3.163,
   and `score_weights.json` was fitted and retired once already (refitting drove two of six
   parameters to their clip bounds and lost 2.76 -> 2.91 A on the benchmark while gaining 65%
   on the dev objective). CALIBRATION IS CLOSED, NOT REDIRECTED.
2. THE SENTENCE THE REPORT WILL USE for lane P's curve, from lane T, and it replaces the word
   "contraction" everywhere in my own notes and in the published S28 page: averaging does not
   contract the structure; by the exact identity d(C)^2 = <d_k^2> - s^2 the emitted distance is
   the members' RMS distance less their own spread, so the distortion is a SEPARATION-DEPENDENT
   SHEAR -- a large shrink where distances are short and the spread is comparable to them
   (-23% at the virtual bond) fading to nothing where distances are long -- multiplied by the
   pool's inherited long-range over-extension (the posterior's own -0.048 to -0.589 A across
   separations, S25 L1), which the fading shrink no longer masks; the ratio therefore crosses
   1.0 near |i-j| = 8 and ends above it, and "the 22% contraction" is the left-hand end of that
   curve rather than a property of the operator. T also gives lane P a free check: the shrink
   factor is native-free and must be <= 1 and monotone, the pool-bias factor is ORACLE and must
   carry the entire crossing, and their product must reproduce P's measured profile to floating
   point.
3. THE ONE LIVE EXPERIMENT, assigned to lane M after F1. Same S12 table: a scorer knowing ONLY
   the native's per-separation profile -- 11 to 14 numbers, no pair detail -- reaches 2.299 A
   against 3.048, while the deployable version is 3.163. It is the LOWEST-DIMENSIONAL NAMED
   QUANTITY whose ORACLE version clears the bound's load-bearing threshold of 0.14 by a wide
   margin, and its leading component is compactness, where the record's achievable native-free
   proxies reach 0.24 to 0.37. The new angle over S12: fit the profile RATIO leave-fold-out from
   native-free features chosen for compactness explicitly, rather than fitting the profile.
   Registered prediction (T's, adopted): inside 0.05 A of production, and a result only if it
   clears 0.7x MDE with the fold CI excluding zero AND its cosine clears 0.140 with the shrink
   signature. If it lands near the ORACLE 0.75 A, assumption B2 is falsified and the sprint's
   headline changes from a ceiling to an opening.






## THE SIGN IS THE INFORMATION, MEASURED ON 21 FIELDS AT ONCE (2026-09-20 01:32, S29-L35)
Lane D's field survey confirms the bound AND produces the sprint's cleanest statement of its
thesis. Three results, all on 126 targets with fold CIs.
1. B2 SURVIVES. No field's signed mean cosine clears the 0.140 random reference with a 2-SE
   margin; the best is DISTPOT +0.1128 with a fold CI whose upper end is 0.137, just under it.
   The largest ORACLE gain any of the 21 fields buys through a global step is 0.019 A. Eight of
   the 21 have fold CIs excluding ZERO, so the honest form of B2 is "no field's mean alignment
   reaches the random reference", not "every field is uninformative".
2. THE MAGNITUDE IS THERE AND THE SIGN IS NOT -- THE RESULT THAT ADDS SOMETHING. The PER-TARGET
   |cos| is 0.251 to 0.325 on EVERY ONE of the 21 fields, about TWICE a random direction's
   0.140. These fields are genuinely aligned with the direction to the native on each target;
   what destroys the mean is that the sign is near a coin toss (fraction positive 0.468 to
   0.643). Granting an ORACLE sign and an ORACLE step per target and nothing else:
     EXPAND (pure de-contraction) |cos| 0.296, frac+ 0.468 -> 2.708 A
     RG_LAW 0.294 / 0.579 -> 2.735 | PROJ 0.294 / 0.587 -> 2.740 | ENV 0.305 / 0.571 -> 2.740
     LEG 0.299 / 0.619 -> 2.756 | MSET_500 0.303 / 0.643 -> 2.759 | MEDOID 0.251 / 0.476 -> 2.894
   against production's 3.0483. This reproduces S29-L23's structure on 21 fields instead of one
   and puts the SIGN-ORACLE CEILING at 2.708 A on the point cloud, BELOW lane T's 2.98 figure.
   Assumption B4 (a field good on some targets and reversed on others enters at |rho|(2q-1)) is
   therefore not a technicality: it is the entire gap, 0.021 signed against 0.296 unsigned for
   EXPAND. "The sign is the information" is now measured for a whole operator class, and the
   class includes the shipped projection.
3. AND A THIRD SIDE OF THE CONTRACTION STORY. EXPAND has the LARGEST unsigned alignment of all
   21 fields and a signed mean of -0.021 with 46.8% positive. Production is systematically 22%
   contracted in the bond, yet expanding it toward the pool's mean radius points AWAY from the
   native as often as toward it. Contraction is real and is NOT a globally correctable
   direction -- lane L's Jensen framing, lane D's own "the shipped cost's descent EXPANDS", and
   now this, meeting from three sides.
4. THE ALIGNMENT IS A FAIL18 SET PROPERTY: mean cosine on the 18 against the 108 is +0.348 vs
   +0.046 for CONS_TRIM, +0.316 vs +0.078 for MSET_250, +0.249 vs +0.040 for PROJ. The
   displacement information lives where the error is largest, which is where a per-target sign
   would be worth the most and where the sprint has repeatedly found no way to get it.

## THE BOUND SURVIVES ITS DIRECT ATTACK (2026-09-20 01:28, `s29/results/s29_D_fields.json`; lane D's entry to follow)
Lane D's field survey is the one experiment that could have falsified the sprint's headline, and
it does not. It measured the ORACLE cosine of EIGHTEEN displacement fields the original survey
never covered -- nine channel re-rankings, the m-set differences at m = 5, 50, 150, 250, 500, the
MEDOID minus mean displacement, CONS_TRIM, EXPAND (pure de-contraction) and PROJ (the projection's
own displacement) -- against the random-shape-field reference of 0.1398.
NOT ONE FIELD BEATS THE RANDOM REFERENCE. The ranking, signed mean cosine over 126 targets:
DISTPOT 0.1128 (SE 0.0316), MSET_250 0.1118, MSET_150 0.0946, RG_LAW 0.0933, CONTACT 0.0933,
CONS_TRIM 0.0890, LEG 0.0883, MSET_500 0.0815, CONTACT_LL 0.0786, SS_MATCH 0.0750, ENV 0.0730,
PROJ 0.0694, MSET_5 0.0629, CAGEO 0.0565 ... down to DIS_MEAN 0.0231, MEDOID -0.0097,
MSET_50 -0.0187 and EXPAND -0.0208. `beats_random_reference` is False on every one of them.
AND THE PRICE IS THE POINT: the implied point-cloud RMSD at each field's own BEST (ORACLE) step
runs from 3.0289 for the best field to 3.0482 for the worst, against production's 3.0483. The
single best displacement field this project can construct, stepped by an amount chosen with the
native in hand, is worth 0.019 A on the point cloud.
SO ASSUMPTION B2 -- rho_max <= 0.14 for every field constructible from the present information --
SURVIVES a direct, pre-registered attack with eighteen new candidates, and the bound of S29-L23
stands: no native-free operator over this information reaches below about 3.18 A on the built
chain. Lane D designed this attack on my instruction to try to break the sprint's own headline,
and it is worth recording that the attack was real: PROJ, MEDOID and EXPAND are exactly the three
fields a critic would name, and all three are at or below the reference.

## RECOGNITION IS NOW CLOSED THREE WAYS (2026-09-20 01:25, S29-L33)
The within-realism-band measurement, the sprint's only open recognition question, has landed.
Lane D designed it, ran it, and reports against its own registered criterion.
THE DESIGN WAS VALIDATED FIRST, on the coordinator's relay of lane L's warning: of the three
realism statistics, R1 (CAGEO percentile) and R3 (consensus percentile) ARE compactness-loaded
(rho with Rg +0.460 and +0.404) and their cells are read with that stated; R2 (distance from the
pool's median bond/Rg) is linearly orthogonal at -0.015 and is the clean definition. Lane D adds
the caveat that makes the result honest: R2 is a FOLDED function of Rg, so banding on it removes
the MAGNITUDE of the compactness deviation and leaves its SIGN free -- and the sign is exactly
what S14 and rung 8 identify as the per-target information.
F1 FIRES ON 50 OF 70 CELLS: information orthogonal to realism EXISTS in the library and is
substantial (DIS +0.5397 in-band under the clean R2; max-over-scorers p 0.000 in every arm).
This is NOT the flat null lane L's caution anticipated.
F2 FAILS ON 58 OF 70, AND F2 IS THE ONE THAT MATTERS: conditioning on realism does not REVEAL
hidden accuracy information, it REMOVES ordering that was already there. DIS goes +0.568 global
to +0.508 in-band under R1 and +0.540 under R2, with the fold CI of the paired difference below
zero on every informative scorer. The reading: a large part of what the library knows about
accuracy IS realism, which is the perception-distortion theorem's content, and the part that is
not is the same ordering the record already prices globally.
ONE GENUINE SIGN CHANGE, PRICED AND DECLINED: SS_MATCH (secondary-structure compatibility) is the
only scorer uninformative globally and informative in-band (-0.009 to +0.060, replicated across
two realism definitions and both arms, 5/5 folds, 1.15x MDE). By the squared-skill law it is
worth of order 1e-3 A. Lane D explicitly does not propose an arm on it and calls it a mechanism
result. That is the correct handling.
A METHODOLOGICAL CORRECTION FOR THE RECORD: the jointly-Gaussian closed form rho_SY.R, which lane
L supplied and lane D pre-registered as its prediction, SYSTEMATICALLY OVERSTATES what survives
narrow conditioning here (predicted +0.596 against +0.508 measured for DIS under R1; the width
curve converges BELOW the prediction). On this instrument it is an UPPER BOUND, not an estimate,
and any future use must be labelled so. The width curve itself behaved exactly as lane L
predicted and is the right object to report.
SO RECOGNITION IS CLOSED THREE WAYS, each independent: ACROSS realism levels by the
perception-distortion theorem (no realism measure can prefer the distortion-optimal answer);
WITHIN a realism band by this measurement (conditioning removes ordering rather than revealing
it); and PER TARGET by the incidental-parameter theorem (the sign is not estimable from other
targets' answers). The sprint has no open recognition question left.

## THE UNIFIED FINDING HAS A THEOREM (2026-09-20 01:11, S29-L31): THE SIGN IS AN INCIDENTAL PARAMETER
Lane L identified what the sprint has been circling: the missing per-target sign is an INCIDENTAL
PARAMETER in the sense of Neyman & Scott (Econometrica 16:1-32, 1948) -- one nuisance parameter
per target with a bounded number of observations per target. In their canonical case the MLE of
the variance converges to HALF its true value: inconsistent, not merely inefficient. Mapped here,
s(x) = f(x) + g(j) with f the structural part cross-target training estimates and g(j) the
per-target sign. S14's 0.986 within-target against 0.600 across-target, with an overfitting gap of
-0.0005 and a LEARNING CURVE THAT DECLINES WITH MORE TARGETS, is that theorem measured on this
instrument.
AND THE HALF THAT CLOSES THE LOOP. The standard remedy for an incidental parameter is to ELIMINATE
it (conditional or fixed-effects likelihood), which is exactly the within-group pairwise design
lane D is running -- and ELIMINATING IS NOT ESTIMATING. Conditional likelihood buys a consistent
estimate of the shared f precisely by discarding all information about g(j). So the choice is:
estimate f consistently and discard the sign, or estimate f inconsistently and contaminate it.
Neither returns the sign. THE PER-TARGET SIGN IS NOT ESTIMABLE FROM OTHER TARGETS' ANSWERS AS A
MATTER OF STATISTICAL THEORY, NOT OF MODEL CAPACITY. The literature names exactly two escapes:
replication WITHIN the instance (an independent second source), or a covariate observed at
inference. Three adjacent literatures do NOT transfer and are recorded so S30 does not spend time
on them: phase retrieval (sign information is disproportionately VALUABLE, not RECOVERABLE), PCA
sign ambiguity (a gauge fixed by convention; ours is a real latent state with a right answer), and
one-bit compressed sensing (recovers a direction FROM many sign measurements, the converse).
WHY THIS INSTRUMENT IS HARDER THAN THE PROTEIN CASE, stated cleanly for the first time: the
field's per-target conditioner is MSA depth, and a 13-residue query has no family -- any hits are
the fragment's parent proteins, which is the leakage the fold clustering exists to exclude. The
conditioner is not weak here, it is ABSENT. Recycling is optimisation, not information (it feeds
back the network's own outputs, so no external evidence enters); templates are leakage; confidence
heads are functions of the model's own output and were measured to have no in-band skill at
peptide length.
LANE L PRICED ITS OWN PROPOSAL AND THE PRICE KILLS IT AS A ROUTE. Using lane T's formula at
rho = 0.37, a sign classifier at 60% accuracy buys 0.009 A, at 70% 0.035 A, at 80% 0.080 A, at 90%
0.144 A and at 100% 0.228 A: quadratic in (2q - 1), the squared-skill law again, and the whole
direction ceilings at 2.985 A. Lane L had priced it about 4x too high in its own draft and
corrected it rather than quietly fixing it. CONSEQUENCE, adopted: the per-target sign regression
remains S30's cheapest decisive FIRST measurement and is decisive in both directions, but it is a
MECHANISM measurement and not a route to the charter's target. S30's only plausible route to a
materially better number is A BETTER DISTANCE PRIOR (-2.15 A per unit, S24 L13), one of the three
classes lane T's section 7 names as able to break assumption B2. Those are two different projects
and the report must not let the first stand in for the second.

## THE ARCHITECTURAL CEILING (2026-09-20 01:04, S29-L30, rung 9): 2.76 A WITH THE NATIVE IN HAND
Nobody had ever measured the ORACLE ceiling of the TOP-128 PREFIX, which is the only part of the
pool the deployed quantum stage can see (`core/pipeline.py:758` widens the prefix to 2**n = 128
when the stage is on). Because the set-equality theorem says the realised CVaR tail is ALWAYS a
prefix of the energy order, the best prefix-m average IS the tail's exact reachable set, so its
ORACLE value is a hard ceiling on every arm the deployed quantum architecture can emit.
| operator class (ORACLE, point cloud, 126) | top-75 | TOP-128 | K = 500 |
|---|---|---|---|
| best single member (a perfect ranker) | 2.3062 | 2.1458 | 1.7108 |
| best prefix-m average (THE CVaR TAIL'S REACHABLE SET) | 2.8267 | **2.7605** | 2.6062 |
| convex hull in the shared frame (any weighting) | 1.9975 | 1.8071 | 1.1167 |
| production (deployable incumbent) | 3.0483 | 3.0483 | 3.0483 |
THE STATEMENT THE REPORT WILL LEAD WITH ON THE QUANTUM SIDE: the best structure the deployed
quantum selection can emit, WITH THE NATIVE IN HAND, is 2.7605 A. The charter's 2.5 A target is
therefore unreachable through this architecture even with perfect ORACLE selection -- not because
the search is hard, but because the tail is a prefix and the prefix's ORACLE ceiling is above the
target. The whole quantum stage's ORACLE headroom over production, inside its own field of view
and through its own readout, is -0.2879 A.
WIDENING THE PREFIX HELPS (ORACLE) AND TRANSFERS NOTHING. 75 -> 128 buys -0.0663 on the prefix-m
average (1.84x MDE, 5/5 folds, 50W/0L/76T) and -0.1903 on the hull; but the DEPLOYABLE version,
the prefix chosen leave-fold-out (m per fold 70/108/72/72/72), gives +0.0079 A at 0.30x MDE, NOT
MEASURED. The per-target choice is 99% of the ORACLE gap and -3% of it transfers.
THIS IS THE THIRD INDEPENDENT INSTANCE TONIGHT of one pattern: rung 6 (the typicality axis, ORACLE
best global step exactly 0.0), rung 8 (PC1, ORACLE best global eta exactly 0.0), rung 9 (the
prefix, ORACLE per-target 99% and transfer -3%). It also sharpens S22 L4, which measured the
per-target m transferring 65% WITHIN a target's own pool halves: ACROSS targets it transfers
nothing. Every one of the three says the same thing in a different operator class -- the quantity
the system is missing is per-target, and nothing native-free supplies it.

## THE SPRINT'S QUANTUM RESULT (2026-09-20 00:52, S29-L25; provisional until lane D's attack)
THE SET-EQUALITY THEOREM FAILS ON REAL POOLS, AND THE ESCAPE BUYS NO ACCURACY. Both halves
measured in one job, and together they are this sprint's thesis demonstrated rather than asserted.
MECHANISM (deterministic, per target, no statistics needed; exhaustive over all C(500,2) = 124,750
pairs of the deployed pool; 12 targets): with the objective placed on the tail's own coordinate
average, V(S) = f(mean of the subset), the f-optimal PAIR is NOT the energy-order prefix on 11 of
12 targets (mean objective gap +0.0847), the f-optimal m = 5 subset is not a prefix on 12 of 12
(mean gap +0.1522, clearing lane T's registered 0.10 at m = 5 rather than m = 75), supports reach
ranks 324 and 269 of 500, AND IT IS NOT A DIFFERENT SORT EITHER: the pair of the two lowest
per-state f values equals the exhaustive optimum on only 1 of 12 targets. So the optimum is
reachable neither by sorting E nor by sorting f. It is a genuine subset-selection problem, exactly
as lane T's three-state two-dimensional witness (S29-L17/L18) requires.
THIS IS THE FIRST MEASUREMENT IN THE PROJECT'S RECORD of a formulation in which "which set" is a
real optimisation variable rather than the read-out of a sort, and it is an affirmative answer, on
the deployed pool, to charter section 11 questions 5 and 7. The deployed spine cannot do this: its
tail is provably a prefix (Barkoutsos eq 12), which is why lane T could derive that the whole
deployed quantum stage reduces at the endpoint to choosing one number m.
PRICE (ORACLE, point cloud, n = 12, a PROBE and labelled one): the f-optimal m = 5 subset's
average is +0.4278 A WORSE than production (0.59x MDE, NOT MEASURED at this n, fold CI straddling
zero, 4W/8L, power 0.38, Type-M 1.60). The direction is the one that matters and it is the wrong
one.
CHRONOLOGY (contract rule 27): I registered this exact prediction in integration note 11 at 00:26
-- "the flatness falls (the mechanism works) AND the endpoint gets WORSE, because f is a
marginal-class objective bound by theorem 2" -- and lane B posted the measurement at 00:53. The
prediction preceded the measurement and both halves came in as registered.
WHAT IT MEANS. The flatness WAS fixable: a lift exists, it is derived, it costs the same 2P
evaluations, it is more device-realisable than the deployed objective, and on real pools it
genuinely escapes the prefix. And fixing it did not help, because the scoring function it consumes
is still bound by theorem 2 and the bound of S29-L23. THE FLATNESS WAS NEVER THE BARRIER; THE
INFORMATION IS. That is the sprint's thesis, and it is now demonstrated from the constructive side
as well as the destructive one.


## Integration note 17 (2026-09-20 00:55, after S29-L26): THE ADVERSARY FALSIFIES PART OF THE THEORY, AND REFINES THE M6 CONTROL
Lane D wrote its checkers from the STATEMENTS in THEORY.md, importing neither lane T's nor lane
B's code, so agreement cannot come from a shared implementation error. Three outcomes.
1. THE INFORMATIVENESS IDENTITY HOLDS to a relative 2e-14 on 126/126, and lane D is right that
   this is a CORRECTNESS CHECK on two implementations rather than evidence for the theorem. The
   geometric compression reproduces too: ker(Jc^T) is 48.0% of pair space averaged over the 126
   (T's illustrative 45%).
2. THE SIGN LAW IS FALSIFIED AT ITS OWN REGISTERED BAR, BOTH LIMBS. T registered "beta > 1 on at
   least 2/3 of targets and sign(cos) = -sign(beta - 1) on at least 70%". Measured: beta median
   0.756, above 1 on 13%, sign agreement 49% with a coin-toss CI of [41%, 59%]. Lane D attacked
   the STRONGEST version first (the w*kappa-weighted regression the expectation actually
   implies, a third reference map, the discrete posterior median, and a per-shell breakdown):
   every variant is below 1 and none reaches the bar. The third sub-clause fails in the opposite
   direction to its prediction.
   WHAT THIS KILLS: corollary 2b as a predictive law on this instrument. WHAT IT DOES NOT TOUCH:
   theorem 2's central result, that the expectation contains no term in the native's deviation
   from typical -- a statement about what the marginals CANNOT do, independent of beta's value.
   The bound (S29-L23) rests on that half, not on 2b.
   OPERATIONAL CONSEQUENCE FOR MY OWN CONTRACT RULE 20: the rule (report the implied shrink and
   the native percentile beside any cosine gain) STANDS as prudence, but its stated MECHANISM
   ("shrinking moves beta and therefore the cosine") is now unsupported on this instrument and
   must not be asserted until lane D's shrink grid measures it directly. That job is running.
3. M6, THE FIXED-PROFILE CONTROL, REFINES T's Q1 RATHER THAN CONFIRMING IT. Lane D's independent
   computation of the optimal rank-weight profile lands on T's numbers to the third decimal
   (t* -1.5329 vs -1.534; F* -4.7221 vs -4.7237; m* 30 vs 29; PR 340.4 vs 342.3). But T's
   prediction that the deployed arm's emitted structure matches the fixed profile on >= 120 of
   126 targets FAILS: exact set equality on 3 to 10 targets, structure agreement within the
   built-chain floor on 43/126 and within 0.02 A on 77/126, worst target 0.58 A apart.
   AND YET the endpoint contrast does not fire: the fixed prefix minus the deployed VQE arm is
   -0.0097 A at 0.43x MDE (NOT MEASURED, fold CI straddling zero, power 0.22).
   THE HONEST STATEMENT, which the report will use: the deployed quantum stage is NOT
   bit-reproducible by a fixed target-independent rank profile, but it is STATISTICALLY
   INDISTINGUISHABLE from one at the endpoint. T's reduction is right in effect and too strong
   in detail, and lane D's measurement is the version that goes in the report.


## Integration note 18 (2026-09-20 00:56, after S29-L27 and S29-L28): TWO OF MY OWN ERRORS
1. A FALSE PROVENANCE CLAIM I MADE AND REPEATED TO THE USER, now corrected in the headline block
   above. Lane D certified all three chronology claims from git commit times, as I asked it to
   rather than taking my word. Two stand with 24-minute leads: my own prediction of lane B's
   result (note 11 at 00:27:56 against S29-L25 at 00:51:58) and lane T's rung-8 prediction
   (S29-L11 at 00:07:02 against lane O's S29-L21 at 00:31:19, with a stated failure threshold,
   the strongest provenance in the sprint). The third does NOT: the sign formula postdates the
   numbers it reproduces by 7.5 minutes. Lane D invited lane T to produce an earlier uncommitted
   artefact and amend if one exists; until then the agreement is post hoc.
2. A GATE I SET THAT WOULD HAVE KILLED A LIVE IDEA BY A CONSTANT. I told lane B "recompute the
   flat fraction for the new objective against the deployed 85.5%; if it does not fall
   materially, the idea is dead". Lane B registered BEFORE measuring that it would not fall and
   COULD not fall, with the derivation and a replacement gate: for any state strictly above the
   VaR the same envelope argument that zeroes dCVaR/dp_y also zeroes dR_alpha/dp_y, so the f
   term is flat on EXACTLY the CVaR term's subspace at every lam. Measured: flat fraction 0.8552
   at every lam, identical to the CVaR term's, reproducing lane T's 85.5% independently. My gate
   was a quantity that is constant by construction and reading it as written would have killed
   the idea for the wrong reason.
   THE QUANTITY THAT ACTUALLY MOVED, and it is the right one: under the DEPLOYED readout the
   emitted structure moves along 0 of 511 continuous directions (it changes only through the
   discrete scalar m); under tail-then-aggregate it moves along 74, and the objective is
   non-constant along ALL 74 (overlap of flat-objective with moving-readout directions = 0.0000
   at every lam). The sensitivity is real, not a near-zero derivative: the mean norm of
   dC/dp_y over the tail directions is 55.7 A per unit probability. That is the mechanism
   result, and lane B got it by correcting my gate rather than obeying it.


## Integration note 19 (2026-09-20 01:01): A DEADLOCK I CAUSED, AND THE FIX (governor v2.6)
FIVE JOBS ACROSS FOUR LANES WERE STARVED FOR UP TO 40 MINUTES and I caused it. When I raised the
governor's bands to run the box at the 94 to 95% the user asked for (v2.5: suspend above 94% RAM,
hard 95.5, resume below 92), I did not touch CPU_RESUME, which was 80.0 from S26. The resume
condition is RAM below 92 AND smoothed CPU below 80. With eight jobs deliberately pinning the CPU
at 88 to 96%, THE CPU CLAUSE COULD NEVER FIRE, so any job the governor suspended stayed suspended
for ever. Lane O's ladder shard 3 sat at 60 of 625 cells from 00:19 to 00:59; lane P's shard 0,
lane D's band run, lane B's gradient job and lane X's configuration run were all stopped too, and
none of the four lanes could see it because the processes were alive, registered and simply not
scheduled.
DIAGNOSIS AND FIX. Spotted by comparing shard progress across siblings (540, finished, 60, 420 of
625) rather than by any alarm, which is itself a gap: nothing in the harness reports a job that is
registered, alive and making no progress. Resumed all five by pid, raised CPU_RESUME to 97.0 with
the reason written in the source, and restarted the governor (v2.6). RAM, not CPU, is the real
constraint on this box; gating a resume on low CPU is incompatible with an instruction to run the
CPU hot.
COST: about 40 minutes of five jobs' wall time, no computation lost (every one is checkpointed and
resumed where it stopped). RECORDED HERE because the sprint's standard applies to the coordinator:
this is the second operational error of mine this hour, after the flatness gate that would have
killed a live idea by measuring a constant.
FOLLOW-UP FOR A LATER SPRINT, not now: the governor should log a STARVED warning when a suspended
job's suspension exceeds some multiple of MIN_SUSPEND, and jobrun should surface it.


## Integration note 20 (2026-09-20 01:02, after S29-L29): A MODEL WITHDRAWAL, AND ONE FINDING INSIDE IT
Lane T accepted all three of lane D's findings without qualification, amended THEORY.md IN PLACE
rather than rewriting it, re-ran nothing and recomputed no number of its own to fit. Corollary 2b
is withdrawn at its own registered bar. Theorem 2's central result and corollaries 2a and 2c are
untouched, and section 8's bound rests on the MAGNITUDE of the achievable cosine rather than on
its sign, so the sprint's headline is unaffected.
THE POST-MORTEM IS THE VALUABLE PART, because a withdrawal without a diagnosis is not a result.
The measured beta < 1 substituted into the theory predicts a POSITIVE cosine while the measurement
is negative, so a dropped term is bigger than the kept one. Two candidates, both measurable:
(i) THE LINEARISATION SATURATES. The risk's derivative 2F - 1 saturates at +-1, and lane T measured
    |2F-1| > 0.9 on 3.8% of pairs for one target, 5.5% for another and 53.8% for a third, with the
    worst saturation on a FAIL18 target, which is exactly where the cosine is most negative
    (-0.143). Where it saturates the coefficient stops being a function of the regression slope at
    all and depends on corr(sign(a-b), a), which beta does not measure.
(ii) ASSUMPTION A4 IS AN IDEALISATION, and this is the interesting half: the measured sign implies
    cov(a, n) > cov(b, n) -- THE POOL'S DEVIATION FROM TYPICAL TRACKS THE NATIVE BETTER THAN THE
    POSTERIOR'S MEDIAN MAP DOES, which is consistent with the pool being made of real structures
    and the median map not being one. The system does hold a weak native channel (sequence
    conditioning is worth 0.776 A, S12), and this says it sits in the POOL rather than in the
    posterior's summary of it. That is a pointer worth carrying into the report and the next
    sprint: it is an argument for operators that read the pool's own dispersion rather than the
    posterior's marginals, and it is orthogonal to everything the sprint has closed.
FOLLOW-UP, cheap and available to lane D when its queue clears: recompute the beta law on the
UNSATURATED pairs only (|2F-1| < 0.5). If the sign agreement clears the coin-toss CI there, (i) is
the mechanism and 2b is a statement about the risk's linear regime; if not, (ii) is.
ALSO AMENDED BY THE SAME SATURATION: section 1.4's mid-range over-weighting is 2.7x in the linear
regime and 1.6x where the derivative is saturated, so the honest figure is a target-dependent 1.6
to 2.7x. Direction unchanged. And lane T endorses my M6 phrasing without change and restates its
clause in the file as an ENDPOINT-equivalence control rather than a bit-equivalence one.

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
