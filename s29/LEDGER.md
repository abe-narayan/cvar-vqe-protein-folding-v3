# S29 LEDGER (append-only; `## S29-L<n> -- TITLE (date time, lane)`; number from the tail in one process)

## S29-L0 -- THE CHARTER, THE READING, THE POSITION, AND WAVE 1 (2026-09-19 23:35, coordinator)
Charter: `s29/BRIEF.md` (verbatim). Contract: `s29/S29_CONTRACT.md` (inherits S28's 12 rules and
its two addenda; adds rules 13 to 19: findings engaged per hypothesis, the information test for
every import, the nine questions and ten controls for quantum claims, 12-target probes before
126, the multiplicity count, mechanism beside outcome, the cost-RMSD meter). Running state:
`s29/STATE.md`. Governor v2.5 bands: suspend above 94% RAM, hard kill 95.5%, resume below 92%,
launch below 91% (the user's 94 to 95% instruction; `s26/governor.py` lines 57 to 65). Cap 6.

Reading done by the coordinator before this entry: `s27/REPORT_S28.md`, the S28 ledger tail,
`s27/RESUME_S28.md`, `s27/RETRACTIONS_S28.md`, lane A's damaged expectations,
`docs/STATE_BRIEF_2026-09-12.md` sections 4 and 5, `s26/REPORT.md` VII.2 (ESMFold infeasible on
this box, L13; the distogram already uses ESM-2 650M PCA-32, `core/data.py:886`), `s17`'s ESM
contact-head result (in-band content, length-gated, no shortlist value), `s23/LEDGER.md` L6 (the
scale is unreachable in principle), `s24/LEDGER.md` L2/L3 (score-selected pools have parallel
biases, cos 0.943; the unselected blind source is 31% independent and 0.76 A worse; mixtures sit
on a line). The "more information" routes the charter's second-order questions point at are
therefore already in the record: a larger language model (hardware), a structure-trained
predictor (leakage and hardware), the ESM contact head (measured), a learned residual generator
(S24), calibration of the posterior (S25 L2: worse under the readout), scale (S23 L6).

Position (`s29/STATE.md`): H1, the typicality axis (extrapolate the conditioned answer away
from the sequence-blind answer; the common-mode error as signal), with its strongest objection
already in the record (S24 L3's geometry says the sign is wrong) and a 20-minute ORACLE kill
(lane O rung 6). H0, the fallback that is probably the truth: the ceiling is the absence of
per-target information orthogonal to the typicality axis; S29 measures it decisively by
recognisable-operator class (lane O's ladder), searches for a source of such information (lanes
L, T, M), and builds the CVaR-VQE formulation in which the quantum stage consumes the pool
jointly so that the quantum stage is where any such information would act. H2, the divergent
direction (lane X): the pool is the wrong state space; a configuration-space CVaR-VQE with an
interacting Hamiltonian and a tail-ENSEMBLE readout, the one readout S13 to S21 never used.
Findings engaged: H1 attacks 3, 4, 5, 8, 11; H2 attacks 9, 10, 11 (aggregation) and accepts 8
as the risk; H0 accepts all eleven and measures 8 and 11.

Wave 1 (six lanes, spawned after this entry): L literature (permanent), M data-path map and
convenience choices and the harness audit, T theory (the contraction theorem, the locally
informative class, the spectral condition, CVaR over a structural observable, the reachable
set, the entropy term, what an objective must know), O the ORACLE ceiling ladder by operator
class and the typicality-axis probe, D adversary (permanent, rotating) and the cost-RMSD meter
and the suite, X the divergent state space. No build starts before O's ladder, T's sections 2
and 3, and M's map are posted. Budget: reading and theory 35%, probes 35%, full runs 30%.
Multiplicity: 0 endpoint comparisons at this entry.
Artefacts: this file; `s29/STATE.md`; `s29/briefs/S29{L,M,T,O,D,X}.md`.

## S29-L1 -- TOPIC 1, NATIVE-FREE STRUCTURE QUALITY: NO PUBLISHED QA METHOD HAS EVER BEEN TRAINED OR EVALUATED BELOW 40 TO 50 RESIDUES; THE FIELD'S FOUR SIGNAL CLASSES ARE EACH ABSENT, ALREADY OURS, LEAKAGE, OR MEASURED DEAD HERE; THE ONLY NATIVE-FREE SELECTORS SHOWN TO WORK AT 9 TO 25 aa RANK AN ENSEMBLE THEIR OWN ENERGY GENERATED (PEP-FOLD 2.6 A, APPTEST 1.96 A) OR ARE A CONVERGED FREE ENERGY (LINDORFF-LARSEN); AND AF2's pLDDT HAS NO WITHIN-TARGET SKILL ON 588 PEPTIDES OF 10 TO 40 aa, WHICH IS EXTERNAL CONFIRMATION OF CHARTER FINDING 8 (2026-09-19 23:54, L)

Question (brief topic 1, the coordinator's first need): is there ANY native-free quality signal
for a 9 to 16-mer that the field has shown to work, and on what data? Falsifier for the negative
half: a single published QA method trained or benchmarked at 9 to 16 residues, or any
general-purpose QA signal whose construction survives the length. No experiment; literature
only; 15 primary papers read for their construction, not their abstracts.

ANSWER. Two native-free selectors have peptide-length evidence, neither transfers.
(a) A coarse-grained force field ranking the ensemble IT generated: PEP-FOLD's lowest-sOPEP-
energy cluster centroid, 2.6 A mean cRMSD on 25 NMR peptides of 9 to 25 aa (Maupetit 2009);
APPTEST's XPLOR-NIH restraint energy, 1.96 A on 42 peptides of 9 to 25 aa (Timmons 2021).
Generation and selection share the energy, so the score is a convergence diagnostic over a
self-consistent set. Over a retrieval pool the same class of energy is measured WORSE than a
random subset here (AMBER +0.455, Legacy +0.330, 5/5 folds, S25 L16).
(b) The equilibrium population under converged all-atom MD (Lindorff-Larsen, Science 334:517,
2011; chignolin 10 aa, Trp-cage 20 aa), which is a FREE energy at 1e5 to 1e6 CPU-hours per
peptide; our budget is 6.43 core-equivalents and one AMBER process at a time.

THE NEGATIVE, WHICH IS THE POINT. Every general-purpose QA method excludes this length by
construction or by its own text: ProQ3 "targets that were shorter than 50 residues were filtered
out both from the CASP11 and CAMEO data sets"; VoroMQA's learning set is chains "longer than 99
residues" (12,825 entries); DeepAccNet trains on 50 to 300 residues; QMEANDisCo quotes its own
error bar as ~0.12 at 40 residues vs ~0.05 for larger; GraphQA and EnQA say nothing about short
chains. The four signal classes they reduce to (L_1 section 5): S1 packing/burial -- not a
VARIABLE at 9 to 16 aa (VoroMQA's own Table I separates good from bad by the solvent area share,
0.392 vs 0.470, which at peptide length is set by the length), and its CA forms ENV/HP/EXVOL are
coin tosses or anti on the built chain (S28-L48); S2 sequence-to-local-structure agreement -- IS
this project's axis (the distogram + DIS), anti at pref 0.071, and phi is sequence-blind at this
length (S13, 36.1 vs 36.4 deg); S3 consensus -- score_i = (1/N) sum_{j!=i} GDT-TS_ij
(DAVIS-EMAconsensus, verbatim from Kwon 2021) IS our CONS channel, and the family assumes
independent member errors while our pool is 68% common-mode by exact identity (S23 L9), with
CONS_POOL@chain the most anti-recognising of the 31 scorers (pref 0.056, pool-member contrast
-0.127, S28-L48); S4 single-point physics -- measured worse than random here.

THE EXTERNAL CONFIRMATION OF FINDING 8. McDonald et al., Structure 31:111-119 (2023), 588
peptides of 10 to 40 aa against NMR: "the lowest RMSD structures failed to correlate with lowest
pLDDT ranked structures"; full text, "there was no correlation between the first three ranks
assigned by AF2 and the structure that gave the lowest Ca RMSD" (preprint: 13% of lowest-RMSD
models were rank 0, ranks 1 to 4 took 26/18/24/19%). Taking rank-1 instead of best-of-5 costs
0.2 to 1.1 A, largest on the mixed-secondary-structure classes. That is the field's best
self-confidence signal, from the model that solved the protein case, with no within-target
ranking skill at peptide length -- the same failure this project measures, on a different
instrument, with vastly more information. Charter finding 8 is not a defect of the S27 library;
it is the field's position at this length.

ALSO RECORDED. (i) CASP14's assessors (Kwon 2021): single-model QA overtook consensus for the
first time (BAKER-experimental top-1 GDT-TS loss 8.4; BAKER-ROSETTASERVER top-1 lDDT loss 4.0;
DAVIS-EMAconsensus and GOAP "about average"), and "it is more difficult to train an EMA method
to estimate the superposition-dependent quantity, GDT-TS" -- our objective and our endpoint are
both superposition-dependent, which is the hard case by the field's own verdict. (ii) CASP15's
EMA winners are quasi-single-model (MULTICOM_qa, per-target corr 0.66, ranking loss 0.14): their
skill is bought with an INDEPENDENT reference ensemble, which this instrument does not have
(S24 L2/L3: the only second source is 31% angularly independent and 0.76 A worse). The value of
that family is the independence, not the scoring rule -- a statement about the generator, not
about QA. (iii) Training any of these leave-fold-out on our corpus is a dead end for three
measured reasons, not one: the features are already in `s27/cache`, a supervised combination of
all 31 reaches 0.960 held-out sign accuracy but the SAME rule on the RAND_SIGNED control reaches
0.952 with head-to-head 0.492 (it learned "is this a real trace", S28-L48), and in-band
discrimination here is signal-limited, not sample-limited (S13's flat learning curve).

THE ONE ACTIONABLE ITEM, and it is not an import. The only S4 sub-class with peptide-length
precedent is the free energy, and `docs/FINDINGS.md` section B records that S8 BUILT it --
`strain` = E(built) - E(freely relaxed), `E_free`, `F_qh` = <E> - kT S_qh with S_qh = 1/2 sum
log lambda of a 4 ps 300 K Langevin CA covariance, `F_boltz` = -kT log <exp(-E/kT)>, `width`,
`S_msf`, each through the same builder and relaxation as the native -- and it "completed only 1
of its 24 targets before the box filled up", committed and resumable as `python -m s8.relax
best`, with no claim made from one target. So the one class the literature says works at this
length is the one class this project started and never finished. Caveats stated with it: S_qh at
4 ps is noisy and downward-biased; the underlying channel (AMBER single point) is +0.455 A worse
than random, so the claim would be that the ENTROPY term rescues a channel measured anti; the
cost is one Langevin ensemble per candidate under the one-AMBER-process rule. Proposed as a
12-target pre-registered probe for lane D or O, NOT a build, and not by lane L (I run no
endpoint experiments).

VERDICT. Topic 1 yields no importable native-free scorer. KEPT (4, all as evidence or anchors):
Kwon 2021 (CASP14 EMA assessment); McDonald 2023 (the peptide pLDDT result); PEP-FOLD 2.6 A /
APPTEST 1.96 A as the published native-free anchors at 9 to 25 aa; Lindorff-Larsen 2011 (the
working selector is a free energy). REJECTED (5 families, with the reason each fails HERE):
packing/burial QA (the statistic is not a variable at this length); sequence-agreement QA (it is
our own axis, and it is anti); consensus QA (assumes independent member errors; 68% common-mode
by measurement); single-point physics QA (worse than random on this pool, with controls);
AF2-feature QA (leakage plus hardware, and its base signal fails at peptide length).
Multiplicity: 0 endpoint comparisons; no number in this entry is new measurement, every one is
cited to its sprint or its paper.
Artefacts: `s29/lit/L_1_native_free_qa.md` (15 papers with construction, assumption, check,
information test, verdict; sections 5 to 9 are the synthesis), `s29/lit/L_INDEX.md`.

## S29-L2 -- THE COST-RMSD METER IS BUILT AND THE S28 VALUES FOR THE SHIPPED COST REPRODUCE THROUGH IT (LADDER rho -0.402 CHAIN / -0.182 CA, COSINE -0.034, NATIVE PERCENTILE 0.369, PREFERENCE 0.071 CHAIN / 0.206 CA WITH THE POOL-MEMBER CONTROL AT 0.020 / 0.126); AND A CAVEAT THE FOUR NUMBERS CARRY: THE ANTI-CORRELATION IS A PROPERTY OF THE NEAR-NATIVE HALF OF THE LADDER, NOT OF THE COST OVER ALL STRUCTURES (2026-09-19 23:56, D)

EVERY NUMBER IN THIS ENTRY IS ORACLE (the meter reads the native to score the ladder, to point
the gradient and to place the native in its pool). The meter never tunes anything.
Question (contract rule 19, brief S29D duty 2): does a single command exist that prices any
proposed cost function on the four axes the charter names, and does it reproduce S28's values
for the shipped cost through its own code path?
Falsifier: any of the four S28 anchors not reproduced to the third decimal.

Code `s29/s29_D_cost_audit.py`; tests `tests/test_s29_D.py` (7 pass); results
`s29/results/s29_D_cost_audit_{DIS_ca,DIS_chain-s28rows,DIS_SURR_ca}.json`; ladder cache
`s29/results/s29_D_ladder_structs/` (126 npz, 1.1 MB, built in 72 s). Commit 7c8547dd. No job
needed: the cache build is 72 s at 0.3 GB and one meter run is 13 s.

WHAT IT DOES. `python s29/s29_D_cost_audit.py meter --f <cost> [--basis ca|chain|chain-s28rows]`
for `<cost>` either a name from the S27 channel library (DIS = the shipped lookup, DIS_SURR =
the linear-interpolation surrogate S~, the 15 CA scorers, the 16 backbone scorers) or a python
callable `module:function` with the contract `f(W, ctx) -> (m,)`, lower is better, over a stack
of CA clouds. THE CONTEXT IS NaN-POISONED BY CONSTRUCTION (`ctx.nat_ca`, `ctx.oracle_rr` are
NaN): a cost that reads a native quantity returns NaN or raises, and the meter REFUSES it,
naming the failure. It reports, on the 126 targets with fold-clustered CIs:
 (a) LADDER  the per-target Spearman rho(f, ORACLE RMSD) over three ladders (below);
 (b) COSINE  cos(-grad f at production, the direction to the native), rigid body removed from
     both (S28-L23b's machinery), analytic for S~/RG_LAW/EXVOL and batched central differences
     (h = 1e-3 A; within 3e-4 of the analytic on S~) otherwise, with the 16-draw random-
     direction reference per target;
 (c) PCTILE  the native's percentile in its own 500-member pool under f (`s28_A_objdiag`);
 (d) PREF    pref(ORACLE circ_best vs PROD) with the S28-L36 POOL-MEMBER CONTROL
     (pref(random real pool trace vs PROD) = PROD's own pool percentile) and their paired
     `ST.compare` contrast, plus the RAND_SIGNED contrast C2 registered as its clause 2.
Nine rungs per target, in lane A's frame, cached with every S28 assertion re-run at build time:
PROD (asserted equal to the deployed average to 1e-10), circ_opt (lane A's NATIVE-FREE
recognition optimum, `circ_l1_i80`), RAND_SIGNED[0], GAUSS_MATCHED[0], GAUSS_0.3[0] (C2's
controls, native-free directions and ORACLE scales), sub0, circ_s0, circ_best, NATIVE (ORACLE).
circ_s0 and sub0 are regenerated and asserted against lane A's `per_start[0]` / `per_sub[0]`
(max deviation < 1e-6 on 126/126) and every regenerated rung against C2's stored DIS scores
(< 1e-5 on 126/126), so the cache IS C2's ladder and no pinned artefact was regenerated.

REPRODUCTION OF THE FOUR S28 ANCHORS (the falsifier is silent; all to the third decimal):
| anchor | S28 | through the meter | source |
|---|---|---|---|
| ladder rho, built chain | -0.402 (`s28_C2_chain_summary.json :: ladder_rho :: DIS@chain`, S28-L48) | -0.4023 SE 0.0401 fold [-0.477, -0.322], 11/126 positive | `--f DIS --basis chain-s28rows` |
| ladder rho, CA level | -0.182 (S28-L35 table) | -0.1818 SE 0.0488 fold [-0.308, -0.053], 27/126 positive | `--f DIS --basis ca` |
| gradient cosine at production | -0.034 SE 0.021, random reference 0.140 (S28-L23b) | -0.0339 SE 0.0214 fold [-0.059, -0.014], 56/126 positive, random reference 0.140, Spearman(cos, production RMSD) -0.372 | `--f DIS` (the shipped table is piecewise constant, so the gradient is S~'s, exactly as S28-L23b) |
| native's pool percentile | 0.369, production below the native on 99/126 and below the pool's best on 6/126 (S28-L30, `s28_A_objdiag.json`) | 0.3688 fold [0.309, 0.406], 99/126, 6/126 on S~; 0.3676 [0.306, 0.405] under the shipped lookup | `--f DIS_SURR` / `--f DIS` |
| pref(ORACLE circ_best vs PROD), chain | 0.071 (S28-L48) | 0.0714 fold [0.034, 0.109] | `--basis chain-s28rows` |
| pref, CA level | 0.206 (S28-L35) | 0.2063 fold [0.116, 0.293] | `--basis ca` |
| pool-member control, chain | pct(PROD) 0.020, contrast +0.051, 0.80x MDE (S28-L48/L49) | 0.0201, +0.0513, +0.80x, fold [+0.007, +0.090], 4/5 folds | `--basis chain-s28rows` |
| pool-member control, CA | pct(PROD) 0.126, contrast +0.080, 0.83x MDE (S28-L36/L37) | 0.1265, +0.0799, +0.83x, fold [+0.006, +0.144], 4/5 folds | `--basis ca` |
The head-to-head numbers of S28-L36 reproduce too (the ORACLE structure beats a random pool
member under DIS on 0.626 at CA level and 0.632 on the chain, the native on 0.632 / 0.634).

THE CAVEAT EVERY FUTURE QUOTE OF "rho = -0.40" MUST CARRY (new here, not in S28). The rho
depends on WHICH ladder, and the sign flips between them. For the shipped cost:
| ladder (rungs) | CA level | built chain |
|---|---|---|
| S28 (PROD, sub0, circ_s0, circ_best, NATIVE: production plus four near-native structures) | -0.182 [-0.308, -0.053] | -0.402 [-0.477, -0.322] |
| CHARTER (circ_opt, RAND_SIGNED, PROD, GAUSS_MATCHED, circ_best, NATIVE: the brief's six rungs) | +0.260 [+0.201, +0.321] | -0.092 [-0.200, +0.016] |
| FULL (all nine) | +0.118 [+0.064, +0.166] | -0.236 [-0.325, -0.146] |
The shipped cost ORDERS THE BULK correctly at CA level (+0.26 over the charter's rungs: it
rejects a random signed combination at 3.88 A and a matched Gaussian at 4.12 A, preferring
production to each on 0.968 / 0.952 of targets) and ANTI-ORDERS the near-native half (-0.18 /
-0.40). "Garbage rejection reads as skill" is the project's own standing warning
(`decoy-bank-not-a-pool-proxy`); the S28 ladder is the one that excludes garbage, which is why
it is the one that matters, and the charter's six rungs are NOT a stricter test than S28's --
they are a weaker one. Any lane proposing a cost must beat the shipped cost on the S28 ladder,
and a positive CHARTER rho with a negative S28 rho is the failure mode to expect, not a pass.
Two further readings, both ORACLE: (i) on the chain the CHARTER rho is -0.478 on FAIL18 against
-0.028 on the 108, the same regime split as the cosine (-0.143 / -0.016) and as S28-L23b's
Spearman(cos, RMSD) = -0.372; (ii) the shipped cost prefers lane A's own native-free optimum
(circ_opt, mean 3.385 A, 0.34 A WORSE than production) to production on 0.786 of targets and
scores it below 0.969 of real pool traces: the objective's minimiser is a structure no real
trace resembles, which is S28-L18b's global statement read as a preference.

Verdict: the meter is built, reproduces every S28 anchor it was asked to reproduce, and is
open for business. Every lane's proposed objective goes through it before an endpoint run:
post the cost as `module:function` (native-free, batched) and quote all four numbers with
their ladders named. The meter is a GATE, not a result: passing it is necessary, never
sufficient, and nothing in it is deployable.
Multiplicity: 0 endpoint comparisons (this entry runs none; `s29/STATE.md` count unchanged).
Artefacts: `s29/s29_D_cost_audit.py`, `tests/test_s29_D.py`,
`s29/results/s29_D_cost_audit_DIS_ca.json`, `..._DIS_chain-s28rows.json`, `..._DIS_SURR_ca.json`,
`s29/results/s29_D_ladder_structs/`.

## S29-L3 -- ADVERSARY CHECK OF PREREG_S29_O (the ORACLE ceiling ladder and the typicality-axis probe): THE LADDER IS SOUND AND ITS GATES ARE REAL; RUNG 6 HAS NO CONTROL IN THE OPERATOR'S SPACE AND F6a's SECOND CLAUSE COMPARES A MEAN AGAINST A SINGLE-DRAW MAGNITUDE; THREE CONTROLS REQUIRED BEFORE ANY POSITIVE IS READ (2026-09-19 23:56, D)

Question (brief S29D duty 3): is each falsifier in `s29/PREREG_S29_O.md` falsifiable, and does
each control match the operator's space? Checked against the charter's sections 14 and 16, the
contract's rules 6, 7, 12, 16, 17, and the memory entries `control-must-match-the-operators-space`,
`zero-information-control-must-be-plausible`, `grid-oracles-are-order-statistics`.

WHAT IS RIGHT (stated first, because most of it is).
1. The labelling is exact: one deployable arm (`lfo_LIB75`), everything else ORACLE, said in
   section 0 and repeated per rung. The ORACLE ladder is a ceiling table by operator class and
   is not read as a result.
2. The reproduction gates are real gates, not decorations: production's 3.048338 to 1e-9 per
   target against the S28 frame; rung 1 against S28-L1b's 2.3062 / 1.7108; rung 2 at m = 75
   equal to production to 1e-9; LIB75 and BPRIME rebuilt from S24's own stable RNGs and gated
   against `s24/results/biasalign.json` and `qmatch.json` to 1e-6 PER TARGET, with the run
   stopping on failure. That is the right way to reuse an artefact whose structures were never
   stored, and it is the strongest part of the prereg.
3. Ties: rung 2's DIS order is `np.lexsort((key, DIS))` with the S27 stable key (rule 12), and
   rung 1 reports tied `rr` values instead of taking an argmin silently.
4. Order statistics are priced where they arise (rungs 2, 3, 4 and the per-target t*), with
   `best_of_k_within` AND the split-half transfer, which is the construction that nulls itself
   (`grid-oracles-are-order-statistics`: only split-half transfer arms survived that audit).
5. The chain comparator is production RE-PROJECTED IN THE SAME JOB (`prod`), so the S28-L18
   branch-flip floor (12/126 above 0.02 A, one at 0.5) cannot enter the one contrast that
   matters. The prereg says so and cites the reason.
6. The regime clause is the S28-L40 construction (a FAIL18 sign is a regime claim only at
   random-18 p < 0.05), which is the right bar and is pre-registered rather than invented after.
7. F6b is falsifiable and correctly ordered: the point cloud gates, the built chain decides,
   and a chain effect below 0.7x MDE is NOT A RESULT.

WHAT MUST CHANGE BEFORE ANY RUNG-6 POSITIVE IS READ (three items; all cheap).
(a) THE STEP HAS NO CONTROL IN ITS OWN SPACE. The operator is "displace production by t|u| along
    u"; its only comparator is production. S28-L39/L40 measured exactly this class and found the
    shipped objective's own descent step INDISTINGUISHABLE from a random direction of the same
    size (0.0x to 0.3x MDE at every e) -- the step "beat production" question is not the
    question; "beat a displacement of the same size" is. Required beside `lfo_LIB75`, on the
    same targets and the same basis: (i) RANDOM DIRECTION, >= 8 draws per target, rigid-body
    removed, scaled to the SAME per-target RMS displacement as the LFO step, the MEAN of the
    draws as the control (S28-L39's construction; a best-of-8 is an order statistic and is not
    a control); (ii) SCALE-ONLY, X = C dilated about its centroid to the same RMS displacement.
    (ii) is the decisive one on this axis: production is a 22%-contracted trace (S25; the
    meter's own head-to-head has the shipped cost preferring production to 87% of real traces
    on contraction alone, S28-L36), the blind average is contracted too, and u = C - B_sup is a
    difference of two contracted structures, so "extrapolating along u" and "de-contracting"
    are confounded until measured apart. Without (i) and (ii) a positive at t > 0 cannot be
    told from "any displacement of this size helps", and the entry would have to say so.
(b) F6a's SECOND CLAUSE COMPARES QUANTITIES IN DIFFERENT SPACES. "the mean cos is at or below
    the measured mean |cos| of the random fields (about 0.14) with the fold CI including that
    value": 0.14 is the typical magnitude of ONE random draw's cosine; the statistic being
    tested is a MEAN OVER 126, whose random-reference distribution has SE about 0.14/sqrt(126)
    = 0.012. As written the clause can essentially never fire (a mean near zero has a CI
    nowhere near 0.14), so F6a reduces to its first clause. The first clause is the right test
    -- but its reference must be MEASURED, not assumed: compute the mean over 126 of the
    random fields' SIGNED cosines with its own fold CI and compare the observed mean against
    THAT (`mde-is-per-comparison-not-per-instrument`; and when an analytic null and a measured
    null disagree, the measurement is the null). Recommended wording: F6a fires if the fold CI
    of mean cos(u, v) contains the random fields' signed-mean reference. Keep the mean |cos|
    reported as the per-draw scale, and stop using it as a bar.
(c) THE COSINE'S ZERO-INFORMATION CONTROL IS NOT PLAUSIBLE IN THE OPERATOR'S CLASS. A Gaussian
    shape field is this space's "uniform on the torus" (`zero-information-control-must-be-plausible`:
    a uniform control is a WORSE measure, not an uninformative one). u is a difference of two
    pool-averaged protein structures: low-frequency, dominated by the contraction mode, nothing
    like isotropic. The matched null costs one extra draw: u_null = B1_sup - B2_sup from TWO
    independent blind LIB75 draws (the same `SD.stable_rng` construction at a second draw
    index), superposed identically, rigid body removed, scaled to |u|; cos(u_null, v) is then
    the cosine of a direction with u's geometry and NO conditioned-vs-blind content. If
    cos(u, v) is not above cos(u_null, v) with the fold CI excluding it, the typicality axis
    carries nothing beyond the geometry of pool-average differences. I regard this as the
    single most informative addition to the prereg, and it is 20 minutes.
FOUR SMALLER POINTS.
(d) t* at a grid endpoint: T = {-1.0 .. +2.0} is fixed and pre-registered (good). If the
    leave-fold-out argmin lands at +2.0 the grid is truncating the operator and the entry must
    say "at the grid edge", not quote t* as an estimate.
(e) `lfo_LIB75` is DEPLOYABLE in the project's sense (a scalar fitted on other targets'
    natives), and its NaN-poison test is registered. Keep the poison at the HELD-OUT target
    only, as written; poisoning the training folds would test nothing.
(f) Rung 5's gate ("within 0.05 A of S10-5's 1.802 / 0.953; if BETTER the convention differed"):
    a better ORACLE hull number is a better solver, never a result, and the entry should say
    which convention differed rather than quoting the improvement.
(g) Multiplicity: 4 endpoint comparisons + 2 stratum tests is the right count and is declared.
    With the sprint at 0 endpoint comparisons before this lane, a rung-6 positive at 1x MDE
    would be priced against max-over-4; say so in the entry rather than at review time.

Verdict: PREREG ACCEPTED FOR THE ORACLE LADDER (rungs 1 to 5 and 7) AS WRITTEN. RUNG 6 IS
ACCEPTED FOR THE MEASUREMENT AND HELD FOR ANY POSITIVE: no `lfo_LIB75` result is read as a
result until (a) the random-direction and scale-only controls and (c) the blind-difference
cosine null are in the same entry, and F6a is decided on (b)'s measured signed-mean reference.
If the registered prior holds and rung 6 dies, none of this matters and the entry says so in
one line; these three are the cost of the positive being believable if it does not.
Multiplicity: this check runs no comparison. Artefacts: `s29/PREREG_S29_O.md` (read at
2026-09-19 23:56), `s27/LEDGER.md` S28-L39/L40 (the step-vs-random construction), `s24/LEDGER.md` L2/L3.

## S29-L4 -- ADVERSARY CHECK OF PREREG_S29_X (the configuration-space CVaR-VQE probe): FALSIFIERS ARE FALSIFIABLE AND THE CONTROL SET IS THE BEST IN THE RECORD; FOUR HOLES -- THE ORACLE-BEST CHIMERA IS AN 8^S ORDER STATISTIC WITH NO MATCHED PARENT CONTROL, D1's PRIOR IS UNFALSIFIABLE AS WORDED, THE R2-R1 CONTRAST IS SET-MATCHED BUT NOT WEIGHT-MATCHED, AND 12 TARGETS CANNOT CARRY A 0.7x-MDE GO (2026-09-19 23:56, D)

Question (brief S29D duty 3): is each falsifier in `s29/PREREG_S29_X.md` falsifiable, and does
each control match the operator's space? Checked against charter sections 11 (the nine
questions and the ten controls), 14 and 16, contract rules 9, 10, 13, 15, 16, 17, and the
memory entries `concentration-is-wrong-when-discrimination-binds` (the control discipline),
`grid-oracles-are-order-statistics`, `exhaustive-enumeration-closes-the-search-half`.

WHAT IS RIGHT.
1. The closure table (section 1) is the model of what contract rule 10 asks: every closed line
   cited with its ledger entry and a stated reason the closure may not apply. The set-equality
   theorem is ACCEPTED AS BINDING up front, and the entry states in advance that the quantum
   stage's whole possible contribution beyond m is R2 - R1 on the same set. That is the honest
   framing S28 had to be argued into.
2. The controls cover nine of charter section 11's ten (classical equivalent GIBBS, diagonalised
   GS, permuted PERM, random/untrained, matched budget SA at 2^q, untrained circuit, simpler
   ansatz = the product-state restriction, order statistic = the exact top-m at the VQE's
   realised m, second seed). The tenth, "a classical equivalent at matched ENTROPY", is
   included and is the right matching variable.
3. The space is exactly enumerable (8^S <= 262,144), so every classical control is EXACT rather
   than sampled: no search-quality confound can hide in this probe. That is the single best
   design choice in the prereg.
4. Priors are registered in the failing direction for every falsifier (P1 worse by 0.1 to 0.4 A;
   P2 to P5 null), and the probe is 12 targets before the 126 (rule 16) with a GO rule stated.
5. Nothing is tuned on RMSD (section 8), Gamma comes from a native-free median-gap rule, and
   T = 1 is fixed by the posterior's own units with the S25 L2 reason for not calibrating it.
6. The cost is offered to this meter as `s29.s29_X_config:cost_nll` (contract rule 19). Note
   that H_diag is defined on TORSIONS/built chains, so it meters at `--basis ca` on the rebuilt
   CA clouds; the meter will accept it as a callable and the ladder rho it returns is the
   number that decides whether the objective is worth a VQE at all.

FOUR HOLES.
(a) D1's ORACLE BEST IS AN ORDER STATISTIC OVER 8^S AND ITS CONTROL IS NOT MATCHED. "chimera
    ORACLE best vs the DIS top-8's best member": the left side is a minimum over up to 262,144
    structures, the right a minimum over 8. A minimum over more things is smaller whatever the
    space contains (`grid-oracles-are-order-statistics`; S23's -0.077 A "oracle" was 101%
    accounted for by its own best-of-K null). The matched control is a RANDOM RECOMBINATION
    SPACE of the same size and the same parents: the same 8 members, the same S segments, but
    the segment-to-member assignment drawn from a null that destroys structural compatibility
    -- e.g. each parent's segment torsions independently permuted ACROSS segment positions
    (a "scrambled-chimera" space of identical cardinality and identical marginal fragment
    content), ORACLE-minimised the same way. If the scrambled space's ORACLE best matches the
    chimera space's, the 0.3 to 0.8 A is the order statistic, not recombination. Cheap: the
    same enumeration code with a permuted index map.
(b) D1's PRIOR IS NOT FALSIFIABLE AS WORDED. "if it is worth < 0.1 A the space is not richer
    ... and the lane's remaining arms are formalities" names no action: formalities still get
    run and still get quoted. Register the branch: below 0.1 A (against the matched control of
    (a), not against the top-8), the lane reports D1 and STOPS, or states explicitly why the
    remaining arms are worth the compute given the space is no richer than its parents.
(c) P3 IS SET-MATCHED BUT NOT WEIGHT-MATCHED. R2 - R1 on the same tail set is the right
    contrast, but a difference between a uniform average and a p-weighted average over the same
    set is, mechanically, a difference in EFFECTIVE SET SIZE: the weighted average has lower
    participation ratio and is therefore LESS averaged and LESS contracted. Two quantities must
    be printed beside P3 or the sign cannot be read: the weights' participation ratio (1/sum
    w^2) and the emitted structure's Rg or mean virtual bond. `averaging-space-beats-the-objective`
    and S28-L36 both say contraction is the axis on which these comparisons move. The matched
    control is a RANDOM-WEIGHT arm on the same set at the same participation ratio (Dirichlet
    weights fitted to the realised PR, mean of >= 8 draws): if R2 - R1 is reproduced by random
    weights of the same concentration, the quantum stage contributed concentration, not
    information. This is the one control charter section 11 does not name and this comparison
    needs.
(d) 12 TARGETS CANNOT CARRY THE GO RULE THE WAY IT IS WRITTEN. "GO to the 126 iff effect <=
    -0.7 x MDE" on n = 12: the MDE there is computed from 12 paired differences, so 0.7x MDE is
    roughly a 1.4-sigma effect with power near 0.25; the Type-M factor at that power is above
    2, and the prereg's own `ST.compare` will print it. A GO is a decision to spend compute,
    not a claim, so a lenient bar is defensible -- but the entry must print power and Type-M
    beside the GO and must say "GO, not a result" in the same sentence, and the 126-target run
    must re-register its own falsifier rather than inheriting the probe's. Also: with folds
    among 12 targets the fold-clustered CI is over at most 5 clusters of 2 to 3 targets and
    should be quoted as descriptive only at this n (`ST.compare` will supply it; the verdict
    line will read NOT MEASURED for anything smaller than a large effect, which is correct).
TWO SMALLER POINTS.
(e) The PERM control permutes the posterior's rows within sequence separation: good (it
    preserves the separation-distance marginal, which is where most of the distogram's
    information about a peptide sits). State the seed and report the permuted arm's H_diag
    spectrum beside the real one, so "PERM did not reproduce it" is not confounded by PERM
    having a flatter spectrum (the S28-L8b/L11 lesson about same-spectrum controls).
(f) R3 caps at the top-512 configurations by probability. At q = 18 that is 0.2% of the space;
    the entry must report the mass captured, or a "full-state" readout is a top-512 readout.

Verdict: PREREG ACCEPTED, PROBE MAY RUN. D1 and P1 to P5 are all falsifiable and their priors
are registered in the failing direction. No number from this lane is read as a positive until
(a)'s scrambled-chimera control accompanies D1, (c)'s participation-ratio control accompanies
P3, and any GO carries its power, its Type-M factor and the words "GO, not a result". (b) and
(f) are wording. The lane's cost function goes through the meter (S29-L2) before any endpoint
run, as rule 19 requires; the meter's S28-ladder rho is the number to look at, not the charter
ladder's (see the caveat in S29-L2).
Multiplicity: this check runs no comparison. Artefacts: `s29/PREREG_S29_X.md` (read at
2026-09-19 23:56), `s27/LEDGER.md` S28-L21/L41/L43, `s21/LEDGER.md` L14/L17/L18.

## S29-L5 -- SUITE STATUS AT LANE D's FIRST GATE: 17 LIGHT FILES, 381 TESTS, 378 PASSED / 3 SKIPPED / 0 FAILED IN 115 s AT 0.93 GB PEAK; THE HEAVY FILES AND THE TWO AMBER FILES WAIT FOR A QUIET WINDOW (2026-09-19 23:56, D)

Question (brief S29D duty 3): is the suite green after the first S29 commits that touch
`tests/` and `s29/`?
`python s26/jobrun.py --agent S29D --tag TEST --name s29D_pytest_light --est-ram 1.2 -- python
-m pytest <17 files> -q -p no:cacheprovider`: exit 0, wall 115.2 s, PEAK RSS 0.93 GB
(`s26/jobs_done/s29D_pytest_light.json`, log `s26/logs/s29D_pytest_light.log`). 381 tests:
378 passed, 3 skipped (the `VERIFY_SLOW` opt-ins in `tests/test_equivalence.py`), 0 failed,
0 errors. Files: the eight light core files (cvar data energy equivalence geometry instrument
project quantum), the seven S28 lane files, and the two S29 lane files that exist
(`tests/test_s29_D.py` 7 pass, `tests/test_s29_O.py`). The box was at 61% RAM with no other
job running; the governor did not touch it.
DEFERRED, as in S28 (S28-L5): `tests/test_pipeline.py` (forks two workers, tree 1.46 GB),
`tests/test_integration.py` with VERIFY_SLOW=1 inside the command, and the two AMBER files, one
TEST job each -- they go in the quiet window the coordinator announces, not beside a lane's
CPU job. Nothing outside `s29/` and `tests/` has changed on this branch since the S28 close,
so the light-file result carries over from S28-L45's full-tree run (440 passed / 2 skipped)
for every file neither lane has touched.
Standing offer to every lane: if you leave an operator untested I write `tests/test_s29_<lane>.py`
for it (identity parameters reproduce the deployed operator; NaN-poison). Lane O's file exists
and is in this run; lane X's is registered in its prereg and is not yet on disk.
Multiplicity: no comparison. Artefacts: `s26/jobs_done/s29D_pytest_light.json`,
`s26/logs/s29D_pytest_light.log`.

## S29-L6 -- THE METER'S FIRST CUSTOMER: LANE X's PAIR LOG-SCORE IS THE FIRST COST IN THE RECORD THAT DOES NOT ANTI-ORDER THE NEAR-NATIVE LADDER (+0.200 rho ABOVE THE SHIPPED COST, FOLD CI [+0.100, +0.343], 5/5 FOLDS) -- AND IT STILL DOES NOT RECOGNISE NATIVENESS (NATIVE AT THE 37.8th PERCENTILE OF ITS OWN POOL, PRODUCTION PREFERRED TO THE NATIVE ON 62% OF TARGETS); ITS GRADIENT IS UNDEFINED ON 114/126 (2026-09-20 00:01, D)

EVERY NUMBER HERE IS ORACLE. Nothing is deployable, nothing is tuned, and this is a GATE
reading, not a result: it says what lane X's Hamiltonian can and cannot be expected to do
BEFORE its probe spends compute (contract rule 19; my S29-L2's standing offer).
Cost metered: `s29.s29_X_config:cost_nll` (commit 08b5153c), the PAIR half of lane X's H_diag
-- the 17-bin distogram posterior's negative log score summed over pairs, eps 1e-4. The
Ramachandran half is a function of torsions and is not visible from a CA cloud, so the meter
sees the pair term only; lane X's own docstring says so and this entry repeats it. Basis CA,
n = 126, 18.6 s. Artefact `s29/results/s29_D_cost_audit_X_cost_nll_ca.json`; the comparator is
`s29_D_cost_audit_DIS_ca.json` (S29-L2); the contrast block `s29/results/s29_D_x_vs_dis_fmt.txt`.

THE FOUR NUMBERS (shipped cost beside it, same 126 targets, same rungs, same code path):
| | X pair log-score | shipped DIS |
|---|---|---|
| (a) ladder rho, S28 rungs (the binding ladder, pre-specified in S29-L2) | +0.018 [-0.069, +0.123] | -0.182 [-0.308, -0.053] |
| (a) ladder rho, CHARTER rungs | +0.294 [+0.209, +0.343] | +0.260 [+0.201, +0.321] |
| (a) ladder rho, FULL nine rungs | +0.192 [+0.107, +0.266] | +0.118 [+0.064, +0.166] |
| (b) gradient cosine at production | +0.110 on 12 of 126 targets -- NOT A MEASUREMENT (below) | -0.034 on 126 |
| (c) native's percentile in its own pool | 0.378 [0.311, 0.429] | 0.368 [0.306, 0.405] |
| (c) production scores below the native on | 78/126 | 98/126 |
| (d) pref(ORACLE circ_best vs PROD) | 0.341 [0.268, 0.423] | 0.206 [0.113, 0.293] |
| (d) pref(NATIVE vs PROD) | 0.381 | 0.222 |
| (d) pool-member control pct(PROD) | 0.249 | 0.126 |
| (d) contrast vs the pool member | +0.092, 0.88x MDE, fold [+0.060, +0.128], 5/5 | +0.080, 0.83x, fold [+0.006, +0.144], 4/5 |
| (d) head-to-head, the NATIVE below a random pool member | 0.622 | 0.632 |

THE ONE MEASURED DIFFERENCE, `ST.fmt` verbatim (paired per target; POSITIVE = X better; the
word WORSE in the block is `ST.compare`'s RMSD convention and is not the reading here, and for
the same reason its W/L columns are reversed: 57 targets where X's rho is HIGHER, 30 lower,
39 tied):
```
  ladder rho (S28 rungs, CA level): X cost_nll - shipped DIS (POSITIVE = X better; ST's WORSE label is its RMSD convention)
    a 0.0179 (med -0.0526)   b -0.1818 (med -0.2000)   n=126
    effect +0.1998   median +0.0000   SE 0.0490   MDE 0.1374   effect/MDE +1.45
    iid  CI95 [+0.1077, +0.2961]
    fold CI95 [+0.0995, +0.3431]   folds same sign 5/5   per-fold 0:+0.153 1:+0.472 2:+0.205 3:+0.161 4:+0.054
    30W/57L/39T   worst degradation +2.0000 (1MF6)   p90 +0.9737   power 0.98  Type-M 1.01
    concentration: drop-top10 +0.2846 vs uniform-effect null p10/p50/p90 +0.2188/+0.2809/+0.3430 -> pctile 0.530
    VERDICT: WORSE
```
On the CHARTER ladder the difference is 0.43x MDE (NOT MEASURED) and on the FULL nine 0.93x
(NOT MEASURED): the gain is specifically on the near-native half, which is the half that
matters and the one S29-L2 pre-specified as binding. The preference difference is +0.135
[+0.044, +0.252], 1.33x, 4/5 folds; the native-percentile difference is +0.011, 0.25x (NOT
MEASURED) -- the two costs place the native in the same place in the pool.

READING (ORACLE diagnostic; three sentences the lane should carry into its probe).
1. LANE X's COST IS NOT ADVERSARIAL WHERE THE SHIPPED COST IS. The shipped cost's ladder rho
   on the near-native rungs is -0.182 with the fold CI below zero (S29-L2); lane X's is +0.018
   with the fold CI straddling zero, and the paired difference +0.200 clears 1.45x its own MDE
   with 5/5 folds and power 0.98. This is the first cost in this project's record to be
   UNINFORMATIVE rather than ANTI-INFORMATIVE on the ladder from production to the native.
   Mechanism, stated because rule 18 asks: the shipped cost is a BAYES RISK (an expected L1
   distance error under a ~2x over-confident posterior, S25 L2), whose minimiser is a
   contracted structure; the pair log-score is a PROPER SCORING RULE of the same posterior,
   whose minimiser is not driven to contract. Same information, different functional, and the
   functional is worth +0.200 of ladder rho.
2. IT DOES NOT SOLVE RECOGNITION, AND CHARTER FINDING 8 SURVIVES IT. The native sits at the
   37.8th percentile of its own 500-member pool under lane X's cost (the shipped cost: 36.8th;
   difference 0.25x MDE, NOT MEASURED), production scores below the native on 78/126 targets,
   and the cost prefers production to a 0.29 A ORACLE structure on 66% of targets. The
   pool-member contrast (+0.092) is in the Type-M zone at 0.88x MDE, exactly where the shipped
   cost's sits. Lane X's prereg accepted finding 8 as its risk; this is that risk, priced,
   before the probe runs. An objective that ORDERS the ladder better but still ranks the native
   at the 38th percentile cannot, by itself, move the endpoint: it can only stop making things
   worse.
3. THE GRADIENT IS NOT DEFINED FOR THIS COST AND THE +0.110 MUST NOT BE QUOTED. `cost_nll` is a
   bin lookup (`np.digitize` into 17 bins), so 99.2% of its finite-difference gradient
   components are exactly zero at h = 1e-3 A and the cosine is NaN on 114 of 126 targets. The
   +0.110 is computed on the 12 targets where some pair distance happened to sit within 1e-3 A
   of a bin edge -- a selected subsample of the targets nearest a discontinuity, not a
   measurement of the cost's local behaviour. If a local reading is ever wanted for this cost
   it needs A2's smoothed finite difference (h = 0.5 A, labelled "smoothed FD" as S28-L23b did
   for the four step-function channels), and the meter takes `--fd-h 0.5`. This costs lane X
   nothing: it optimises over a discrete space of at most 262,144 configurations and never
   takes a gradient in coordinates.
WHAT THIS DOES NOT SAY. It does not say lane X's probe will help the endpoint (its own prereg's
registered prior is that P1 is WORSE by 0.1 to 0.4 A on the built chain, and I agree with that
prior). It does not license the Ramachandran half, which was not metered. It does not measure
the chimera space, only the cost on S28's nine rungs. And the +0.200 is a rho difference, not
an Angstrom: `averaging-space-beats-the-objective` priced the whole objective channel at 0.171 A
through a fixed readout, and `better-matrix-worse-ranking` is on the record for exactly this
shape of result (five interventions improved a proxy and not the ranking).
Multiplicity: 4 ORACLE diagnostic contrasts in this entry (three ladders and the percentile,
plus the preference), 0 endpoint comparisons; the S28-ladder contrast was pre-specified as the
binding one in S29-L2 before this cost existed, so it is not a best-of-three pick, and the
other two are reported because they are the same statistic on other rungs. The sprint's endpoint
comparison count is unchanged (`s29/STATE.md`: 0).
Verdict: PASSES THE METER AS A BETTER-BEHAVED COST THAN THE INCUMBENT ON THE BINDING LADDER, AND
FAILS IT ON RECOGNITION. Lane X's probe may run on its own pre-registered falsifiers (S29-L4);
this entry is the gate reading rule 19 requires, and its second reading is the one to carry:
the objective is now uninformative rather than adversarial, and uninformative does not fold a
peptide.
Artefacts: `s29/results/s29_D_cost_audit_X_cost_nll_ca.json`, `s29_D_cost_audit_DIS_ca.json`,
`s29/results/s29_D_x_vs_dis_fmt.txt`, code `s29/s29_D_cost_audit.py`, `s29/s29_X_config.py :: cost_nll`.

## S29-L7 -- THEORY SECTION 2: AN OBJECTIVE BUILT FROM THE PER-PAIR MARGINALS IS LOCALLY INFORMATIVE ONLY TO SECOND ORDER AND ONLY ABOUT THE POOL; THE COSINE'S SIGN IS MINUS THE SIGN OF (beta - 1) WHERE beta IS THE POSTERIOR'S OVER-DEVIATION FROM TYPICAL; AND A POSITIVE COSINE IS PURCHASABLE WITH NO INFORMATION BY SHRINKING THE TARGET MAP TOWARD TYPICALITY (2026-09-20 00:01, T)

Question (`s29/briefs/S29T.md` item 2; charter finding 5 and 8): given only the marginals
p_ij, which objectives f(C) have positive expected cosine between their steepest descent at the
production cloud and the direction to the native, under the record's error model (68% common
mode, S23 L9; every predictor emits the typical peptide, S18/S19; score-selected sources have
parallel biases, cos 0.943, S24 L3)? DERIVATION, no experiment; `s29/THEORY.md` section 2.

THE EXACT IDENTITY (2.5), which is the whole of "the gradient is blind" in one scalar product:
with g = grad_D Phi the objective's per-pair coefficients at c, r = Jc u the first-order per-pair
signal toward the native and u = t - c,
    cos(-grad f(c), u) = -<g, r> / (||Jc^T g|| ||u||).
An objective is locally informative iff its per-pair force coefficients are negatively correlated
with production's own signed per-pair error against the native. Nothing about its functional form,
normalisation, temperature, readout or ansatz enters. Two corollaries: (C1) at a true minimiser
grad f = 0 and the cosine is 0/0 -- the shipped cosine is read off a residual (|grad S| 0.0488
A^-1 RMS per atom against |u| 3.048 A, S28-L23b); (C2) grad f = Jc^T g annihilates ker(Jc^T),
dimension P - (3N-6) = 25 of 55 at N = 12, so 45% of pair space is invisible to ANY marginal
objective at the gradient level and a "make the target map embeddable" repair buys exactly zero
there.

THEOREM 2 (assumptions A1 marginal class, A2 proper risk so phi'(median) = 0, A3 coherence of the
pool's and the posterior's deviations from typical, A4 nothing in the system sees the native's
deviation from typical): with d(c) = tau + a, m = tau + b, d(t) = tau + n,
    E[<g, r>] = -sum_ij w kappa var(a) (1 - beta),   beta = cov(a,b)/var(a),
so E[cos] contains NO term in n, is second order in the small quantities, and is exactly zero at
beta = 1 -- when the posterior's median map and the production structure deviate from typical in
the same way and by the same amount, which is what "the average already sits at the posterior's
per-pair median" (S28-L26b) says. Corollary 2a: no function of the marginals, separable or not,
can have an expected cosine whose size is set by n. Corollary 2b: the SIGN is -sign(beta - 1);
the measured -0.034 (and -0.143 on FAIL18) says beta > 1, which is expected because m is an
unconstrained per-pair object and is generally not a realisable distance matrix while c is an
average of real windows. Corollary 2c: non-separability buys nothing (the extra freedom cannot
introduce a term in n, and its component in ker(Jc^T) is annihilated).
Magnitude check (2.8): E[cos] ~ cos_random sqrt(P(1-rho)/2) sigma_a/sigma_n = 0.036 at P = 55,
rho = 0.943 (S24 L3, used as a proxy and labelled one), sigma_a/sigma_n = 0.2 (S12 coord_null:
sequence conditioning is worth 0.776 A of a 3.989 A blind pipeline), against the measured 0.034.
The magnitude match uses two proxies and is not a fit; the sign is a separate statement.

THE WARNING FOR RULE 19, which is the operational output. Shrinking the objective's target map
toward typicality, m -> tau + s(m - tau), drives beta -> s beta and turns the cosine POSITIVE with
zero information added, while moving the emitted structure toward the typical map (production is
already 22% contracted, S23 L1). Meter number 2 is therefore gameable on its own; any candidate
objective that gains cosine must report the implied shrink, or be read beside the native
percentile, which the shrink moves the wrong way.

WHAT A NEW SOURCE MUST BREAK (section 2.5's table): a second differently biased pool, the pool's
dispersion and a joint over the same information all attack A3 or nothing, and A3 is worth a
second-order quantity (at rho = 0.65, the UNSELECTED blind library, |E[cos]| ~ 0.12, still at the
random reference 0.140 -- and that source is 0.76 A worse, S24 L2). Only a channel that sees n
changes the ORDER of the answer: a physics term evaluated on the emitted structure (outside M, so
outside the theorem; measured worse than random, S25 L16) or a learned residual (first order by
construction; S19's error-coherence tax). This is charter finding 8 and S24's prior-derivative
lever (-2.15 A per unit) derived rather than observed.

PREDICTION, three clauses, all ORACLE diagnostics on existing artefacts, checkable by LANE D in
minutes with `s27/results/s28_A2_cosine_rows.jsonl` + `s12/cache/disto_*.npz` + lane O's blind-pool
mean map: (1) the identity holds to relative error < 1e-6 on 126/126; (2) beta > 1 on at least 2/3
of targets and sign(cos_DIS) = -sign(beta - 1) on at least 70%, exceptions concentrated where
|beta - 1| is smallest -- falsified if median beta <= 1 or the sign agreement is inside a coin-toss
CI; (3) re-scoring with the shrunk target map at s in (1.0, 0.75, 0.5) raises the cosine
monotonically, crossing zero near s = 1/median(beta), while the native percentile worsens over the
same grid.
Multiplicity: 0 endpoint comparisons; 3 registered predictions, none yet measured.
Artefacts: `s29/THEORY.md` section 2 (commit d4b52a6c). No job, no native read by me.

## S29-L8 -- TOPIC 2, CORRELATED ERROR IN ENSEMBLES: S23 L9's IDENTITY IS THE KROGH-VEDELSBY AMBIGUITY DECOMPOSITION, SO FINDING 11 IS A LAW AND NOT A DEFECT; THE (1 - 1/M) COEFFICIENT ON THE COVARIANCE MEANS MORE MEMBERS IS WORTH <= ~0.008 A (ARITHMETIC ON S23's OWN NUMBERS); AND EVERY LITERATURE METHOD THAT BREAKS CORRELATED ERROR NEEDS TRAINABLE MEMBERS, A CONTROL WITH A KNOWN MEAN, SAMPLES OF THE TRUTH, OR A KNOWN BIAS RATIO -- THIS INSTRUMENT HAS NONE, WHICH IS ALSO THE MISSING INGREDIENT IN H1 (2026-09-20 00:02, L)

Question (brief topic 2): what does the literature offer for ensembles under common-mode error,
and does any of it break the 68%? Falsifier for the negative half: one method whose stated
inputs this instrument has. No experiment; literature only.

THE LAW IS ALREADY OURS. Brown, Wyatt & Tino, JMLR 6:1621-1650 (2005), eq (10), give Krogh &
Vedelsby (1995): (f_ens - t)^2 = sum_i c_i (f_i - t)^2 - sum_i c_i (f_i - f_ens)^2, the second
term being the AMBIGUITY. With uniform weights and the project's notation (w_k = t + e_k,
c the average, d_k = w_k - c) that is exactly S23 L9's identity mean_k|e_k|^2 = |ebar|^2 +
mean_k|d_k|^2 = 160.36 + 63.82, verified to 2.7e-14 A on 126/126
(`s23/errdecomp.py`, `s23/results/errdecomp.json`). The project's "common-mode / idiosyncratic"
split is the field's "ensemble error / ambiguity" split, re-derived independently. Finding 11 is
therefore a statement about where this pool sits on a universal law, not a pathology, and it is
citable as such.

THE PREDICTION THE EXPECTATION FORM ADDS. Ueda & Nakano (1996), as eq (9) of the same paper:
E{(fbar - t)^2} = bias^2 + (1/M) var + (1 - 1/M) covar. The coefficient on covar does not decay.
With rho = covar/var = 0.676 the reducible factor 1/M + (1 - 1/M) rho runs 1.000 (M=1), 0.689
(25), 0.680 (75), 0.677 (500), 0.676 (limit). Between the shipped top-75 and the whole 500-member
pool it moves 0.0037. DERIVED, and labelled as arithmetic on S23's aggregate numbers rather than
a new measurement: under member exchangeability |c_inf - t|^2 = 160.36 - (63.82*75/74)/75 =
159.50, so an INFINITE pool of the same kind returns 3.0483 * sqrt(159.50/160.36) = 3.040 A on
the point cloud -- more members is worth about 0.008 A, and residual correlation among the d_k
makes it smaller, so read it as an upper bound. Caveat stated in the note: the same conversion
maps the member scale to 3.604 A against the recorded typical-member 3.7037 A, so the
aggregate-to-RMSD map is approximate at ~3%; an exact per-target version is five minutes in lane
O off `s23/results/errdecomp.json` and I run no experiments. This explains three empirical
results at once: S17 (widening K makes the answer worse -- the averaging gain from extra members
is ~0.4% of a member's variance, so displacing good members from the shortlist dominates it),
S24's pool union at +0.0022 A, and `operator-consumes-set-mean`.

WHAT EACH METHOD NEEDS, AND WHY WE HAVE NONE OF IT. (a) Negative correlation learning (Liu & Yao
1999; the diversity error e_i^div = (1/M) sum (1/2)(f_i - t)^2 - kappa (1/M) sum (1/2)(f_i -
fbar)^2, eq 17, with the proven bound lambda_upper = M/(M-1), gamma_upper = M^2/(2(M-1)^2),
eq 39) needs members that are ESTIMATORS BEING TRAINED; ours are retrieved real windows with no
parameter to push. (b) Control variates (Glynn & Szechtman 2002): "suppose that there exists a
random variable Y ... for which EY is KNOWN", lambda* = cov(X,C)/var C, a Hilbert-space
projection onto the span of ZERO-MEAN controls. The known mean here is the native; an uncentred
control does not reduce variance, it moves the estimator by an unknown amount. This is the formal
statement of "a bias shared by every member is invisible to any within-pool statistic". (c)
Multifidelity Monte Carlo (Peherstorfer, Willcox & Gunzburger, SIAM Rev 60:550-591, 2018;
alpha_i* = rho_i sigma_hi/sigma_i, eq 3.12; the efficiency condition sqrt(1-rho_1^2) + sum
sqrt(c_lo/c_hi)(rho_i^2 - rho_{i+1}^2) < 1, eq 3.16) is UNBIASED precisely because it is anchored
by m0 samples of the HIGH-FIDELITY model: it reduces variance around an unbiased anchor and never
removes a biased model's bias. We have no anchor. (d) Boosting/residual fitting needs a residual
learnable from inference-time features with errors decorrelated from the first model's; S24 closed
it, and `error-coherence-decides-correctors` gives the exact violated condition (at identical
0.688 sign accuracy, coherent mistakes emit +0.31 A and i.i.d. mistakes -0.14 A, because a
corrector trained on the predictor's own features inherits its error structure). (e) Recycling
(Jumper et al., Nature 596:583, 2021: 4 iterations, the pair and single representations and
predicted CB fed back, final loss applied at each) is a TRAINING-TIME property -- it would require
a structure-conditioned distogram, which is leakage and unaffordable.

DIRECTLY ON H1. The two-source contrast is Richardson's family, and Richardson cancels the leading
error term only when the two evaluations' biases differ by a KNOWN factor along a SHARED
direction. S24 L2/L3 measured that geometry: bias cosine 0.647 between the blind library source
and the incumbent (31% angularly independent), q = 1.231 (the blind source 0.76 A WORSE),
score-selected mixtures on a line at cos 0.943. A known ratio along a shared direction is exactly
what that is not, and the literature names the failure mode the coordinator's own objection
predicts: the non-parallel 31%, which belongs to the WORSE source, is what an extrapolation
amplifies. The literature does not rescue H1 and does not kill it either; it says the probe must
be priced as a ONE-PARAMETER leave-fold-out fit with its control in the operator's own space.
Separately, Abe et al. (arXiv:2302.00704, ensemble risk = R_avg - Jensen gap, eq 3; ~600
ensembles) is the external form of "diversity-maximising selection: dead": diversity
interventions harm ensembles whose members are already good, and even free diversity carries an
opportunity cost because the best members predict nearly identically.

ONE MORE, WORTH THE SPRINT'S ATTENTION. AF2's pLDDT fits lDDT-Ca = 0.997 pLDDT - 1.17 at Pearson
r = 0.76 across a broad quality range (Jumper 2021) and has NO within-target ranking skill on 588
peptides of 10-40 aa (S29-L1). A signal can be strongly calibrated globally and worthless in-band:
the project's "in-band is the only ranking metric" confirmed on the field's best confidence score.

VERDICT. KEPT (4, as framing or confirmation): the ambiguity decomposition (it is ours); the
bias-variance-covariance form (it supplies the M-law and the 0.008 A bound); Abe et al. on
predictive diversity; the global-vs-in-band contrast from 2.7. REJECTED (5 families, each with
the input it needs and we lack): negative correlation learning (trainable members); control
variates (a control with a known mean); multifidelity (samples of the truth); Richardson-style
two-source extrapolation as stated (a known bias ratio along a shared direction); boosting /
residual correction (decorrelated, feature-learnable residuals). The one place the family could
act here is a PARAMETERISED GENERATOR, where eq (17)'s kappa and eq (39)'s bound would apply --
lane X's configuration space is the only S29 direction with that property.
Multiplicity: 0 endpoint comparisons. The only new number is the derived 0.008 A bound, labelled
as arithmetic on existing artefacts with its caveat.
Artefacts: `s29/lit/L_2_correlated_error.md`; `s29/lit/L_INDEX.md`; inputs `s23/results/errdecomp.json`,
`s27/results/chain_rows.jsonl` (the 3.0483 anchor).

## S29-L9 -- THE DATA-PATH MAP, THE CONVENIENCE-CHOICE LIST, AND THE HARNESS AUDIT: NINE CHECKS, ALL PASS, ONE DECLARED NON-BIT-EXACTNESS THAT MOVES NOTHING; AND THE DISTOGRAM MEMORISES ITS TRAINING PEPTIDES BY 8x (2026-09-20 00:04, M)

Question (charter section 1 items 5 to 7, and section 9's "evaluation-harness failure -- yes, check
this too"): reconstruct the full data path and state where information is created, transformed,
compressed and destroyed; list every choice made for convenience rather than for a stated scientific
reason; and check the evaluation harness itself. Falsifier, pre-stated in the brief: "the harness is
sound" may be written only if every check passes; anything that fails is reported as a failure.

**Deliverable 1, `s29/DATAPATH.md`** -- 13 stages, every function cited by file and line, every
object's shape and unit, the compression factor or entropy count at each step. Measured counts:
183 features per pair (42 base + 141 ESM), 372,881 MLP parameters, per-fold training sets of
6,609 to 6,640 chains of which **90.5% are the fold-SHARED fragment bank** (`fold_fragments` removes
0 to 15 of 6,003 per fold, so leave-fold-out constrains 9.5% of the training data), universes of
7,016 to 39,410 windows (median 17,088), log2 C(17088,500) = 3,252 bits of retrieval choice resolved
by a key whose rank correlation with true RMSD over the universe is +0.066 (`s8/generate.py`
RESULT 2), log2 C(500,75) = 300.6 bits of set choice, and 3n-6 = 21 to 42 cloud degrees of freedom
collapsing to 2n-2 torsions. Two statements every lane should carry: (i) **the 3.2126 anchor never
passes through the quantum stage** (`s27/run_vqe_chain.py:124-130` is the classical top-75;
production `Config.quantum = False`, `core/pipeline.py:179`, and the production cache record carries
`quantum: null`, `n_top: 75`); (ii) the only stage that CREATES three-dimensional coordinates is
retrieval -- everything after it is a selection, an average or a projection of deposited geometry.
Where the 68% common-mode error (S23 L9) enters: stage 3 (the predictor's coherent "typical peptide"
error) and stage 5 (the pool's shared bias); stages 4/6/7 transmit it (S24 L5/L6); stage 10 cannot
remove it, by S23 L9's exact identity. Where it could in principle be observed: only against a
referent that does not share it -- never from a within-pool statistic.

**Deliverable 2, `s29/CONVENIENCE_CHOICES.md`** -- 31 entries, each with file, line, the
alternatives, the ledger line that tested it, or, where nothing tested it, a cheap decisive test.
**17 of 31 are UNTESTED.** The five with the largest reach: C1, the 17 bins and their edges (the
record tested only the CONSUMPTION of those bins -- S25 L6/L12's DEQUANT family is worth 0.003 A
with full leakage -- never the binning itself; the two outer centres, 4.0 and 25.0, are invented by
`CENTRES[0] = edge - 0.5` and `CENTRES[-1] = edge + 2.0`); C2, the soft-bin sigma 0.6 and the
`separation_weights` function that exists and is never called on the deployed path; C6's open half
(every alternative retrieval key ever tested was another native-free QUALITY score, and S17 L12
names the untested class -- a coverage-preserving shortlist -- in its own conclusion); C13, the
top-128 truncation, whose ORACLE ceiling is NOT in the record and is one line to compute; C24, the
medoid frame, the one readout choice with no measurement anywhere and the frame in which S23 L9's
common-mode decomposition is defined. Tested, and therefore not redesign targets: the uniform
readout (S23 L5/L8, S28-L21 -- the p-weighted readouts are +0.24 to +0.33 A WORSE at n = 126),
m = 75 (leave-fold-out, S22 L4/L7), zrank (S25 L17), the ansatz depth and chi (S21 L16), alpha =
0.18 (pinned by a native-free tail-size probe), the L1 functional form (S25 L12), ramah@0.3
(declared post hoc, +0.004 A), the four multi-start conformations (S9-2). **T = 0.5 has no recorded
criterion anywhere in `s25/PREREG_PHYS.md` or the S25 ledger** -- the one harness number with no
stated derivation.

**Deliverable 3, the audit** (`s29/s29_M_harness.py`; job `m_harness_audit3`, exit 0, wall 65.2 s,
peak RSS 0.549 GB, `s26/logs/m_harness_audit3.log`; artefact `s29/results/s29_M_harness_audit.json`,
source sha256 4f10460fce236435):

    [1] ca_rmsd vs an independent 5-line Kabsch written from the definition, 10 targets x 6 pairs
        plus a rigid-motion invariance check: worst |diff| 3.675e-14; a mirrored copy scores
        differently, so reflections are forbidden in both                                    PASS
    [2] native: universe nat_ca vs peptide_db 0.000e+00 A on all 126; vs a FRESH parse of the
        deposited PDB 0.000e+00 A on 10/10                                                   PASS
    [3] benchmark leak: a static grep over the whole import closure hits only core/data.py's two
        path CONSTANTS and one docstring; an in-process POISON (builtins.open, os.path.exists and
        np.load all raise on the sealed artefacts) live through 6 full target rebuilds and 126
        cold distograms: 0 forbidden reads                                                   PASS
    [4] ST.pinned_folds == the universes == peptide_folds.json on all 126; 10 targets served by
        fold_model(own fold) with their own sequence absent from that fold's training set     PASS
    [5] anchor: chain_rows DIS chain 3.2126 / cloud 3.0483 at n = 126; a FRESH re-projection on 6
        targets reproduces both at worst 0.000e+00 A                                          PASS
    [6] MDE: the literal MDE_K = 2.8016 (`s24/stats_lib.py:56`); the fold CI resamples FOLDS with
        replacement (lines 119-127), 4,000 draws, beside a separate iid bootstrap; reproduced
        SE 0.026133 and MDE 0.073215 = 2.8016 x SE exactly                                    PASS
    [7] 126 targets, 126 distinct pdb ids, 126 distinct sequences, lengths 9..16, 126 universe
        files, pdb-sorted (histogram 9:9 10:11 11:13 12:19 13:23 14:14 15:16 16:21)           PASS
    [8] cold distogram on ALL 126 with the poison live: in-process determinism exactly 0.0e+00,
        but the s12 cache is NOT bit-identical to a fresh recomputation -- max|prob| 1.97e-06,
        max|risk| 1.34e-04 (relative 5.63e-06), max|score| 1.91e-06 against a score sd of ~1.25.
        Consequence MEASURED rather than assumed: the 500-candidate score ORDER differs on 2/126,
        the top-75 SET on 0/126, and the emitted point-cloud RMSD by 0.00e+00 A on 126/126    PASS
    [9] EMPIRICAL leave-fold-out (ORACLE DIAGNOSTIC; reads the native distances, tunes nothing):
        mean NLL of the DEPLOYED own-fold model 3.28777 against 1.21275 for the four models that
        SAW that fold; delta +2.07502, SE 0.15166 = 4.88x MDE, correct sign on all five folds
        (per-fold deployed NLL 3.064 / 3.125 / 3.621 / 3.472 / 3.181)                         PASS

**Verdict: the harness is sound.** The check-8 non-bit-exactness is declared here rather than
buried: an S29 arm whose operator reads the score's ORDER below the top-75 cut -- a full-pool rank
correlation, a reranker's tail, an order-statistic null -- should recompute the posterior instead of
reading the s12 cache, or accept a 2-in-126 chance that its order differs from another lane's. No
result in the record is affected (the SET and the emitted RMSD are identical on 126/126).

**Two corrections to the brief.** `s25/qcand_lib.py` does not exist -- the `Encoding` class is
`s22/qcand_lib.py:120` (an older one is at `s14/vqe_encoding.py:50`). And `catrace_prior.npz` is not
a native store: it holds one key, `nlp`, a (24, 36) pseudo-angle log-density table
(`core/geometry.py:884`); the pinned native is the universe's `nat_ca`.

**The finding that falls out of check 9, and it is for the other lanes.** The models that trained on
a peptide assign its true distance bins **e^2.075 = 8.0x** more probability than the held-out model
that actually ships. The distogram memorises its training peptides heavily (372,881 parameters
against 787 peptides plus 6,003 fragments, dropout deliberately 0, `core/predict.py:335`, where the
docstring records that regularisation improves the distance matrix and WORSENS the ranking). Two
consequences: (i) any in-sample diagnostic computed on the corpus -- calibration, sharpness, MAE, a
fitted residual correction -- describes memorisation and must be recomputed out-of-fold before it
means anything about the deployed model; (ii) the sequence channel the record prices is the
3.288-nat held-out model, not the 1.213-nat in-sample one.

Comparisons (multiplicity): **0 endpoint comparisons.** This entry contains no deployable contrast;
check 9 is an ORACLE diagnostic of the instrument, not of an arm.
Artefacts: `s29/DATAPATH.md`, `s29/CONVENIENCE_CHOICES.md`, `s29/s29_M_harness_audit.md`,
`s29/s29_M_harness.py`, `s29/results/s29_M_harness_audit.json`, `s26/logs/m_harness_audit3.log`.

## S29-L10 -- CONTRACT ADDENDUM 20 IS NOW MECHANICAL IN THE METER: EVERY COSINE CARRIES ITS SHRINK SIGNATURE, AND THE SHIPPED COST's DESCENT DIRECTION EXPANDS (BOND x1.044, Rg x1.025, ONLY 18/126 CONTRACT); PLUS TWO DEFECTS IN MY OWN METER, FOUND AND FIXED WITH REGRESSION TESTS (2026-09-20 00:06, D)

ORACLE diagnostic throughout; nothing deployable; the meter still tunes nothing.
Question: the coordinator's addendum 20 (from lane T's theorem 2, S29-L7) says a cosine gain is
purchasable with zero information by shrinking an objective's target map toward typicality, and
makes lane D veto any cosine gain that does not report (a) the implied shrink, (b) the native
percentile and (c) the emitted geometry. Can the meter make that mechanical rather than a
reviewer's memory?

1. THE SHRINK SIGNATURE IS NOW PRINTED UNDER EVERY COSINE (`s29/s29_D_cost_audit.py`, SHRINK_E
   = 0.3 A). The implied shrink of an arbitrary cost's target map is not computable in general
   (most costs have no explicit target map), but its SIGNATURE is, and it is the quantity
   theorem 2 actually points at: descend 0.3 A along -grad f from production and report what
   happened to the mean virtual CA-CA bond and to Rg, with the count of targets that contract.
   A shrink-bought cosine moves the emitted structure toward the typical map, so its ratios sit
   below 1; information does not have to.
   THE SHIPPED COST's BASELINE, measured (n = 126, the number every candidate is read against):
   descending 0.3 A along -grad S~ multiplies the mean bond by 1.0438 and Rg by 1.0248, and only
   18 of 126 targets contract. So the shipped cost's descent direction EXPANDS a structure that
   is already 22% contracted (S23 L1; production's mean bond 2.96 A against the native's 3.81 A,
   S28-L35). Read with S28-L30 (S~ falls on 126/126 at e = 0.1 while the RMSD rises) this says
   the shipped objective is NOT sitting at a contraction optimum and its blindness is not a
   contraction artefact -- which is worth knowing, because "the cost just likes contracted
   things" was the obvious explanation of the -0.034 and it is now measured false at the
   gradient level. (It remains true of the cost's RANKING: production beats 87% of real traces
   under it, S28-L36.)
   Lane X's pair log-score, for comparison: bond x1.0232, Rg x1.0065, 3 of 12 targets contract
   (on the 12 targets where its gradient exists at all -- see 3 below), so its cosine is not
   shrink-bought either. Both costs are on the expanding side; neither buys its cosine the cheap
   way. Addendum 20's veto has nothing to fire on yet, and now it will fire automatically when
   it does.
2. DEFECT 1 IN MY OWN METER, FIXED: an all-NaN axis crashed the renderer. Metering CAGEO (a
   step-function channel) produced a zero gradient on EVERY target, `_summ` returned
   `ci95_fold = None`, and `render` raised `TypeError: 'NoneType' object is not subscriptable` --
   after the run had completed, so the 126-target result was lost to a print statement. Fixed:
   every summary key is always present, the cosine block prints "UNDEFINED on k of n targets,
   100% zero gradient components, re-run with --fd-h 0.5 for A2's smoothed finite difference",
   and the preference rows render with or without a CI. Regression test
   `tests/test_s29_D.py::test_all_nan_axis_does_not_crash_the_render`.
3. DEFECT 2, FIXED IN THE SAME PLACE: a PARTIALLY defined cosine was being reported as if it
   were a measurement. Lane X's cost has a defined gradient on 12 of 126 targets (99.2% of
   components are exactly zero at h = 1e-3 A); the meter printed "+0.110" with a fold CI and no
   warning, and S29-L6 had to add the caveat by hand. The meter now prints, in the header line
   itself, "UNDEFINED on 114 of 126 targets (zero gradient); the mean below is on the survivors
   and is a SELECTED SUBSAMPLE, not a measurement" -- the survivors are exactly the targets
   whose pair distances sit within h of a bin edge, which is a selection on proximity to a
   discontinuity, not a random subsample. No number of S29-L6 changes; its caveat is now
   emitted by the code instead of remembered by me.
4. `Spearman(cos, production RMSD)` is now computed on the finite subset (it was NaN whenever
   any cosine was NaN); for the shipped cost it is unchanged at -0.372 (S28-L23b).
Suite: `tests/test_s29_D.py` 9 pass (7 + the two regressions). The three meter artefacts and
lane X's were regenerated so every stored summary carries the signature; every number in
S29-L2 and S29-L6 is unchanged by the fix (the shipped cost's four anchors re-verified by the
same test file after the patch).
Multiplicity: no comparison in this entry. Artefacts: `s29/s29_D_cost_audit.py`,
`tests/test_s29_D.py`, `s29/results/s29_D_cost_audit_{DIS_ca,DIS_SURR_ca,DIS_chain-s28rows,X_cost_nll_ca,CAGEO_ca}.json`.
Verdict: addendum 20 is mechanical. Any lane reporting a cosine gain gets the shrink signature
printed beside it whether or not it asks, and a gain with ratios below 1 and an unmoved native
percentile will be vetoed on sight, as the coordinator's rule says.

## S29-L11 -- THEORY SECTION 3, THE SPECTRAL CONDITION: GRADIENT VARIANCE IS r_stable/D^2 AND NOTHING ELSE, WHICH PREDICTS S28's GAUSSIAN (7%), S28-B2's kNN (46x) AND LANE D's J* = 85.7 (88) WITH NO FREE PARAMETER; MEASURED ON 12 REAL POOLS, CENTERING RAISES lambda_2/lambda_1 FROM 0.138 TO 0.465 AND MAKES THE DECAY WORSE, NOT FLAT (-2.30 PER QUBIT AGAINST -1.83); AND A CENTERED HAMILTONIAN CONSUMED BY ANY p-READOUT IS SELF-CANCELLING BY THE POLE SYMMETRY (2026-09-20 00:06, T)

Question (`s29/briefs/S29T.md` item 3; charter finding 7; the coordinator's integration note 2,
which spawned lane B on the strength of a spectral number): (a) where does lambda_2/lambda_1 =
0.11 to 0.14 come from, (b) what is the gradient variance of <psi|A|psi> as a function of A's
spectrum, (c) what are the spectra of the double-centered graph A_c and the signed agreement
matrix G = Delta Delta^T on real pools and is the decay then flat in n, (d) what IS the ground
state of diag(E) - J A_c. Derivation plus one governed property job; no native, no RMSD, nothing
deployable. `s29/THEORY.md` section 3.

THE JOB. `s29/s29_T_spectra.py` -> `s29/results/s29_T_spectra.json`, `s29_T_spectra_rows.jsonl`
(72 spectra cells), `s29_T_grad_rows.jsonl` (216 gradient cells); `s26/jobs_done/s29T_spectra.json`
exit 0, wall 100.2 s, peak RSS 0.337 GB; 12 targets (S27 T11's set), n = 4..9, 120 theta draws,
seed 1009, depth 3, the S28 lane B law exactly. ANCHOR: my hop-only variances for A reproduce
`s27/results/s28_B_train.json :: hop_only|J1` EXACTLY at n = 4..8 (4.0612e-03, 7.2844e-04,
2.8167e-04, 5.8983e-05, 3.9483e-05) and to 7% at n = 9 (4.4466e-06 vs 4.1567e-06; S28's n = 9
register orders the 500 candidates differently under the same graph).

(a) lambda_2/lambda_1 ~= f_1 <d^2>/(2 sigma^2) where f_1 is the share of the pool's shape
variance in its top mode: measured f_1 = 0.324, <d^2>/2sigma^2 = 0.534, predicted 0.173 against
measured 0.138 (median over 12, range [0.104, 0.211]; S28-L8b's 0.138 reproduced). The ratio is a
property of the POOL (about 6 effective shape modes, the leading one a third of the variance), not
of the kernel: any smooth decreasing kernel gives the same to first order, and widening sigma
flattens the spectrum only by making the graph a constant plus noise.

(b) THE LAW. dF/dtheta_k = <psi|[A, Gt_k]|psi> exactly (Gt_k real antisymmetric, Gt_k^2 = -I/4),
so for a 2-design state Var[dF/dtheta] ~= tr(A_0^2)/D^2 = r_stable(A)/D^2 at unit spectral norm,
r_stable = ||A||_F^2/||A||_2^2. The variance depends on A ONLY through its Frobenius norm, which
is basis-independent -- so a diagonal matrix with the same spectrum must decay identically, which
is S28-L11 caveat (b) derived rather than observed. THREE CHECKS, NO FREE PARAMETER: S28's
Gaussian graph, r_stable 1.036 -> 3.95e-6 against measured 4.16e-6 (7%); at n = 4, 4.36e-3 against
4.06e-3 (7%). S28-B2's kNN graph, for which the entry recorded a mechanism-free "30 to 60x": a
k-regular graph has r_stable = M/k, so Var = 500/(10 * 512^2) = 1.9e-4 = 46x the Gaussian's.
Lane D's J* = 85.7 (S28-L11 item 1): J* = D sqrt(Var_diag/r_stable) = 512 * 0.171 = 88.
DESIGN RULE: an off-diagonal term is gradient-visible at the deployed width iff its stable rank
grows with the register (parity at n = 9 needs r_stable ~ 16 D ~ 8000). No dense kernel can,
centered or not; a k-regular graph reaches M/k; and ANY Gram matrix of structural deviations is
capped by its rank, r_stable <= 3 N_res - 6 <= 42 on this instrument.

(c) MEASURED, AND IT CORRECTS THE READING IN `s29/STATE.md` INTEGRATION NOTE 2. At n = 9, median
over 12 pools: lambda_2/lambda_1 = 0.138 (A) / 0.465 (A_c) / 0.634 (G); r_stable = 1.036 / 1.591 /
1.675; top-eigenvector uniform overlap 0.969 / 0.000 / 0.007. So the note is right that centering
removes the degeneracy in the lambda_2/lambda_1 sense and that S28 closed ONE similarity measure
rather than the class. BUT THE TRAINABILITY HALF OF THE READING IS WRONG: the hop-only gradient
variance decays at -2.305 per qubit for A_c and -1.900 for G against A's -1.830 (predicted by
(3.4): -2.51, -1.96, -2.02), and at n = 9 A_c is 3.599e-6 and G is 5.051e-6 against A's 4.447e-6.
Centering makes the decay WORSE, not flat, because r_stable rises only 1.04 -> 1.6; a matrix with
lambda_2/lambda_1 = 0.63 and a fast tail still has ||A||_F^2 = O(1) at unit spectral norm. For
lane B: J* falls from 88 to 71, a 19% reduction, not an order of magnitude, and the Gram rank cap
means no re-weighting fixes it. A build justified by "the spectrum is no longer degenerate" is
justified by the wrong number; what is genuinely new in A_c and G is the MEANING of the ground
state, (d).

(d) THE GROUND STATE IS THE POOL'S PRINCIPAL CONTRAST, AND A p-READOUT CANCELS IT. <v|G|v> =
|sum_i v_i delta_i|^2/N_res, so the top eigenvector is exactly the signed combination of members
whose deviations from the pool mean add to the largest displacement: positive on one pole of the
pool's first shape mode, negative on the other, bimodal in probability. Every deployed readout
(R1 tail, R2 p-weighted, R3 p-top-m) is a function of p_i = psi_i^2 and is blind to the sign, so
sum_i v_i^2 W_i = cbar + delta (mass on the + pole minus mass on the - pole) + O(eps), and
centering makes that imbalance small BY CONSTRUCTION. A centered off-diagonal Hamiltonian
consumed by any probability readout is self-cancelling to first order; the residual imbalance is
whatever diag(E) contributes, which is no new information (section 2). Lane B must therefore use
a SIGNED readout (S28 lane A's, refuted for accuracy under the shipped OBJECTIVE, not as a
readout) or break the pole symmetry explicitly. And under the signed readout the family collapses
to one parameter: the ground state emits production +- eta PC1(pool). Not the S27 section-6
consistency mechanism (that selects the pool's MODE and moves the set mean toward the centre; this
selects the extremes of a contrast and moves it away), but the same law decides both -- the
terminal operator consumes the retained set's mean, so everything rests on the SIGN, which
section 2 says the marginals do not supply.

PREDICTIONS (four, all checkable in minutes): (1) LANE B -- on A_c and G the hop-only slope is
-2.3 +- 0.3 and -1.9 +- 0.3, NOT flat, and both n = 9 variances are within 1.5x of the Gaussian's;
falsified if either slope is above -1.3 or either n = 9 variance exceeds 3e-5; reproducible in
100 s with `python s29/s29_T_spectra.py --grad`. (2) LANE B or D -- for any unit-spectral-norm
observable, Var[dF/dtheta_0] = r_stable/D^2 within 2x at n >= 7, J* = D sqrt(0.0305/r_stable);
for G at n = 9, Var 6.4e-6 and J* 71; falsified by any observable at n >= 7 departing by over 3x.
(3) LANE B -- with a p-readout the J -> infinity ground state of diag(E) - J A_c emits the pool
mean to within the projection floor on >= 80% of targets. (4) LANE O, minutes -- the whole
family's ORACLE ceiling is min_eta mean RMSD(c + eta PC1(pool), native) at the best global eta
chosen leave-fold-out: predicted UNDER 0.15 A better than production, per-target sign priced by
`best_of_k_within`; if it is above 0.30 A the prediction fails and the B build is worth much more
than I think.
Multiplicity: 0 endpoint comparisons; 4 registered predictions; the 18 (matrix, n) gradient cells
are a property measurement, not a contrast.
Artefacts: `s29/THEORY.md` section 3 (commit 70d45dec); `s29/results/s29_T_spectra.json`,
`s29_T_spectra_rows.jsonl`, `s29_T_grad_rows.jsonl`; `s26/jobs_done/s29T_spectra.json`;
`s26/logs/s29T_spectra.log`.

## S29-L12 -- TOPIC 3, DECISION THEORY: THE PERCEPTION-DISTORTION THEOREM (BLAU & MICHAELI 2018) MAKES CHARTER FINDING 8 A NECESSITY -- FOR *ANY* DISTORTION MEASURE THE DISTORTION-OPTIMAL ESTIMATOR'S OUTPUT DISTRIBUTION MUST DIVERGE FROM THE REAL ONE, MOST STEEPLY AT LOW DISTORTION, SO EVERY REALISM-TYPE SCORER (WHICH IS EVERY SCORER IN THE LIBRARY) MUST DISPREFER THE RMSD-OPTIMAL ANSWER; THE 25.8% CONTRACTION IS JENSEN'S INEQUALITY; AND POST-HOC CALIBRATION IS CLOSED BY ARGMAX/MEDIAN INVARIANCE (2026-09-20 00:07, L)

Question (brief topic 3): what does decision theory say about structure point estimates (L1/L2 in
coordinate vs distance space, medoid vs mean under multimodality), about proper scoring rules for
distograms and about calibrating an over-confident distance posterior -- and is there a result
that the Bayes estimator of an over-confident posterior is contracted? No experiment; literature
only. Falsifier for the framing claim: a construction that attains both minimal distortion and
the correct output distribution under a non-invertible degradation.

THE RESULT. Blau Y, Michaeli T, "The Perception-Distortion Tradeoff", CVPR 2018 (arXiv:1711.06077).
Distortion is E[delta(X, Xhat)] for any per-sample measure; perceptual quality is a divergence
d(p_X, p_Xhat) between the estimator's OUTPUT DISTRIBUTION and the distribution of real signals
(with total variation, d is the best achievable probability of telling an output from a real
signal). P(D) = min d(p_X, p_Xhat) s.t. E[Delta] <= D (their eq 11). THEOREM 3, verbatim: "If
d(p,q) of (4) is convex in its second argument, then the perception-distortion function P(D) of
(11) is 1) monotonically non-increasing; 2) convex", and "Theorem 3 requires no assumptions on
the distortion measure Delta. This implies that a tradeoff between perceptual quality and
distortion exists for any distortion measure". Convexity puts the steepest part of the tradeoff
at the LOW-DISTORTION end, which is where this project operates. THEOREM 1: if the degradation is
non-invertible and the distortion-minimising estimator is unique, that distortion measure is not
stably distribution-preserving. The mechanism, verbatim, is ours: "the average of valid images is
not necessarily a valid image, so that the MMSE estimate frequently 'falls off' the natural image
manifold".

WHY THIS IS FINDING 8. Our distortion measure is CA-RMSD (the endpoint). Every native-free scorer
in the S27 library is, up to a monotone map, a log-density of real structures -- a realism
measure, i.e. an estimate of the d-axis. Theorem 3 then says the RMSD-optimal answer is precisely
the answer a realism scorer is most confident is fake. S28-L48 is that statement measured: 20 of
31 scorers prefer the projected production average to a 0.25 A ORACLE structure, and CAGEO
"prefers any real trace to a contracted one" (S28-L36). The scorers are not broken; they read the
correct axis, and the endpoint rewards the other one. The project's own two operator families are
the two ends of P(D): the coordinate average is the low-distortion off-manifold end (contracted
25.8%, cloud 3.0483) and the medoid/pool member is the on-manifold higher-distortion end (which
is why realism prefers it, and why consensus is capped at the pool's mode). The built chain's
projection is a MEASURED point on the curve: +0.164 A to return to ideal geometry (3.0483 cloud
-> 3.2126 chain, `s27/results/chain_rows.jsonl`), so the coordinator's wave-2 probe P1 is a
measurement of P(D)'s local slope and should be framed that way -- interpretable whichever way it
comes out.

WHAT IT DOES NOT LICENSE. The theorem constrains an estimator's output DISTRIBUTION, not the
per-target ranking of a fixed candidate set: ranking two structures of MATCHED realism by accuracy
is not forbidden. So the one version of finding 8's question that is not answered "no" in advance
is IN-BAND-BY-REALISM ranking -- match plausibility first, then rank. The pool-member control
(S28-L36/L37) is a crude form of exactly that, and the same paper supplies the converse that makes
the control mandatory: "we could achieve perfect perceptual quality by randomly drawing natural
images that have nothing to do with the original ground-truth images. In this case the distortion
would be quite large."

THE CONTRACTION IS JENSEN, NOT A BIAS. For any X, Y: ||E[X] - E[Y]|| = ||E[X-Y]|| <= E[||X-Y||],
strict unless X-Y is a.s. a fixed direction. So a coordinate average shrinks every interatomic
distance, by more when the pool disagrees more. The measured 25.8% backbone contraction
(`averaging-space-beats-the-objective`) needs no distogram-bias explanation and cannot be
reweighted away -- only projected (the +0.164 A above) or avoided by not averaging in coordinate
space. This is also the formal sibling of the shrink-toward-typicality loophole behind contract
rule 20: both improve a distance-space objective while degrading realism.

CALIBRATION IS CLOSED FOR A MECHANISM REASON. Guo et al., ICML 2017 (arXiv:1706.04599): temperature
scaling divides logits by one scalar fitted by NLL and "does not alter which class receives the
highest predicted probability (argmax) or change classification accuracy". For a per-pair L1 Bayes
readout it is also MEDIAN-invariant whenever the per-pair posterior is symmetric about its centre,
so a pure temperature change can act only through (i) the 17-bin asymmetry and (ii) the
discretisation (S25 L7: 17 distinct per-pair target values, gaps to 4 A). Both are artefacts, so
S25 L2's "the posterior is ~2x over-confident and calibrating it makes RMSD WORSE" is now explained
rather than merely observed. Note for the record: this ANSWERS the brief's question "is there a
result that the Bayes estimator of an over-confident posterior is contracted?" -- NO, not from
over-confidence per se. A symmetric width error does not move the median or the mean. The
contraction in this system comes from AVERAGING (Jensen) and from shrinkage toward the prior's
centre (the typicality axis), not from miscalibrated width.

PROPER SCORING RULES. Gneiting & Raftery, JASA 102(477):359-378 (2007): S is proper if honest
reporting minimises expected score. The distogram is trained by 17-bin cross-entropy, which IS the
logarithmic score and IS strictly proper -- the training objective is already decision-theoretically
correct, and S25's 51 arms of re-reading gave corr(progress-toward-truth, endpoint) = +0.054. What
is worth importing is the diagnostic: CRPS is the distance-sensitive generalisation of absolute
error to predictive distributions, where "the logarithmic score assigns harsh penalties regardless
of closeness". Any future comparison of two posteriors should use CRPS, never MAE -- which the
project reached empirically twice (`prior-mae-prices-selected-rmsd` r = 0.19;
`error-shape-not-mae-decides-ranking`).

ROUTE NOTE. AlphaFold1 (Senior et al., Nature 577:706, 2020) minimised a potential fitted to the
distance histogram; AlphaFold2 (Jumper et al., Nature 596:583, 2021) abandoned the distogram
potential for a structure module and predicts its own accuracy (lDDT-Ca = 0.997 pLDDT - 1.17,
r = 0.76). Neither post-hoc calibrates the distogram. The field's answer to an over-confident
per-pair posterior was to STOP TAKING A PER-PAIR POINT ESTIMATE. The project's analogue of
"consume the posterior jointly" is lane X's configuration space.

VERDICT. KEPT (5, as framing/diagnostic/route, none importable as an operator): the classical
Bayes estimators plus the Jensen contraction; Blau & Michaeli Theorems 1 and 3; their converse
(realism alone buys nothing = the pool-member control's theoretical form); Gneiting & Raftery's
CRPS as the posterior-comparison metric; the AlphaFold route note. REJECTED (2 levers): post-hoc
temperature calibration of the distogram (argmax- and median-invariant, no information added;
independently closed by S25 L2) and changing the training scoring rule (the log score is already
strictly proper). Multiplicity: 0 endpoint comparisons; no new measurement in this entry.
Artefacts: `s29/lit/L_3_decision_theory.md`; `s29/lit/L_INDEX.md`.

## S29-L13 -- TOPIC 4, QUANTUM: THE SET-EQUALITY THEOREM IS EQ (12) OF BARKOUTSOS ET AL. (CVaR IS DEFINED ON SORTED SAMPLES, FOR DIAGONAL HAMILTONIANS), SO IT IS A DEFINITION AND NOT A DISCOVERY; NON-CLASSICALITY NEEDS BOTH NON-COMMUTING TERMS *AND* A TARGET THAT IS NOT AN EIGENVECTOR (GIBBS/THERMAL, QUANTUM BOLTZMANN MACHINES); CVaR's GLOBAL OPTIMUM IS THE FLAT SET {overlap >= alpha}, A NEW CHEAP FALSIFIER FOR OUR OWN QUANTUM ARMS; AND arXiv:2312.09121 MAKES "IT TRAINS" EVIDENCE *FOR* CLASSICAL SIMULABILITY (2026-09-20 00:13, L)

Question (brief topic 4; charter section 6 group 1 and section 11): under what condition does the
quantum object differ from its classical counterpart, what Hamiltonian structure would make the
variational state not reproducible by a sort or an eigensolver, and which papers show a
variational state doing something an eigensolver on the same H cannot? No experiment; literature
only. Nine papers read for their constructions.

1. CVaR IS CLASSICAL BY DEFINITION ON A DIAGONAL H. Barkoutsos, Nannicini, Robert, Tavernelli,
Woerner, Quantum 4:256 (2020), arXiv:1907.04769: "CVaR_alpha(X) = E[X | X <= F_X^{-1}(alpha)]"
(eq 11) and, for samples sorted nondecreasing, CVaR_alpha = (1/ceil(alpha K)) sum_{k=0}^{ceil(alpha
K)} H_k (eq 12); "the limit alpha -> 0 corresponds to the minimum, and alpha = 1 corresponds to
the expected value". The abstract scopes the method: "In the case of classical optimization
problems, which yield diagonal Hamiltonians, we argue that aggregating the samples in a different
way than the expected value is more natural." The project's set-equality theorem (S25: the CVaR
tail's support is always a subset of an initial prefix of the energy order; 2,592 adversarial
cells, 0 violations) is therefore eq (12) restated, not a property of our implementation. RECORD
CORRECTION OF EMPHASIS, not of fact: cite Barkoutsos eq (12) beside it. On a diagonal H the only
quantum content is the sampling distribution the trial state induces.

2. A NEW CHEAP FALSIFIER FOR OUR OWN ARMS. The same paper: "for any problem (1) and parameters
theta* such that |psi(theta*)> has overlap rho > 0 with the ground state, theta* is a global
minimum of CVaR_alpha(X(theta*)) for alpha <= rho". So CVaR_alpha's global-minimiser set is
{theta : overlap with the best candidate >= alpha} -- large and flat. With an averaging readout
the objective is indifferent to precisely the freedom the readout consumes (which OTHER
candidates populate the tail), which is a candidate mechanism for the sprint's central puzzle:
the optimiser reaches the optimum on 126/126 while the emitted structure does not move (S28-L18b,
S28-L26b). I have not seen this stated in the record. CONSEQUENCE, offered to lane D as a meter
clause and to lane T to check: any accuracy change attributed to CVaR optimisation must be shown
NOT to be a tie-break within that flat set. Also recorded: the empirical CVaR estimator's variance
is O(1/(K alpha^2)), so the standard error grows as 1/alpha and matching the expectation's
accuracy needs K/alpha samples -- not binding on an exact statevector, but the right price to
quote if any arm ever samples. And Proposition 5.1: local minima of the expectation and of CVaR do
not map to each other (their two-qubit example has a constant expectation and an informative
CVaR_0.5 = sin^2(theta/2)) -- the paper's real argument for CVaR is LANDSCAPE repair, which is not
this project's binding problem since our objective already falls on 126/126.

3. THE CONDITION FOR NON-CLASSICALITY. The classical counterpart of the current formulation is a
SORT (H diagonal, CVaR reads a prefix). It becomes an EIGENSOLVER as soon as H is non-diagonal but
the target is its ground state. So a formulation is non-classical in the relevant sense only if
BOTH (C1) the Hamiltonian's terms do not commute, so the eigenbasis is not the computational basis,
AND (C2) the prepared object is NOT an eigenvector, so diagonalisation is also not the counterpart.
(C2) is satisfied by a thermal/Gibbs state e^{-beta H}/Z, by a state whose role is to be a sampling
distribution, and by any FREE-energy objective (energy minus entropy). The project has never used
(C2): S28's non-diagonal Hamiltonian satisfied (C1) with a near-rank-one (degenerate) off-diagonal
(S28-L8b/L11), and everything else has been diagonal. Suggestive from our own record: S25 found the
trained state sits 0.902 nats / 45% of its mass from its own analytic Gibbs optimum and that RMSD
tracks READOUT ENTROPY (rho -0.74) rather than alpha or T -- the entropy is where the output lives.
TWO CAVEATS SO NOBODY OVER-READS IT: for a diagonal H the Gibbs state is a classical Boltzmann
distribution over candidates, and S21 exhaustively enumerated the 2^n latent on 75/126 targets, so
(C2) alone buys nothing; and see item 4.

4. THE FRAMING RESULT, AND IT CUTS AGAINST US. Cerezo, Larocca, Garcia-Martin, Diaz, Braccia,
Fontana, Rudolph, Bermejo, Ijaz, Thanasilp, Anschuetz, Holmes, "Does provable absence of barren
plateaus imply classical simulability?", Nat Commun 16:7907 (2025), arXiv:2312.09121, verbatim:
"many commonly used models whose loss landscapes avoid barren plateaus can also admit classical
simulation, provided that one can collect some classical data from quantum devices during an
initial data acquisition phase ... barren plateaus result from a curse of dimensionality, and ...
current approaches for solving them end up encoding the problem into some small, classically
simulable, subspaces." The project's one genuine quantum positive -- the optimiser trains, beats
best-of-200 from the untrained circuit, closes 78-89% of the free-energy gap, no barren plateau at
any measured width (S13/S25) -- sits squarely in the regime this paper warns about. It does NOT
retract that result (a correctly scoped trainability claim; rule 9 already forbids treating the
simulator as a cause). It changes how any S29 quantum claim must be positioned: "our circuit
trains well" is now evidence FOR classical simulability, so the charter's first control (a
classical equivalent) is the central test rather than a formality, and the reachable claim is "a
specific quantity the matched classical control does not reproduce", not "classically impossible".

5. THE CHEAPEST OPEN QUANTUM QUESTION IN THE PROJECT. Larocca, Czarnik, Sharma, Muraleedharan,
Coles, Cerezo, Quantum 6:824 (2022), arXiv:2105.14377, diagnose trainability and reachability
through the DYNAMICAL LIE ALGEBRA g generated by the ansatz's generators under commutation (the
reachable set is exp(g)); Ragone et al., Nat Commun 15 (2024), arXiv:2309.09342, give "an exact
expression for the variance of the loss function of sufficiently deep parametrized quantum
circuits", resolving "a standing conjecture about a connection between loss concentration and the
dimension of the Lie algebra of the circuit's generators" -- the dependence is INVERSE in dim(g)
(sourcing caveat: I read the abstract and surrounding text, not the displayed equation; anyone
using the exact P_g(rho) P_g(O)/dim(g) form should pull it from the paper). `docs/STATE_BRIEF` 5.7
item 4 lists the deployed ansatz's DLA as NOT MEASURED. It is classical linear algebra on the
generators, costs minutes, is not an endpoint experiment, and answers charter question 5 ("what
the ansatz can represent, and what it provably cannot"). RECOMMENDED to lane T or X; I do not run
it.

6. THE REST, BRIEFLY. Cerezo et al., Nat Commun 12:1791 (2021): global observables give
exponentially vanishing gradients even at shallow depth while local ones give at worst polynomial
decay at depth O(log n) -- assumes blocks forming local 2-designs, which S13 says we are nowhere
near, so NOTED and not actionable (rule 9). Amin, Andriyash, Rolfe, Kulchytskyy, Melko, Phys Rev X
8:021050 (2018), quantum Boltzmann machines: the model distribution is that of a thermal state of
a TRANSVERSE-FIELD Hamiltonian, which is not a classical Boltzmann distribution over the same
variables -- the cleanest existing instance of (C1)+(C2) together, and what the charter's QBM
bullet points at. QAOA's mixer B = sum_i X_i is the structural answer to "what does an off-diagonal
term mean": it is the only thing that moves amplitude BETWEEN candidates, so S28's rank-one
similarity graph was a mixer that moves amplitude along one direction only -- the failure was the
degeneracy, not the idea. ADAPT-VQE (Grimsley, Economou, Barnes, Mayhall, Nat Commun 10:3007, 2019,
arXiv:1812.11173) "grows it systematically one operator at a time", selecting by the energy
gradient at zero angle, i.e. the commutator expectation |<psi|[H, A_i]|psi>| over an operator pool;
it is ansatz construction, not a source of non-classicality, and it is the principled way to build
a circuit matched to a structured non-diagonal H, with the resulting DLA reported.

VERDICT. Topic 4 yields NO RMSD lever. It yields: a definition (cite Barkoutsos eq 12 for the
set-equality theorem), a falsifier (CVaR's flat optimum set), a cheap unmeasured diagnostic (the
ansatz's DLA), a condition ((C1) and (C2) together), and a framing constraint (arXiv:2312.09121).
KEPT 7 (Barkoutsos definition and properties; Larocca; Ragone with caveat; Cerezo 2025; Amin QBM;
the QAOA mixer reading; ADAPT conditionally). NOTED 1 (Cerezo 2021, assumption not met here).
Multiplicity: 0 endpoint comparisons; no new measurement in this entry.
Artefacts: `s29/lit/L_4_quantum.md`; `s29/lit/L_INDEX.md`.

## S29-L14 -- TOPIC 5, THE PUBLISHED PEPTIDE CEILING: NO PAPER STATES ONE; THE PUBLISHED NUMBERS AT 9 TO 25 aa ARE 1.96 TO 2.6 A BUT NONE IS LIKE-FOR-LIKE (CURATED NMR SETS, BEST-OF-N REPORTING, AND SELECTION INSIDE AN ENSEMBLE THE SAME ENERGY GENERATED); SPLIT THE COMPARISON -- OUR GENERATION (ORACLE 1.71 TO 2.31 A) IS INSIDE THE PUBLISHED BAND AND OUR SELECTION IS THE 0.9 TO 1.5 A GAP; AND THE PEPTIDE LITERATURE'S CENTRAL ARCHITECTURAL CLAIM (RETRIEVE BY PREDICTED LOCAL CONFORMATION) IS ALREADY CLOSED HERE, WHICH I NEARLY RE-IMPORTED (2026-09-20 00:16, L)

Question (brief topic 5): what accuracy does the field reach at 9 to 16 residues, which inputs
are legitimately available here, what is the published ceiling for native-free peptide
prediction at this length, and how does 3.21 A compare? No experiment; literature only.

THERE IS NO PUBLISHED CEILING. No paper reports an upper bound for native-free peptide
prediction at this length; the literature reports method scores on small curated sets. The
nearest thing to a ceiling statement is McDonald et al. 2023's negative result (S29-L1): the best
confidence signal in the field cannot rank within a peptide's own five models. That is a ceiling
on SELECTION, and it agrees with this project's own measurement of the same thing. If S29
produces a measured bound on native-free selection for 9-16-mers with controls, that is a
contribution to the field, not only to the project.

THE NUMBERS, WITH WHAT EACH ONE IS. PEP-FOLD1 (Maupetit, Derreumaux, Tuffery, NAR 37:W498, 2009):
"averaged on 25 peptides and five runs, the PEP-FOLD LEC reproduces the NMR structure at 2.6 A
cRMSD", where the LEC is the lowest-sOPEP-energy cluster centroid -- a genuine native-free single
answer, 9-25 aa. APPTEST (Timmons & Hewage, Brief Bioinform 22:bbab308, 2021): 1.96 A on 42
peptides of 9-25 aa, selected by lowest XPLOR-NIH energy / CYANA target function (PEP-FOLD 2.05,
PEPstrMOD 4.66). PEP-FOLD3 (Lamiable et al., NAR 44:W449, 2016): near-native in the top five
scored models for 80% of 56 targets of 25-52 aa -- BEST-OF-5, and out of our length range.
PEP-FOLD4 (Rey et al., NAR 51:W432, 2023): no RMSD table; "the three methods performed similarly
on a total 17 peptides" against trRosetta and AF2, its edge being pH-dependent and poly-charged
cases, and it notes "TrRosetta and AlphaFold2 failed on two peptides of 10 and 17 amino acids
which are described as beta-hairpins experimentally". AF2 on peptides (McDonald et al., Structure
31:111, 2023): best-of-5 by class 2.2 (disulfide-rich), 2.3 (helical membrane), 2.9 (beta-hairpin),
4.4 (mixed soluble), 4.5 A (helical soluble), rank-1 being 0.2-1.1 A worse. MD (Lindorff-Larsen
et al., Science 334:517, 2011): folds chignolin (10 aa) and Trp-cage (20 aa) to the native, at
100 us to 1 ms per system. Ours: 3.2126 A built chain, 126 targets of 9-16 aa, one deployable
answer, no best-of-N (`s27/results/chain_rows.jsonl :: DIS`; cloud 3.0483).

WHY NONE OF IT IS LIKE-FOR-LIKE, three specific ways. (1) COMPOSITION: the published sets are
curated NMR peptides in solution with regular secondary structure, and PEP-FOLD states its own
scope excludes membrane-bound, ligand-bound and metal-stabilised peptides; our instrument is 126
identity-clustered PDB targets whose hard stratum is chemically identifiable -- S12 measured that
56% of FAIL18 is steric-zipper amyloid segments and lasso peptides, which no linear-window
retrieval can represent. (2) REPORTING: PEP-FOLD3's 80% and McDonald's per-class figures are
best-of-N; our matched quantities are the ORACLE top-75 ceiling 2.31 A and pool-best 1.71 A
(charter finding 10), which sit INSIDE or below the published band. (3) REGIME: every published
method selects within an ensemble its own energy or restraints generated (sOPEP ranks sOPEP-
assembled fragments; APPTEST ranks by violation of the restraints it folded under; AF2 ranks by
its own pLDDT head), which S29-L1 established is the only regime where native-free selection is
shown to work at this length. We rank 500 real windows retrieved from other proteins with an
independently constructed objective.

SO THE COMPARISON SPLITS, and this is the useful output. GENERATION: 1.71-2.31 A ORACLE is
competitive with the published field at this length; generation is not the project's problem.
SELECTION: the 0.9-1.5 A between 3.2126 and those ceilings is the whole gap, and the field's own
best method has no in-band skill here either. Reporting 3.21 against 1.96 without the three
caveats would be misleading in the project's own disfavour and should not be done.

WHAT IS LEGITIMATELY AVAILABLE. A torsion head: NO (S13, phi carries no sequence signal at this
length, 36.1 vs 36.4 deg blind; direct build 4.151 A). A structure-trained predictor (AF2,
ESMFold, OmegaFold, trRosetta): NO (leakage against the natives; ESMFold infeasible on this box,
S26 VII.2). sOPEP or another peptide-tuned coarse-grained field: NO VALUE (Legacy's functional
class, +0.330 A worse than a random subset as a ranker, S25 L16). Converged MD free energy: NO
(10^5-10^6 CPU-hours per peptide). Chemical-shift torsion restraints: CLOSED (54/126 coverage;
ORACLE-perfect torsions still 2.021 A). A fold-from-restraints terminal operator: YES IN
PRINCIPLE and it is architectural, not a QA import -- `core/project.py` already does L-BFGS over
(phi,psi) with multi-start, and what has never been run is that optimiser against the DISTOGRAM
rather than against a retrieved coordinate average (S12 Part 4); flagged to lanes M and X, priced
by nobody yet.

A NEAR-MISS OF MY OWN, RECORDED AS THE BRIEF REQUIRES. I first wrote up PEP-FOLD's central
architectural claim -- retrieve fragments by PREDICTED LOCAL CONFORMATION rather than by BLOSUM
sum -- as the strongest unexploited lead in the peptide literature, on the strength of S12's
diagnostic (2.284 -> 1.640 A pool-best on FAIL18, near-native recall 16/500 -> 92/500, robust to
30% prediction error). It is CLOSED, in two places. S13 built the predictor: "a leave-fold-out
4-state torsion-bin predictor, 0.690 accuracy overall, 0.517 on FAIL18 (majority baseline 0.562)
-- it failed as a retrieval key" (`s13/BRIEF.md`). And the 22-key screen over the whole window
universe plus 24 chain arms on 126 paired targets (`docs/FINDINGS.md`) tested `disto`, `dconf`,
`dshort` (a secondary-structure key, the closest available analogue of a structural alphabet) and
fusions: the distogram key gives the largest pool-DISTRIBUTION move in the project (best-window
percentile 39.3 -> 16.5 median, sub-2 A count 54.5 -> 101.5, pool mean 4.453 -> 3.605) while
making the pool worse where it matters (pool-best 1.711 -> 2.161, distance-ORACLE ceiling 1.994 ->
2.542) "because the distogram key concentrates the pool on one predicted structural type: many
near-native candidates on the targets it gets right, none at all on the targets it does not"; on
the chain no key beat BLOSUM (incumbent 3.454, best alternative fuse_be 3.447, CI [-0.100,
+0.083]), and S26's report lists the retrieval key as closed. The mechanism that kills it is the
typicality axis again. The S12 line that tempted me is a FAIL18 pool-best diagnostic and the
chain result that closes it lives in a different file, so anyone reading `s12/lit_FINDINGS.md`
without `docs/FINDINGS.md` beside it will re-propose this. Logged so they do not.

VERDICT. KEPT (as anchors and route notes, none importable as an operator): PEP-FOLD's 2.6 A and
APPTEST's 1.96 A as the published native-free anchors at 9-25 aa; McDonald 2023 as the field's
only ceiling-like statement (and it is a ceiling on selection); AlphaFold1's minimise-the-posterior
route and APPTEST's fold-from-restraints operator as architectural notes for lanes M and X.
REJECTED: the structural-alphabet retrieval key (closed here twice, with the mechanism);
sOPEP and peptide-tuned coarse-grained fields (Legacy's class, anti here); torsion heads (closed);
structure-trained predictors (leakage and hardware); MD free energy (cost).
Multiplicity: 0 endpoint comparisons; no new measurement in this entry.
Artefacts: `s29/lit/L_5_peptide_ceiling.md`; `s29/lit/L_INDEX.md`.

## S29-L15 -- THE COORDINATOR'S Q1 AND Q2: (Q1) THE DEPLOYED CVaR IS EXACTLY CONSTANT ON 85.5% OF THE SIMPLEX DIRECTIONS, THE ENTROPY TERM SETS THOSE TO UNIFORM, AND THE WHOLE QUANTUM STAGE REDUCES TO A TARGET-INDEPENDENT RANK-WEIGHT PROFILE WHOSE ONLY ENDPOINT CHANNEL IS m -- WHICH IS S25 L17 DERIVED AND EXPLAINS "THE OPTIMISER REACHES THE OPTIMUM AND THE STRUCTURE DOES NOT MOVE"; (Q2) A ZERO-DIAGONAL COUPLING'S THERMAL STATE IS, AT SECOND ORDER, A CLASSICAL REWEIGHTING BY THE SQUARED-SIMILARITY DEGREE -- THE CONTROL S28 ALREADY RAN -- AND THE ONLY NON-COMMUTING OPERATOR WITH EXTENSIVE STABLE RANK IS A LOCAL MIXER, SO THE CELL EXISTS ONLY IN CONFIGURATION SPACE (2026-09-20 00:22, T)

Both questions taken ahead of my remaining sections at the coordinator's request; `s29/THEORY.md`
sections Q1 and Q2 (commit 1f082628). Derivations plus one property computation; no native, no
RMSD, nothing deployable.

Q1, THE FLAT SET. Barkoutsos's statement (p(x*) >= alpha => p is a global minimiser of CVaR_alpha)
is the special case at the optimum and is not the operative one here: for a Haar-random real state
P(p(x*) >= 0.18) ~ (1-alpha)^((D-3)/2) = e^-46 at D = 512, and the deployed state has p_max ~
0.0023. The operative statement needs no overlap assumption and holds at EVERY point: by the
envelope theorem dCVaR/dp(x) = (E_x - q)/alpha on the strict tail and EXACTLY ZERO for every state
above the VaR, so at the realised m = 74.1 of D = 512 the CVaR term is blind to D - m - 1 = 437 of
the 511 simplex directions, 85.5%.
WHAT THE ENTROPY TERM SELECTS, EXACTLY. The project's CVaR is the LOWER tail, whose
Rockafellar-Uryasev form is a MAXIMUM; the bracket is linear in p and concave in t and the simplex
is compact convex, so Sion's minimax theorem gives F* = max_t [t - T log sum_x exp((t - E_x)_+/
(alpha T))] with p*(x) proportional to exp((t* - E_x)_+/(alpha T)). The code asserts the alpha = 1
limit (the Gibbs free energy and the Boltzmann law, S25 section 6.3, to 6 decimals). So p* is
EXACTLY UNIFORM above the VaR -- the 437 flat directions, resolved by flattening -- and rises
exponentially below it with scale alpha*T (per-rank factor exp(dE/(alpha T)) = 1.078, total 8.4x
across the prefix). At the deployed cell (n = 9, alpha 0.18, T 0.5, rank ladder): uniform F
-4.5394, m 92; EXACT OPTIMUM F -4.7237, m 29, prefix 29, PR 342.3, H 8.819 bits, t* -1.534;
DEPLOYED CIRCUIT F -4.5610 (sd 0.022), m 74.07, H 8.836 bits, PR 407 (seed 0, 126 targets). The
optimum is a near-uniform state with a modest prefix enhancement -- "uniform ABOVE the VaR,
enhanced below it", not "uniform ON a prefix" -- and the circuit is under-trained toward it.
THE REDUCTION. E = zrank(score[top[:D]]) is the standardised rank ladder on every target to 1.18%
of range (S25 L17), so p* depends on (alpha, T) and nothing else: the quantum stage is a FIXED
weight profile over RANKS applied to each target's own sorted list, and the emitted structure's
only target-specific input is which candidate the distogram put at which rank. Corroboration,
native-free: the realised m has sd 6.74 and corr(m, chain length) = +0.026 at seed 0 (71.2, 7.94,
-0.104 at seed 1) over 126 targets, against a between-cell spread of 29 to 92 as (alpha, T) moves;
and the tail is the classical top-m prefix to 1.1e-13 on 4,914 cells (S28-L21). Hence the full
mechanism for "the optimiser reduces the objective on 126/126 and the structure does not move":
85.5% of the directions are flat and the entropy sets them uniform; the non-flat directions are
the tail's internal weights; the readout consumes only the tail SET; the set is a prefix fixed by
the ordering; the ONE scalar left is where the prefix cuts, and m is a function of (alpha, T)
alone. THE DEPLOYED CVaR-VQE IS EQUIVALENT AT THE ENDPOINT TO CHOOSING ONE NUMBER m, and the
m-ladder has been priced three times (S22, S27 T5, `s28_B_mladder.json`).
TWO CLAUSES FOR LANE D's METER. (M5) FLAT FRACTION: every proposed objective reports the fraction
of readout-relevant directions along which it is exactly constant at production (deployed CVaR:
437/511 = 85.5%, and after the readout everything except one scalar); a candidate is worth a build
only if that fraction is materially lower AND the non-flat directions are ones the readout
consumes. (M6) THE FIXED-PROFILE CONTROL, strictly stronger than "a classical equivalent": replace
the entire quantum stage by the target-independent profile p*(alpha, T) applied to the target's
own rank order -- no circuit, no optimiser, no per-target computation. PREDICTION: the deployed
arm's emitted structure equals that control to within the built-chain input floor on at least
120/126 targets. Any formulation claiming the quantum stage contributes must break this control.

Q2, THE NON-COMMUTING FREE-ENERGY CELL. DERIVED NO for the candidate-index encoding with a
pool-geometry coupling, on a new second-order result rather than on the earlier nulls. For
H = diag(E) - J A with A_xx = 0 (every graph in the record), the Duhamel expansion of
<x|exp(-H/T)|x> has its FIRST-ORDER term vanish identically, so p_x proportional to
exp(-E_x/T)[1 + J^2 sum_y A_xy^2 g(E_x, E_y; T) + O(J^3)]: the entire quantum content of the
thermal state at leading order is a classical reweighting of the Boltzmann law by the candidate's
squared-similarity degree, a native-free scalar computable without a circuit -- and that is what
S28's degree-matched RAND control held fixed, with F3/F4 silent on 72 contrasts (S28-L41). Three
conditions the cell would need: (a) INFORMATION, the coupling must break theorem 2's (A4), which
no pool-geometry operator does -- BINDING; (b) TRAINABILITY, by my variance law the coupling is
gradient-visible at D = 512 only if its stable rank grows with the register (dense kernels ~1.6;
any Gram of structural deviations capped at 3 N_res - 6 <= 42), and at J ~ 70 to 90 the coupling
IS the energy and the object is an eigenvector problem again, violating (C2); (c) READOUT, a
p-based readout is blind to the sign structure that is the centered coupling's only new content
(section 3d). Plus a fourth from the record: S25 L15 measured that this readout cannot resolve a
45%-of-mass distributional difference, and a J^2 A^2 correction is far smaller.
WHAT WOULD HAVE TO BE TRUE INSTEAD, AND THE SMALLEST QUALIFYING FORMULATION. The obstruction in
(b) points at exactly one operator class: a sum of LOCAL Pauli terms. For a transverse field
M = sum_q X_q, ||M||_F^2 = n D and ||M||_2 = n, so r_stable = D/n, Var = 1/(nD), slope EXACTLY -1
per qubit, and J* = sqrt(n D Var_diag) = 11.8 at n = 9 -- an O(10) coupling instead of O(90). A
local mixer is only meaningful when basis states have local structure, i.e. in a
CONFIGURATION-SPACE encoding (lane X's), not in a candidate-index encoding where X_q flips a bit
of an arbitrary label. Smallest qualifying formulation: basis state = a per-residue configuration
assignment; H = H_diag(1- and 2-body posterior) + Gamma sum_q X_q; the prepared object a
free-energy or thermal state, never an eigenvector; readout = the CVaR tail's coordinate average.
CHEAPEST FALSIFIER, in order: (1) Gamma = 0 vs Gamma > 0 must differ on the SAMPLED distribution
by TV > 0.45, because below that the readout provably cannot resolve it (S25 L15); (2) the
classical counterpart must be named correctly -- for a local mixer it is a classical thermal
sampler or simulated annealing over the same configuration space at matched evaluations, NOT an
eigensolver; (3) the endpoint contrast against that sampler on the built chain. If (1) fails the
cell is empty for this instrument and no build follows. What the cell still cannot do: by (a) it
creates no information about the native's deviation from typical, so its honest upside is charter
section 17's "classically irreproducible contribution at unchanged RMSD", unless the configuration
space's own posterior carries more than the pool's marginals -- lane X's premise, measured by lane
D's meter, not by me.
Multiplicity: 0 endpoint comparisons; 3 new registered predictions (M6, the Q2 falsifier chain
clause (1), and the -1 per qubit mixer slope).
Artefacts: `s29/THEORY.md` sections Q1 and Q2; `s29/s29_T_reach.py`;
`s29/results/s29_T_reach.json` (the alpha = 1 Gibbs assertion and the (alpha, T) grid; the
per-target free-energy-gap block is queued as job s29T_reach behind the launch cap and is not
quoted above -- every number quoted comes from the rank-ladder computation or from
`s27/results/s28_B_rows.jsonl`).

## S29-L16 -- TOPIC 6, RANKING INSIDE A MATCHED-REALISM BAND: NOBODY HAS DONE IT (THE CLOSEST METHODOLOGICAL PAPER DOCUMENTS THE CONFOUND AND DECLINES THE FIX), SO LANE D's MEASUREMENT IS NOVEL; THE DESIGN IS PUBLISHED AND HUMAN-VALIDATED (PIRM 2018) AND ITS OWN RESULT IS A WARNING -- THE REALISM INDEX CORRELATES 0.83 *BETWEEN* BANDS AND UNRELIABLY *WITHIN* THEM; AND THE QUANTITY BEING MEASURED HAS A CLOSED FORM, THE PARTIAL CORRELATION rho_SY.R, WHICH IS EXACTLY ZERO IFF A SCORER's LINK TO ACCURACY IS FULLY MEDIATED BY REALISM (2026-09-20 00:23, L)

Question (coordinator's topic 6, following S29-L12): supply lane D with the literature on
building a matched-realism band and on what to expect inside it -- (a) band-restricted evaluation
in the perception-distortion literature, (b) equal-energy-shell ranking in statistical physics and
energy-matched decoys in QA, (c) whether any published QA evaluation has conditioned on a realism
statistic before correlating with GDT/RMSD, (d) a bound on surviving in-band signal analogous to
the Ueda-Nakano coefficient. No experiment; literature only.

(a) THE DESIGN IS PUBLISHED, AND ITS RESULT IS A WARNING. Blau, Mechrez, Timofte, Michaeli,
Zelnik-Manor, "The 2018 PIRM Challenge on Perceptual Image Super-resolution", ECCV 2018 Workshops,
arXiv:1809.07517, is exactly this construction in transpose (fix distortion, rank by realism). Its
rationale is ours: methods at different distortion levels "cannot be compared or ranked using
these common metrics". "the perception-distortion plane was divided into three regions by setting
thresholds on the RMSE values (regions 1/2/3 were defined by RMSE <= 11.5/12.5/16 respectively).
In each region, the goal was to obtain the best mean perceptual quality." The realism axis
PI = (1/2)((10 - Ma) + NIQE) is built from NO-REFERENCE measures, i.e. native-free, the property
we need; it was validated by 35 human raters scoring "how realistic the image looked" WITHOUT
sight of the ground truth. THE RESULT LANE D MUST PRE-REGISTER AGAINST, verbatim: "while the PI is
well correlated with the human-opinion-scores on a coarse scale (in between regions), it is not
always well-correlated with these scores on a finer scale (rankings within the regions)" --
Spearman 0.83 across bands, unreliable within them. That is this project's own global-vs-in-band
split (distogram +0.653 global, +0.091 in-band) reproduced independently in vision on this exact
design. Expect a SMALL in-band effect and power for it. Also verbatim, and it is our contracted
average in another field: "the outputs of EDSR, a state-of-the-art algorithm in terms of
distortion, are mostly voted as 'definitely fake'. This is due to the aggressive averaging causing
blurriness as a consequence of optimizing for distortion." And PIRM observed "the tradeoff appears
to be stronger in the low distortion regime", which is Theorem 3's convexity measured -- we
operate at that end, so the band must be NARROW to mean anything, which is the experiment's
central design tension.

(b) THE PHYSICS GIVES THE CLEANEST STATEMENT, AND IT IS WHY THE EXPERIMENT IS INTERPRETABLE EITHER
WAY. Inside an exactly constant-energy shell the Boltzmann weight is constant, so the energy
induces ZERO ordering within the shell and everything in-band is carried by the density of states.
Generalised: banding on a statistic R destroys R's own discriminating power by construction, so an
in-band ranking experiment measures EXACTLY what is orthogonal to R. A positive is therefore
direct evidence of information orthogonal to realism -- the quantity `s29/STATE.md`'s H0 says the
system lacks; a null says the library carries nothing beyond realism. Umbrella sampling (Torrie &
Valleau 1977) with WHAM (Kumar et al. 1992) supplies the discipline for an imposed selection:
unbias with the known bias, or confine every claim to the band -- a within-band mean is NOT an
estimate of a pool mean. On "energy-matched decoys": the QA literature has restricted sets
(CASP11's best150 = the 150 best by consensus; sel20 = 20 maximally different models; VoroMQA's
BZQ15 = models from three strong servers) but none is realism-matched, and their own lesson is a
warning. Olechnovic & Venclovas (2017), verbatim: "for sel20 sets, every method that is based
solely on analyzing geometric features and applying statistical potentials (GOAP, DOOP, dDFIRE and
all the VoroMQA variations) achieved worse results than the best-performing composite methods
incorporating evolutionary information"; and a trivial HHpred-agreement score (TM-score to a
homology model, i.e. agreement with an INDEPENDENT predictor) matched ProQ2 on sel20 while being
"much worse than all the other tested QA scores for best150 and BZQ15". WHICH METHOD WINS DEPENDS
ON HOW THE SUBSET WAS BUILT, and pure-geometry scorers swap places with information-carrying ones
between subsets: the band's construction is not neutral, must be pre-registered with its
rationale, and a result on one band definition does not transfer to another. Corroboration from
our own field: ANDIS (Yu et al., Bioinformatics 2019) frames it as "native recognition emphasizes
the differences of overall structure quality between native and decoy structures, while decoy
discrimination generally focuses on the backbone differences among decoy structures ... The
potential's abilities of native recognition and decoy discrimination cannot be optimized
simultaneously with the same parameter sets" -- the protein-potential community arriving at
Blau & Michaeli's conclusion empirically, without the theorem.

(c) NOBODY HAS DONE IT, STATED EXPLICITLY AS THE COORDINATOR ASKED. I found no published QA
evaluation that conditions on a native-free realism statistic and then measures a score's
correlation with GDT or RMSD inside the band. The closest and most careful paper is Hamelryck et
al., "Artefacts and biases affecting the evaluation of scoring functions on decoy sets for protein
structure prediction" (PMC2677743): it establishes that for "139 out of 149 of the decoy sets
considered" the native is trivially discriminable through individual energy terms (improper
torsions, vdW clashes) rather than fold quality, that MD decoys violate i.i.d. with correlations
that "arise primarily as a consequence of differences _between_ the five trajectories", and that
near-native enrichment inflates performance -- and it explicitly does NOT implement matching or
conditioning on a confounder, nor analyse a narrow quality range; its mitigations concern
effective sample size. So the field's most careful methodological paper identified the confound
and stopped short of the fix. LANE D's MEASUREMENT IS NOVEL, NOT DERIVATIVE. Standard caveat on a
negative literature claim: this is a targeted search across the QA, decoy-evaluation and
perception-distortion literatures, not a proof of absence.

(d) THE BOUND, AND IT IS EXACT. Let Y be accuracy, R the realism statistic defining the band, S a
scorer. For jointly Gaussian (S, Y, R) the conditional covariance of (S,Y) given R = r is
Sigma_SY - Sigma_SR Sigma_RR^{-1} Sigma_RY, which does NOT depend on r, so the in-band correlation
is the partial correlation, constant across the band:
    rho_SY.R = (rho_SY - rho_SR rho_RY) / sqrt((1 - rho_SR^2)(1 - rho_RY^2)).
Its ZERO is the whole experiment: rho_SY.R = 0 exactly when rho_SY = rho_SR rho_RY, i.e. when a
scorer's entire association with accuracy is MEDIATED by realism. The matched-realism band is
precisely a test of whether any scorer carries accuracy information beyond its realism content.
THREE OPERATIONAL CONSEQUENCES FOR D. (i) The answer is PREDICTABLE BEFORE ANY ENDPOINT RUN from
three native-free-computable correlations (rho_SY needs the native and is ORACLE, so it is a
diagnostic, but rho_SR is not) -- computing them first IS the pre-registration. (ii) Band WIDTH
interpolates: conditioning on R in [a,b] is not conditioning on R = r, so in-band skill runs
monotonically from the global rho_SY (wide) to rho_SY.R (thin). Report skill as a CURVE over band
widths and extrapolate rather than picking one width: it answers the multiplicity objection with
one pre-registered curve instead of k bands, and converts the width-versus-power tension into a
measured trend instead of one underpowered point. (iii) The classical name for the attenuation is
range restriction (Thorndike's case II); a raw within-band correlation must never be compared to a
raw global correlation without stating the selection. THE PROJECT HAS MET THIS LOGIC TWICE
ALREADY: `shared-referent-floor` (two quantities measured against a common reference correlate by
construction; measuring the floor "turned '2/3 sequence-independent' into '1/5'") is this
correction applied as a diagnostic, and `decorrelated-errors-exist-but-are-unusable` reports
truth-partialled error correlations of 0.04-0.26 with fusion worth +0.004-0.011 "because the gain
goes as the SQUARE of the weaker channel's skill". The squared-skill law is the honest prior for
the Angstrom value: a small partial correlation buys nearly nothing, so a positive here would be
most valuable as a MECHANISM result (information orthogonal to realism exists) rather than as an
accuracy result. I quote no invented effect size.

VERDICT. KEPT 7, NOTED 1; no importable operator -- the deliverable is a design, a bound, a prior
and a warning. The five lines for lane D: the design is published and human-validated (PIRM); the
measurement is novel (nobody has conditioned on realism); the target quantity has a closed form
(rho_SY.R, zero iff fully mediated) and is predictable in advance; report skill versus band width
as a curve; and expect the existing library to underperform, because every scorer in it was fitted
for the BETWEEN-band task (ANDIS), so a null from the current library would NOT close the in-band
question -- it would say the library is the wrong instrument for it, which is a weaker and
different claim that must be stated before the run.
Multiplicity: 0 endpoint comparisons; no new measurement in this entry.
Artefacts: `s29/lit/L_6_matched_realism.md`; `s29/lit/L_INDEX.md`.

## S29-L17 -- THEORY SECTIONS 1 AND 4: (1) CONTRACTION IS AN EXACT VARIANCE IDENTITY, d(C)^2 = <d_k^2> - s^2, SO IT IS THE POOL'S OWN 32% AND A POSTERIOR CALIBRATION CANNOT TOUCH IT OR THE MINIMISER (ONLY THE SEPARATION-BAND WEIGHTS, WHICH THE 2x OVER-CONFIDENCE DISTORTS BY 2.7x); (4) UNDER A CVaR OVER THE TAIL'S AVERAGE STRUCTURE THE SET-EQUALITY THEOREM FAILS -- THREE STATES, TWO DIMENSIONS, THE OPTIMAL SUBSET IS THE TWO WORST STATES BY EVERY PER-STATE CRITERION (2026-09-20 00:26, T)

`s29/THEORY.md` sections 1 and 4 (this commit). Derivations; no job, no native, nothing
deployable. Section 6 was answered inside Q1 (S29-L15) and is not repeated.

SECTION 1, THE CONTRACTION THEOREM. (1.1) The L1 Bayes risk is minimised at the per-pair median
and its local form is S(C) = const + (1/2)||D(C) - m||^2 in the w*kappa metric, kappa = 2 p(median):
the objective is a weighted squared distance from the median map, and the median map is generally
not realisable, so the structural minimiser is its projection. (1.2) THE EXACT IDENTITY, not
Jensen: for a coordinate average, d_ij(C)^2 = mean_k d_ij(W_k)^2 - s_ij^2 with s_ij^2 the members'
pair-vector spread, and Rg(C)^2 = mean_k Rg(W_k)^2 - Delta^2 with Delta^2 exactly S23 L9's
IDIOSYNCRATIC term. Contraction is the variance decomposition of the pool: the operator cannot
take the 68% and cannot avoid paying the 32% in contraction. The measured 22% bond contraction
requires s^2/<d^2> = 0.396; the envelope check gives Rg = sqrt(6.80^2 - 63.82/N) = 6.41 against a
measured 6.21 (3% high; the 63.82 is a per-target sum over targets of varying length). The
separation dependence is S23 L1's "averaging smooths" with no extra assumption: s^2 is the
superposition residual and does not grow with |i-j| while d^2 does, so the fractional contraction
falls from 22% at the bond to 6% at the envelope. (1.3) The average beats the typical member by
the pool's dispersion and loses 2<b,c> + |c|^2 to the contraction; the members that beat it are
exactly those inside a ball of radius |b + c| about the median map, and that ball holds an eighth
of the pool (pct(PROD) = 0.126, S28-L36). (1.4) A 2x OVER-CONFIDENT POSTERIOR: the posterior
enters twice and the two entries behave differently. Through the MEDIAN MAP -- unchanged by a
width error (z_mean -0.052, S25 L1) -- so the minimiser is unchanged to first order, which derives
S25 L2's null and S25 L7's retraction. Through the METRIC w*kappa ~ 1/sigma^2 -- and a UNIFORM
over-confidence multiplies the whole metric by a constant and changes no minimiser at all. Only
the heterogeneity matters: z_sd is 1.23 / 2.05 / 1.87 / 1.28 across separations 2-2 / 4-5 / 6-8 /
9-15, so the shipped objective over-weights mid-range pairs by about 2.7x against a calibrated
one. CORRECTION TO THE STANDING PROGRAMME (`s27/REPORT_S28.md` section 12 item 2, "calibrate the
posterior and re-read the meter"): calibration cannot move the median map, hence not the
minimiser, and cannot touch the contraction, which is the pool's dispersion; the only live version
of that item is a SEPARATION-BAND RE-WEIGHTING with one or two parameters.
PREDICTION (12 targets, minutes, lane D or M; inputs all cached): dS/dlambda at production =
sum w d (2F(d) - 1) is NEGATIVE on >= 75% of targets (the objective wants production EXPANDED),
while the native-optimal scale is BELOW 1 on 73/126 (S23 L9) -- so the sign agreement on the scale
axis is below a coin toss, the scale-axis instance of theorem 2; and widening the posterior to
z_sd = 1 leaves that sign unchanged on >= 90% of targets while changing the magnitude by 0.4 to
0.7x. Falsified if clause 1 holds on under half the targets or if calibration flips the sign on
more than 25%.

SECTION 4, CVaR OVER A STRUCTURAL OBSERVABLE. Two inequivalent lifts: (a) scalarise then tail
(a per-state structural loss; everything transfers, including set-equality; this is S27's 80
diagonal energies) and (b) TAIL THEN AGGREGATE, R_alpha(p) = sum_x lambda_x(p) W_x, the tail's own
coordinate average, with the objective f(R_alpha(p)). For (b): the envelope argument gives
dR/dp_y = (W_y - W_{x_q})/alpha on the strict tail, so dF/dp_y = <grad f(R), W_y - W_{x_q}>/alpha
and the parameter-shift chain costs the SAME 2P evaluations as the deployed objective. Non-smooth
exactly where the scalar CVaR is, plus f's own bin kinks, plus -- if the tail's ORDER is allowed to
depend on R -- order permutations, which is why the order should stay a per-state scalar; on each
tail-set cell R is linear in p, so a convex f makes F convex there and the landscape is a union of
convex cells indexed by the SET. Shot noise: (b) needs a quantile plus a MEAN structure over the
tail (variance Cov(W|tail)/(alpha S)), not the full 2^n-outcome distribution the deployed gradient
needs -- (b) is the MORE device-realisable of the two.
THE SET-EQUALITY THEOREM FAILS FOR (b). Smallest witness: three states, two dimensions.
W1 = (+1,0), W2 = (-1,0), W3 = (0,0.1), target (0,0), f = |R|, tail size 2. Per-state loss puts W3
first (0.1 vs 1), so every per-state prefix of size 2 is {3,1} or {3,2} with f = 0.5025, while the
optimal set {1,2} has f = 0 -- the optimal subset is the two WORST states by every per-state
criterion, because their errors cancel. V(S) = f(mean_{i in S} W_i) is neither additive nor
monotone, so no greedy per-state rule is optimal and the gap is not measure-zero. THIS IS THE ONE
PLACE IN THE PROJECT WHERE THE THEOREM THAT MAKES THE DEPLOYED SPINE CLASSICAL GENUINELY BREAKS:
under (b) "which set" is a real optimisation variable rather than a read-out of the sort, which
answers charter questions 5 and 7 in the affirmative IN PRINCIPLE. Three binding caveats: escaping
the theorem creates no information (f is still bound by theorem 2, and if f is the distogram risk
the best-cancelling set is the one whose average sits at the median map -- the contracted average
again); the classical counterpart is not absent, it is greedy plus local search over subsets, and
rule 15 requires it; and the landscape becomes combinatorial, which is the right shape for a
variational method and also the shape where classical local search is hard to beat at D = 512.
PREDICTION (native-free, minutes, lane D or X): a greedy-plus-local-search 75-subset minimising
f(mean of the subset) beats the DIS top-75 on the OBJECTIVE by at least 0.10 (falsified under
0.02), AND its coordinate average is WORSE than production on the built chain by at least 0.1 A
(ladder rho -0.40, S29-L2; the minimiser is away from the native, S28-L18b). So escaping the
theorem is necessary and NOT sufficient -- and if the second clause fails, i.e. the
aggregate-optimal subset is BETTER than production, that is the sprint's first genuine opening and
should go to 126 immediately.
Multiplicity: 0 endpoint comparisons; 2 registered predictions (the scale-derivative triple and
the subset-gap pair).
Artefacts: `s29/THEORY.md` sections 1 and 4. Numbers quoted are from S23 L1/L9, S25 L1/L2/L7,
S28-L18b/L36 and S29-L2, each named in the text.

## S29-L18 -- THEORY SECTIONS 5 AND 7, AND THE LANE CLOSES ITS BRIEF: (5) THE DLA AT n = 3 MAKES S26's UNEXPLAINED n = 6 / n = 9 OBSTRUCTION A RULE (3 | n), THE DEPLOYED n = 9 DEPTH-3 ALGEBRA IS su(256), HALF OF so(512) -- AND THE ALGEBRA IS NOT THE OBSTRUCTION, THE 27 PARAMETERS ARE; THE ENCODING MAKES CLASSICAL PREFIXES RANK-2 AND GENERIC 75-SUBSETS RANK-16, AN ANSATZ-SIDE REASON FOR THE PREFIX RESULT INDEPENDENT OF THE CVaR CLIP; (7) EVERY CANDIDATE INFORMATION SOURCE RANKED BY cov(source, n), AND ONLY THREE CLASSES ARE FIRST ORDER (2026-09-20 00:30, T)

`s29/THEORY.md` sections 5 and 7 (this commit); job `s26/jobs_done/` not required for the light
phases, which ran in-process (`s29/results/s29_T_reach.json`: `dla`, `subset_cut_rank`,
`entropy`); the per-target free-energy-gap phase is queued as `s29T_reach` behind the launch cap
and nothing above depends on it. Property measurements only: no native, no RMSD, no pool in
section 5.

SECTION 5, THE REACHABLE SET. (a) THE DLA, exact Pauli-set closure with an independent dense
nested-commutator cross-check at every n <= 5 cell (12/12 agree): n = 3 gives 3 / 8 / 16 / 28 at
L = 1..4 with dim so(8) = 28, so n = 3 needs depth 4; n = 4 and n = 5 are full at depth 2; n = 6 is
510 / 1023 / 2016 (S26 reproduced). With S26's n = 7..11 this makes a RULE over n = 3..11: the
obstruction occurs exactly when 3 | n, and in those cases the closure is 2 dim su(2^(n-2)) at
L = 2 and dim su(2^(n-1)) at L = 3, climbing a chain of proper subalgebras one rung per layer.
S26 labelled the n = 6 / n = 9 pattern HYPOTHESIS, not explained; the mechanism (the three-step
CNOT chain plus ring leaving the early conjugation rounds inside the stabiliser of a complex
structure when the ring length is divisible by 3) is labelled CONJECTURE, and its next test n = 12
is not runnable at s26/q_dla.py's CAP = 2^21 (dim su(2^11) = 4,194,303).
SCOPE CORRECTION TO `s26/REPORT.md` V.9: "the algebra at the deployed cell is maximal" is true at
n = 7, L = 3 (8128 = so(128)) and FALSE at the S27/S28 register n = 9, L = 3, where it is 65535 =
dim su(256), exactly half of dim so(512) = 130816. The plateau conclusion is unchanged; the
sentence needs its n.
(b) AND THE ALGEBRA IS NOT THE OBSTRUCTION: su(N) acts transitively on the unit sphere of C^N =
R^(2N), so even su(256) generates a group whose orbit of |0...0> is all of S^511. What is
unreachable is everything outside the image of a smooth map R^27 -> S^511: dimension <= 27 of 511,
measure zero, with a generic target's best squared overlap of order P/D = 0.05. S28's measured cap
of 0.80 to 0.88 against the J = 3 ground state therefore says that state is very far from generic
(PR 301 against the family's 465), not that the family is nearly universal. Any argument that
reaches for controllability here answers the wrong question.
(c) THE CUT-RANK ASYMMETRY, measured (n = 9, 200 random subsets per size): a random 75-subset's
uniform superposition has max-over-cuts Schmidt rank 16 -- the register's generic maximum -- while
the PREFIX {0..74} has rank 2. The depth-3 chain gives bond dimension 8 (16 with the ring). So the
deployed encoding makes exactly the classical prefixes cheap and every other set expensive: an
ANSATZ-SIDE mechanism for "the circuit reproduces the classical top-m", independent of the CVaR
clip the set-equality theorem argues from -- two separate reasons for one fact. The gap is a factor
2 at n = 9 and grows with the register (8 against 64 at n = 12), which binds on any formulation
that widens the register to make an off-diagonal term trainable. A design variable nobody has
used: `s22/qcand_lib.py :: Encoding` takes a free label permutation, set to the identity in
deployment; it decides which FAMILY of sets is cheap, and it comes with its own control (a random
relabelling).
(d) THE 0.85 CAP IS AMPLITUDE CONCENTRATION, NOT SIGNS (S28-L41: the sign-aligned |psi_vqe| is
representable at 0.986). PREDICTION, from dimension counting (1 - F ~ 0.151 x 27/P): F = 0.887 at
L = 4, 0.925 at L = 6, 0.943 at L = 8. FALSIFIER: F at L = 6 outside [0.88, 0.96] on the median of
4 targets; the exponential alternative predicts 0.98 and is separated by the same run. About ten
minutes by re-running `s27/s28_B_represent.py` at layers in {4, 6, 8} on 4 of its 12 targets (one
line: the module reads B.LAYERS). What it buys is design guidance, not accuracy -- S28 already
showed reaching the coherent basin moves the endpoint by nothing.

SECTION 7, WHAT AN OBJECTIVE MUST KNOW. The meter's four numbers are not equally informative:
theorem 2 makes the GRADIENT COSINE second order and purchasable with no information (section 2.4's
shrink; section 1.5's scale axis), while the NATIVE PERCENTILE and the PREFERENCE both require
cov(channel, n) != 0 and cannot be gamed that way. Meter numbers 3 and 4 are the binding ones and
2 must always be read beside 3. Every candidate source is then judged by one question -- does it
carry cov(., n) given the distogram -- and the record's ranking is: (1) A BETTER DISTANCE PRIOR
(breaks A4 directly; -2.15 A per unit, the only steep lever, S24; blocked by hardware and leakage,
not by theory); (2) A LEARNED RESIDUAL (first order by construction; the binding constraint is
S19's error coherence, +0.31 A at 0.688 sign accuracy when the mistakes are coherent); (3) A
PHYSICS TERM ON THE EMITTED STRUCTURE (outside class M, so theorem 2 does not bound it; excluded by
measurement in its single-point form, S25 L16, but its FREE-energy form is unmeasured here and is
the only class with peptide-length precedent, S29-L1, with S8's F_qh machinery built and never
finished); (4) A SECOND POOL (breaks A3 only: second order; S24 L2/L3 prices it at cos 0.943 once
selected, 31% independent only when 0.76 A worse); (5) AN ESM-ATTENTION MAP (redundant with a
distogram that already consumes ESM-2 650M PCA-32; the cheap test is the partial correlation with n
given DIS, not an endpoint run); (6) THE POOL'S DISPERSION (a second moment: predicts the error's
MAGNITUDE, never its SIGN, and both binding meter numbers need the sign); (7) A JOINT over the same
marginals (section 2 C2: the repair lands partly in ker(Jc^T), which the gradient annihilates
exactly -- 45% of pair space at N = 12 -- so it buys zero locally by construction).
THE READING: rows 4 to 7 are what a quantum formulation naturally reaches for and theorem 2 prices
them all at second order; rows 1 to 3 are the only first-order classes and two are outside this
sprint's reach. Charter finding 8 is not a gap in the scorer library, it is a corollary of the
information the system holds. The reachable results are therefore (i) mechanism, now in closed form
(theorem 2, the stable-rank law, the Q1 reduction), (ii) the ONE structural opening of section 4.5
(the set-equality theorem genuinely fails under a CVaR over the tail's average structure, so "which
set" becomes a real optimisation variable) -- with the order of work stated: find an f that clears
meter numbers 3 and 4 FIRST, then build the CVaR that consumes it, because building the structural
CVaR around the shipped cost reproduces S28 with more machinery -- and (iii) the free-energy class,
whose natural quantum home is exactly Q2's configuration-space cell with a local mixer.

LANE T's BRIEF IS COMPLETE: sections 2, 3, Q1, Q2, 1, 4, 5, 7 in `s29/THEORY.md` (section 6 is
inside Q1), each with a measurable prediction and a named lane. Multiplicity: 0 endpoint
comparisons by this lane in the whole sprint; 13 registered predictions, 0 measured by me.
Artefacts: `s29/THEORY.md`; `s29/s29_T_spectra.py` + `s29/results/s29_T_spectra{.json,_rows.jsonl}`
+ `s29_T_grad_rows.jsonl` (job s29T_spectra, exit 0, 100.2 s, 0.337 GB);
`s29/s29_T_reach.py` + `s29/results/s29_T_reach.json`.

## S29-L19 -- TOPIC 7, TRAINING AN IN-BAND RANKER: THE RECORD ALREADY ANSWERS IT TWICE (S14's 0.986 WITHIN / 0.600 ACROSS AGAINST 0.638 NEEDED, CAPACITY SATURATED BY A LINEAR MODEL; S12's FLAT LEARNING CURVE AGAINST A LOUD LEAKED CONTROL), THE LITERATURE's PAIRWISE GAINS ARE ALL MEASURED *BETWEEN* BANDS (RMSD 0-12 A), AND -- THE URGENT ITEM -- THE IN-BAND AXIS IS COMPACTNESS (+0.909 WITH THE NATIVE's z-SCORED Rg), SO A COMPACTNESS-LOADED REALISM BAND WOULD PRODUCE A SELF-FULFILLING NULL; D MUST MEASURE rho(R, Rg) BEFORE FIXING THE BAND (2026-09-20 00:30, L)

Question (coordinator's topic 7, commissioned before D's in-band result lands): (a) what does
learning-to-rank say about within-group ranking when the grouping variable is the dominant signal
and has been conditioned away; (b) has anyone trained a structure scorer on same-target PAIRS
rather than absolute labels, and what did it buy; (c) what are the cheap surrogates for
conformational entropy at peptide length and how well do they approximate the converged answer;
(d) given the record, what is the honest prior for a trained in-band ranker? No experiment;
literature plus the record.

**THE URGENT ITEM FIRST, because it affects lane D's design and not just its interpretation.**
S14's flip diagnostic (`in-band-ordering-is-per-target`) found that per-target in-band skill
correlates **+0.909 with the native's z-scored RADIUS OF GYRATION** and **+0.951 with
rho(contacts, RMSD)**. The in-band ordering axis IS compactness. Every realism statistic available
in the S27 library is compactness-like: Legacy is explicitly a compactness model (top-75 Rg
-0.758 A), and RG_LAW, RG_UNIV, EXVOL, CAGEO, CONTACT and ENV are all functions of packing. So if
D bands on a compactness-loaded realism statistic R, the band removes the exact axis that carries
whatever in-band skill exists, and by S29-L16's partial-correlation bound rho_SY.R is crushed
toward zero BY CONSTRUCTION OF THE BAND. The resulting null would be SELF-FULFILLING and
uninterpretable as evidence about in-band signal. THE FIX IS CHEAP AND MUST PRECEDE THE RUN:
measure rho(R, Rg) for the chosen band statistic and report it beside the result. Three cases:
(i) strongly compactness-loaded -- the design is self-defeating and R must be changed or
length-residualised; (ii) weakly loaded -- valid, and rho(R, Rg) is exactly the number needed to
interpret the effect size; (iii) orthogonal -- clean, and the experiment is then testing for a
SECOND in-band axis beyond the one S14 found, which is the most interesting version and should be
stated as the hypothesis. Constructive corollary already in the record: native-free compactness
proxies reach residualised r = 0.244 (pool mean Rg), 0.313 (distogram-predicted Rg), 0.365
(incumbent emitted Rg) against the ORACLE's 0.909, all three CIs excluding zero, i.e. "roughly 61%
sign accuracy, against a project where fusion gain goes as the SQUARE of the weaker channel's
skill".

(d) THE RECORD ALREADY ANSWERS TOPIC 7, TWICE. S14: in-sample 0.986, same-target-held-out 0.986
with an overfitting gap of -0.0005 [-0.0011, +0.0000], cross-target 0.600 [0.549, 0.655] -- a
transfer gap of -0.3859, **770x the overfitting gap**, against the 0.638 required for 2.0 A
through a top-100 operator. "A LINEAR 4,125-parameter pair potential already saturates the
within-target problem at 0.986, so nonlinearity, capacity, equivariance, graph or geometric
architectures cannot improve on it and none of them touches transfer, where the entire loss sits."
The learning curve DECLINES with more targets (m=1 at -0.348 beats m=12 at -0.162) because
averaging across targets cancels a per-target axis whose sign flips. S12: a set-transformer over
the full signed deviation map -- the last untested architecture class -- has a FLAT learning curve
(3.043 / 3.049 / 3.046 / 3.036 / 3.026 at n = 8/16/32/64/75) while the identical harness with a
LEAKED label is loud at n=8 (2.534). Signal-limited, not sample-limited. HONEST PRIOR: a trained
in-band ranker will train beautifully within a target and transfer at ~0.600 against a 0.638
requirement, and no change of loss (pairwise, listwise, lambda-weighted, contrastive) or
architecture will move it, because the measured bottleneck is the per-target SIGN of the axis,
which inference cannot observe. The literature in (a) and (b) offers losses and architectures; the
record says the constraint is neither. IF D's RESULT IS NULL, the follow-up should NOT be a better
in-band ranker but the per-target SIGN question, which is the one direction S14 named as open:
"the only direction with leverage is anything that supplies the PER-TARGET SIGN of the in-band
axis at inference: a conditioning signal, not a better objective."

(a) LTR SUPPLIES THE OPERATION, NOT THE INTENT. RankNet (Burges et al., ICML 2005) models
P(i beats j) = sigma(s_i - s_j) within a query with a cross-entropy loss; LambdaRank weights each
pair's gradient by the metric change from swapping it; ListNet/ListMLE replace pairs with a
Plackett-Luce likelihood. THE PROPERTY THAT MATTERS, in one line: if s(x) = f(x) + g(group), then
within a group s_i - s_j = f(x_i) - f(x_j) and the per-group confound CANCELS EXACTLY. So a
matched-realism band plus a pairwise objective is formally a MATCHED CASE-CONTROL DESIGN with
CONDITIONAL LOGISTIC REGRESSION as its estimator -- that is the correct statistical name for what
the sprint is about to do, and it carries the right warnings. NEGATIVE LITERATURE RESULT: nothing
in LTR treats the grouping variable as a NUISANCE TO BE REMOVED; it treats it as a QUERY to be
served, where cross-query incomparability is not a defect but the task definition. The algebra is
identical, the inferential claim is not: in LTR, good within-query order with no cross-query
calibration is a success; here it may still be useless, because the deployed operator must choose
among candidates the band discarded. The confound framing belongs to matched-design statistics and
causal inference, not LTR. POWER CAUTION: pairwise training on n items per group yields O(n^2)
pairs that are NOT independent; the effective sample size for a cross-target claim is TARGETS, not
pairs. Any in-band pairwise result must be powered on targets. ONE UNEXPLOITED REFINEMENT:
LambdaRank's metric weighting -- weight each in-band pair by |Delta RMSD| instead of equally (the
project's S28-L48 pairwise learner used unweighted sign pairs). NOTED, not KEPT: it is a loss
refinement, not an information change, so contract rule 14 excludes it as an import on its own.

(b) YES, IT EXISTS, AND THE GAINS ARE BETWEEN-BAND. Jing X, Dong Q et al., "Sorting protein decoys
by machine-learning-to-rank" (PMC4987638): SVMrank, linear kernel, and "training instances are
pairs of decoys from the same proteins" -- within-target pairs, our design. Features are
knowledge-based potentials (DFIRE, DOPE, GOAP, RWplus, ...) plus QA-program outputs. On 3DRobot:
classification 0.30 wmPMCC / 0.68 AUC, regression 0.88 / 0.94, **pairwise ranking 0.95 / 0.95**.
THE ASSUMPTION THAT BREAKS THE TRANSFER: 3DRobot spans **RMSD 0 to 12 A** -- the widest possible
between-band comparison -- and the paper does not analyse cross-target transfer at all, which is
exactly the quantity S14 measured here and found to collapse. ProQ4's Siamese rank loss is the QA
literature's own pairwise variant (R-per-model 0.56 vs ProQ3D's 0.61; its edge is in first-rank
loss, the selection metric). The Frontiers in Bioinformatics 2023 paper (PMC10616882) designs a
quality-dependent epsilon-insensitive loss L(y,y') = max(0, |y - y'| - eps(y)) with eps narrowing
as quality rises (0.45 at scores 0-0.1, 0.01 above 0.8) because "the model should not try too hard
to accurately fit poor-quality decoys where we do not need good accuracy anyhow" -- overall
correlation 0.75 -> 0.90 in ablation, but it reports NO disaggregated performance on the
high-quality subset it was designed for, independently corroborating S29-L16's "nobody reports
in-band skill". THE INFORMATION TEST: a pairwise model sees the SAME features with a different
loss -- nothing new. And the project has already run it: S28-L48's nested pairwise-logistic over
all 31 scorers reached 0.960 held-out sign accuracy against a 0.502 null, while the SAME rule on
RAND_SIGNED reached 0.952 with head-to-head 0.492. It learned "is this a real protein trace" --
the realism axis, which is exactly what the band removes. USEFUL PREDICTION FOR D: if the band is
correctly constructed, re-running that pairwise learner inside it should collapse the 0.960 toward
chance. That is a cheap validity check ON THE BAND ITSELF, before anything is concluded from a
null.

(c) THE FREE-ENERGY CLASS: BIASED, BADLY CONVERGENT, AND PROBABLY NOT ORTHOGONAL. Quasi-harmonic
entropy (Karplus & Kushick 1981) diagonalises the coordinate covariance and sums harmonic-oscillator
entropies. Numata, Wan & Knapp (2007), "Conformational entropy of biomolecules: beyond the
quasi-harmonic approximation": the largest classical eigenvalues "tend to be more anharmonic and
show statistical dependence beyond correlation", and a k-nearest-neighbour correction "calculates a
tighter UPPER LIMIT to entropy than the quasi-harmonic approximation" -- so S_qh is an UPPER BOUND
when converged and F_qh = <E> - T S_qh is correspondingly an UNDERestimate of the free energy. The
literature also reports "severe convergence problems" and a single-minimum assumption "too crude
for flexible structures"; normal-mode entropy is "very computationally expensive" and usually
truncated. RECONCILIATION FOR S8's STAGE, labelled as my reading of two sources: `docs/FINDINGS.md`
calls S8's S_qh "noisy and downward-biased" at 4 ps / 100 snapshots over up to 42 CA dof, which is
NOT in conflict with "upper bound" -- the functional form is an upper bound once converged, while
finite sampling underestimates the covariance and pushes S_qh down. The two biases have OPPOSITE
signs and different per-candidate magnitudes, so neither cancels in a within-target ranking; a
4 ps F_qh ranking is not merely noisy. THE POINT NOBODY HERE HAS STATED, flagged as a hypothesis:
F = E - TS ranks differently from E only if S varies across candidates non-collinearly with E, and
what S tracks for a 9-16-mer is BASIN WIDTH -- compact well-packed conformations have narrower
basins, extended ones wider. So the entropy term is plausibly COMPACTNESS-LIKE, i.e. on the same
axis as Legacy, RG_LAW, RG_UNIV, EXVOL and CAGEO -- not the orthogonal information source S29-L1
hoped for. Both of the sprint's live ideas (the matched-realism band and the free-energy stage)
would then be acting on the SAME axis, in opposite directions. RECOMMENDATION: gate the AMBER
spend on a three-target diagnostic answering one question -- is the entropy term orthogonal to
compactness? Correlate S8's `width` / `S_msf` channels against Rg and the existing compactness
channels. If not orthogonal, the class closes cheaply; if orthogonal, that orthogonality is itself
the result and justifies the full cost. THE COORDINATOR's CAVEAT, RESTATED AS ASKED: any
free-energy arm needs the ENTROPY term to rescue a channel (AMBER single-point) measured **+0.455
A worse than a random subset** as a ranker, 5/5 folds (S25 L16) -- a strong claim for a term
estimated from 4 ps of Langevin dynamics with two opposed biases.

VERDICT. KEPT 5 (the LTR framing as a matched case-control design; the high-quality-loss paper as
corroboration of the novelty claim; the quasi-harmonic bias direction; S14 and S12 as topic 7's
honest prior). NOTED 1 (LambdaRank metric weighting -- a loss refinement, excluded by rule 14).
RECORDED 4 (LTR does not address the confound case; my two-bias reconciliation; my
entropy-is-compactness-like hypothesis; the self-fulfilling-null warning). REJECTED 2 (the
decoy-LTR and ProQ4 pairwise imports -- same features, different loss, gains measured between
bands). No importable operator. Multiplicity: 0 endpoint comparisons; no new measurement in this
entry; every number is cited to its sprint or its paper.
Artefacts: `s29/lit/L_7_inband_training.md`; `s29/lit/L_INDEX.md`.

## S29-L20 -- RUNG 6, THE TYPICALITY-AXIS PROBE (H1): THE FALSIFIER FIRES ON BOTH CLAUSES AND THE KILL IS STRONGER THAN THE ONE REGISTERED -- THE ORACLE COSINE BETWEEN THE CONDITIONED-MINUS-BLIND AXIS AND THE NATIVE DIRECTION IS -0.058 (BELOW THE RANDOM-FIELD |cos| 0.144 AT 2.19x MDE, INDISTINGUISHABLE FROM A RANDOM SIGNED DIRECTION AT 0.73x), THE ORACLE BEST GLOBAL STEP IS EXACTLY t = 0 (ONE SCALAR WITH THE NATIVE IN HAND CANNOT BEAT DOING NOTHING), AND THE LEAVE-FOLD-OUT STEP IS BIT-IDENTICAL TO PRODUCTION ON 126/126 TARGETS ON THE BUILT CHAIN; ON FAIL18 THE AXIS POINTS AWAY FROM THE NATIVE HARDER (-0.317 vs -0.015, RANDOM-18 p 0.0006), THE OPPOSITE SIGN TO THE REGISTERED PRIOR (2026-09-20 00:30, O)

H1 (`s29/STATE.md`) IS DEAD ON THE MEASUREMENT IT NAMED. Pre-registration `s29/PREREG_S29_O.md`
sections 4 and 5 (written and committed at 1e9bb035 before the first number). EVERY NUMBER IN THIS
ENTRY IS ORACLE (it reads `nat_ca`) EXCEPT the leave-fold-out arm, which is the lane's one
DEPLOYABLE contrast and is labelled so in its own block.

Question. Per target, in the production cloud's own frame: u = (shipped top-75 average) minus
(sequence-blind average), v = (native) minus (shipped average), both with the 6-dim rigid-body
tangent space removed. Does u point toward the native, and does one scalar step along it help?

Falsifiers, registered. F6a: the fold CI of the mean ORACLE cos(u, v) includes the random
reference, or the mean cos is at or below the random fields' mean |cos| with the fold CI including
it. F6b: the leave-fold-out step does not clear 0.7x its own MDE with the fold CI excluding zero
in the improving direction. Registered prior (S24 L2/L3's geometry, the coordinator's own stated
objection): cos near zero or negative, the step null. STATE.md's alternative prior: cos +0.2 to
+0.4 on the 108, higher on FAIL18.

The blind clouds, rebuilt from S24's OWN stable RNGs and GATED against S24's artefacts (not new
draws): LIB75 = `s24/biasalign.py`'s C1, a uniform 75-window draw from the target's leakage-safe
universe (no retrieval, no score: the "typical peptide of this length"), ORACLE mean 3.8055 A;
BPRIME = `s24/qmatch.py`'s B', the top-75 by the shipped score of 2,000 uniform windows with the
BLOSUM pool excluded, ORACLE mean 3.1128 A. Gate: each cloud's RMSD, its native-frame bias cosine
with production and production's own RMSD must equal `s24/results/biasalign.json` /
`s24/results/qmatch.json` to 1e-6 on EVERY target; max deviation over 126 x 2 gates = 0.0
(the rebuild is exact, `s29/results/s29_O_cloud_rows.jsonl :: axis.<def>.gate`). Production's
cloud reproduces S28-L1b's frame to 2.1e-14 and the pinned 3.048338 over 126.

### F6a: the cosine. FIRES.
ORACLE mean cos(u, v) over 126, PRIMARY definition LIB75: **-0.0584** (SE 0.0324, median -0.0804);
positive on 55/126; above its own target's random-field |cos| on 37/126. The random reference is
16 Gaussian shape fields per target in the same 3n-6 space: signed mean +0.0075, mean |cos| 0.1443
(analytic sqrt(2/(pi(3n-6))) 0.1413 -- measured and analytic agree, `control-must-match-the-operators-space`).
```
  cos(u,v) vs random-field signed cos [LIB75]
    a -0.0584 (med -0.0804)   b 0.0075 (med 0.0067)   n=126
    effect -0.0659   median -0.0830   SE 0.0324   MDE 0.0908   effect/MDE -0.73
    iid  CI95 [-0.1272, -0.0036]
    fold CI95 [-0.0989, -0.0320]   folds same sign 5/5   per-fold 0:-0.014 1:-0.082 2:-0.124 3:-0.078 4:-0.039
    72W/54L/0T   worst degradation +0.7355 (8TXS)   p90 +0.3854   power 0.53  Type-M 1.37
    concentration: drop-top10 -0.0135 vs uniform-effect null p10/p50/p90 -0.0560/-0.0137/+0.0298 -> pctile 0.502
    VERDICT: NOT MEASURED (|effect| 0.0659 <= its own MDE 0.0908, 0.73x)
  cos(u,v) vs random-field |cos| [LIB75]
    a -0.0584 (med -0.0804)   b 0.1443 (med 0.1423)   n=126
    effect -0.2027   median -0.2024   SE 0.0330   MDE 0.0924   effect/MDE -2.19
    iid  CI95 [-0.2669, -0.1389]
    fold CI95 [-0.2355, -0.1663]   folds same sign 5/5   per-fold 0:-0.142 1:-0.216 2:-0.256 3:-0.223 4:-0.183
    89W/37L/0T   worst degradation +0.6383 (8TXS)   p90 +0.2696   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.1492 vs uniform-effect null p10/p50/p90 -0.1925/-0.1480/-0.1044 -> pctile 0.489
    VERDICT: BETTER
```
Read: against the random fields' SIGNED cosine the axis is indistinguishable (0.73x MDE, NOT
MEASURED); against their |cos| it is WORSE at 2.19x MDE with the fold CI [-0.2355, -0.1663] entirely
below zero and 5/5 folds. Either clause of F6a fires. The axis carries no more direction toward
the native than a random shape field, and its signed mean is on the wrong side of zero.
Secondary definition BPRIME: cos -0.0087 (SE 0.0309), same picture (vs signed 0.12x, vs |cos| -1.81x,
5/5 folds), and |u| is only 2.348 A against LIB75's 8.949 A (S24 L3 measured B's bias cosine with
production at 0.943: there is almost no axis there to step along).

### F6b: the step. FIRES, and harder than registered.
The ORACLE mean curve over t (t = -1 is the blind average, t = 0 production, t > 0 the
extrapolation S24 never ran), LIB75:
```
t      -1.0   -0.8   -0.6   -0.4   -0.2    0.0   +0.2   +0.4   +0.6   +0.8   +1.0   +1.5   +2.0
RMSD 3.806  3.572  3.367  3.201  3.089  3.048  3.107  3.256  3.475  3.744  4.047  4.915  5.875
```
**The ORACLE best GLOBAL step is t = -0.0, i.e. production itself, ORACLE mean 3.0483 = the
production anchor exactly.** Interpolating toward the blind answer is monotonically worse (S24 L2's
mixture curve, reproduced from a different construction); extrapolating away from it is worse
faster. The curve's minimum over 31 grid values IS the do-nothing point. That is a stronger
statement than the registered falsifier: not "a deployable step fails", but "no step exists,
with the native in hand, at one scalar for all targets".
Every leave-fold-out choice is therefore t = 0 on all five folds (ties 1 each), the DEPLOYABLE
structure is bit-identical to production on 126/126 targets (max |deviation| 0.0), and both
contrasts are exact zeros:
```
  lfo_LIB75 vs prod (POINT CLOUD, DEPLOYABLE step)
    a 3.0483 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  lfo_LIB75 vs prod (BUILT CHAIN, DEPLOYABLE step)
    a 3.2105 (med 2.9661)   b 3.2105 (med 2.9661)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
```
(The 126 ties are the correct reading of a NOT MEASURED verdict here: the operator emitted
production. Chain price +0.1622 for both arms, identically, because they are the same cloud.)
Secondary BPRIME: leave-fold-out t per fold {'0': -0.0, '1': -0.0, '2': -0.2, '3': -0.1, '4': -0.0}; the arm moves 48/126 clouds and is
**+0.0042 A WORSE on the point cloud (0.45x MDE) and +0.0037 A worse on the built chain (0.29x MDE)**,
NOT MEASURED in either basis, fold CI straddling zero, 1/5 folds same sign:
```
  lfo_BPRIME vs prod (BUILT CHAIN, DEPLOYABLE step)
    a 3.2142 (med 3.0222)   b 3.2105 (med 2.9661)   n=126
    effect +0.0037   median +0.0000   SE 0.0045   MDE 0.0126   effect/MDE +0.29
    iid  CI95 [-0.0047, +0.0128]
    fold CI95 [-0.0040, +0.0145]   folds same sign 1/5   per-fold 0:+0.000 1:+0.000 2:+0.025 3:-0.007 4:+0.000
    25W/23L/78T   worst degradation +0.2693 (7BX2)   p90 +0.0204   power 0.13  Type-M 2.97
    concentration: drop-top10 +0.0110 vs uniform-effect null p10/p50/p90 +0.0054/+0.0106/+0.0165 -> pctile 0.545
    VERDICT: NOT MEASURED (|effect| 0.0037 <= its own MDE 0.0126, 0.29x)
```

### The mechanism, measured beside the outcome (contract rule 18)
Per target the optimal step is t* = (|v|/|u|) cos(u, v) for a locally linear RMSD, so the sign of
t* must follow the sign of cos. It does: **rho(cos(u,v), t*) = +0.810** (LIB75; +0.865 for BPRIME).
And the signs cancel across targets -- 54% of targets want t* < 0 (toward the typical answer),
36% want t* > 0, 10% want exactly zero -- which is precisely why the global step is zero.
The per-target ORACLE step is worth 2.7422 A (-0.3061 vs production) and is an order statistic over 31
grid values whose transferable part is nil: choosing the global t on half the targets and applying
it to the other half is worth **+0.0056 A [+0.0000, +0.0346]** over 2,000 splits, i.e. zero. This is
S23 L6's scale result again from a new direction (a per-target correction that is real, ORACLE,
and provably not a property of anything native-free), and S16's "every arm chose do nothing" with
the arm's own optimum now located rather than assumed.

### The regime question, answered with the sign REVERSED from the prior
STATE.md registered "cos 0.2 to 0.4 on the 108, higher on FAIL18". Measured (ORACLE): the 108 sit
at -0.0153 and FAIL18 at **-0.3171** -- the axis points AWAY from the native three times harder where
the distogram fails. Against the registered random-18 null (20,000 random 18-subsets of the 126,
one-sided in the observed direction) that IS a set property: **p = 0.0006** [null 2.5/97.5
-0.2118, +0.0984], and it survives the Bonferroni bar over this lane's 2 stratum tests. BPRIME's
FAIL18 cos -0.0959 vs +0.0058 on the 108 is NOT a set property (p 0.123) and is not quoted as one.
The mechanism is already in the record: on FAIL18 the sequence-BLIND pipeline BEATS the shipped one
(5.425 vs 6.019, `sequence-conditioning-hurts-the-failures`, S12 `coord_null`), so on those targets
production is FURTHER from the native than the typical answer is and u = production - blind points
outward by construction. H1's premise -- "the conditioned answer lies between the typical answer
and the native" -- is false exactly where it was hoped to be most true.

### Verdict
**H1 IS FALSIFIED on both registered clauses, on the PRIMARY and the SECONDARY blind definition,
on the point cloud and on the built chain.** Nothing about the typicality axis is deployable and
nothing about it is even ORACLE-usable at one scalar. The coordinator's own stated objection
(S24 L2/L3's geometry) is confirmed: the axis has a component along the shared error and a large
component that is minus the blind cloud's orthogonal error. The leading hypothesis becomes H0.
What is NOT closed by this entry: a per-target step exists (-0.3061 A ORACLE) and its SIGN is the
whole content -- any future operator claiming this axis must supply the per-target sign natively,
and this lane's random-18-nulled measurement says the sign correlates with the regime (FAIL18 vs
the 108), which is the one thing the record says is not predictable native-free (S22 L7, S23 L7,
S28-L39's three routers).

Comparisons this entry: 4 deployable endpoint contrasts (2 blind definitions x {point cloud,
built chain}); 2 stratum tests, each against the random-18 null; every other number ORACLE
diagnostic. Multiplicity: all four are NOT MEASURED, two of them exact zeros, so no max-over-K
bar applies.
Artefacts: `s29/results/s29_O_cloud_rows.jsonl` (126 rows, gates included),
`s29/results/s29_O_lfo.json` (the curves, the fold choices, every ST block),
`s29/results/s29_O_chain_rows.jsonl` (prod / lfo_LIB75 / lfo_BPRIME, 378 projections),
`s29/results/s29_O_structs/<pdb>.npz`; code `s29/s29_O_ladder.py`, tests `tests/test_s29_O.py`
(11 pass, including the leave-fold-out NaN-poison and the rigid-body removal); prereg
`s29/PREREG_S29_O.md`; jobs `s26/jobs_done/s29O_{probe1,probe2,cloud126,lfo,chainA}.json`
(probe 1KWE,1KZ2,1LB7 peak RSS 0.308 GB; the 126-target cloud 300 s at peak RSS 0.321 GB; chain
group A 378 projections, peak RSS 0.108 GB).
Chain anchor, reported not gated: production re-projected in this job is **3.2105** over 126 against
S28-L26b's 3.2071 (mean difference +0.00347, 18/126 targets differ by more than 0.02 A, 3 by more
than 0.1, max 0.176) and S27's DIS 3.2126 (-0.00209, max 0.517). That is S28-L18's branch-flip floor
of the multi-start projection under a 1e-14 input perturbation, it is identical on both sides of
every contrast here (same code path, same job), and it is why S27's rows are not used as the
comparator.

## S29-L21 -- RUNG 8, LANE T's ONE-PARAMETER FAMILY (S29-L11 PREDICTION 4): T's PREDICTION IS CONFIRMED AT ITS FLOOR -- THE ORACLE BEST GLOBAL eta ALONG THE POOL's FIRST SHAPE MODE IS EXACTLY ZERO (0.000 A, NOT "UNDER 0.15"), THE BEST PER-TARGET eta IS POSITIVE ON 52% OF TARGETS SO THE SIGN IS A COIN FLIP EXACTLY AS T DERIVED, THE ONE-GLOBAL-SIGN CEILINGS ARE -0.214 / -0.246 A AND STILL PER-TARGET IN MAGNITUDE, AND THE LEAVE-FOLD-OUT eta IS +0.0071 A WORSE THAN PRODUCTION (0.82x MDE); LANE B's CENTERED-HAMILTONIAN BUILD HAS AN ACHIEVABLE CEILING OF ZERO THROUGH THIS READOUT (2026-09-20 00:30, O)

Question, set by the coordinator from lane T's S29-L11 prediction 4: a centered or agreement-matrix
Hamiltonian consumed by a signed readout collapses to X(eta) = production + eta*PC1, PC1 the pool's
first shape mode. What is that family's ceiling? T predicted under 0.15 A better than production
and flagged above 0.30 A as "lane B's build is worth much more than the theory says".
PC1 IS NATIVE-FREE (`s29/s29_O_ladder.py :: pool_pc1`): the top right-singular vector of the
top-75 members' deviations from the production average, in the average's own medoid frame, with
the rigid-body component removed, scaled so that ||PC1|| = sqrt(n) and eta is in ANGSTROM of
point-cloud RMSD displacement. Its SIGN is fixed by a deterministic native-free convention (the
largest-|component| entry positive), which is exactly T's point: the marginals cannot supply the
sign. eta grid {-3.0, ..., +3.0} step 0.1 (61 values), fixed in the module before the first number.
EVERY eta BELOW IS ORACLE except the leave-fold-out arm.

The pool's own spread along PC1 is 1.068 A (sd of the members' projections, a real and large
direction) and PC1 carries 36.3% of the members' deviation variance. So the family is not
degenerate -- it is informative-looking and empty.

| arm | mean (point cloud) | vs production 3.0483 |
|---|---|---|
| ORACLE best GLOBAL eta (one scalar, 126 targets) | **3.0483** at eta = +0.0 | **+0.0000** |
| ORACLE best per-target eta (an ORDER STATISTIC over 61) | 2.5940 | -0.4543 |
| ORACLE best per-target \|eta\| with the sign forced POSITIVE | 2.8341 | -0.2142 |
| ORACLE best per-target \|eta\| with the sign forced NEGATIVE | 2.8028 | -0.2455 |
| leave-fold-out eta (DEPLOYABLE; eta per fold {'0': -0.1, '1': 0.0, '2': 0.0, '3': 0.1, '4': 0.1}) | 3.0554 | **+0.0071** |

```
  ORACLE PC1 global eta vs prod (POINT CLOUD)
    a 3.0483 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  ORACLE PC1 per-target eta vs prod (POINT CLOUD, an order statistic)
    a 2.5940 (med 2.3715)   b 3.0483 (med 2.8373)   n=126
    effect -0.4543   median -0.1659   SE 0.0527   MDE 0.1476   effect/MDE -3.08
    iid  CI95 [-0.5588, -0.3539]
    fold CI95 [-0.5319, -0.4102]   folds same sign 5/5   per-fold 0:-0.405 1:-0.613 2:-0.409 3:-0.425 4:-0.434
    122W/0L/4T   worst degradation +0.0000 (1D0W)   p90 -0.0071   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.3285 vs uniform-effect null p10/p50/p90 -0.3981/-0.3292/-0.2663 -> pctile 0.505
    VERDICT: BETTER
```
The deployable arm, the only non-ORACLE row:
```
  PC1 one-parameter family, LEAVE-FOLD-OUT eta vs prod (POINT CLOUD, DEPLOYABLE)
    a 3.0554 (med 2.8384)   b 3.0483 (med 2.8373)   n=126
    effect +0.0071   median +0.0000   SE 0.0031   MDE 0.0086   effect/MDE +0.82
    iid  CI95 [+0.0012, +0.0130]
    fold CI95 [+0.0011, +0.0125]   folds same sign 3/5   per-fold 0:+0.016 1:+0.000 2:+0.000 3:+0.006 4:+0.012
    34W/44L/48T   worst degradation +0.0924 (2MK7)   p90 +0.0626   power 0.64  Type-M 1.26
    concentration: drop-top10 +0.0124 vs uniform-effect null p10/p50/p90 +0.0084/+0.0123/+0.0162 -> pctile 0.514
    VERDICT: NOT MEASURED (|effect| 0.0071 <= its own MDE 0.0086, 0.82x)
```

Reading. (1) **The achievable ceiling of the one-parameter family is 0.0000 A.** The mean curve's
minimum over 61 values of eta is at eta = 0: with the native in hand, one global step along the
pool's first shape mode cannot beat doing nothing. T predicted "under 0.15"; the answer is the
floor of that interval, and the 0.30 A flag does not fire on any arm a deployable operator could
imitate. (2) **The sign is a coin flip, measured**: the ORACLE best eta is positive on 52% of
targets (median |eta| 0.9 A), so the two one-global-sign ceilings (-0.2142 and -0.2455) are the price
of guessing it, and both of those still choose the MAGNITUDE per target. (3) The per-target eta is
worth -0.4543 A ORACLE and is an order statistic: `best_of_k_within` reads 190% accounted by its
across-target null at k_eff 36.1, and the honest transfer measurement -- the leave-fold-out arm in
the last row -- is **+0.0071 A, i.e. worse than production**, 0.82x MDE, NOT MEASURED, 3/5 folds,
48/126 targets untouched. Nothing about eta transfers between targets.
(4) Mechanism beside outcome: this is the same shape as rung 6 (S29-L20) and as S23 L6's
per-target scale -- a real per-target direction whose sign and magnitude are properties of the
(pool, native) PAIR and carry no native-free content. Two independent one-parameter families
(the typicality axis and the pool's first shape mode) both have an ORACLE global optimum of
EXACTLY do-nothing on this instrument.

Verdict for lane B (the coordinator's gate): **a centered / agreement-matrix Hamiltonian consumed
by a signed readout that collapses to this family has an achievable ceiling of 0.000 A through the
deployed readout, and its ORACLE per-target ceiling (2.594 A) is entirely sign-and-magnitude
information that no arm in this project's record can supply.** T's prediction 4 stands as derived.
Built chain NOT taken: the coordinator's condition was "only if it clears 0.15 A"; the achievable
arm clears nothing and is on the wrong side of zero, so 126 projections are not spent on it.

Comparisons this entry: 1 deployable endpoint contrast (leave-fold-out eta vs production, point
cloud), NOT MEASURED; the rest ORACLE diagnostics. Artefacts:
`s29/results/s29_O_pc1_rows.jsonl` (126 rows: the per-target curve, PC1's sd and variance share),
`s29/results/s29_O_pc1.json` (every ST block and the fold choices); code
`s29/s29_O_ladder.py :: pool_pc1, pc1_row, analyse_pc1`; job `s26/jobs_done/s29O_pc1.json`
(331 s, peak RSS 0.266 GB). The rung is an addendum to `s29/PREREG_S29_O.md` (section 3 as
amended by the coordinator's instruction of 2026-09-20; the eta grid, the sign convention and both
one-global-sign arms were fixed in code before the first number, and this entry is the first
number).

## S29-L22 -- LANE P, THE PROJECTION-STAGE MECHANISM (the 126-arm run is still in flight; these six do not depend on it): PRODUCTION'S CLOUD REPRODUCES BIT-EXACTLY ON 126/126 AND ITS BUILT CHAIN ON 6/6 AT EXACTLY 0.0; THE +0.164 A PROJECTION PRICE TRACKS THE CLOUD'S CONTRACTION AT SPEARMAN +0.604 (+0.481 PARTIALLED); THE AVERAGING DISTORTION IS NOT A CONTRACTION BUT A MONOTONE SHAPE DISTORTION IN SEQUENCE SEPARATION (0.773 AT |i-j|=1, CROSSING 1.00 NEAR 8, 1.10 AT 13); CORRECTING THAT PROFILE IS REFUTED AT +0.578 A (2.83x MDE, 5/5 FOLDS) WITH AN ORACLE-FITTED PROFILE NO BETTER (+0.582) AND AN EXACT r==1 CONTROL; AND S23 L2/L3 REPRODUCE ON THESE CLOUDS (2026-09-20 00:37, P)

Question (`s29/briefs/S29P.md`, wave-2 probe P1; prereg `s29/PREREG_S29_P.md`, written 00:10
before any number, addenda 1-4). Production's point cloud is 3.0483 A and its built chain is
3.2126: the projection COSTS +0.1643 A, the largest single-stage loss in the record and the only
stage no S28 lane touched. The cloud's mean adjacent CA-CA distance is 2.9614 A against the ideal
3.80, so the least-squares fit of a rigid-length chain to it is a BIASED fit. Does removing the
geometric inconsistency BEFORE the projection change the built chain? This entry posts the six
results that do NOT depend on the 126 x 26 projection run (still queued/running behind the
cap-8 limit); the arm verdicts and the ORACLE s-curve follow in a second entry.

WHY THIS IS NOT S23 L6 (contract rule 10; stated in the prereg before any number).
`s23/LEDGER.md` L6(d), lines 165-176: "s* is not a property of the target's fold. It is fitting
the particular (pool, reference-structure) pair ... unreachable IN PRINCIPLE from native-free
information ... because the thing it depends on is the answer." Four differences:
(1) L6's s* is read off the NATIVE; g = 3.80 / (the cloud's own mean adjacent CA-CA distance) is
    read off the cloud and one covalent constant.
(2) L6 measured scale ON THE CLOUD, where scale IS the answer (Kabsch does not fit scale;
    `s23/LEDGER.md` line 40 gives the closed form s* = Ct/Cc). This lane measures ON THE BUILT
    CHAIN, whose CA-CA distance is 3.80 by construction, so an input rescale CANNOT rescale the
    output -- only change which ideal-geometry chain is nearest.
(3) L6's lever multiplied the emitted structure; g multiplies the INPUT to a non-convex
    multi-start optimisation with a branch degeneracy (S26 L88; S28-L18 / S28-L27b).
(4) L6's own L1 (lines 14-19) names the option it did not test: "Only changing the space you
    average in (torsions) or repairing geometry afterwards can". Repairing it BEFORE the
    projection was never run.
L6 STANDS and this entry does not touch it; every ORACLE quantity below is labelled ORACLE.

(a) THE CODE PATH, AND IT IS EXACT. `s29_P_scale.py :: production_cloud` reproduces
`s27/s28_B_prodcheck.py :: project_production` (channels -> DIS energy -> lexsort tie-key top-75
-> `readout_uniform`). All 126 clouds reproduce `s27/results/chain_rows.jsonl :: DIS`'s
`rmsd_cloud` with difference EXACTLY 0.0 (`s29_P_factors.json :: n_cloud_mismatch` 0). The
pre-registered PROBE GATE -- arm PROD reproduces the production BUILT CHAIN to < 1e-9 on the six
registered targets 1A13, 1CS9, 2LNG, 2NB7, 9BFL, 9BAF -- PASSES at max |diff| **0.000e+00 on
6/6** (`s29_P_probe.json`). Job `s26/jobs_done/s29P_probe6.json` exit 0, wall 1243.9 s (about 9
min governor-suspended), **peak RSS 0.296 GB**, 156 cells, 7.95 s/projection under 8-job
contention. So both sides of every contrast in this lane share one code path and one input, and
the S28-L18 / S28-L27b CROSS-PATH floor does not apply; the INPUT-PERTURBATION floor is measured
separately (arm FLOOR, second entry).

(b) THE PREMISE HOLDS: THE PRICE IS NOT A FLAT TOLL. Projection price per target (built chain
minus point cloud, `chain_rows.jsonl :: DIS`): mean **+0.1643**, sd 0.2014, median +0.0958,
positive on 111/126, range -0.347 to +0.652. Spearman with the native-free contraction factor
g: **+0.604**; partialling out the cloud's own RMSD and the chain length, **+0.481**. (s_SPAN
+0.378, s_ISO +0.339, length -0.142, cloud RMSD +0.511, ORACLE s* -0.313.) The projection costs
most exactly where the average is most geometrically inconsistent with the space it is projected
into. This does not establish that removing the inconsistency removes the price.

(c) [ORACLE] S23 L2 AND L3 REPRODUCE ON THESE CLOUDS, AND THE NATIVE-FREE FACTORS POINT THE
OTHER WAY. L6's closed form s* = sum(svd(Cc^T Tc)) / ||Cc||^2 on the 126 cached clouds
(`s29_P_oracle_cloud_sstar.json`; no projection): s* mean **0.9481**, sd 0.2565,
range 0.277 to 1.924, **42.9% want EXPANSION** -- S23 L2 reported
53/126 = 42.1% wanting expansion against 73 wanting contraction; the split, the spread and the
bimodality all reproduce. Spearman(s*, g) **-0.117**, (s*, s_SPAN) -0.033, (s*, s_ISO) +0.000:
inside the +-0.11 band S23 L3 measured for eight other native-free candidates, and the only one
outside it has the WRONG SIGN. Mechanically: g is an expansion (mean 1.352, above 1 on
125/126) while the RMSD-optimal cloud scale is a contraction on 57% of targets. On the POINT-CLOUD
basis BOND must be harmful; the lane's question survives only because the built chain's scale is
fixed by covalent geometry, which is exactly the distinction above.

(d) THE AVERAGING DISTORTION IS NOT A CONTRACTION. IT IS A MONOTONE SHAPE DISTORTION IN SEQUENCE
SEPARATION. `s29_P_contraction_profile.json`, all 126 clouds, every pair pooled at each
separation:

    |i-j|                       1      2      3      4      5      6      7      8      9     10     11     12     13
    cloud/native [ORACLE]   0.773  0.824  0.830  0.870  0.921  0.948  0.968  0.996  1.029  1.047  1.065  1.083  1.102
      median ratio          0.784  0.857  0.886  0.912  0.960  0.975  0.981  0.997  1.009  1.013  1.029  1.038  1.048
    cloud/posterior (NF)    0.776  0.821  0.809  0.846  0.902  0.913  0.924  0.955  0.965  0.967  0.971  1.001  0.952

The curve is monotone, starts **23% SHORT** at the virtual bond, **CROSSES 1.00 near |i-j| = 8**
and ends slightly LONG (1.102 at 13 by ratio of means, 1.048 by median ratio). The NATIVE-FREE
profile, taken against the distogram posterior's own L1-Bayes median map, has the same shape and
the same crossing region, so this is not an ORACLE-only statement. Cloud mean Rg 6.2061 against
the native's 6.6009: S23 L1's two numbers (bond 22% short, envelope 6% short) are the two ENDS
of this curve, and the middle is where the sign changes.
Consequences: BOND's g = 1/0.773 = 1.29 sets the bond exactly right and inflates every
separation beyond 8 by ~29% on top of distances that are ALREADY long, so BOND is predicted
harmful -- the registered null prior's mechanism, measured rather than argued. The best SINGLE
scalar in least squares must sit near the middle of the curve, ~1.05 to 1.08; s_ISO is 1.077 and
s_SPAN is 1.102. This is also the direct measurement of the Jensen mechanism lane L derives in
S29-L12: averaging shrinks a distance by more when the pool disagrees more about it, and the
pool disagrees most about local geometry relative to that distance's size.
METHOD NOTE (my own error, caught before use): my first pass used the MEAN OF PER-PAIR RATIOS and
read 1.37 at sep 13. That statistic is dominated by small denominators (a native hairpin puts two
residues 13 apart close in space). The ratio of means and the median ratio agree at 1.10 and 1.05
and are the honest numbers; 1.37 is withdrawn and was never used. The shape conclusion is
identical under all three statistics.

(e) CORRECTING THAT PROFILE IS REFUTED, AND AN ORACLE-FITTED PROFILE IS NO BETTER.
Falsifier for this sub-arm, registered by construction: the correction should help if the
distortion is what costs the accuracy. Operator: divide each cloud's distance map by the
leave-fold-out population profile r(|i-j|), re-embed by classical MDS (double-centring, top-3
eigenvectors), score. POINT-CLOUD BASIS, diagnostic; the built chain is not measured for this arm.
Means: production cloud **3.0483**; corrected NATIVE-FREE (posterior-fitted) **3.6266**;
corrected [ORACLE] (fitted on the OTHER FOLDS' real natives) **3.6307**.

  P sep-profile NATIVE-FREE - production (POINT CLOUD, diagnostic)
    a 3.6266 (med 3.3807)   b 3.0483 (med 2.8373)   n=126
    effect +0.5782   median +0.4351   SE 0.0730   MDE 0.2045   effect/MDE +2.83
    iid  CI95 [+0.4390, +0.7259]
    fold CI95 [+0.4165, +0.7230]   folds same sign 5/5   per-fold 0:+0.410 1:+0.549 2:+0.793 3:+0.357 4:+0.732
    27W/99L/0T   worst degradation +3.1148 (1S9Z)   p90 +1.8985   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.6746 vs uniform-effect null p10/p50/p90 +0.5784/+0.6714/+0.7745 -> pctile 0.515
    VERDICT: WORSE

  P sep-profile [ORACLE, leave-fold-out native-fitted] - production (POINT CLOUD, diagnostic)
    a 3.6307 (med 3.3521)   b 3.0483 (med 2.8373)   n=126
    effect +0.5824   median +0.4689   SE 0.0813   MDE 0.2278   effect/MDE +2.56
    iid  CI95 [+0.4257, +0.7424]
    fold CI95 [+0.4343, +0.7425]   folds same sign 5/5   per-fold 0:+0.586 1:+0.846 2:+0.719 3:+0.478 4:+0.344
    37W/89L/0T   worst degradation +2.8826 (2MJQ)   p90 +1.9322   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.6975 vs uniform-effect null p10/p50/p90 +0.5874/+0.6940/+0.8079 -> pctile 0.512
    VERDICT: WORSE

THE CONTROL IS EXACT, which is what makes this readable: with r == 1 the same MDS returns the
production cloud at **0.00e+00 A**, so the +0.578 is the correction and not the re-embedding.
And the ORACLE arm -- the best population profile that exists -- is as bad as the native-free one
(+0.5824 vs +0.5782), which rules out "the profile was mis-estimated".
Reading: undistorting the distance map destroys more than the distortion costs, because the
coordinate average's value is the COHERENCE of its errors (68% common-mode, S23 L9) and a
per-pair correction fitted to a population profile breaks that coherence while replacing it with
nothing. SCOPE: this closes the POPULATION-PROFILE class of distance-space corrections. It says
nothing about a per-target profile, which S23 L6 already places out of reach in principle. A
fortiori it strengthens this lane's registered null prior for the scalar arms: if a 13-parameter
profile with an oracle fit cannot help, one scalar drawn from the same curve will not.

(f) WHAT IS STILL IN FLIGHT. The 126 x 26 projection run (arms PROD, BOND, SPAN, ISO, CTRL-INV,
CTRL-GLOBAL, CTRL-LAM, BOND-LAMFIX, 8 x CTRL-RAND, FLOOR, then the 9-point ORACLE s-grid), as 4
concurrent shards with a per-(arm, target) checkpoint. It carries the registered falsifiers
F-P1/F-P2/F-P3, the branch-flip floor on this path, the FAIL18 / 108 split with a random-18
null, and -- prereg addendum 4 -- the WIDER MULTI-START readout MS-OBJ/MS-MEAN/MS-ORACLE, since
every arm's chain is a feasible point of the shipped projection problem and `obj0` is recorded
for each. Prediction registered: MS-OBJ lowers production's own objective and NOT the RMSD.

VERDICT. (a) instrument property, EXACT. (b) the lane's premise HOLDS. (c) S23 L2/L3 REPRODUCED
and the native-free factors are orthogonal to the ORACLE scale. (d) NEW MECHANISM: the distortion
is a monotone shape error in sequence separation, not a contraction. (e) REFUTED for the
population-profile class, with an ORACLE-fitted control that is no better and an exact r == 1
control. Registered null prior for the scalar arms: unchanged and strengthened.
Multiplicity: **2 endpoint-style comparisons** in this entry, both on the POINT CLOUD and both
declared diagnostic (the sep-profile arms); **0 built-chain endpoint comparisons** -- those are
in the second entry, where the lane's budget is K = 6.
Artefacts: `s29/PREREG_S29_P.md`, `s29/s29_P_scale.py`, `tests/test_s29_P.py` (20 pass),
`s29/results/s29_P_factors.json`, `s29_P_probe.json`, `s29_P_rows_probe.jsonl`,
`s29_P_oracle_cloud_sstar.json`, `s29_P_contraction_profile.json`, `s29_P_sepprofile_cloud.json`;
jobs `s26/jobs_done/s29P_factors.json`, `s26/jobs_done/s29P_probe6.json`;
findings `s29/s29_P_FINDINGS.md` P1-P10.

## S29-L23 -- THEORY SECTION 8, THE ACHIEVABLE NATIVE-FREE BOUND: EVERY OPERATOR IS A DISPLACEMENT AND ITS VALUE IS ONE COSINE, RMSD = RMSD_prod sqrt(1 - rho^2), SO 2.5 A COSTS rho = 0.628 AGAINST A MEASURED 0.04 AND A RANDOM-FIELD 0.140; THE BOUND ASSEMBLES TO >= 3.18 A WITH 2.98 A AS THE PERFECT-SIGN CEILING ON THE BEST FIELD ANYONE HAS BUILT, AND THE CEILING IS THEREFORE NOT THE POOL, THE CIRCUIT, THE OPTIMISER OR THE READOUT BUT THE ABSENCE OF A NATIVE-FREE VECTOR WITH A COSINE ABOVE 0.14 (2026-09-20 00:43, T)

Question (the coordinator, 00:35): given exactly what the system has -- 500 real windows with a
68% common-mode error, the 17-bin posterior's marginals, and no channel carrying covariance with
the native's deviation from typical -- what is the best expected BUILT-CHAIN RMSD attainable by
ANY native-free operator? Derivation from pieces already on the board; `s29/THEORY.md` section 8;
12 tests in `tests/test_s29_T.py` including the identity and the price table.

THE FORM. Write any operator's output as A = c + u with u native-free. Exactly,
|A - t|^2 = |e|^2(1 - 2 rho s + s^2) with e = t - c, rho = cos(u, e), s = |u|/|e|; minimising over
the step gives s = rho and
    RMSD_achievable = RMSD_prod * sqrt(1 - rho_max^2),
rho_max the best EXPECTED cosine of a native-free displacement field with the direction to the
native. Selection, re-weighting, a gradient step, a basin average, a signed readout, a projection
and a quantum tail are all displacements, so this collapses every operator class -- and every rung
of lane O's ladder -- onto one number. Inverted (RMSD_prod = 3.2126): 3.00 A needs rho 0.358,
2.50 needs 0.628, 2.31 needs 0.695, 1.71 needs 0.847. The random-shape-field reference is 0.140
(S28-L23b), so the charter's 2.5 A is 4.5 random directions. A CONFLATION TO AVOID: S14's "0.638
needed" is an in-band SPEARMAN for a 2.0 A target; 0.628 here is a COSINE in R^3N for a 2.5 A
target -- two different quantities and two different targets, never to be quoted as one threshold.

WHAT THE RECORD MEASURES FOR rho (all ORACLE diagnostics). Shipped objective's descent -0.034 (SE
0.021; -0.143 on FAIL18, S28-L23b); every other S27 channel with a gradient +0.034 to -0.006; the
typicality axis -0.058 (S29-L20); the pool's first shape mode PC1 per-target |rho| = 0.37 to 0.39
(inverted from lane O's ORACLE one-global-sign ceilings -0.214 / -0.246 A) with the sign correct
on 52% of targets, hence a signed rho of 0.015 and a leave-fold-out arm +0.0071 A WORSE than
production (S29-L21). The sign enters as gain = RMSD(1 - sqrt(1 - rho^2(2q-1)^2)): that formula
reproduces lane O's PC1 triple (magnitude, 52%, +0.0071) exactly, and it was derived before those
numbers were read, which is the bound's one end-to-end check.

THE BOUND. Every native-free field ever built here (|rho| <= 0.04) -> >= 3.210 A, a 0.003 A gain.
A field as strong as a RANDOM one WITH the sign right on every target (0.140) -> >= 3.181 A. The
best structured field known, PC1, with a PERFECT ORACLE per-target sign (0.37) -> >= 2.98 A. The
charter's 2.5 A needs 0.628. Three independent ceilings agree and do not add, because each is a
displacement already priced by its own cosine: enlarging or re-weighting the set is <= 0.008 A by
the Ueda-Nakano coefficient (S29-L8) and is capped in principle by Krogh-Vedelsby (any selection
that raises the members' mean quality destroys ambiguity, and S27 section 6 measured the two
cancelling with the wrong sign, +0.286 A for CONS); in-band ordering transfers at 0.600 across
targets with capacity saturated by a linear model (S14, S12, S29-L19); and the scale lever is
unreachable in principle because s* is a function of ebar alone (S23 L6/L9).
CENTRAL STATEMENT: **no native-free operator over the present information reaches below about
3.18 A on the built chain, the honest central estimate is 3.21 A (production), and the nearest
thing to a margin is 2.98 A, which requires a per-target SIGN for the pool's own principal mode
that is measured at chance.** The gap to 2.5 A is not a search gap, not expressivity (a
27-parameter family holds a 0.25 A structure on every target, S28-L26b), not aggregation (lane O's
ladder: the hull of the shipped top-75 contains a 2.00 A point and the pool's hull a 1.12 A point)
and not optimisation. It is one number.

ASSUMPTIONS. (B1) native-free; (B2) rho_max <= 0.14 for every field constructible from the present
information -- THE LOAD-BEARING ONE, derived for the marginal class by theorem 2 and measured
outside it (31 scorers S28-L48, 38 signals S12, eight routers, the typicality axis S29-L20, PC1
S29-L21); (B3) rigid-body components projected out and the identity read to first order in
s = |u|/|e|, where at s <= 0.14 the neglected term is under 2% of the gain; (B4) rho_max is an
expectation over the 126, so a field that is excellent on some targets and reversed on others
enters at |rho|(2q-1) -- buying the per-target sign is worth exactly as much as buying the
direction. A NEW SOURCE MUST BREAK (B2) AND ONLY (B2), and by section 7 the only classes that can
are a better distance prior (which moves e itself, slope -2.15 A per unit), a learned residual with
errors decorrelated from the predictor's (blocked by error coherence, S19), and a physics term on
the emitted structure in its FREE-energy form (outside the marginal class, so theorem 2 does not
bound it; its single-point form is measured worse than random, S25 L16). A second pool, the pool's
dispersion, an attention map and a joint all leave (B2) intact, which is why they price at 0.002
to 0.03 A whenever they are run.

THE FALSIFIER, ONE MEASUREMENT: exhibit a native-free displacement field whose mean cosine with
t - c over 126 targets exceeds 0.140 with the fold-clustered CI excluding it. It is meter number 2
and it runs in minutes, now with the mandatory shrink signature (S29-L10) so a cosine bought by
contracting toward typicality is rejected automatically. Secondary: (i) any operator whose
BUILT-CHAIN mean beats 3.181 A with a fold CI excluding production, since the identity says that
requires rho > 0.14 whether or not the cosine was measured; (ii) a per-target sign predictor for
PC1 with held-out accuracy above 0.58, which at |rho| = 0.37 buys 0.02 A and at 0.75 buys 0.09 A --
small, but it would prove the sign is learnable, the one clause of (B4) with no derivation behind
it.

FOR THE SPRINT. This is the charter's "decisive measurement of what imposes the ceiling" in the
form the charter asked for, and it ranks the remaining work: (i) lane D's band experiment, the one
route to a cosine never measured under a matched-realism control; (ii) section 4.5's structural
CVaR, the first formulation whose classical counterpart genuinely goes away, which by this bound
will not move the RMSD unless (i) succeeds first; (iii) the free-energy class. If (i) closes
negative, S29's headline is this bound with its table, and the next sprint's question is not an
architecture but whether the distance prior can be improved -- the only lever the record has ever
measured with a steep slope.
Multiplicity: 0 endpoint comparisons; 3 registered falsifiers. Every rho quoted is ORACLE and
named as such; nothing here tunes or selects anything.
Artefacts: `s29/THEORY.md` section 8 (commit 53c42bd4); `tests/test_s29_T.py` (12 pass);
inputs are S28-L23b, S29-L8, S29-L19, S29-L20, S29-L21, lane O's `s29/results/s29_O_cloud_rows.jsonl`
(the ORACLE ladder: hull_top75 1.9975, hull_pool 1.1167, best1_pool 1.7108, bestm 2.6062) and
S23 L6/L9.

## S29-L24 -- SECTION 9, ON LANE P's SEPARATION PROFILE (S29-L22): THE CURVE IS THE PRODUCT OF THE POOL'S INHERITED LONG-RANGE OVER-EXTENSION AND A FADING AVERAGING SHRINK, SO "CONTRACTION" IS THE LEFT-HAND END OF A SHEAR AND NOT A PROPERTY OF THE OPERATOR; AN ORACLE PROFILE CORRECTION CANNOT HELP BECAUSE THE PROFILE SUBSPACE IS MOSTLY UNREALISABLE; AND I WITHDRAW MY OWN SECTION 1.4 BUILD RECOMMENDATION -- SEPARATION-BAND RE-WEIGHTING IS A GENUINELY DIFFERENT OPERATION FROM P's RESCALE AND IS NEVERTHELESS CLOSED, MEASURED LEAVE-FOLD-OUT AT +0.010 A IN S12 AND RETIRED FROM THE SHIPPED CODE ONCE ALREADY (2026-09-20 00:43, T)

Question (the coordinator, 00:45): does the separation-band re-weighting of the OBJECTIVE that my
section 1.4 recommended survive lane P's refutation of the post-hoc cloud rescale, or is it the
same operation seen twice; and what is the one-sentence mechanism the report should use instead of
"contraction"? Derivation plus a record search; no job; `s29/THEORY.md` section 9.

THE MECHANISM. Section 1.2's identity compares the average to its own MEMBERS; lane P measured the
ratio to the NATIVE, and that ratio factorises exactly:
    d_ij(C)/d_ij(t) = [mean_k d_ij(W_k)^2]^(1/2)/d_ij(t)  x  sqrt(1 - s_ij^2/mean_k d_ij(W_k)^2),
i.e. POOL BIAS at that separation TIMES AVERAGING SHRINK at that separation. The shrink is <= 1
always and rises toward 1 with separation, because s^2 is the members' superposition residual and
does not grow with |i-j| while d^2 does. The pool bias is not an averaging effect at all: the pool
is selected by a posterior that S25 L1 measured as predicting distances systematically TOO LONG,
increasingly so with separation (signed error -0.048 A at 2-2 growing monotonically to -0.589 A at
9-15). The two factors move oppositely in separation and trade places near |i-j| = 8.
THE SENTENCE FOR THE REPORT: "Averaging does not contract the structure: by the exact identity
d(C)^2 = <d_k^2> - s^2 the emitted distance is the members' RMS distance less their own spread, so
the distortion is a SEPARATION-DEPENDENT SHEAR -- a large shrink where distances are short and the
spread is comparable to them (-23% at the virtual bond), fading to nothing where distances are
long -- multiplied by the pool's inherited over-extension at long range, which the fading shrink no
longer masks; the ratio therefore crosses 1.0 near |i-j| = 8 and ends above it, and 'the 22%
contraction' is the left-hand end of that curve, not a property of the operator."
A ZERO-COST CHECK FOR LANE P, and it is an identity so it must pass exactly: compute the two
factors separately on the clouds P already has. The shrink factor d_ij(C)/[mean_k d_ij(W_k)^2]^(1/2)
is NATIVE-FREE and must be <= 1 at every separation and monotonically rising; the pool-bias factor
is ORACLE and must carry the ENTIRE crossing; their product must reproduce P's measured profile to
floating point. If the crossing appears in the shrink factor, the identity is wrong and section 1.2
with it.
WHY P's ORACLE-FITTED PROFILE CORRECTION COULD NOT HELP, DERIVED: a per-separation multiplicative
profile is a vector in pair space, which has P coordinates against a realisable set of dimension
3N - 6 (55 against 30 at N = 12, section 2 C2). The only exactly realisable profile direction is
the UNIFORM one -- a global scale -- and that is closed twice over (+0.0002 A by nested CV, and
unreachable in principle per target because s* depends on ebar, S23 L2/L6/L9). Every differential
component is mostly outside the tangent space of the realisable set and is mostly annihilated when
the corrected map is returned to a structure. P's null is what the factorisation plus C2 predicts.

IS THE OBJECTIVE RE-WEIGHTING THE SAME OPERATION? NO. It changes the SCORE of each of the 500
windows, hence which 75 are retained, hence the emitted average; the displacement is a difference
of real windows and is realisable by construction, so no post-hoc profile edit can imitate it and
"the same operation seen twice" is not the right closure. BUT IT IS CLOSED ANYWAY, BY A MEASUREMENT
THAT PREDATES THIS SPRINT. `s12/obj_FINDINGS.md` section 4, all 126 targets, point cloud: shipped
Bayes risk 3.048; 1/sd-weighted L1 3.047; **LFO per-shell weighted L1 3.058** (a separation-band
re-weighting fitted leave-fold-out: +0.010 A, WORSE); shell-profile only, deployable, 3.163. And
the mechanism was built and retired once already in the project's own history: `score_weights.json`
is "intentionally absent" because "refitting drives two of the six parameters to their clip bounds
and loses on the benchmark (top-1 2.76 -> 2.91 A) while gaining 65% on the dev objective"
(`core/predict.py:377-383`) -- the canonical shape of optimising the objective and losing the
structure.
I THEREFORE WITHDRAW THE BUILD RECOMMENDATION IN MY SECTION 1.4 (the withdrawal is marked in the
file at the point of the original sentence). The derivation stands -- a calibration can act on the
Bayes-risk minimiser through the separation-band weights and through nothing else, because a width
error leaves the median map unchanged -- but the conclusion was wrong to call that a live route.
CALIBRATION IS CLOSED, NOT REDIRECTED: its only channel is measured at +0.010 A leave-fold-out, was
retired once before, and is bounded by theorem 2 in any case (a re-weighted objective is still in
class M with g -> gamma g, so its expected informativeness is sum gamma w kappa var(a)(1 - beta), a
re-weighting of second-order terms, capped by section 8 at the 0.03 A a random-magnitude cosine
buys). No lane should spend a build on it.

THE ONE THING IN THAT TABLE WORTH THE SPRINT'S ATTENTION. In the same S12 screen, ORACLE
shell-profile only (n-2 numbers) = 2.299 against the shipped 3.048: a scorer that knows ONLY the
native's per-separation distance profile -- 11 to 14 numbers, no pair detail -- selects a top-75
whose average is 0.75 A better than production, and the complement in the same file is that the
prediction's entire pair-specific detail is worth 0.155 A ("the shipped distogram is, operationally,
a predictor of a 12-to-14-number curve"). That is the lowest-dimensional named instance of what
would break section 8's assumption (B2). Its price is already set by the record: the profile's
leading component is COMPACTNESS, the in-band axis is compactness at oracle +0.909 with the
native's z-scored Rg (S29-L19), the achievable native-free proxies reach 0.24 to 0.37 against that
0.909, and the DEPLOYABLE version of the same scorer is 3.163, i.e. worse than production -- so the
gap 3.163 -> 2.299 is exactly the value of the native information in those numbers.
CHEAPEST FALSIFIER (lane M or P, minutes, no new code path): fit leave-fold-out a predictor of the
per-separation profile ratio from native-free features only (pool profile, posterior profile,
length, predicted secondary-structure content), re-score the 500 windows with the corrected
profile, take the top-75, average, project. REGISTERED PREDICTION from section 8: the
leave-fold-out arm lands INSIDE 0.05 A of production, and it is a result only if it clears 0.7x MDE
with the fold CI excluding zero AND its cosine clears 0.140 with the shrink signature attached. If
it lands near the ORACLE 0.75 A, section 8's (B2) is falsified and that is the sprint's result.
Multiplicity: 0 endpoint comparisons; 1 registered prediction; 1 WITHDRAWAL of my own earlier
recommendation.
Artefacts: `s29/THEORY.md` section 9 (commit 591fab47); sources S29-L22 (lane P), S25 L1,
S23 L2/L6/L9, `s12/obj_FINDINGS.md` sections 3 and 4, `core/predict.py:377-395`, S29-L19.

## S29-L25 -- THE SET-EQUALITY THEOREM FAILS ON REAL POOLS, MEASURED: THE f-OPTIMAL 2-SUBSET IS NOT THE ENERGY-ORDER PREFIX ON 11/12 TARGETS AND NOT THE TWO BEST BY THE PER-STATE CRITERION ON 11/12 EITHER (EXHAUSTIVE OVER ALL 124,750 PAIRS OF THE DEPLOYED 500-POOL), THE f-OPTIMAL m = 5 SUBSET IS NON-PREFIX ON 12/12 WITH A MEAN OBJECTIVE GAP OF +0.152, AND THE ESCAPE BUYS NOTHING: THOSE SUBSETS' AVERAGES ARE WORSE THAN THE PREFIX ON 10/12 AND WORSE THAN PRODUCTION ON 8/12 (NOT MEASURED AT n = 12, DIRECTION ONLY) -- LANE T's S29-L17 SECTION 4.5 CONFIRMED ON THIS INSTRUMENT, BOTH CLAUSES (2026-09-20 00:51, B)
Question (`s29/PREREG_S29_B.md` ADDENDUM 2 section B2.6 falsifier F5c; the coordinator's
instruction of 00:35; lane T's **S29-L17 section 4** and **S29-L18**, whose three-state, two-
dimensional witness shows the set-equality theorem MUST fail under a CVaR over the tail's AVERAGE
structure): does it fail on a REAL pool of this instrument, and does escaping it buy accuracy?
Falsifier as registered (T's, adopted verbatim): the objective gap of the aggregate-optimal subset
over the energy prefix is under 0.02 (it is not), and -- the clause that decides whether this is an
opening -- the aggregate-optimal subset's coordinate average is BETTER than production, in which
case it goes to 126 immediately.

**THE OBJECT.** f = the shipped distogram risk (`s27/s28_A_amp.py :: Surrogate`, the piecewise-
linear reading of the shipped risk table) evaluated on the coordinate AVERAGE of a subset, in the
deployed readout's own frame (`Frame(W, DIS top-75)`, every window superposed on the top-75
medoid). V(S) = f(mean_{i in S} W_i) is neither additive nor monotone, so no per-state sort can be
optimal for it; the deployed CVaR-VQE's tail, by contrast, is provably a PREFIX of the E order
(S24/S25, and S29-L13: it is Barkoutsos eq (12), a definition). Selection is native-free
throughout; the RMSD columns are ORACLE and are attached afterwards.

**PROVENANCE.** `s29/s29_B_compat.py` + `s29/s29_B_tta.py :: subset_target` (18 tests,
`tests/test_s29_B.py`, all pass); job `s29B_tta_subset` (`s26/jobs_done/s29B_tta_subset.json`
exit 0, wall 92 s, peak RSS 0.34 GB); rows `s29/results/s29_B_tta_subset_rows.jsonl`, 12 rows =
S27 T11's 12 trainability targets. The pair search is EXHAUSTIVE over all C(500, 2) = 124,750
pairs of the deployed 500-candidate pool; the m = 5 search is greedy plus swap local search to
convergence. Ties are averaged over the tied argmin set, never read off array order (memory
`tie-breaking-leaks-the-pool-order`).

**RESULT 1, THE MECHANISM -- deterministic, per target, no statistics needed
(`s29_B_tta_subset_rows.jsonl`):**

    pdb    best pair ranks   f_pair   f_prefix    gap     greedy m=5 ranks        f_greedy  f_pref5    gap
    1A13   [0, 10]           0.7402   0.7503    +0.0101  [0, 1, 161, 162, 163]    0.7257   0.7549   +0.0292
    1I6Y   [1, 33]           1.4692   1.5311    +0.0619  [0, 1, 2, 29, 33]        1.4008   1.5737   +0.1729
    1M02   [24, 37]          1.6091   1.9231    +0.3140  [37, 52, 53, 73, 130]    1.5927   1.8614   +0.2688
    2BFI   [6, 38]           0.6047   0.6158    +0.0111  [2, 4, 38, 61, 118]      0.5975   0.6155   +0.0180
    2LWS   [0, 5]            1.6476   1.8429    +0.1952  [0, 5, 6, 8, 52]         1.6574   1.9611   +0.3037
    2MP9   [0, 12]           2.1534   2.2874    +0.1340  [44, 229, 234, 324, 438] 2.0912  2.3700   +0.2788
    2P5H   [0, 4]            1.3083   1.3516    +0.0432  [0, 1, 4, 5, 15]         1.3081   1.3247   +0.0166
    5Z5W   [4, 118]          0.7098   0.7114    +0.0015  [1, 5, 26, 118, 217]     0.7006   0.7116   +0.0110
    6MBM   [15, 82]          1.0228   1.0374    +0.0147  [0, 15, 82, 86, 97]      1.0141   1.0352   +0.0211
    7JGX   [0, 1]            1.8154   1.8154    +0.0000  [0, 1, 2, 7, 11]         1.8406   1.8989   +0.0582
    8HVS   [1, 2]            1.3878   1.4390    +0.0512  [10, 41, 92, 185, 269]   1.3202   1.5338   +0.2135
    9KAR   [0, 3]            2.1736   2.3526    +0.1790  [52, 57, 84, 127, 206]   2.2120   2.6461   +0.4341

- The f-optimal PAIR is **not** the energy-order prefix {rank 0, rank 1} on **11 of 12** targets
  (7JGX is the one where it happens to coincide, gap exactly 0.0000). Mean objective gap +0.0847.
- The f-optimal m = 5 subset is **not** a prefix on **12 of 12**; mean gap +0.1522, which clears
  lane T's registered 0.10 at m = 5 rather than at m = 75. Supports reach rank 324 (2MP9) and 269
  (8HVS) of 500.
- **It is not a different sort either.** The pair of the two lowest-f_single members (f as a
  PER-STATE score, the other obvious classical rule) equals the exhaustive optimum on only
  **1 of 12** targets. So the optimum is not reachable by sorting E, and not by sorting f: it is a
  genuine subset-selection problem, which is exactly what lane T's three-state witness says.
- This is the first measurement in this project's record of a formulation in which "which set" is a
  real optimisation variable rather than a read-out of a sort. On the charter's section 11
  questions 5 and 7 it is an affirmative IN PRINCIPLE, measured on the deployed pool.

**RESULT 2, THE PRICE -- the escape buys no accuracy (ORACLE, point cloud, n = 12, paired through
`s24.stats_lib.compare` with `pinned_folds`; production = DIS top-75 uniform, mean 3.2529):**

  greedy m=5 (f-optimal subset) - production
    a 3.6808 (med 3.2640)   b 3.2529 (med 2.6820)   n=12
    effect +0.4278   median +0.0813   SE 0.2567   MDE 0.7193   effect/MDE +0.59
    iid  CI95 [+0.0154, +0.9677]
    fold CI95 [-0.0785, +0.8049]   folds same sign 3/5   per-fold 0:+1.092 1:-0.078 2:+0.065 3:-0.150 4:+0.527
    4W/8L/0T   worst degradation +2.7500 (2MP9)   p90 +1.3450   power 0.38  Type-M 1.60
    VERDICT: NOT MEASURED (|effect| 0.4278 <= its own MDE 0.7193, 0.59x)

  greedy m=5 - the m=5 ENERGY PREFIX
    a 3.6808 (med 3.2640)   b 3.4293 (med 2.8810)   n=12
    effect +0.2515   median +0.0688   SE 0.1739   MDE 0.4873   effect/MDE +0.52
    iid  CI95 [-0.0007, +0.6295]
    fold CI95 [+0.0386, +0.6086]   folds same sign 5/5   per-fold 0:+0.841 1:+0.023 2:+0.085 3:+0.025 4:+0.063
    2W/10L/0T   worst degradation +2.0548 (2MP9)   p90 +0.5000   power 0.30  Type-M 1.80
    VERDICT: NOT MEASURED (|effect| 0.2515 <= its own MDE 0.4873, 0.52x)

  f-optimal PAIR - the energy prefix pair
    a 3.7370 (med 3.4339)   b 3.4236 (med 3.0863)   n=12
    effect +0.3134   median +0.0980   SE 0.1473   MDE 0.4126   effect/MDE +0.76
    iid  CI95 [+0.0527, +0.6012]
    fold CI95 [+0.1456, +0.3883]   folds same sign 4/5   per-fold 0:+0.421 1:-0.208 2:+0.263 3:+0.365 4:+0.363
    3W/8L/1T   worst degradation +1.3620 (2LWS)   p90 +1.1016   power 0.57  Type-M 1.33
    VERDICT: NOT MEASURED (|effect| 0.3134 <= its own MDE 0.4126, 0.76x)

**How to read those three blocks, stated plainly.** At n = 12 NONE of them is a result: every
effect is inside its own MDE and the entry does not claim otherwise (contract rule 6; S29-L4(d)).
What IS readable is the DIRECTION and the W/L, which are descriptive: the f-optimal m = 5 subset
loses to the m = 5 energy prefix on 10 of 12 targets with 5/5 folds the same sign, and loses to
production on 8 of 12; the f-optimal pair loses to the prefix pair on 8 of 12. **T's F5c second
clause therefore does NOT fail**: the aggregate-optimal subset is not better than production, so
this is not the sprint's opening, and nothing here goes to 126. The worst cell is 2MP9, where the
objective moves the support from the prefix to ranks 44 to 324 and the ORACLE RMSD goes from 1.188
(production) to 3.938: the S29-L2 anti-correlation (ladder rho -0.182 CA / -0.402 chain) doing
exactly what it says it does, on a target where production is already good.

**THE TWO SENTENCES THIS ENTRY IS FOR.** (1) Escaping the set-equality theorem is achievable on
this instrument and is achieved: the optimum is non-prefix on 11/12 and 12/12 and is unreachable
by any sort of E or of f. (2) Escaping it creates no INFORMATION -- f is a marginal-class objective
(lane T's theorem 2, S29-L7; contract addendum 1 rule 21) and the sets it prefers are worse. The
combination is lane T's "necessary and not sufficient" measured rather than cited, and it prices
the endpoint arm of measurement 5 before it is run: the VQE cannot find a better set than an
exhaustive/greedy classical search over the same objective, and that search is worse than
production.

**Charter section 11, the questions this answers.** (7) Does the quantum output contain information
unavailable classically? For THIS objective, no: the classical greedy-plus-local-search is the
correct counterpart (rule 15, S29-L17) and it is a strict upper bound on what the circuit could
find. (8) Can a classical control reproduce the effect? Yes, and better. The affirmative on (5) and
(7) "in principle" stands as a statement about the FORMULATION, not about accuracy, and the entry
says so in the same sentence.

Multiplicity: 4 ORACLE diagnostic comparisons on 12 targets, 0 endpoint comparisons. No native
chose any parameter; f, the subsets and the frame are native-free.
Artefacts: `s29/results/s29_B_tta_subset_rows.jsonl`; `s26/jobs_done/s29B_tta_subset.json`;
`s26/logs/s29B_tta_subset.log`; code `s29/s29_B_tta.py`, tests `tests/test_s29_B.py` (18 pass);
prereg `s29/PREREG_S29_B.md` addendum 2 (commit d3ccf6b6, before the job).
