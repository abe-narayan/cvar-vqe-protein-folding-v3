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
