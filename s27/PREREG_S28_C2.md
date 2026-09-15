# S28 LANE C2 PRE-REGISTRATION -- THE RECOGNITION AUDIT: WHICH NATIVE-FREE SCORER, IF ANY, PREFERS THE ORACLE STRUCTURES THE CIRCUIT FAMILY CAN EXPRESS?

Written 2026-09-14 22:35 before any C2 number was read. Brief `s27/briefs/S28C2.md`; contract
`s27/S28_CONTRACT.md` with addendum 1. Record read: S28-L1b (ORACLE ceiling 0.288 A best-of-5,
0.36 to 0.42 single start, random 27-subspace 0.608; S~ prefers the ORACLE structure on 25/126,
S on 26/126), S28-L13 (like-for-like: quote the median), S28-L17, S28-L18b/L20 (recognition
arms, point cloud), S28-L26b/L27b (built chain: every lam arm worse; the objective is the failing
leg; ORACLE emitted 0.252), S28-L23b/L30 (the objective's gradient at production is blind; the
native scores worse than the average under DIS on 99/126, 36.9th percentile of its own pool).
EVERYTHING HERE IS AN ORACLE DIAGNOSTIC: every structure scored except production and the two
native-free controls' DIRECTIONS is chosen against the native, and every sentence with a number
from one of them says ORACLE. Nothing here is deployable; the falsifier is for the CLOSURE claim.

## 1. The new angle, stated up front
S27 measured 18 native-free channels on POOL MEMBERS (`s27/results/pool_rows.jsonl`) and found
none that helps the average (S27 L2/L3); S28-L1b/L26b measured ONE objective (the shipped Bayes
risk and its smooth surrogate) on the ORACLE structure. No native-free scorer in the S27 library
has been evaluated on OFF-POOL structures, the signed combinations the amplitude family can
express, against production. That is the question: if no scorer in the library prefers a 0.3 A
ORACLE structure to the 3.05 A production average, recognition is closed for the whole library
on this family in one measurement, and the report can say what an objective would have to know.

## 2. Inputs, per target (126 targets, pinned folds; `s27/s28_C2_recog_audit.py`)
Frame: lane A's `Frame(cand.W, top)` (every pool window posed on the DIS-top-75 medoid,
S28-L1b: reproduces the deployed average to 5.7e-14). Structures, all (n, 3) CA clouds:
- PROD: the production average (`s27/results/s28_A_structs/<pdb>.npz :: prod`), asserted equal
  to `s24.d_harness.readout_uniform` on the DIS top-75 to 1e-10.
- ORACLE circ_best: `s28_A_structs :: oracle_circ` (best of 5 starts, S28-L1b; an order
  statistic, S28-L13).
- ORACLE circ_s0: the single-start seed-0 ORACLE optimum, regenerated with lane A's
  `oracle_circuit_ceiling(starts=1)`; its RMSD is asserted equal to `s28_A_oracle_rows.jsonl ::
  arms.oracle_circ.per_start[0]` to 1e-6.
- ORACLE sub0: the random-27-subspace least-squares structure (subspace 0), regenerated with
  lane A's `oracle_subspace_ls(subspace_matrix(pdb, 0, k))`; RMSD asserted equal to
  `per_sub[0]` to 1e-6.
- ORACLE native-equivalent: `s28_A_structs :: oracle_aff500` (the affine-500 least squares,
  RMSD ~1e-7 to the native by dimension counting, S28-L1b); it is the native's CA trace to
  floating point and is labelled NATIVE(aff500).
The ladder is PROD (3.05) > sub0 (0.61) > circ_s0 (0.36 to 0.42) > circ_best (0.29) > NATIVE
(0.00), all ORACLE except PROD.
Controls (the two the brief names, plus one scale variant):
- (i) RAND_SIGNED: a random signed AFFINE combination of the 500 posed windows (w = B z with a
  fresh stable-seeded B (500 x 27) and z ~ N(0, 1), normalised to sum w = 1), moved along its
  displacement from PROD, C = (1 - s) PROD + s C_w, with s chosen by bisection so that
  CA-RMSD(C, PROD) equals CA-RMSD(circ_best, PROD) on that target (the ORACLE structure's own
  distance from production; the SCALE is ORACLE, the DIRECTION is native-free). It stays in the
  affine family. 4 draws per target, the fraction averaged over draws. Note: a random member of
  the circuit family cannot be placed at a chosen distance without an optimisation (untrained
  circuit states sit at ~87 A, S28-L26b), so the matched family is the affine span that
  contains the circuit family; stated here, not discovered later.
- (ii) GAUSS_0.3: PROD plus isotropic Gaussian noise scaled so that CA-RMSD(C, PROD) equals the
  ORACLE structure's RMSD to the native on that target (the brief's "0.3 A"); 4 draws.
- (ii-b) GAUSS_MATCHED: the same at the ORACLE structure's distance from PROD (the scale of
  (i)); 4 draws.
- (iii) the FAIL18 / 108 split (ORACLE label) of every fraction.

## 3. Scorers (native-free functions of one CA trace, or of one built chain, and the pool)
CA level, through a one-Context adapter over `s27/ham_lib.py` (the Context's "pool" W is the
stack of the structures to score; the universe fits and the distogram are the target's own):
DIS (`s12.instrument.shipped_score`), DIS_MEAN, CONTACT_LL, DISTPOT, CONTACT, ENV, HP, RG_LAW,
RG_UNIV, EXVOL, CAGEO, SS_MATCH (12), plus three pool-relative adapters: CONS_POOL (mean CA-RMSD
to the 500 pool windows), DMAP_CONS_POOL (L1 to the pool's median distance map), POOLGO_POOL
(the pool's contact frequencies). Test: on a pool member, every adapter reproduces the
`s27/cache` pool-row value (CONS_POOL to the (k-1)/k self-exclusion identity) to 1e-9.
Built chain (each structure projected by `s12.instrument.project`, its phi/psi): RAMA, DSSPHB,
ELEC, TORS_CONS_POOL, LEG (the 11 terms and the total). AMB: only if one AMBER process fits the
box and it is the only AMBER job (contract rule 8); deferred by default and said so.
Lower is better for every scorer (the S27 convention).

## 4. Statistics, per scorer and per structure X in the ladder and the controls
- pref(X) = fraction of the 126 targets with score(X) < score(PROD), with the Wilson 95% CI and
  a fold-clustered bootstrap CI (4,000 draws over the 5 folds).
- The recognition contrast: pref(circ_best) - pref(RAND_SIGNED) as a paired per-target
  difference of indicators (`ST.compare` on the 0/1 vectors, fold CI); the same against
  GAUSS_MATCHED and GAUSS_0.3.
- Ladder ordering: per target the Spearman rho between the scorer's values and the ORACLE RMSD
  over the five ladder points (PROD, sub0, circ_s0, circ_best, NATIVE); mean rho, fold CI, and
  the fraction of targets with rho > 0. A scorer that recognises orders the ladder (rho -> +1).
- The linear combination: a pairwise logistic (no intercept, L2, sign-augmented pairs
  (Delta, 0) and (-Delta, 1) with Delta = standardised score(circ_best) - score(PROD) over the
  CA scorers), alpha chosen by inner leave-one-fold-out, held-out sign accuracy over 126 targets;
  null: 200 draws of random sign flips of Delta per target through the identical nested
  procedure. Reported also for the chain scorers when they land.
- Multiplicity: 15 CA scorers (and 16 chain scorers) are tested; the best single is priced as an
  order statistic with a max-over-scorers sign-flip null (500 draws), as in S28-L6.
- Strata: every pref and contrast on FAIL18 (n = 18) and the 108.

## 5. Falsifier (for the closure claim) and registered prior
"Recognition of the ORACLE structures is closed on the S27 library" is FALSIFIED if any single
scorer, or the nested linear combination, has pref(circ_best) with the fold-clustered CI
excluding 0.5 on the side above 0.5 AND beats control (i): pref(circ_best) - pref(RAND_SIGNED)
with the fold CI above zero. A scorer that passes the first clause and fails the second
"prefers any signed combination", not the ORACLE one. A scorer whose pref(circ_best) has the
fold CI below 0.5 is ANTI-recognition on this family (as CONS was on the pool, S27 L3).
Registered prior: DIS prefers the ORACLE structure on a minority (25 to 26 of 126, S28-L1b; the
native itself is at the 36.9th percentile of its pool, S28-L30); DIS_MEAN and CONTACT_LL (the
same posterior) the same; the geometric channels (EXVOL, CAGEO, RG_LAW, RG_UNIV) prefer PROD
because the contracted average has no clashes and a typical Rg while a signed combination is
off the manifold; the consistency channels prefer PROD (the average IS the pool's mode); the
statistical potentials (DISTPOT, CONTACT, ENV, HP) are the open cell, prior null. Nothing clears
both clauses. If one does, it is the candidate objective for a lane A2 run, not a result.
Second seed and reversed order: the ORACLE regeneration is deterministic per target; a positive
is re-run with a second seed of the controls (4 more draws) and a second seed of the
sign-flip null before it is believed.

## 6. Cost, memory, checkpoints
CA audit: one Context per target (universe pair distances in float32 as S27), the single-start
ORACLE regeneration (300 Adam iterations on the 9-qubit statevector, about 0.5 s), the subspace
LS (exact), 4 + 4 + 4 control draws: a few seconds per target; peak RSS expected under 0.6 GB
(S27's `run_pool` shape); one-target probe under jobrun first, peak quoted. Built-chain audit:
8 projections per target at 3 to 6 s each, about 1 to 1.5 h uncontended, checkpointed per
target in `s27/results/s28_C2_chain_rows.jsonl`, resumable. Results:
`s27/results/s28_C2_ca_rows.jsonl`, `s28_C2_ca_summary.json`, `s28_C2_chain_rows.jsonl`,
`s28_C2_chain_summary.json`. Tests `tests/test_s28_C2.py`: the adapters reproduce pool-row
values on a pool member; the controls land at the stated RMSD (to 1e-6); the ladder order of
ORACLE RMSDs; the pairwise logistic recovers a planted preference and not a random one.

## 7. Rule-0 forks not taken
- Scoring the ORACLE structures with a learned scorer trained on them: that would be an ORACLE
  objective; not run.
- Re-optimising the circuit under a scorer that passes: that is lane A2's job if any passes.
- AMBER: deferred (box at 93 to 96% from the user's load all evening; one AMBER process is not
  guaranteed to fit); recorded in "what I did not do" if it stays deferred.

## ADDENDUM 1 (2026-09-14 22:50, written after the CA-level `analyse_ca` printed and before any
number of the additions below was read)
The CA audit's linear combination prefers the ORACLE structure on a held-out 0.968 but also
"prefers" the random signed control on 0.887 and the matched Gaussian control on 0.895: a rule
that is mostly ANTI-PRODUCTION (it recognises the contracted average, not the ORACLE
structure). To separate the two readings, two statistics are added, for every scorer and for
the linear combination: (a) HEAD-TO-HEAD, pref(circ_best over RAND_SIGNED) = the fraction of
targets on which the scorer scores the ORACLE structure below the random signed combination at
the same distance from production (production does not enter; this is recognition AMONG signed
combinations), with the fold CI, and the same against GAUSS_MATCHED; (b) the paired contrast
pref_lc(circ_best) - pref_lc(control) for the linear combination, as already registered for the
single scorers. The second seed of the controls (draws 4 to 7, `s28_C2_ca_rows_seed2.jsonl`)
is run for every scorer, as section 5 requires for a positive (CAGEO clears both registered
clauses at seed 0: pref 0.611, fold CI [0.551, 0.691]; vs RAND_SIGNED +0.190 [+0.136, +0.252];
its max-over-15 sign-flip null p is 0.072, at the 95th percentile). The falsifier is unchanged;
these are the readings the ledger entry must carry beside it.
