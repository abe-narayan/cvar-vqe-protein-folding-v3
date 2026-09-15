# S28 LANE A PRE-REGISTRATION -- THE AMPLITUDE READOUT (a CVaR-VQE read as a SIGNED combination)

Written 2026-09-14 before any endpoint number was read. Never edited after a result exists;
addenda are appended. Brief: `s27/briefs/S28A.md`. Contract: `s27/S28_CONTRACT.md`.

## 0. The question, in two halves that are never in one sentence without the label

Every operator the project ships is convex (non-negative weights summing to one). The ORACLE
bound ladder (`docs/FINDINGS.md` S10-5, all ORACLE) puts the uniform top-75 average at 2.954 A
emitted, the convex hull of the K = 500 pool at 0.853, and the AFFINE span of the same 500
windows at 0.064. The real amplitudes of the deployed RY/CNOT circuit are signed for free.

- EXPRESSIVITY (ORACLE): how near the native can the family
  C(theta) = sum_i psi_i(theta) W_i / sum_i psi_i(theta)   (i over the 500 real candidates)
  reach, when theta is chosen with the native in hand? This is the family's ceiling. It is a
  diagnostic, labelled ORACLE at every appearance, never a result.
- RECOGNITION (deployable): does the NATIVE-FREE objective
  F_lam(theta) = CVaR_alpha(E; p_theta) - T H(p_theta) + lam * S(C(theta))
  find those structures? lam = 0 is the deployed selector's objective; lam > 0 adds the
  shipped distogram Bayes risk of the EMITTED structure's own distance map.

## 1. What is held fixed

| | value | source |
|---|---|---|
| targets | the 126 dev targets, `sorted(pdb)`, all of them | `s25.phys_lib.targets()` |
| folds | the 5 pinned folds | `s24.stats_lib.pinned_folds` |
| pool | the shipped K = 500 BLOSUM pool, `W` (500, n, 3) | `s24.d_harness.Candidates.from_universe` via `s27.run_pool.channels_for` |
| Hamiltonian | `E = Encoding(zrank(DIS)).E` over 512 basis states; DIS = the shipped score from `channels_for(pdb)[1]["DIS"]` (asserted bit-equal to the S25 cache); 12 padding states at the pad energy | `s22.qcand_lib.Encoding`, `s25.phys_lib.zrank` |
| circuit | `core.quantum.StatevectorCircuit(n=9, layers=3)`, 27 parameters, real amplitudes, exact statevector | frozen production module |
| CVaR, entropy | alpha = 0.18, T = 0.5 (the S25 / S27 suite settings for the 9-qubit pool, `s27/PREREG.md` section 1); `core.quantum.free_energy` gives the exact parameter-shift gradient of CVaR - T H | |
| optimiser | Adam, lr 0.15, 80 iterations, theta0 ~ N(0, 0.6) from `np.random.default_rng(seed)` -- the deployed loop of `core.quantum.run_cvar_vqe`, rewritten only to add the lam term; at lam = 0 it must reproduce `run_cvar_vqe`'s `p` to 1e-12 (a soundness gate, asserted in the tests and reported per target) | |
| production arm | DIS top-75 (ties by the S27 stable key, `s27.run_pool.topm`) -> `s24.d_harness.readout_uniform` -> `s12.instrument.project`; anchors 3.048338 point cloud, 3.2126 built chain (`s27/results/chain_rows.jsonl :: DIS`) must reproduce | |
| basis | BUILT CHAIN (`s12.instrument.project`, ramah 0.3, multi-start) is the reporting basis; the point cloud is an intermediate and is reported beside it, never across from it | contract rule 5 |
| statistics | `s24.stats_lib.compare` (paired, SE, MDE = 2.8016 SE, iid CI beside fold CI, W/L/T, concentration null); decide on the fold CI; every contrast printed with `ST.fmt` verbatim | contract rule 6 |
| seeds | seed 0 primary; any positive is re-run at seed 1 and with the target order reversed | |

Nothing native enters any deployable arm. `nat_ca` is read only inside functions named
`oracle_*`, and the NaN-poison test (`tests/test_s28_A.py`) asserts every deployable output is
bit-identical when `nat_ca` and `oracle_rr` are replaced by NaN.

## 2. The readout, exactly

1. FRAME. The deployed average superposes the retained set on its own medoid
   (`s12.instrument.coordinate_average`). Here the retained set is the DIS top-75 and its medoid
   `b` is `I.medoid(I.pairwise_rmsd(W[top75]))`; ALL 500 windows are posed on `W[top75][b]` with
   `I.superpose_batch`. This is native-free (DIS chooses the set, the set chooses its medoid)
   and it makes the uniform-on-top-75 readout equal to `readout_uniform(cand, top75)` (asserted
   in the tests to 1e-9, the max deviation reported).
2. WEIGHTS. psi = `circ.state(theta)` (512 real amplitudes). The 12 padding states are EXCLUDED
   from the readout; their probability mass is reported. Affine weights over the 500 real
   candidates: w_i = psi_i / sum_j psi_j, so sum_i w_i = 1 for either sign of the denominator.
   C(theta) = sum_i w_i W_i. If |sum_j psi_j| < 1e-9 the readout is undefined and the arm
   records a failure for that target (reported, never silently repaired).
3. OBJECTIVE TERM. S(C) is the shipped Bayes-risk of C's own |i-j| >= 2 distance map,
   `s12.instrument.shipped_score(dg, pair_dists(C))`. The shipped lookup is piecewise constant on
   a 0.05 A grid and has zero gradient almost everywhere, so the OPTIMISED term is S~(C): the
   same risk table (`dg["risk"]`, `dg["grid"]`) read by LINEAR interpolation between grid points,
   clipped to the grid ends. S~ equals S at every grid point and differs from it by at most one
   bin's variation elsewhere (max |risk step| 0.08 on 1A13, `s12/cache/disto_1A13.npz`). BOTH
   S (shipped) and S~ (surrogate) are reported on every emitted structure.
4. GRADIENT. For RY generators the state itself obeys dpsi/dtheta_k = psi(theta + pi e_k) / 2
   exactly (one shifted simulation per parameter, batched). Chain rule through
   dC/dpsi_i = (W_i - C) / sum_j psi_j and dS~/dC (linear-interpolation slopes through the
   pair distances). The CVaR - T H part uses `core.quantum.free_energy`'s exact parameter-shift
   gradient. The sum is the exact gradient of F_lam; tests check it against central finite
   differences (cosine > 0.999, max abs error < 1e-6 at h = 1e-5, on synthetic pools).

## 3. Arms (per target; seed 0 unless stated)

(1) PRODUCTION: DIS top-75 uniform average -> built chain. The comparator for everything.

(2) CIRCUIT, lam = 0: the deployed CVaR-VQE state, read as signed weights. Does the plain
    trained state, read this way, already differ from the average?

(3) CIRCUIT, lam in {0.3, 1, 3}, 80 iterations (the deployed loop). A grid of three: the lam is
    chosen leave-fold-out (nested; the held-out effect is what is reported), and the per-target
    minimum over the grid is priced with `ST.best_of_k_within` (split-half transfer, k_eff).
    Secondary: the same three lam at 400 iterations, as a convergence check (the F trace at
    0/20/40/60/80/.../400 is stored) and as three further grid cells, priced the same way.

(4) CLASSICAL CONTROLS at matched evaluation budget (80 objective+gradient evaluations,
    L-BFGS-B `maxfun=80`) and at convergence (`maxfun=2000`), both reported. Every control
    replaces ONLY the map theta -> psi; p = psi^2 / |psi|^2, w = psi / sum psi, C and F_lam are
    the identical code path (the "matched-objective" form, PRIMARY for falsifier F2):
    (4a-500) psi = z, z in R^500 (no padding);
    (4a-75)  psi = z on the DIS top-75, 0 elsewhere (75 numbers);
    (4b)     psi = B z, B a (500, 27) Gaussian matrix from a stable per-target RNG, z in R^27;
             8 subspaces, the MEAN over subspaces is the control (best-of-8 reported apart,
             against `best_of_k_null`);
    (4c)     psi = sqrt(softmax(z)), z in R^500: non-negative amplitudes, so w is on the SIMPLEX
             and p = softmax(z).
    Inits: "rand" (z ~ N(0, 1), the analogue of the circuit's random theta) is PRIMARY; "prod"
    (z = the production weights: uniform on the top-75, zero or -20 in the logit elsewhere) is
    SECONDARY and answers a different question (does the objective ever improve on production
    from production?).
    The brief's literal form, S~-only over affine weights w (no CVaR, no entropy), is ALSO run
    for 4a-500, 4a-75, 4b, 4c ("S-only"), SECONDARY: it conflates the family with the objective.

(5) ZERO-INFORMATION CONTROL: the UNTRAINED circuit, theta ~ N(0, 0.6) (the same law), read as
    signed weights; 16 draws; the MEAN over draws is the control; best-of-16 is an ORACLE
    selection and is reported against `ST.best_of_k_null` of the draws.

(6) ORACLE EXPRESSIVITY (diagnostic, never a result): minimise the point-cloud CA-RMSD of
    C(theta) to `nat_ca` over theta, Adam lr 0.05, 300 iterations, 5 starts (theta0 ~ N(0, 0.6),
    seeds 0..4), best of 5 (a best-of-5 order statistic, reported with the mean over starts);
    the best is projected -> the family's ORACLE emitted ceiling. Beside it, in the SAME frame,
    the ORACLE least-squares ceilings of the control families: affine-500, affine-75 and the
    random 27-dim affine subspace (mean over the 8 subspaces), all closed-form, all ORACLE.

## 4. Built chain: which arms are projected (3.4 s per projection, S27's rate)

Job 1 (point cloud + every optimisation, all arms, per-target checkpoints). Job 2 (built chain,
UNCONDITIONAL, the primary list): production; circuit lam 0, 0.3, 1, 3 at 80 iterations;
4a-500 rand matched-budget at lam 0.3, 1, 3; 4a-75 rand matched-budget at lam 0.3, 1, 3; one
untrained draw (draw 0); the ORACLE circuit optimum; the ORACLE affine-500 least squares.
That is 14 projections per target. Job 3 (built chain, CONDITIONAL, decided by a rule fixed
here): if ANY circuit lam arm is at or below production on the built chain (mean over 126),
project 4b (subspaces 0 and 1) and 4c at the three lam, and the 400-iteration circuit arms;
otherwise those arms are reported on the point cloud only and the reason is this sentence.

## 5. Falsifiers and registered priors

F1 "THE FORMULATION MOVES ACCURACY". Fires if some circuit lam arm (lam in {0.3, 1, 3}, 80
   iterations, seed 0) beats PRODUCTION on the BUILT CHAIN: effect < 0, |effect| > its own MDE,
   fold-clustered CI excluding zero, 5/5 folds the same sign; then replicated at seed 1 and with
   the target order reversed, inside the seed-0 CI; and, because lam is a grid, the nested
   leave-fold-out held-out effect must also clear its MDE and `best_of_k_within` must NOT read
   "NOT A SIGNAL". Registered prior: DOES NOT FIRE. Arm (2), lam = 0, is expected WORSE than
   production by a large margin (the trained state carries 8.8 of 9 bits, S27 T9: near-uniform
   |psi| with unstructured signs, so the affine combination is a near-random signed sum of all
   500 windows with an O(1) denominator). Arms (3) are expected worse than production too: S15
   measured direct optimisation of a structure against the distogram objective at 3.321 A and
   S8-9 put the native at the objective's 37th percentile; the open question is whether the
   circuit's 27-dimensional family regularises that failure, and the prior is that it does not.

F2 "THE CIRCUIT MATTERS". Fires only if F1 fires AND that arm beats (4b) (the random 27-dim
   subspace, mean of 8, matched objective, matched budget) beyond its MDE with the fold CI
   excluding zero, AND beats (4a-500 rand, matched budget) beyond its MDE. Registered prior:
   does not fire (F1 is expected not to fire). If F1 fires and F2 does not, the result is a
   property of the objective, not of the circuit, and is reported as such.

F3 EXPRESSIVITY (ORACLE, a measurement with a prior, not a result). Prior: the family's ORACLE
   point-cloud ceiling lies between the affine-500 span (0.035 raw in S10-5's joint-transform
   convention; re-measured here in the fixed frame) and the uniform top-75 average (2.815 raw),
   and BELOW the ORACLE best single member of the top-75 (2.306): a signed 27-parameter family
   should beat selection. If the ceiling is ABOVE 2.306 the family is less expressive than
   picking one window, and the whole formulation is closed by its ceiling. If the ceiling is
   low but F1 does not fire, the finding is "the structures exist in the family and the
   objective does not recognise them", and section 6 says what the objective would need.

F4 (4a) PRIOR. Unconstrained signed weights minimising the same objective are expected WORSE
   than production (S15's 3.321 and S8-9's 37th percentile; and an unconstrained affine
   combination can leave the manifold of peptide-like structures entirely: its Rg and
   contraction will say so).

F5 (5) PRIOR. The untrained circuit read as signed weights is expected far worse than
   production and far worse than a random 75-subset (3.4209): it is a random signed sum.

## 6. Diagnostics (per arm, per target; all in the rows)

Sign structure of the optimum (fraction of negative w_i; negative mass sum |w_i| over w_i < 0),
effective number of members 1 / sum w_i^2, the denominator sum psi_i, padding mass, the emitted
structure's Rg and mean virtual CA-CA bond against the pool's and production's (contraction),
S (shipped) and S~ (surrogate) of the emitted structure and of production's average, the F
trace, seconds. Aggregate: per-target scatter against production, FAIL18 vs the 108 others
(`I.FAIL18`), the ORACLE ceilings side by side with the recognition arms (labelled), the
correlation across targets between the ORACLE ceiling and the recognition arm's RMSD, and
what the objective would have to know: the RMSD of the arm's optimum against the RMSD of the
ORACLE optimum, and S~ at the ORACLE optimum against S~ at the arm's optimum (does the
objective even prefer the ORACLE structure? if S~(oracle) > S~(arm) on most targets the
objective is wrong, not the optimiser).

## 7. Statistics and what counts

`ST.compare(arm_chain, production_chain, folds, names=pdbs, label=...)` for every arm; the
built chain is the basis on both sides; the point cloud is a second table with the same
contrasts. A result must clear its own MDE and have the fold CI exclude zero; 0.7 to 1.3x MDE
is the Type-M zone and is not a result; below 0.7x is underpowered, never "trend". Any
minimum over the lam grid or over subspaces or over draws is an order statistic and is priced
(`best_of_k_within`, `best_of_k_null`); only the nested leave-fold-out and split-half numbers
are quoted as achievable. Ties are broken by the S27 stable key, never by array order.

## 8. Resources

Point-cloud job: about 20 to 30 s per target (all optimisations are on a 512-amplitude
statevector or in <= 500 dimensions), about 1 h; built-chain job: 14 x 3.4 s = 48 s per
target, about 1.7 h. Memory estimate 0.6 GB (S27's chain job peaked at 0.321 GB,
`s26/jobs_done/s27_chain.json`); a one-target probe under `s26/jobrun.py --agent S28A` is run
first and its peak RSS quoted. Checkpoints: one JSONL row per (arm, target) appended
atomically after each target; resumable. Tag CPU. No AMBER, no ESM (`s12/cache/disto_*.npz`
exists for all 126 targets, so `I.distogram` never touches ESM).

## 9. Rule-0 forks (the alternative not taken, and why)

(a) Frame: a joint transform solved with the weights (S10-5's `cf` convention) would give a
    lower ORACLE ceiling; the fixed DIS-medoid frame is the deployed average's frame and the
    only one in which the uniform-on-set check holds exactly. Taken: fixed frame.
(b) Padding: treating the 12 padding amplitudes as zero-windows would bias every readout;
    excluding them and reporting their mass is the honest choice. Taken: exclude.
(c) The surrogate S~: the shipped lookup has zero gradient a.e.; the interpolated table is the
    same table. Both values are reported on every structure. Taken: S~ optimised, S reported.
(d) alpha, T: the production 7-qubit LFO table (alpha 1.0 on 3/5 folds, T 0.3) is a different
    register (top-128); the 9-qubit-over-500 setting is what S27 T9 and the brief name. Taken:
    0.18 / 0.5.
(e) The classical controls' objective: the brief's S-only form conflates family and objective;
    the matched-objective form isolates the family. Both are run; matched-objective is primary
    for F2, S-only is reported as the brief's literal control.
(f) A gradient-free optimiser for the circuit (Nelder-Mead) would dodge the surrogate; the
    exact gradient is available and is what the spine uses. Taken: exact gradient, Adam.
(g) An amplitude readout normalised by sum |psi_i| (signed but not affine) would break
    translation invariance of the combination. Taken: affine.

## 10. What I will not do

No AMBER. No touching `core/`. No parameter chosen on a native quantity. No number reported
without its artefact path (`s27/results/s28_A_rows.jsonl`, `s27/results/s28_A_summary.json`).

## ADDENDUM 1 (2026-09-14, before any code was run or any endpoint number read)

Section 2.3 said the surrogate S~ is "clipped to the grid ends". Corrected before implementation:
BEYOND the grid ends (d < 2.0 A or d > 39.95 A) S~ continues LINEARLY with the slope of the
last grid interval. Reason: the true per-pair Bayes risk E_post|t - d| has slope exactly
-w_p / +w_p once d lies beyond all posterior mass, so the linear extension is the faithful
risk, and a clipped surrogate would give the optimiser zero gradient whenever an affine
combination blows up (small denominator), leaving it no way back. The SHIPPED S is still
reported exactly as shipped (clipped lookup) on every emitted structure. Nothing else changes.

## ADDENDUM 2 (2026-09-14, coordinator steer received before any code was run)

1. Order of work: the ORACLE expressivity ceiling (section 3, arm 6; built chain of the best
   start) is run on all 126 targets FIRST, as its own governed job, and posted as the first
   ledger entry after the pre-registration, labelled ORACLE. Recognition follows.
2. No cosmetic variants: the 400-iteration circuit cells are NOT grid cells. They are run at
   lam = 1 only, as a convergence DIAGNOSTIC (F trace), reported in the findings, never as an
   arm. The lam grid stays {0.3, 1, 3} at 80 iterations. The "prod"-init classical controls
   stay (cheap, point cloud only unless section 4's rule fires).
3. Every arm that clears 0.7x MDE against production on the point cloud is projected to the
   built chain before it is called anything (this is stricter than section 4's list and
   supersedes it where they differ); the ledger verdict is written on the built chain.
4. Any positive is posted promptly with, in the same entry: (a) whether it persists through
   `s12.instrument.project` and whether it holds on the 108 non-FAIL18 targets; (b) what is
   new relative to S26/S27 with the ledger line it goes beyond; (c) the three-way split
   (objective optimum vs what the circuit reached vs the emitted structure); and nothing is
   built on it until lane D posts STANDS.

## ADDENDUM 3 (2026-09-14 20:05, after the ORACLE phase, BEFORE the recognition run; answers lane D's S28-L1 caveats)

(b) The target-order reversal is not a replication: every target is seeded on its own, so the
    reversed run can only detect state leaking between targets and is reported as a
    determinism check. The replication of any positive is the second seed (seed 1) alone.
(e) Undefined readouts (|sum psi| < 1e-9): a target whose arm readout is undefined is scored,
    in that arm's paired comparison, at the target's ZERO-INFORMATION value (the mean point-cloud
    RMSD of the 16 untrained draws for that target, and on the chain the projection of draw 0),
    never dropped and never repaired; the count per arm is reported in the summary. If the count
    is zero the rule is moot and says so.
(a) The recognition rows record, for every circuit arm, the three parts of F (CVaR, T H, lam S~)
    at theta0 and at the optimum, so the grid's degeneracy (k_eff) can be judged from the
    parts, not assumed; if lam S~ dominates at every lam the grid is one cell and no nested
    choice is quoted.
(c) The inner (nested) choice of lam is on the BUILT CHAIN, per held-out fold, the chosen lam
    printed per fold; if it differs across folds, k_eff and the split-half transfer are the
    quoted numbers.
(d) Arm (2), lam = 0: the signs are fixed by the initialisation trajectory (CVaR - T H is blind
    to the sign of psi); any difference from the average is a statement about random signs.

## ADDENDUM 4 -- A2, THE OBJECTIVE'S LOCAL BEHAVIOUR AT THE PRODUCTION POINT (2026-09-14 21:30, written after S28-L18b (point cloud) and BEFORE any A2 number; brief `s27/briefs/S28A2.md`)

New angle (contract rule 10): S15 measured the objective's global optimum and S8-9 its ranking
of the native among pool members; S28-L18b measured its optimum over signed families. Nobody
has measured its LOCAL behaviour at the pipeline's own output. Code `s27/s28_A2_local.py`,
tests `tests/test_s28_A2.py`, results `s27/results/s28_A2_*.json`.

A2.1 ORACLE DIAGNOSTIC (labelled in every sentence). At the production point cloud C0 (the DIS
top-75 uniform average in its medoid frame, `s28_A_amp.Frame` + uniform readout, anchored
3.048338): g = dS~/dC, the analytic gradient of the shipped Bayes risk read by linear
interpolation (`s28_A_amp.Surrogate.value_grad`, tested against central finite differences),
and the ORACLE direction u = (native Kabsch-aligned onto C0) - C0. Both -g and u have their
rigid-body components removed (least-squares projection onto the 6-dimensional space of
translations and infinitesimal rotations about C0's centroid; tested: a pure rigid field
projects to zero). Reported per target: cos(-g_perp, u_perp); the distribution over 126
(mean, median, fraction positive); the FAIL18 / 108 split with `ST.compare` of the cosine
against zero is NOT a paired contrast, so the split is reported as means with SE and a sign
test. The same cosine for the S27 CA channels: RG_LAW and EXVOL with analytic gradients;
DISTPOT, CONTACT, ENV, CAGEO are step functions of the coordinates (histogram lookups) with
zero gradient almost everywhere, so for them a SMOOTHED central-difference gradient at h = 0.5 A
per coordinate is used and labelled "smoothed FD (step function)"; the distogram's cosine is
placed among them. Prior: cosine near zero or negative on most targets (the distogram's errors
are common-mode with the pool's, S23 L9); if positive on the 108 and negative on FAIL18, that is
the sequence-conditioning sign flip seen from the objective's side, and is said so.

A2.2 DEPLOYABLE STEP LADDER (native-free; never sees u). C(e) = C0 - e * g / rms(g) for
e in {0.1, 0.3, 1.0} A of RMS displacement (rigid-body components removed from g first, so the
step is pure shape), then the built chain `s12.instrument.project`, paired against production
(re-projected in the same job, S28-L18/L19). Controls, matched in the operator's space:
(i) a RANDOM DIRECTION of the same RMS displacement (Gaussian field, rigid-body removed,
normalised; 8 draws from a stable per-target RNG; the MEAN over draws is the control and the
best-of-8 is an order statistic reported against `best_of_k_null`); on the point cloud all 8
draws, on the built chain draws 0 and 1 per e (16 projections per target is the budget's
limit; stated); (ii) the SAME step restricted to the CIRCUIT FAMILY: theta_P = the native-free
argmin over theta of rms(C(theta) - C0) (Adam, 300 iterations, 5 starts, the family's nearest
point to production; its residual rms is reported, since the family need not contain C0), then
one steepest-descent step on S~(C(theta)) in theta, scaled by a scalar line search so the
C-space RMS displacement from C(theta_P) is e; paired against C(theta_P) itself (its own
projected baseline) AND against production. Falsifier: some e beats production on the built
chain beyond its MDE with the fold CI excluding zero on 5/5 folds AND beats the random-direction
mean beyond its MDE. Prior: the step degrades at every e; the random direction degrades about
equally; the circuit's one-step arm degrades from its own baseline. Any positive is replicated
at seed 1 (the random control's draws) and answered by lane D before anything is built on it.
Diagnostics: |g| (RMS), the S~ decrease per e, Rg and virtual bond of C(e), FAIL18 / 108 split.
NaN-poison: the deployable ladder is run with `nat_ca` NaN and its emitted clouds must be
bit-identical (test).

Order of work (coordinator): addendum (this), code and tests, A2.1 on 126 targets as its own
governed job and posted as an ORACLE DIAGNOSTIC entry, A2.2 on the point cloud, then its
built chain only AFTER the primary chain job has landed and its verdict is posted.
