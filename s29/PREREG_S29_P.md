# PREREG S29 LANE P -- THE PROJECTION PRICE: IS A GEOMETRICALLY INCONSISTENT CLOUD A BIASED FIT?

Written 2026-09-20 00:10 Pacific, BEFORE any number of this lane was computed. Brief
`s29/briefs/S29P.md`; contract `s29/S29_CONTRACT.md` (rules 7, 10, 12, 13, 16, 17, 18, 19 and
addendum 1 rules 20, 21); charter `s29/BRIEF.md` sections 2, 9, 10, 14.

---

## 1. THE QUESTION

Production emits a point cloud C (the DIS top-75 uniform coordinate average in the retained
set's medoid frame) at mean CA-RMSD **3.0483 A**, and then projects it onto the nearest
ideal-geometry chain, emitting **3.2126 A**. The projection **COSTS +0.159 A**, the largest
single-stage loss in the record, and it is the one stage no S28 lane touched.

The cloud is **geometrically inconsistent with its own output space**: its mean adjacent
CA-CA distance is 2.961 A against the ideal 3.812 (22% short, `s23/LEDGER.md` L1 lines 5-8),
while its radius of gyration is only 6% short. Averaging smooths: it shortens local bonds and
largely preserves global extent. The projection

    minimise over (phi, psi)   RMSD( build(phi, psi), C )  +  lam * ramah(phi, psi),  lam = 0.3
    (`s12/instrument.py :: project` -> `core/project.py :: lam_path`, multi-start, maxiter 300)

fits a chain whose bonds are FIXED at 3.80 A to a target whose bonds are 2.96 A. A least-squares
fit of a rigid-length object to a systematically contracted target is a **biased** fit: the
residual cannot be made isotropic, and the optimiser trades global placement against local
over-length. **Does removing the geometric inconsistency BEFORE the projection change the built
chain?**

## 2. WHY THIS IS NOT S23 L6 (contract rule 10)

`s23/LEDGER.md` **L6** (lines 143-186) closed the per-target RMSD-optimal scale s*:

- L6(b), lines 147-151: the population-optimal GLOBAL scalar captures 1.7% of the ceiling,
  in-sample, and is CLOSED.
- L6(d), lines 165-176: *"s\* is not a property of the target's fold. It is fitting the
  particular (pool, reference-structure) pair"* ... *"That makes the scale correction
  unreachable IN PRINCIPLE from native-free information, not merely unreachable at n=126 ...
  it provably carries no native-free content, because the thing it depends on is the answer."*
- The closed form quoted at `s23/LEDGER.md` line 40: Kabsch rotation is scale-invariant, so
  `s* = Ct/Cc` -- a ratio whose numerator Ct is **a contraction with the native**.

Four statements distinguish this lane's quantity from L6's:

1. **L6's s\* is read off the native; this lane's g is read off the cloud's own bonds.**
   g = 3.80 / (mean adjacent CA-CA distance of C) involves no reference structure of any kind.
   L6's impossibility argument is *"the thing it depends on is the answer"*; g depends only on
   C and on a covalent constant. The argument therefore does not reach it.
2. **L6 measured the scale on the CLOUD, where scale is the answer; this lane measures it on
   the BUILT CHAIN, where scale cannot be the answer.** The emitted chain is built by
   `core.project.build_ca_exact` at ideal covalent geometry: its CA-CA distance is 3.80 A by
   construction, whatever s was applied to the input. **Rescaling the cloud cannot rescale the
   output.** It can only change WHICH ideal-geometry chain is nearest, i.e. the shape and the
   branch. L6's mechanism (a scale error propagating into Kabsch RMSD because Kabsch does not
   fit scale) is structurally unavailable here.
3. **L6's lever was terminal; this lane's is upstream of a non-linear operator.** s* multiplies
   the emitted structure. g multiplies the INPUT to a multi-start non-convex optimisation whose
   branch degeneracy is worth 0.08 A to a perfect chooser (S26 L88) and whose branch flips under
   a 1e-13 input perturbation (S28-L18, S28-L27b). A monotone rescale of the input is not a
   monotone transformation of the output.
4. **L6's own diagnosis names this lane's quantity as the one it did NOT test.** `s23/LEDGER.md`
   L1 lines 14-19: *"A single global scale cannot fix that -- it trades one error for the other
   ... Only changing the space you average in (torsions) or repairing geometry afterwards can."*
   Repairing the geometry BEFORE the projection is the third option, and it was never run.

What this lane therefore does NOT claim and will not claim: that a native-free predictor of s*
exists. L6 stands. ORACLE-SCALE below is L6's s* re-measured on the built chain, labelled ORACLE
in every line, purely as the ceiling of this operator family.

## 3. WHICH CHARTER FINDINGS THIS ATTACKS (contract rule 13)

Attacks **none of 1-11 directly**; it attacks a stage the findings do not name. The hypothesis
is not about recognition, objectives, Hamiltonians or pools: it is that **+0.159 A of the
shipped 3.2126 is spent by a biased fit at the last stage, and part of it is native-free
recoverable**.

Accepted as binding: 1, 2, 3, 4, 5, 6, 7, 9 (nothing here touches the objective or the circuit).
**Finding 8 (recognition)**: this lane needs no recogniser -- g, s_SPAN and s_ISO are determined,
not selected, so there is no ranking step and no recognition failure to inherit. **Finding 11
(68% common-mode error)**: the common-mode error is exactly WHY the cloud is contracted (the mean
of points on a curved manifold lies inside it, `s23` L4), so this lane does not try to cancel the
common-mode error -- it tries to stop the projection from being MISLED by its most visible
signature. If it works, it works because the contraction is common-mode, which is the one way
finding 11 could be an asset.

## 4. THE DATA PATH (contract non-negotiable: both sides of every contrast on ONE code path, S28-L43)

Every arm in this lane takes the SAME production cloud and the SAME projection call. No arm is
read from a stored row of another job.

- `s29/s29_P_scale.py :: production_cloud(pdb)` reproduces `s27/s28_B_prodcheck.py ::
  project_production` exactly: `RP.channels_for(pdb)` -> `RP.rng_for(pdb,"tiekey")` ->
  `RV.energy_for("DIS", ...)` -> `np.lexsort((key, E))[:75]` -> `H.readout_uniform`. That recipe
  reproduced S27's `chain_rows.jsonl :: DIS` **bit-exactly on 126/126** (`s27/results/
  s28_B_prodcheck.json`: `n_identical_chain` 126, `max_d_chain` 0.0, mean 3.212625219723354).
- The cloud is cached once per target to `s29/results/s29_P_clouds/<pdb>.npz` (float64, lossless)
  and asserted to reproduce `rmsd_cloud` from `chain_rows.jsonl :: DIS` with **difference exactly
  0.0** before any arm runs. A target that fails that assertion is dropped and named.
- Every arm is `C_s = mu + s * (C - mu)`, mu = C.mean(0), followed by
  `I.project(C_s, seq, fold)` with the shipped defaults (lam 0.3, multi=True, maxiter=300),
  scored by `I.ca_rmsd(ca, nat_ca)`. **s = 1.0 is the exact float identity on C**, so arm PROD is
  the production path with zero numerical difference.
- The native is read ONLY by the scorer of the emitted chain, and by ORACLE-SCALE, which is
  labelled ORACLE everywhere it appears.

## 5. THE ARMS

s below is the factor applied about the cloud centroid. Seven registered arms (the brief's),
plus three mechanism controls declared here and counted in the multiplicity budget.

| # | arm | s | native-free? | role |
|---|---|---|---|---|
| 1 | **PROD** | 1.0 | yes | the gate and the comparator |
| 2 | **BOND** | g = 3.80 / mean_i \|C_{i+1} - C_i\| | yes | candidate |
| 3 | **SPAN** | Rg_post / Rg(C) | yes | candidate |
| 4 | **ISO** | sum d_ij m_ij / sum d_ij^2 over posterior pairs | yes | candidate |
| 5 | **CTRL-RAND** | 8 derangements of the realised g across targets | control | matched magnitude, zero per-target information |
| 6 | **CTRL-INV** | 1/g | control | same magnitude, wrong direction |
| 7 | **ORACLE-SCALE** | per-target argmin over an s-grid of the BUILT CHAIN's RMSD | **ORACLE** | the ceiling of the whole family |
| 8 | **CTRL-GLOBAL** | mean(g) over all 126 (a constant; g is native-free so no fold split is needed) | yes | the zero-per-target-information version of BOND at the same central magnitude |
| 9 | **CTRL-LAM** | 1.0, projected with lam = 0.3 * mean(g) | yes | isolates the effective-lambda confound (see 5.1) |
| 10 | **BOND-LAMFIX** | g, projected with lam = 0.3 * g | yes | the lambda-neutral version of BOND |
| F | **FLOOR** | 1.0, on C perturbed by 1e-13 relative | -- | the branch-flip floor measured on THIS code path |

**Definitions, fixed here.**
- `m_ij` is the posterior's **L1-Bayes (median) distance** for pair (i,j): `grid[argmin(risk_ij)]`
  from `s12.instrument.distogram(pdb)`, whose `risk` row is the L1 Bayes risk on the grid
  2..40 by 0.05. The posterior supplies pairs with |i-j| >= 2 only.
- `Rg_post^2 = (1/N^2) * sum_{i<j} m_ij^2` with the |i-j| = 1 terms completed at the ideal
  3.80 A (the classical identity Rg^2 = (1/2N^2) sum_ij d_ij^2 needs every pair; the completion
  constant is covalent, not native). `Rg(C)` is computed from C's coordinates directly.
  A **restricted variant** s_SPANr = sqrt( sum_{|i-j|>=2} m_ij^2 / sum_{|i-j|>=2} d_ij^2 ), which
  injects no 3.80, is reported as a DIAGNOSTIC number (not an arm, not projected).
- ISO uses the posterior's own pair set (|i-j| >= 2) and no completion, so SPAN and ISO differ in
  the functional (RMS ratio vs regression slope) and, deliberately, SPAN alone carries the bond
  constant. By Cauchy-Schwarz s_ISO <= s_SPANr always; their gap measures how well the cloud's
  map is PROPORTIONAL to the posterior's.
- CTRL-RAND's derangements are drawn with `np.random.default_rng(sha256("s29P|ctrlrand|<draw>"))`
  and verified to fix no target (no target receives its own g). The multiset of factors is
  exactly the multiset of realised g, so the control is matched in magnitude by construction
  ("a control must match the operator's space", the project's most repeated error).
- ORACLE-SCALE's grid: s in {0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30}
  (1.00 reused from PROD; 9 new projections per target). The grid is priced as the order
  statistic it is (`ST.best_of_k_within`, k = 10) and reported ONLY as a ceiling. The full
  curve (mean built-chain RMSD as a function of s) is the lane's main mechanism deliverable.

### 5.1 The confound this lane must not fall into
Rescaling C by s scales the RMSD term of the projection objective by roughly s while leaving
`lam * ramah` unchanged, so **BOND silently weakens the torsion prior by about 1/g ~ 0.78**. If
BOND moves the answer, the move could be the geometry OR the effective lambda. CTRL-LAM (lambda
changed, geometry unchanged) and BOND-LAMFIX (both changed so the ratio is restored) complete a
2x2 that separates them. These two arms are **mechanism diagnostics only**: they are read only
if BOND is non-null, and they can never themselves be the positive.

## 6. THE FALSIFIER (registered)

**F-P1 (primary).** BOND beats PROD on the **BUILT CHAIN**, n = 126, with

  (a) |effect| >= 1.0x its own MDE, and
  (b) the **fold-clustered CI excluding zero**, and
  (c) **5/5 folds same sign**, and
  (d) BOND beats **CTRL-RAND's mean of 8** by >= 1.0x that contrast's MDE.

All four clauses must hold. SPAN and ISO carry the identical falsifier (F-P2, F-P3).

**Multiplicity (contract rule 17).** Three deployable candidate arms are tested against PROD;
with CTRL-GLOBAL there are four deployable contrasts, plus 2 mechanism contrasts, plus the
ORACLE grid (priced as an order statistic, never a result). **The lane's endpoint comparison
count is 4.** A positive at 1.0-1.3x MDE among four is in the Type-M zone AND under a
max-over-4 null: it will be reported with the max-over-K null stated and called SUGGESTIVE, not
demonstrated, and a second seed / second tie-key is owed before any endpoint claim.

**The branch-flip floor gate.** S28-L27b measured a mean floor near 0.006 A with a 0.5 A tail on
2LNG for contrasts crossing code paths. Here both sides share one path and one input, so the
floor formally does not apply -- but the rescale IS an input change, so arm FLOOR measures the
real thing on this path: how many of the 126 move by more than 0.02 A when C is perturbed by
1e-13 relative. **Any mean effect below the measured FLOOR mean will be reported as inside the
floor, whatever its CI says.** This number is printed beside every small effect in the entry.

## 7. THE REGISTERED PRIOR

**Null: 0.0 to 0.3x MDE for BOND, SPAN and ISO on the built chain.** Reasons, stated so the prior
is falsifiable:

1. The projection already enforces ideal geometry. The RMSD term's minimiser over (phi, psi) is
   dominated by the cloud's SHAPE; a similarity transform of the target changes the objective's
   level far more than its argmin's shape.
2. `s23` L1's trade: BOND's g ~ 1.28 fixes the 22% bond deficit but inflates the envelope by 28%
   against a cloud whose envelope is only 6% short -- it over-expands globally by ~20%. SPAN's
   s ~ 1.06 fixes the envelope and leaves the bonds 22% short. Neither is the right object;
   **the contraction is not a uniform shrink** and one scalar cannot undo a non-uniform
   distortion. This is a mechanism prediction, not a hedge: it predicts BOND and SPAN move the
   answer in OPPOSITE directions if they move it at all, and that the ORACLE-SCALE argmin sits
   between 1.00 and g, nearer 1.00.
3. Lane T's theorem 2 (S29-L7, addendum 1 rule 21): the posterior's median map is an
   unconstrained per-pair object that over-deviates from typical (beta > 1), so SPAN and ISO
   should over-expand relative to the truth; and neither is a gradient claim, so addendum rule 20
   applies only through its clause (c), which this prereg adopts (emitted bond and Rg reported).

**Registered mechanism prediction, separate from the endpoint prior.** The **projection price**
(chain minus cloud, per arm) **FALLS under BOND**, toward the +0.026 A that S28-L39 measured for
a 1 A de-contracting displacement (`s27/LEDGER.md` S28-L39, the `chain_projection_price` line:
prod +0.159, step_e1 +0.026, rand1_e1 +0.059). If the price falls and the RMSD does not, the
reading is stated in advance: **the +0.159 A is not a recoverable loss but the price of the
geometry constraint, and paying it earlier costs the same.** That reading is the lane's deliverable
if the prior holds, together with the ORACLE-SCALE ceiling.

**If the prior holds I will say so in one line and report the ORACLE-SCALE ceiling as the number
that matters** (brief, line 49).

## 8. DIAGNOSTICS REGISTERED IN THE SAME ENTRY

1. The realised **g** per target: distribution (mean, sd, deciles), Spearman with n (length),
   with production's built-chain RMSD, with production's cloud RMSD, and with the ORACLE s*.
2. s_SPAN, s_ISO, s_SPANr distributions; their correlation with g and with each other.
3. The **projection price per arm** (mean chain RMSD minus mean cloud RMSD, and the per-target
   paired version).
4. Emitted **mean virtual CA-CA bond and Rg** of the cloud before and after each rescale, and of
   the built chain (addendum 1 rule 20(c)).
5. **FAIL18 vs the 108**, each through `ST.compare`, with a **random-18 null** (1000 random
   18-subsets, the observed FAIL18 effect's percentile) before any stratum claim is made.
6. The **branch-flip floor** on this path (arm FLOOR): mean |diff|, max, count above 0.02 A.
7. The **s-curve**: mean built-chain RMSD against s over the ORACLE grid, and the histogram of
   the per-target argmin.
8. Wall time and peak RSS of every job.

## 9. THE METER (contract rule 19)

**No new cost function is defined or optimised by this lane.** Every arm uses the SHIPPED
projection objective, unchanged, including its lam and its multi-start. g, s_SPAN and s_ISO are
scalars that calibrate the projection's INPUT; they do not rank, score or select any structure,
so there is nothing for the ladder Spearman, the gradient cosine, the native percentile or the
ORACLE preference to be computed over. **If and only if BOND, SPAN or ISO produces an endpoint
positive**, the native-free discrepancy that motivates the winning arm is exposed as a cost over
clouds (`s29.s29_P_scale:cost_isoresid`, the least-squares residual of the cloud's map against
the posterior median map after the optimal scalar) and run through
`python s29/s29_D_cost_audit.py meter --f s29.s29_P_scale:cost_isoresid --basis ca`
BEFORE the claim is written into the ledger. Lane D is invited to attack in either case.

## 10. THE PROBE AND ITS GATE (contract rule 16)

Six targets, registered here before they are run:
**1A13** (n=14, the record's reference target), **1CS9** (n=9, the shortest; S28-L18's +0.0186
chain difference), **2LNG** (n=13, the -0.513 A branch-flip tail of S28-L27b), **2NB7** (n=14,
FAIL18), **9BFL** (n=10, -0.160 in S28-L27b), **9BAF** (n=16, the longest).

**GATE: arm PROD must reproduce `chain_rows.jsonl :: DIS`'s `rmsd_chain` on all six to
< 1e-9 A** (the recipe is bit-exact on 126/126 in `s28_B_prodcheck.json`, so the expected
difference is 0.0). If any target misses the gate, the 126-target run does not start and the
failure is reported instead. The probe runs under
`python s26/jobrun.py --agent S29P --tag CPU --name s29P_probe6 --est-ram 0.5` and its peak RSS
is quoted in the ledger entry. **The probe is never quoted as evidence for the instrument.**

## 11. TESTS (`tests/test_s29_P.py`), all written before the 126 run

1. s = 1.0 rescale returns the input array bit-for-bit, and PROD's projection equals a direct
   `I.project(C, ...)` call bit-for-bit.
2. The Rg identity: Rg from coordinates equals Rg from the pair-distance identity to < 1e-12 on
   random clouds and on a real cached cloud.
3. BOND: the rescaled cloud's mean adjacent CA-CA distance equals 3.80 to < 1e-12.
4. ISO: the analytic slope equals a 1-D numerical minimiser of the same least-squares objective.
5. SPAN: Rg(rescaled) equals Rg_post to < 1e-12.
6. CTRL-INV = 1/g exactly; CTRL-RAND's derangements fix no target and preserve the multiset.
7. **NaN-poison**: every deployable scale function is called with a context whose `nat_ca` is NaN
   and must return a finite scale; a cloud containing NaN must propagate NaN (never silently
   return 1.0).
8. The posterior median is the argmin of the stored L1 risk row (independent recomputation from
   `prob` and `CENTRES`).

## 12. WHAT WOULD MAKE ME RETRACT

If arm PROD fails the 1e-9 gate; if any arm's scale function is found to read a native quantity;
if the cached cloud fails the exact `rmsd_cloud` assertion on any target; if a positive does not
survive the CTRL-RAND mean, the CTRL-GLOBAL constant, the FLOOR, the max-over-4 null and a second
tie-key seed. Each of these is checked and reported, not assumed.

---

## ADDENDUM 1 (2026-09-20 00:12, coordinator instruction; execution only, no scientific change)

The 126-target run is **sharded into 4 concurrent jobs by target range** (pinned pdb-sorted
order, contiguous quarters), each with its own rows file
`s29/results/s29_P_rows_shard<k>.jsonl` and its own **per-(arm, target) checkpoint**; the
analysis step concatenates the four and asserts 126 x (number of arms) complete rows with no
duplicate (pdb, arm) key. Sharding changes nothing numerically: each (pdb, arm) projection is a
deterministic function of the cached cloud, the scalar s and lam, computed in a single-threaded
process (`OMP/MKL/OPENBLAS_NUM_THREADS = 1`), and the PROD arm's bit-exact reproduction gate is
re-asserted **per target inside the 126 run**, not only on the probe. The 6-target probe and its
1e-9 gate run first, unsharded, exactly as pre-registered in section 10.

## ADDENDUM 2 (2026-09-20 00:21; execution only, no scientific change)

Inside each shard the work runs in two phases: **primary** (PROD, BOND, SPAN, ISO, CTRL-INV,
CTRL-GLOBAL, CTRL-LAM, BOND-LAMFIX, the 8 CTRL-RAND draws and FLOOR -- every cell the
pre-registered falisifer F-P1/F-P2/F-P3 needs) and then **oracle** (the 9 remaining points of
the ORACLE s-grid, which is a labelled ceiling and not part of any falsifier). The box is
CPU-saturated (8 to 9 concurrent single-threaded jobs against 6.43 core-equivalents, the
governor suspending on the 94% CPU ceiling), so ordering the work this way lets the endpoint
verdict be read before the ceiling grid completes. Each cell is a deterministic function of the
cached cloud, its scalar and lam, so phase order changes no number. The analysis reports the
primary-complete target count and the grid-complete target count separately and refuses to form
the ORACLE contrast until the grid is complete on every analysed target.

## ADDENDUM 3 (2026-09-20 00:23, before any 126-target number; a CONSERVATIVE amendment)

CTRL-INV (s = 1/g) was registered in section 5 as a *direction control*. It is, however,
computed from the cloud alone and is therefore **deployable**, exactly like BOND, SPAN, ISO and
CTRL-GLOBAL. If it were to beat PROD it would be a native-free operator, not a control result,
and it must pay the same multiplicity price. It is therefore moved into the max-over-K set,
making **K = 5** deployable contrasts rather than 4. This raises the bar and lowers nothing:
a positive on CTRL-INV is additionally a POST-HOC finding (it was not a registered candidate),
and the entry will say so, will price it against the max-over-5 sign-flip null, will require the
contraction check (emitted Rg against the native's) and will owe a second tie-key seed before
any endpoint claim. The registered falsifiers F-P1/F-P2/F-P3 for BOND, SPAN and ISO are
unchanged.

## ADDENDUM 4 (2026-09-20 00:28, before any 126-target number)

**A derived, native-free, deployable readout that costs no extra projection: THE WIDER
MULTI-START.** Every arm's emitted structure is an ideal-geometry chain, hence a FEASIBLE POINT
of the shipped projection problem (minimise RMSD(chain, C) + 0.3 * ramah over ideal-geometry
chains, where C is PRODUCTION's own unmodified cloud). Each row therefore also records
`obj0` = that shipped objective evaluated for the arm's chain against production's cloud. This
makes the whole arm set a wider multi-start of `core.project.fit_multi`, whose four generic
starts demonstrably under-optimise (its own docstring; S26 L88 prices the branch degeneracy at
0.08 A to a perfect chooser). Three readouts, all on the built chain:

- **MS-OBJ** (deployable, native-free): the chain of the arm with the lowest `obj0`, ties
  averaged over the argmin set with `ST.argmin_tied` (never broken by array order).
- **MS-MEAN** (zero-information control): the mean outcome over the same arm set, i.e. the
  expectation of picking an arm at random, matched in selection budget.
- **MS-ORACLE** [ORACLE]: the minimum RMSD over the same arm set, the ceiling of selection.

**Registered prediction, written before the numbers:** MS-OBJ strictly lowers the quantity
production minimises (by construction, whenever any arm beats PROD's `obj0`) and **does NOT
lower the built-chain RMSD** -- it is the projection-stage instance of "optimise the shipped
objective harder, get a worse structure" (charter finding 3; S21's budget trap; S28-L26b).
If MS-OBJ instead lowers RMSD beyond 1.0x MDE with the fold CI excluding zero, that is a
deployable improvement to production's own solver and it is reported as the lane's primary
result, with the second tie-key seed and lane D's attack owed before any claim. MS-OBJ is
POST-HOC relative to sections 5 and 6 and is labelled so; it is counted in the lane's
multiplicity as one further contrast (K = 6).

## ADDENDUM 5 (2026-09-20 00:50; CHRONOLOGY, so the write-up cannot be read as hindsight)

Lane T's achievable native-free bound (**S29-L23**, 00:43) prices every native-free operator as a
displacement worth exactly one cosine, `RMSD = RMSD_prod * sqrt(1 - rho^2)`, and reports every
field the project has built at or below **rho = 0.04** against the **0.358** that 3.00 A requires
and the **0.628** that 2.50 A requires, with a random-shape-field reference of 0.140. Every arm
of this lane -- BOND, SPAN, ISO, CTRL-INV, CTRL-GLOBAL -- is a displacement of exactly that kind,
so the bound predicts all of them null or harmful.

**The chronology matters and is recorded here, not argued later.** This lane's falsifiers
F-P1/F-P2/F-P3 and its registered null prior (0.0 to 0.3x MDE) were written at **00:10**, before
S29-L23 existed, and were derived from S23 L1's non-uniform contraction and the geometry of the
projection -- not from T's bound. S29-L22 result 4 (the monotone shape distortion, 00:37) then
independently predicted BOND harmful by mechanism, again before the bound was published.

**Consequences for the second ledger entry, fixed now:**
1. The falsifiers are honoured **exactly as registered**. The prior is not revised, the MDE
   thresholds are not moved, and no arm is dropped because a bound now predicts it null.
2. The entry is written as a **PRE-REGISTERED CONFIRMATION, not a discovery**: it says plainly
   that two independent arguments (T's cosine bound, S29-L23; this lane's own separation profile,
   S29-L22 result 4) predicted the outcome in advance, and it gives the dates.
3. **If any arm nevertheless clears its falsifier**, that is evidence AGAINST the bound and is
   reported as such -- loudly, with the second tie-key seed run before any claim -- rather than
   explained away. A bound that forbids a measured effect is the bound's problem.
4. The lane's value in that case is the same either way: the projection was the one stage the
   sprint had not audited, and it is now audited on the reporting basis with matched controls.
