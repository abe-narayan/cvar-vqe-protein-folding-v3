# s29_P_FINDINGS (Sprint 29, lane P: the projection price)

Standing format (S12 to S28): DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN,
every number with its artefact path, then "what damaged my own expectations" and "what I did not
do and why". Pre-registration `s29/PREREG_S29_P.md` (written 2026-09-20 00:10, before any
number). Code `s29/s29_P_scale.py`, tests `tests/test_s29_P.py`, results `s29/results/s29_P_*`.
Ledger entries `## S29-L<n> -- ... (date time, P)` in `s29/LEDGER.md`.

The lane's object: production's point cloud is 3.0483 A and its built chain is 3.2126 A, so the
last stage of the pipeline COSTS +0.1643 A (measured here on all 126, not quoted from the
record). That stage is the only one no S28 lane touched. The cloud is geometrically inconsistent
with the space it is projected into. The lane asks whether removing that inconsistency before
the projection changes the built chain, and what the ceiling of the whole rescale family is.

---

## DEMONSTRATED

P1. **The production point cloud reproduces bit-exactly in this lane's process on 126/126.**
`s29/results/s29_P_factors.json` (`n_cloud_mismatch` 0): every target's `rmsd_cloud` recomputed
through `s29_P_scale.py :: production_cloud` (the `s27/s28_B_prodcheck.py :: project_production`
recipe) differs from `s27/results/chain_rows.jsonl :: DIS` by **exactly 0.0**. Job
`s26/jobs_done/s29P_factors.json` exit 0, 20.0 s, peak RSS 0.319 GB. Consequence: every contrast
in this lane has both sides on one code path and one input, so the S28-L18 / S28-L27b
cross-path floor does not enter; the *input-perturbation* floor is measured separately (arm
FLOOR).

P2. **S23 L1's geometry reproduces on the current production cloud, to four decimals.**
Mean adjacent CA-CA distance **2.9614 A** (S23 L1: 2.961), radius of gyration **6.2061 A**
(S23 L1: 6.206), against the ideal covalent 3.80 and a native Rg near 6.60. The cloud is 22%
short in the bond and 6% short in the envelope: **the contraction is not a uniform shrink**, and
that is a fact about the current pipeline, not only about the S23 one.

P3. **The bond-consistency factor g is far from a global constant.**
g = 3.80 / (mean adjacent CA-CA distance), native-free: **mean 1.352, sd 0.318, range 0.998 to
2.141**, deciles 1.002 / 1.026 / 1.342 / 1.611 / 1.758 (`s29_P_factors.json :: per_target`).
A third of the panel wants essentially no rescale and a tenth wants more than 1.75x. This is
why the matched-random control (a derangement of the realised g across targets) is a real
control here and not a near-copy of BOND.

P4. **The posterior's implied envelope is larger than the cloud's, but only by 10%.**
s_SPAN (the factor that makes the cloud's Rg equal the Rg implied by the distogram posterior's
L1-Bayes median map) has **mean 1.102, sd 0.084**; s_ISO (the least-squares slope of the cloud's
map onto the posterior median map) **1.077, sd 0.073**; the restricted variant s_SPANr 1.099.
**s_SPAN and s_ISO agree at Spearman +0.989** but correlate with g only at +0.63 and +0.56:
the two native-free rescale signals are NOT the same quantity, and they differ by a factor of
about 3.5 in magnitude (1.10 vs 1.35).

*(sections below completed when the 126-target run lands)*

## ORACLE DIAGNOSTIC

P5. **[ORACLE] The RMSD-optimal cloud scale points the OTHER WAY from every native-free
rescale, and S23 L2/L3 reproduce exactly on this lane's clouds.**
`s29/results/s29_P_oracle_cloud_sstar.json`, S23 L6's closed form s* = sum(svd(Cc^T Tc)) /
||Cc||^2 evaluated on the 126 cached clouds against their natives (ORACLE; computed with no
projection): **s\* mean 0.9481, sd 0.2565, range 0.277 to 1.924, and 42.9% of targets want
EXPANSION** -- S23 L2 reported 53/126 = 42.1% wanting expansion against 73/126 wanting
contraction. The bimodality, the spread and the split all reproduce.
**The native-free factors are orthogonal to it**: Spearman(s\*, g) **-0.117**, (s\*, s_SPAN)
**-0.033**, (s\*, s_ISO) **+0.000** -- inside the +-0.11 band S23 L3 measured for eight other
native-free candidates, and the sign of the only one outside it is NEGATIVE.
**Mechanistically decisive, and registered before the 126-target run reports:** g is an
expansion on 125 of 126 targets (mean 1.352) while the RMSD-optimal cloud scale is a
contraction on 57% of them (mean 0.948). On the POINT-CLOUD basis BOND must therefore be
harmful. The lane's question survives that only because the BUILT CHAIN's scale is fixed by
covalent geometry, so the cloud-basis argument does not transfer -- which is exactly the
distinction the prereg's section 2 rests on, now with a number attached.

P6. **[ORACLE for the price, native-free for g] The lane's premise HOLDS: the projection costs
more exactly where the cloud is more contracted.**
Projection price (built chain minus point cloud, per target, from `s27/results/chain_rows.jsonl
:: DIS`): **mean +0.1643 A, sd 0.2014, median +0.0958, positive on 111/126, range -0.347 to
+0.652**. Spearman with the native-free contraction factor g: **+0.604**; partialling out the
cloud's own RMSD and the chain length, **+0.481**. (Against s_SPAN +0.378, s_ISO +0.339, length
-0.142, cloud RMSD +0.511, ORACLE s* -0.313.)
So the +0.164 A is not a flat toll: it is concentrated on the targets whose averages are most
geometrically inconsistent with the ideal-geometry space they are projected into. That is the
mechanism the lane hypothesises, measured before the intervention. It does NOT establish that
removing the inconsistency removes the price -- the correlation is equally consistent with
"contracted clouds are the ones that disagree most, and disagreement is what the projection
cannot repair" (S29-L12: the contraction is Jensen's inequality, and Jensen's shrink is not
uniform, so it cannot be undone by one scalar). The 126-target run separates those two.

P7. **The lane sits on lane L's perception-distortion curve, and the framing was theirs before
it was mine.** S29-L12 (lane L, topic 3) reads the same +0.164 A as "a MEASURED point on P(D)"
and says the probe "should be framed that way -- interpretable whichever way it comes out":
the coordinate average is the low-distortion off-manifold end (Jensen contraction, 22% in the
bond) and the ideal-geometry chain is the on-manifold end. Under Blau & Michaeli's Theorem 3
the +0.164 A is the local slope of returning to the manifold. Two consequences this lane adopts:
(a) a rescale can only move ALONG the curve, it cannot move the curve, so a scalar that lowers
RMSD by making the emitted chain more compact is buying distortion with realism and must be
caught -- hence the contraction check (emitted Rg against the native's) in the analysis;
(b) the Jensen shrink is NOT uniform (bond 22% short, envelope 6%), so one scalar cannot undo
it, which is the mechanical form of this lane's registered null prior.

P8. **A second, free question the same run answers: is production's projection even SOLVED?**
Every arm emits an ideal-geometry chain, which is a feasible point of the shipped projection
problem for production's own cloud. Recording the shipped objective `obj0` for each arm turns
the arm set into a native-free WIDER MULTI-START of `core.project.fit_multi`, whose docstring
already concedes that its four generic starts under-optimise and which S26 L88 prices at 0.08 A
to a perfect chooser. PREREG addendum 4 registers MS-OBJ (native-free pick), MS-MEAN (the
zero-information pick) and MS-ORACLE, and registers the prediction that MS-OBJ lowers the
objective and NOT the RMSD -- the projection-stage instance of charter finding 3.

P9. **[native-free profile confirmed against an ORACLE one] The averaging distortion is NOT a
contraction. It is a MONOTONE SHAPE DISTORTION IN SEQUENCE SEPARATION, and that is why no single
scalar can undo it.** `s29/results/s29_P_contraction_profile.json`, all 126 cached clouds,
pooled over every pair at each separation (ratio of means; the median ratio beside it):

    |i-j|      1      2      3      4      5      6      7      8      9     10     11     12     13
    cloud/native   .773   .824   .830   .870   .921   .948   .968   .996  1.029  1.047  1.065  1.083  1.102
    median ratio   .784   .857   .886   .912   .960   .975   .982   .997  1.009  1.013  1.029  1.038  1.048
    cloud/posterior-median (NATIVE-FREE)  .776  .821  .809  .846  .902  .913  .924  .955  .965  .967  .971 1.001  .952

The curve is monotone, starts 23% SHORT at the virtual bond and CROSSES 1.00 at |i-j| ~ 8,
ending slightly LONG. The native-free profile (against the distogram posterior's own median map)
has the same shape and the same crossing region, so this is not an ORACLE-only statement.
Cloud mean Rg 6.2061 against the native's 6.6009 (6% short): S23 L1 reproduced.

**Consequences, all of them mechanical:**
(a) S23 L1's "22% short in the bond, 6% short in the envelope" is the two ENDS of this curve;
    the middle is where the distortion changes sign.
(b) BOND's g = 1/0.773 = 1.29 sets the bond exactly right and inflates every separation beyond
    8 by ~29% on top of distances that are ALREADY long. BOND is predicted harmful, and its
    harm should grow with chain length -- this is the registered null prior's mechanism, now
    measured rather than argued.
(c) The best SINGLE scalar in least squares must sit near the middle of the curve, around
    1.05 to 1.08. s_ISO is 1.077 and s_SPAN is 1.102, so ISO and SPAN are near-optimal members
    of a family whose optimum is near 1.0 -- which is exactly why the registered prior for them
    is null rather than harmful.
(d) It is the direct measurement of the Jensen mechanism S29-L12 names: averaging shrinks a
    distance by more when the pool disagrees more about it, and the pool disagrees most about
    LOCAL geometry (adjacent-residue placement) relative to that distance's size.

**METHOD NOTE, because I got this wrong first.** My first pass used the MEAN OF PER-PAIR RATIOS
and reported the sep-13 ratio as 1.37. That statistic is dominated by small denominators (a
native hairpin puts two residues 13 apart very close in space, and the ratio explodes). The
ratio of means and the median ratio agree with each other at 1.10 and 1.05 and are the honest
numbers; the 1.37 is withdrawn before it was ever used. The qualitative conclusion -- monotone,
crossing 1 near |i-j| = 8 -- is identical under all three statistics.

## HYPOTHESIS

## REFUTED

P10. **The shape distortion of P9 is REAL but NOT CORRECTABLE -- and an ORACLE-fitted profile is
no better than a native-free one.** `s29/results/s29_P_sepprofile_cloud.json`. Operator: divide
each cloud's distance map by the population separation profile r(|i-j|), re-embed by classical
MDS (double-centring, top-3 eigenvectors), score against the native. **POINT-CLOUD BASIS, a
diagnostic; the built chain is not measured for this arm.**

    production cloud                                                        3.0483
    separation-profile corrected, NATIVE-FREE (leave-fold-out, posterior-fitted)   3.6266
    separation-profile corrected, [ORACLE]   (leave-fold-out, native-fitted)       3.6307

    P sep-profile NATIVE-FREE - production (POINT CLOUD, diagnostic)
      effect +0.5782  median +0.4351  SE 0.0730  MDE 0.2045  effect/MDE +2.83
      fold CI95 [+0.4165, +0.7230]  folds same sign 5/5  27W/99L  VERDICT: WORSE
    P sep-profile [ORACLE] - production (POINT CLOUD, diagnostic)
      effect +0.5824  median +0.4689  SE 0.0813  MDE 0.2278  effect/MDE +2.56
      fold CI95 [+0.4343, +0.7397]  folds same sign 5/5  37W/89L  VERDICT: WORSE

**THE CONTROL IS EXACT, which is what makes this readable.** With r == 1 the same MDS returns
the production cloud at **0.00e+00 A** on the test target, so the +0.58 A is the correction and
not the re-embedding. And the ORACLE arm -- the profile fitted on the OTHER FOLDS' real natives,
i.e. the best population profile that exists -- is **as bad as the native-free one (+0.582 vs
+0.578)**, which rules out "the profile was mis-estimated" as the explanation.

**Reading.** Undistorting the distance map target-by-target destroys more than the distortion
costs: the coordinate average's value is the COHERENCE of its errors (68% common-mode, S23 L9),
and a per-pair correction fitted to a population profile breaks that coherence while replacing
it with nothing. **Scope:** this closes the POPULATION-PROFILE class of distance-space
corrections; it says nothing about a per-target profile, which S23 L6 already places out of
reach in principle. It also strengthens this lane's registered null prior a fortiori: if the
full 13-parameter profile with an oracle fit cannot help, a single scalar drawn from the same
curve is not going to.


## OPEN

P11. **The arm run, cut to the falsifiers (coordinator's decision, 2026-09-20 00:45).**
The box is CPU-bound on 6.43 core-equivalents, so extra jobs split the same throughput rather
than adding any; the coordinator declined to reassign slots from lane O's ceiling ladder and
directed lane P to run only what carries F-P1/F-P2/F-P3, in one slot. Implemented as a `--phase`
flag and a single consolidated job `s29P_primary`.

    PRIMARY  (running)  15 cells/target = 1890: PROD, BOND, SPAN, ISO, CTRL-INV, CTRL-GLOBAL,
                        CTRL-RAND0..7, FLOOR
    DEFERRED            11 cells/target = 1386: the 9-point ORACLE s-grid, CTRL-LAM,
                        BOND-LAMFIX, and with them MS-OBJ / MS-MEAN / MS-ORACLE

**Arithmetic correction, stated rather than hidden:** the primary set is 58% of the 3276 cells,
not the quarter the instruction estimated. The difference is CTRL-RAND. F-P1 clause (d) requires
BOND to beat the **mean of 8** matched-magnitude draws, so reducing the draws would weaken a
pre-registered falsifier to save compute; all 8 are kept. CTRL-GLOBAL is kept because it is
native-free and therefore sits in the max-over-K set beside CTRL-INV. CTRL-LAM and BOND-LAMFIX
are deferred rather than dropped because the prereg reads them ONLY if BOND is non-null.
`analyse` gates MS-OBJ on grid completeness, so the deferral is mechanical, not remembered.

**A cost I inflicted and should own:** consolidating the four shards, my first process filter was
too broad and terminated the running shard as well as the three queued launchers, which returned
lane P's one active slot to the queue. No computation was lost -- all 52 completed cells are
checkpointed and the resume now reads every shard rows file -- but the consolidated job is
queued behind 8 others rather than running, and that delay is mine, not the scheduler's.

**RESUME (idempotent; no cell is ever recomputed):**

    python s26/jobrun.py --agent S29P --tag CPU --name s29P_primary --est-ram 0.6 --         python s29/s29_P_scale.py run --shard 0 --nshards 1 --phase primary
    python s29/s29_P_scale.py analyse          # primary- and grid-complete counts reported apart

When slots free, the deferred half is the same command with `--phase oracle`.

P12. **METHODOLOGICAL NOTE (the 126-target run was NOT spent, and should not be spent by a
future sprint either): separation-band re-weighting of the selection objective is a GENUINELY
DIFFERENT operator from this lane's rescale -- and it is closed anyway, by a measurement that
predates S29.**

Two independent statements, and they must not be collapsed:

*(i) The operators are distinguishable -- my result 5 does NOT close this route.*
`s29/s29_P_bandsel.py`, falsifier registered in its docstring before any number (the operators
differ if a band weighting moves the emitted |i-j| = 1 ratio by more than 0.02). Replace
production's UNIFORM mean over pairs in `s12.instrument.shipped_score` with a weighted mean whose
weights depend only on |i - j|, re-select the top 75 with the same tie key and the same m,
re-average. The `uniform` arm reproduces production's cloud RMSD at **max |diff| 0.000e+00**, so
the re-scoring path IS the production path. On a 4-target probe, `local_only` (all weight on
|i-j| <= 3) moves the ratio from production's 0.885 to **0.802 -- 0.083, four times the
threshold** -- with only **46%** of the top-75 shared with production's. THE FALSIFIER FIRES.
The reason is structural: band re-weighting acts UPSTREAM of the averaging and changes WHICH
structures are averaged, while P10's rescale acts DOWNSTREAM of it on an already-averaged map.
Lane T reached the same conclusion independently and withdrew its own section 1.4 build
recommendation on those grounds (**S29-L24**).

*(ii) And it is dead on the merits, measured leave-fold-out before this sprint existed.*
`s12/obj_FINDINGS.md` section 4, all 126 targets, coordinate-average arm:

    shipped Bayes risk                    3.048   (production)
    LFO per-shell weighted L1             3.058   WORSE by +0.010
    LFO per-shell recalibrated L1         3.262
    separation-component only             3.148
    shell-profile only (DEPLOYABLE)       3.163
    per-shell standardised L1             3.382
    RANDOM score                          3.445

Every per-shell / per-band re-weighting of the objective that has been tried loses, and the one
fitted weight set (`score_weights.json`) was fitted and RETIRED from the shipped code once
already -- refitting drove two of six parameters to their clip bounds and lost 2.76 to 2.91 A on
the benchmark while gaining 65% on the dev objective. **So the route is distinguishable from
mine and separately closed.** The 126-target confirmation of (i) was dropped from the queue by
the coordinator rather than run, correctly: it would have confirmed a distinction that a 4-target
probe already establishes qualitatively, in order to characterise an operator already measured
dead. `s29/s29_P_bandsel.py` is kept, with its falsifier, so the measurement is one command away
if anyone doubts (i).

**This note exists so the route is not proposed a third time.** It has now been proposed twice
(S12's arms, lane T's section 1.4) and closed twice, on different grounds each time.

## WHAT DAMAGED MY OWN EXPECTATIONS

## WHAT I DID NOT DO AND WHY
