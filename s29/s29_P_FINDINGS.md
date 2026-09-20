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

## HYPOTHESIS

## REFUTED

## OPEN

## WHAT DAMAGED MY OWN EXPECTATIONS

## WHAT I DID NOT DO AND WHY
