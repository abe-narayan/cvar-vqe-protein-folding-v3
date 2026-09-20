#!/usr/bin/env python
"""s29/s29_P_ledger1.py -- lane P's FIRST ledger entry, appended in ONE process.

Numbers the entry from the ledger's tail, stamps it from the system clock read in this same
process (`check-the-clock-before-stamping`), regenerates every ST.fmt block verbatim from the
stored arrays, and appends. Append-only: it never edits an existing entry.
"""
from __future__ import annotations

import datetime
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST            # noqa: E402

LEDGER = os.path.join(HERE, "LEDGER.md")
RES = os.path.join(HERE, "results")


def next_number():
    txt = open(LEDGER, encoding="utf-8").read()
    nums = [int(m) for m in re.findall(r"^## S29-L(\d+)", txt, re.M)]
    return max(nums) + 1


def main():
    n = next_number()
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    sp = json.load(open(os.path.join(RES, "s29_P_sepprofile_cloud.json"), encoding="utf-8"))
    order = sp["order"]
    fac = json.load(open(os.path.join(RES, "s29_P_factors.json"), encoding="utf-8"))
    folds = np.array([fac["per_target"][p]["fold"] for p in order])
    base = np.array(sp["base"]); nf = np.array(sp["nf"]); orc = np.array(sp["orc"])
    b1 = ST.fmt(ST.compare(nf, base, folds=folds, names=order,
                           label="P sep-profile NATIVE-FREE - production (POINT CLOUD, diagnostic)"))
    b2 = ST.fmt(ST.compare(orc, base, folds=folds, names=order,
                           label="P sep-profile [ORACLE, leave-fold-out native-fitted] - production (POINT CLOUD, diagnostic)"))

    prof = json.load(open(os.path.join(RES, "s29_P_contraction_profile.json"), encoding="utf-8"))
    probe = json.load(open(os.path.join(RES, "s29_P_probe.json"), encoding="utf-8"))
    sstar = json.load(open(os.path.join(RES, "s29_P_oracle_cloud_sstar.json"), encoding="utf-8"))
    s = np.array([sstar[p] for p in order])
    g = np.array([fac["per_target"][p]["g"] for p in order])

    def row(k, vals, fmt="%6.3f"):
        return "    %-22s %s" % (k, " ".join(fmt % v for v in vals))

    body = f"""
## S29-L{n} -- LANE P, THE PROJECTION-STAGE MECHANISM (the 126-arm run is still in flight; these six do not depend on it): PRODUCTION'S CLOUD REPRODUCES BIT-EXACTLY ON 126/126 AND ITS BUILT CHAIN ON 6/6 AT EXACTLY 0.0; THE +0.164 A PROJECTION PRICE TRACKS THE CLOUD'S CONTRACTION AT SPEARMAN +0.604 (+0.481 PARTIALLED); THE AVERAGING DISTORTION IS NOT A CONTRACTION BUT A MONOTONE SHAPE DISTORTION IN SEQUENCE SEPARATION (0.773 AT |i-j|=1, CROSSING 1.00 NEAR 8, 1.10 AT 13); CORRECTING THAT PROFILE IS REFUTED AT +0.578 A (2.83x MDE, 5/5 FOLDS) WITH AN ORACLE-FITTED PROFILE NO BETTER (+0.582) AND AN EXACT r==1 CONTROL; AND S23 L2/L3 REPRODUCE ON THESE CLOUDS ({now}, P)

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
(`s29_P_oracle_cloud_sstar.json`; no projection): s* mean **{s.mean():.4f}**, sd {s.std(ddof=1):.4f},
range {s.min():.3f} to {s.max():.3f}, **{100*float((s>1).mean()):.1f}% want EXPANSION** -- S23 L2 reported
53/126 = 42.1% wanting expansion against 73 wanting contraction; the split, the spread and the
bimodality all reproduce. Spearman(s*, g) **-0.117**, (s*, s_SPAN) -0.033, (s*, s_ISO) +0.000:
inside the +-0.11 band S23 L3 measured for eight other native-free candidates, and the only one
outside it has the WRONG SIGN. Mechanically: g is an expansion (mean {g.mean():.3f}, above 1 on
125/126) while the RMSD-optimal cloud scale is a contraction on 57% of targets. On the POINT-CLOUD
basis BOND must be harmful; the lane's question survives only because the built chain's scale is
fixed by covalent geometry, which is exactly the distinction above.

(d) THE AVERAGING DISTORTION IS NOT A CONTRACTION. IT IS A MONOTONE SHAPE DISTORTION IN SEQUENCE
SEPARATION. `s29_P_contraction_profile.json`, all 126 clouds, every pair pooled at each
separation:

{row("|i-j|", prof["sep"], "%6d")}
{row("cloud/native [ORACLE]", prof["rom_nat"])}
{row("  median ratio", prof["med_nat"])}
{row("cloud/posterior (NF)", prof["rom_post"])}

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

{b1}

{b2}

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
"""
    with open(LEDGER, "a", encoding="utf-8") as fh:
        fh.write(body)
    print("appended S29-L%d at %s" % (n, now))


if __name__ == "__main__":
    main()
