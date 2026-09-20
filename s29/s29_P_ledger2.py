#!/usr/bin/env python
"""s29/s29_P_ledger2.py -- lane P's SECOND ledger entry, appended in ONE process.

Numbers from the ledger tail, stamps from the clock read in this same process, and takes every
ST.fmt block VERBATIM from `s29/results/s29_P_summary.json :: contrasts` so the entry cannot
drift from the artefact.
"""
from __future__ import annotations

import datetime
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST            # noqa: E402

LEDGER = os.path.join(HERE, "LEDGER.md")
SUMM = os.path.join(HERE, "results", "s29_P_summary.json")


def main():
    txt = open(LEDGER, encoding="utf-8").read()
    n = max(int(m) for m in re.findall(r"^## S29-L(\d+)", txt, re.M)) + 1
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    S = json.load(open(SUMM, encoding="utf-8"))
    C = S["contrasts"]
    def f(k):
        """ST.fmt on a contrast restored from JSON. json.dump stringifies the integer fold keys
        of `per_fold`, and ST.fmt formats them with %d, so they are coerced back on the way in."""
        o = dict(C[k])
        if o.get("per_fold"):
            o["per_fold"] = {int(q): v for q, v in o["per_fold"].items()}
        return ST.fmt(o)
    M, P = S["means"], S["price"]
    fl = S["floor"]

    body = f"""
## S29-L{n} -- LANE P, THE VERDICT ON THE PROJECTION STAGE (126/126, BUILT CHAIN): EVERY REGISTERED FALSIFIER FAILS TO FIRE AND THE HEADLINE IS NOT "BOND IS WORSE" BUT "BOND IS INDISTINGUISHABLE FROM A MATCHED RANDOM DISPLACEMENT" -- BOND - CTRL-RAND(mean of 8) IS -0.0087 A AT 0.03x MDE WITH THE FOLD CI SPANNING ZERO, WHILE BOTH ARE ABOUT +0.73 A WORSE THAN PRODUCTION; SPAN +0.1220 (1.65x, 5/5) AND ISO +0.0737 (1.25x, 4/5) ARE ALSO WORSE; AND THE MECHANISM IS CONFIRMED WHILE THE HYPOTHESIS IS REFUTED -- MAKING THE CLOUD GEOMETRICALLY CONSISTENT GENUINELY IMPROVES THE FIT (RESIDUAL 0.8135 -> 0.6191 A) AND MAKES THE STRUCTURE WORSE ({now}, P)

Question (`s29/briefs/S29P.md`, wave-2 P1; prereg `s29/PREREG_S29_P.md` written 00:10 with six
addenda). Production's cloud is 3.0483 and its built chain 3.2126: the projection costs
+0.1643 A, the largest single-stage loss in the record and the only stage no S28 lane touched.
The cloud's mean adjacent CA-CA distance is 2.9614 against the ideal 3.80, so the least-squares
fit of a rigid-length chain to it is a BIASED fit. Does removing that inconsistency BEFORE the
projection change the built chain?

NOT S23 L6 (contract rule 10), as registered before any number: L6(d) lines 165-176 closes the
RMSD-optimal scale s\\* because it is read off the NATIVE and is a property of the (pool, native)
pair. g = 3.80 / (the cloud's own mean adjacent CA-CA distance) is read off the cloud and one
covalent constant; and the emitted chain's CA-CA distance is 3.80 BY CONSTRUCTION, so an input
rescale cannot rescale the output, only change which ideal chain is nearest. L6 stands untouched.

CHRONOLOGY, certified from git and pinned by prereg addendum 5 BEFORE these numbers existed,
because a "the bound predicted this" claim is worthless if the reader cannot check the order:
**00:10** F-P1/F-P2/F-P3 and the null prior registered (commit f721ab3b); **00:37** S29-L22
result 4 predicts BOND harmful BY MECHANISM from the separation profile; **00:43** lane T's
cosine bound S29-L23 is published. The mechanism prediction PRECEDES the bound by six minutes
and the falsifiers precede both by half an hour.

THE GATE, re-asserted inside the 126-target run and not merely on the probe: arm PROD reproduces
`s27/results/chain_rows.jsonl :: DIS` with **max |diff| 0.000e+00, identical on 126/126**. Both
sides of every contrast below share one code path and one input.

---

### 1. THE HEADLINE: BOND IS WHAT A RANDOM DISPLACEMENT OF THE SAME SIZE IS WORTH

The interesting statement is NOT that the bond-length correction is worse than production. It is
that it is worth exactly what an equally large displacement carrying NO information is worth --
and that is only sayable because the eight matched-magnitude draws survived the compute cut.
CTRL-RAND is a derangement of the realised g across targets: identical marginal distribution of
factors, zero per-target information, no target receiving its own.

{f("BOND|RAND")}

**0.03x MDE, the fold CI spans zero, 2/5 folds same sign, 60W/66L.** Against production both are
far worse and by the same amount:

{f("BOND")}

{f("CTRL-RAND(mean8)")}

So the whole of BOND's effect is its MAGNITUDE and none of it is its DIRECTION. The bond-length
rule identifies how far to move and carries no usable information about where.

### 2. THE OTHER ARMS

{f("SPAN")}

{f("ISO")}

{f("CTRL-GLOBAL")}

{f("CTRL-INV")}

Means (built chain): PROD **{M['PROD']:.4f}**, ISO {M['ISO']:.4f}, SPAN {M['SPAN']:.4f},
CTRL-INV {M['CTRL-INV']:.4f}, CTRL-GLOBAL {M['CTRL-GLOBAL']:.4f}, BOND {M['BOND']:.4f},
CTRL-RAND mean-of-8 {M['CTRL-RAND(mean8)']:.4f}.

**FALSIFIER VERDICT: F-P1, F-P2 and F-P3 ALL FAIL TO FIRE, and not narrowly.** Every arm is
worse than production or not measurable, and clause (d) -- beat CTRL-RAND's mean by >= 1.0x MDE --
fails at 0.03x. No arm clears its bar, so addendum 5's clause 3 (report a clearing arm as evidence
AGAINST lane T's bound) is not invoked. It remains the standing rule for anyone who reruns this.

**Where my own registered prior was WRONG, stated plainly.** I registered "null, 0.0 to 0.3x MDE"
for BOND, SPAN and ISO. The DIRECTION was right (no arm helps) but the MAGNITUDE was badly wrong
for BOND, which is +2.24x MDE WORSE, not null. The prior that was right is the later, sharper one:
S29-L22 result 4 predicted BOND harmful by mechanism at 00:37, and that is what happened. A null
prior and a harm prediction are different claims and I am not entitled to credit for the second
by having written the first.

### 3. THE MECHANISM IS CONFIRMED AND THE HYPOTHESIS IS REFUTED -- THE MOST USEFUL LINE IN THIS ENTRY

The lane's premise was that a geometrically inconsistent cloud is a biased fit. **It is, and
fixing it works, and the structure still gets worse.** The projection residual (native-free: the
emitted chain's CA-RMSD to the cloud it was fitted to):

    arm            fit_resid   chain RMSD
    PROD             0.8135       3.2126
    BOND             0.6191       3.9348     <- fits its own cloud 24% BETTER, 0.72 A worse
    SPAN             0.7280       3.3346
    ISO              0.7511       3.2863
    CTRL-INV         1.2455       3.3797
    CTRL-GLOBAL      0.8951       3.8705

Making the cloud's bonds ideal lets the ideal-geometry chain sit **24% closer to it**
(0.8135 -> 0.6191 A). The fit really was biased and the bias really is removable. The structure
is nonetheless 0.72 A worse, because the object the fit is now faithful to is a worse object: by
S29-L22 result 4 the distortion is a MONOTONE SHEAR in sequence separation (0.773 at |i-j| = 1,
crossing 1.00 near 8, 1.10 at 13), so setting the bond right inflates every long-range distance
that was already too long. The emitted Rg confirms it, against a native mean of 6.6009:
PROD 6.4820 (1.011x native), BOND **8.1353 (1.272x)**, CTRL-GLOBAL 8.2216, SPAN 6.9974,
ISO 6.8630, CTRL-INV 5.5133 (0.861x). **BOND does not contract or expand the answer toward the
native; it inflates it by 27%.** The contraction check was registered to catch a winner that was
merely a compact blob; it here catches a loser that is merely an inflated one.
The projection price per arm falls exactly as S28-L39 predicted it would for a de-contracted
cloud -- PROD {P['PROD']:+.4f}, BOND {P['BOND']:+.4f}, CTRL-GLOBAL {P['CTRL-GLOBAL']:+.4f} --
**and the price falling is not the price being saved.** It falls because the cloud got worse, not
because the chain got better. That settles the reading registered in advance in the prereg: the
+0.164 A is not a recoverable loss but the cost of the geometry constraint, and paying it earlier
costs more, not less.

### 4. THE BRANCH-FLIP FLOOR, MEASURED ON THIS PATH, AND WHY IT IS THE RIGHT REFERENCE

Arm FLOOR re-projects the SAME cloud perturbed at 1e-13 relative -- zero geometric change, a pure
redraw of the multi-start's branch choice:

    mean |diff| {fl['mean']:.4f} A   max {fl['max']:.4f} A (2LNG)   above 0.02 A on {fl['n_above_0p02']}/126   mean signed -0.0077

The tail is half an angstrom on one target in 126 from a perturbation at the last bit of the
input. Every deployable arm is therefore also contrasted against FLOOR rather than against
production's particular branch, and the conclusions do not move: BOND {C['BOND|FLOOR']['effect']:+.4f}
({C['BOND|FLOOR']['effect_over_mde']:+.2f}x), SPAN {C['SPAN|FLOOR']['effect']:+.4f}
({C['SPAN|FLOOR']['effect_over_mde']:+.2f}x), ISO {C['ISO|FLOOR']['effect']:+.4f}
({C['ISO|FLOOR']['effect_over_mde']:+.2f}x), CTRL-INV {C['CTRL-INV|FLOOR']['effect']:+.4f}
({C['CTRL-INV|FLOOR']['effect_over_mde']:+.2f}x, still NOT MEASURED). The signed FLOOR mean is
-0.0077 A, so production's own branch choice is very slightly WORSE than a random redraw of it --
inside its own noise, but the sign is a free reminder that the multi-start is not solved.

### 5. STRATA, CONCENTRATION, MULTIPLICITY

FAIL18 vs the 108, each through `ST.compare`, with a 2000-draw RANDOM-18 null so no stratum claim
rests on the 18 being special: BOND FAIL18 +1.1116 (0.74x, NOT MEASURED) against other108 +0.6573
(2.32x, 5/5, WORSE); the FAIL18 effect sits at **percentile 0.902** of random 18-subsets
(null p10/p50/p90 +0.3775/+0.7225/+1.1043). SPAN: +0.1984 / +0.1092, percentile 0.888.
ISO: +0.1285 / +0.0646, percentile 0.859. **All three are inside the random-18 band at the 90th
percentile or below, so there is NO regime claim here** -- the harm is carried by the 108, and
FAIL18 is merely noisier. This is the opposite of the pattern S28 kept finding and it is reported
because it did not go my way to report it.

Concentration: BOND's median (+0.2502) is far below its mean (+0.7222), which is the free early
warning -- but the valid test disagrees with the warning. Drop-top10 +0.8796 against a
uniform-effect null of p10/p50/p90 +0.7300/+0.8737/+1.0268 puts it at **percentile 0.522**: the
effect is a uniform harm with a heavy right tail, not a few targets. Median and mean are both
printed above, as the rule requires.

Multiplicity: **5 deployable endpoint contrasts** (BOND, SPAN, ISO, CTRL-GLOBAL, CTRL-INV;
CTRL-INV was moved into this set by addendum 3 because it is native-free, which raises the bar).
Max |effect|/MDE observed **3.52** against a joint sign-flip max-over-K null with p50 0.47,
p90 0.79, p95 0.90, **p = 0.000**. The effects are real; they are real HARM.

### 6. WHAT IS DEFERRED, AND HOW IT MUST BE READ WHEN IT RUNS

Deferred under the coordinator's one-slot decision, 1386 cells: the 9-point ORACLE s-grid, the
2x2 effective-lambda controls (CTRL-LAM, BOND-LAMFIX -- the prereg reads them only if BOND is
non-null, and BOND is harmful rather than null, so they are informative but not load-bearing),
and MS-OBJ/MS-MEAN/MS-ORACLE. `analyse` gates all of them mechanically and reports
grid-complete 0/126, so nothing here is quoted from a partial grid.
**Binding instruction for whoever runs it** (prereg addendum 6, from S29-L31): the grid's
per-target argmin is a Neyman-Scott INCIDENTAL PARAMETER -- one nuisance parameter per target,
not estimable from other targets' answers as a matter of theory rather than model capacity. It
may NOT be written up as "headroom" or as a target for a future predictor; it must be reported
with `ST.best_of_k_within`'s split-half transfer beside it as the measured price.

### 7. VERDICT

**REFUTED, and the projection stage is closed as a source of native-free accuracy.** Every
pre-registered falsifier fails. The stage's +0.164 A is real, it is concentrated where the cloud
is most distorted (S29-L22 result 2, Spearman +0.604), the distortion is genuinely removable at
the level of the fit (residual 0.8135 -> 0.6191), and removing it costs +0.72 A because the
distortion is a shear and a scalar is the wrong instrument for it. The bond-length rule is worth
precisely a matched random displacement (-0.0087 A, 0.03x MDE).
This is the FOURTH independent instance tonight of one pattern and is written as such rather than
as a discovery: a real, large, per-target structure with nothing transferable (S23 L3/L6's scale;
S29-L30's top-128 prefix at ORACLE 2.7605 vs a deployable +0.0079; lane O's rungs; this).
It agrees, through a different construction, with lane D's 21-field survey and with lane T's
bound S29-L23 -- and the agreement is worth something only because the prediction is timestamped
ahead of the bound.

Multiplicity: 5 endpoint comparisons here; the lane's running total is 5 built-chain endpoint
comparisons plus the 2 point-cloud diagnostics of S29-L22.
Artefacts: `s29/PREREG_S29_P.md` (6 addenda), `s29/s29_P_scale.py`, `tests/test_s29_P.py`
(20 pass), `s29/results/s29_P_rows_shard{{0..3}}.jsonl` (1896 rows, 126 targets x 15 arms),
`s29/results/s29_P_summary.json` (every block above verbatim in `text`), `s29_P_factors.json`,
`s29_P_probe.json`, `s29_P_contraction_profile.json`, `s29_P_sepprofile_cloud.json`,
`s29_P_oracle_cloud_sstar.json`, `s29_P_bandsel.json`; jobs `s26/jobs_done/s29P_factors.json`,
`s29P_probe6.json`, `s29P_prim_s{{0..3}}.json`; findings `s29/s29_P_FINDINGS.md` P1-P12.
This entry SUPERSEDES the coordinator's recomputation in the report's section 9.2c.
"""
    with open(LEDGER, "a", encoding="utf-8") as fh:
        fh.write(body)
    print("appended S29-L%d at %s" % (n, now))


if __name__ == "__main__":
    main()
