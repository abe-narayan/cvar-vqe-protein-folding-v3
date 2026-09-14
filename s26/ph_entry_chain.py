"""s26/ph_entry_chain.py -- builds the built-chain steric-reject ledger entry from the report log
(verbatim ST.fmt blocks) and appends it through ph_ledger. One-shot, PH lane."""
import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
log = open(os.path.join(HERE, "results", "ph_reject_report_both.log"), encoding="utf-8").read()


def block(label):
    i = log.find("  " + label)
    assert i >= 0, label
    j = log.find("VERDICT", i)
    j = log.find("\n", j)
    return log[i:j]


WANT = ["built_chain [all, n=126] R@1e4 minus anchor (negative = arm better)",
        "built_chain [all, n=126] R@1e4 minus RANDR@1e4 (matched control)",
        "built_chain [all, n=126] R@1e4 minus PERMR@1e4 (matched control)",
        "built_chain [all, n=126] S@1e4 minus anchor (negative = arm better)",
        "built_chain [all, n=126] S@1e4 minus RANDS@1e4 (matched control)",
        "built_chain [all, n=126] RANDR@1e4 minus anchor (negative = arm better)",
        "built_chain [all, n=126] RANDS@1e4 minus anchor (negative = arm better)",
        "built_chain [all, n=126] R@1e3 minus anchor (negative = arm better)",
        "built_chain [all, n=126] S@1e3 minus anchor (negative = arm better)",
        "built_chain [moved, n=116] R@1e4 minus anchor (negative = arm better)",
        "built_chain [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)"]
blocks = "\n".join("    " + block(w).replace("\n", "\n    ").rstrip() for w in WANT)

rows = []
for T in ("1e3", "1e4", "1e5", "1e6"):
    for arm in ("R", "S"):
        b = block(f"built_chain [all, n=126] {arm}@{T} minus anchor (negative = arm better)")
        eff = re.search(r"effect ([+-]\d\.\d+)", b).group(1)
        med = re.search(r"median ([+-]\d\.\d+)", b).group(1)
        x = re.search(r"effect/MDE ([+-]\d\.\d+)", b).group(1)
        ci = re.search(r"fold CI95 (\[[^\]]+\])", b).group(1)
        fs = re.search(r"folds same sign (\d)/5", b).group(1)
        wl = re.search(r"(\d+W/\d+L/\d+T)", b).group(1)
        v = re.search(r"VERDICT: (.+)", b).group(1).split(" [")[0].split(" (")[0]
        rows.append(f"    {T:<4} {arm}  {eff:>8} (median {med:>8}, {x:>5}x MDE)  fold {ci:<20} {fs}/5  {wl:<12} {v}")
table = "\n".join(rows)

# the anchor deviation from the stored production chain (addendum 2 of the prereg)
r = json.load(open(os.path.join(HERE, "results", "ph_reject_chain.json"), encoding="utf-8"))["rows"]
dev = np.array([x["anchor_vs_prod"] for x in r]); da = np.array([x["anchor"] - x["prod_rmsd_arm"] for x in r])
worst = sorted([(float(v), x["pdb"]) for v, x in zip(dev, r)])[-4:]
anchor_mean = float(np.mean([x["anchor"] for x in r])); prod_mean = float(np.mean([x["prod_rmsd_arm"] for x in r]))

ENTRY = f"""`s26/ph_reject.py chain` then `report`, `s26/results/ph_reject_chain.json` (complete 126/126),
`s26/results/ph_reject_report.json` (both bases), job `s26/jobs_done/ph_reject_chain.json` (exit
0, 11,976 s = 3.3 h, 22.8 projections per target, peak RSS 0.09 GB). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5, addendum 1 (fallback and moved subset) and
addendum 2 (the re-projected anchor). This is the falsifier's last leg and the cross-basis
replication of L43.

BASIS: BUILT CHAIN on both sides (`rmsd_arm`): every retained set is coordinate-averaged and
then projected through `s12.instrument.project` (ramah at 0.3, multi-start, exact gradient,
the production call) and the chain is scored ORACLE against `nat_ca`. The ANCHOR is the
shipped top-75 RE-PROJECTED through the same call, not the stored production chain, so all arms
share one instrument. Its deviation from the stored chain, reported as registered: mean
{anchor_mean:.4f} A against the production {prod_mean:.4f} (the "leaderboard rebuild 3.2126" of the
state brief and L57: the production path round-trips the cloud through float32 before
projecting, `Config.reference_precision`, and this instrument does not); per-target RMSD
difference mean {da.mean():+.4f}, max |{np.abs(da).max():.3f}|; CA deviation between the two chains mean
{dev.mean():.3f} A, above 0.05 A on {int((dev > 0.05).sum())} targets and above 0.5 A on {int((dev > 0.5).sum())}
({", ".join(f"{p} {v:.2f}" for v, p in worst)}), which is the projection's own branch
degeneracy (`core/project.py` docstring: a 1e-13 input difference can route L-BFGS-B into the
other torsion branch, up to 1.6 A). Arms, controls and fallbacks as in L43 (R = reject e_amber
> T and refill to 75, the PRIMARY reading; S = reject, no refill; RANDR/RANDS same count at
random; PERMR/PERMS the same threshold on a permuted energy); the built-chain controls at 1e4
carry 4 draws each (registered in section 3), the other thresholds carry R and S only.
Eleven `ST.fmt` blocks verbatim:

{blocks}

All four thresholds, `[all, n=126]`, arm minus anchor, built chain:

{table}

Threshold sweep (order statistic): R oracle minimum -0.3442, split-half transfer -0.1473 (43%%),
k_eff 3.86; S -0.1345, transfer -0.0535 (40%%): the mildest rung transfers, not a gain (R@1e6
+0.092 at 0.65x MDE, S@1e6 +0.051 at 0.76x, both NOT MEASURED).

READING. (1) The built chain REPLICATES the point cloud in sign on the primary reading: R@1e4
is +0.248 A worse than the re-projected shipped top-75 (fold CI [+0.120, +0.348], 49W/67L/10T)
and +0.157 worse than rejecting the same count at random with a judgment-free refill (fold CI
[+0.048, +0.281]); L43 had +0.228 and +0.167 on the point cloud. Both are Type-M-zone
magnitudes (1.17x and 1.15x MDE), as on the point cloud, and both are tail-carried (median
+0.002 and +0.009 against means +0.248 and +0.157; p90 +1.39; worst +4.20 on 8T61), so the
Adversary's L54 caveats 1 and 2 transfer unchanged to this basis: near zero on the median
target, catastrophic on the minority whose pool has no survivor or whose refill reaches deep.
One difference from L43 that is reported rather than smoothed: fold 1 has the opposite sign on
R (-0.002 against +0.267 / +0.226 / +0.421 / +0.310 on the other four), so R@1e4 is 4/5 folds
on the built chain where it was 5/5 on the point cloud; the fold CI still excludes zero. S@1e4
is +0.104 worse (5/5 folds, 1.11x MDE, Type-M) against +0.108 on the point cloud. (2) The
dose is again monotone: 1e3 +0.563 (1.77x MDE, 5/5 folds, MEASURED), 1e4 +0.248, 1e5 +0.119
(0.78x), 1e6 +0.092 (0.65x); the limit is the anchor. (3) The built chain is noisier than the
point cloud by construction (the projection adds its own branch noise on every arm), so the
matched-control contrasts that were Type-M on the point cloud are UNDERPOWERED here: S vs
RANDS +0.055 at 0.73x MDE (fold CI [-0.0005, +0.106]), R vs PERMR +0.102 at 0.88x, S vs PERMS
+0.050 at 0.59x; RANDR vs anchor +0.091 at 0.75x, RANDS vs anchor +0.049 at 0.99x. The L54
caveat 3 stands: the refill-cost / choice-cost split is a point estimate on both bases. (4)
The moved subset (secondary, n = 116 / 113) agrees: R +0.270 WORSE (Type-M), R vs RANDR +0.169
WORSE (Type-M), S +0.116 WORSE (Type-M).

POWER. Built-chain SEs 0.03 to 0.11 A, MDEs 0.09 to 0.32; at 1e5 and 1e6 R sits at 0.78x and
0.65x its MDE (SE 0.055 and 0.050), so "the mildest threshold is null" is UNDERPOWERED for a
gain below about 0.14 A on this basis and MEASURED against any harm above it. Nothing is
positive; the replication of the harmful direction is this entry (cross-basis, same sign, same
monotone dose, same tail structure).

DISPOSITION. The falsifier of `s26/PREREG_amber_reject.md` section 4 required R or S at 1e4 to
BEAT the anchor and its matched control on the built chain; both are WORSE than the anchor
with fold CIs excluding zero (R 4/5 folds, S 5/5), and R is worse than its matched control.
AMBER as a steric reject filter at a physical threshold, with or without refill, is CLOSED on
both bases in the two forms the record had not measured (physics-set count, refill), with a
harmful sign at 1e3 and 1e4 and an underpowered null at 1e5 and 1e6. The functional lever's
filter form (s19 Q3, s24 D1-C) keeps its closure and gains a sixth and seventh instrument. Why,
from L23: the members the threshold condemns are condemned by the builder's side-chain
placement (96.8%), not by their backbones; removing them removes compact members whose error
cancelled in the average (S23 L5). The remaining physics question on this pool is
`IDEA_rotamer_relief` part B (does the singularity move when the side chains are relieved),
which runs next under the tournament.
"""
out = os.path.join(HERE, "results", "_entry_chain.md")
open(out, "w", encoding="utf-8", newline="\n").write(ENTRY)
title = ("STERIC REJECT, BUILT CHAIN: THE POINT-CLOUD RESULT REPLICATES ACROSS BASES. R@1e4 IS +0.248 A WORSE THAN "
         "THE RE-PROJECTED SHIPPED TOP-75 (fold CI [+0.12, +0.35], 4/5 FOLDS, TYPE-M) AND +0.157 WORSE THAN THE SAME "
         "COUNT REJECTED AT RANDOM; DOSE MONOTONE; THE FILTER IS CLOSED ON BOTH BASES (2026-09-13, PH)")
p = subprocess.run([sys.executable, os.path.join(HERE, "ph_ledger.py"), "--title", title, "--body", out],
                   capture_output=True, text=True)
print(p.stdout.strip(), p.stderr.strip())
os.remove(out)
print(table)
