"""s26/ph_entry_cloud.py -- builds the point-cloud steric-reject ledger entry from the report log,
verbatim ST.fmt blocks, and appends it through ph_ledger (collision-safe). One-shot, PH lane."""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "results", "ph_reject_report_cloud.log")
log = open(LOG, encoding="utf-8").read()


def block(label):
    i = log.find("  " + label)
    assert i >= 0, label
    j = log.find("VERDICT", i)
    j = log.find("\n", j)
    return log[i:j]


WANT = ["point_cloud [all, n=126] R@1e4 minus anchor (negative = arm better)",
        "point_cloud [all, n=126] R@1e4 minus RANDR@1e4 (matched control)",
        "point_cloud [all, n=126] R@1e4 minus PERMR@1e4 (matched control)",
        "point_cloud [all, n=126] S@1e4 minus anchor (negative = arm better)",
        "point_cloud [all, n=126] S@1e4 minus RANDS@1e4 (matched control)",
        "point_cloud [all, n=126] S@1e4 minus PERMS@1e4 (matched control)",
        "point_cloud [all, n=126] RANDR@1e4 minus anchor (negative = arm better)",
        "point_cloud [all, n=126] RANDS@1e4 minus anchor (negative = arm better)",
        "point_cloud [moved, n=116] R@1e4 minus anchor (negative = arm better)",
        "point_cloud [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)",
        "point_cloud [moved, n=113] S@1e4 minus anchor (negative = arm better)",
        "point_cloud [moved, n=113] S@1e4 minus RANDS@1e4 (matched control)"]
blocks = "\n".join("    " + block(w).replace("\n", "\n    ").rstrip() for w in WANT)

rows = []
for T in ("1e3", "1e4", "1e5", "1e6"):
    for arm in ("R", "S"):
        b = block(f"point_cloud [all, n=126] {arm}@{T} minus anchor (negative = arm better)")
        eff = re.search(r"effect ([+-]\d\.\d+)", b).group(1)
        mde = re.search(r"MDE (\d\.\d+)", b).group(1)
        x = re.search(r"effect/MDE ([+-]\d\.\d+)", b).group(1)
        ci = re.search(r"fold CI95 (\[[^\]]+\])", b).group(1)
        wl = re.search(r"(\d+W/\d+L/\d+T)", b).group(1)
        v = re.search(r"VERDICT: (.+)", b).group(1).split(" [")[0].split(" (")[0]
        ctrl = "RANDR" if arm == "R" else "RANDS"
        bc = block(f"point_cloud [all, n=126] {arm}@{T} minus {ctrl}@{T} (matched control)")
        effc = re.search(r"effect ([+-]\d\.\d+)", bc).group(1)
        xc = re.search(r"effect/MDE ([+-]\d\.\d+)", bc).group(1)
        cic = re.search(r"fold CI95 (\[[^\]]+\])", bc).group(1)
        vc = re.search(r"VERDICT: (.+)", bc).group(1).split(" [")[0].split(" (")[0]
        rows.append(f"    {T:<4} {arm}  vs anchor {eff:>8} ({x:>5}x MDE {mde})  fold {ci:<20} {wl:<12} {v:<13}"
                    f" | vs {ctrl} {effc:>8} ({xc:>5}x)  fold {cic:<20} {vc}")
table = "\n".join(rows)

ENTRY = """`s26/ph_reject.py cloud` then `report`, `s26/results/ph_reject_cloud.json` (complete 126/126),
`s26/results/ph_reject_report.json`, log `s26/results/ph_reject_report_cloud.log`; job
`s26/jobs_done/ph_reject_cloud2.json` (exit 0, 170 s for the last 10 cells, peak RSS 0.041 GB;
the first 116 cells came from `ph_reject_cloud`, killed with the governor by the host at about
10:10, L40; the per-target cells made the restart lossless). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5 and addendum 1; the moved-subset secondary and the
R-first reading were fixed by the coordinator at 09:10 (findings section 2.3) before any RMSD
was read.

BASIS: POINT CLOUD on both sides (`I.coordinate_average` of the retained windows, scored ORACLE
against `nat_ca`; the anchor reproduces the production `rmsd_avg` 3.0483 to 1e-6 on every
target). The built-chain twin (`rmsd_arm` basis) is job `ph_reject_chain`, running. Arms: R =
reject e_amber > T from the shipped top-75 and refill from the next-ranked survivors to m = 75
(the PRIMARY reading; it empties only when the whole pool has no survivor: 8 targets at 1e4,
which fall back to the anchor and count as ties); S = reject without refill (empties 11 sets at
1e4, same fallback). Controls matched in the operator's space: RANDR / RANDS reject the SAME
COUNT at random (16 draws, mean); PERMR / PERMS apply the same threshold to a permuted energy
vector (16 permutations). Twelve `ST.fmt` blocks verbatim, primary threshold 1e4:

%(blocks)s

All four thresholds, `[all, n=126]`, arm minus anchor and arm minus its matched random control:

%(table)s

Threshold sweep (an order statistic): `ST.best_of_k_within` over the four thresholds, R: oracle
per-target minimum -0.3368 A, split-half transfer -0.1470 (44%%), k_eff 3.74; S: -0.1177,
transfer -0.0556 (47%%), k_eff 3.72. The transfer is negative because the sweep contains 1e6,
which rejects least and so damages least: "the least harmful threshold transfers" is not a gain.
R@1e6 is +0.066 against the anchor (0.51x MDE, NOT MEASURED) and S@1e6 +0.050 (0.97x).

READING. (1) The falsifier did not clear; it fired the other way. At the primary threshold the
refill arm R is WORSE than the shipped top-75 by +0.228 A [fold +0.147, +0.323], 49W/67L/10T,
5/5 folds (1.17x MDE, Type-M zone: the size is an upper bound, the sign is measured), and WORSE
than rejecting the same number at random and refilling with no energy judgment by +0.167
[+0.070, +0.279], 5/5 folds (1.26x MDE, Type-M). The shrink arm S is worse than the anchor by
+0.108 [+0.065, +0.153], 1.34x MDE, measured; against its random twin +0.071 at exactly 1.00x
MDE, not measured. Restricting to the targets the operator actually moved (the secondary,
n = 116 / 113) changes nothing: R +0.247 and S +0.120, both WORSE; R vs RANDR +0.176, WORSE.
(2) The dose is monotone: 1e3 (54.5 of 75 rejected) +0.532 A; 1e4 (40.1) +0.228; 1e5 (30.3)
+0.093 (0.66x MDE); 1e6 (23.3) +0.066 (0.51x). The less the reject does, the less it costs, and
its limit is the anchor: S16 L27's k-ladder shape and S23 L5's mechanism (the members a physics
score removes carry error that cancels in the average) through a third operator. (3) The
matched controls say where the harm comes from. Rejecting the same COUNT at random and
refilling (RANDR) costs +0.061 at 1e4 (0.57x MDE, not measured), so refilling from ranks 76 to
147 is nearly free and the energy's CHOICE of what to reject costs the other +0.167. Rejecting
the same count at random without refill (RANDS) costs +0.037 (1.47x MDE, measured): shrinking
a 75-set to 35 at random costs 0.04 A and the energy's choice of which 40 costs 0.07 more. The
permuted-energy controls sit between random and real (PERMR +0.122, PERMS +0.033) because a
permuted single point rejects the same energies from other candidates. (4) Why: the census
(L23) found 96.8%% of the members the threshold condemns are condemned by a side-chain contact
placed by the deterministic builder, and the retained set at 1e4 is +0.254 A more expanded in
Rg than the anchor; the reject keeps expanded members and discards compact ones. S25's "AMBER
selects expansion" (+1.10 A Rg on its own top-75) reproduced as a filter.

POWER. Every arm-vs-anchor contrast at 1e3 and 1e4 is above its MDE with the fold CI excluding
zero and 5/5 folds. At 1e5 and 1e6 R is 0.66x and 0.51x its MDE (SE 0.050 and 0.047, so a gain
of 0.13 A or more would have been seen): "the mildest threshold is null" is UNDERPOWERED for a
gain below 0.13 A and MEASURED against any harm above it. Nothing here is positive, so no seed
replication is due; the cross-basis replication is the built-chain run, and the direction must
hold there or the point-cloud result is a basis artefact (S16 showed the two bases can disagree
about a repair). DISPOSITION, pending the built chain: the steric reject at a physical
threshold, with refill, is closed on the point cloud in the two forms the record had not
measured (physics-set count, refill), with a sign, on all 126 and on the moved subset; the
filter form of the functional lever (s24 D1-C, s19 Q3) stays closed and gains its sixth
instrument.
""" % {"blocks": blocks, "table": table}

out = os.path.join(HERE, "results", "_entry_cloud.md")
open(out, "w", encoding="utf-8", newline="\n").write(ENTRY)
title = ("STERIC REJECT, POINT CLOUD: THE FALSIFIER FIRED THE OTHER WAY. AT 1e4 THE REFILL ARM IS "
         "+0.228 A WORSE THAN THE SHIPPED TOP-75 AND +0.167 WORSE THAN REJECTING THE SAME COUNT AT "
         "RANDOM (5/5 FOLDS); THE DOSE IS MONOTONE AND ITS LIMIT IS THE ANCHOR (2026-09-13, PH)")
r = subprocess.run([sys.executable, os.path.join(HERE, "ph_ledger.py"), "--title", title, "--body", out],
                   capture_output=True, text=True)
print(r.stdout.strip(), r.stderr.strip())
os.remove(out)
print(table)
