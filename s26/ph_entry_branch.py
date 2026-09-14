"""s26/ph_entry_branch.py -- ledger entry for branch_select from the report log and artefact. One-shot."""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
log = open(os.path.join(HERE, "logs", "ph_branch_report.log"), encoding="utf-8").read()


def block(label):
    i = log.find("  " + label)
    assert i >= 0, label
    j = log.find("VERDICT", i)
    j = log.find("\n", j)
    return log[i:j]


WANT = ["[all, n=126] pick_e1 minus production (built chain)",
        "[all, n=126] pick_e1 minus random pick (built chain)",
        "[all, n=126] random pick minus production (built chain)",
        "[all, n=126] pick_e0 minus production (built chain)",
        "[all, n=126] pick_e1 minus production (relaxed chain)"]
blocks = "\n".join("    " + block(w).replace("\n", "\n    ").rstrip() for w in WANT)
d = json.load(open(os.path.join(HERE, "results", "ph_branch_report.json"), encoding="utf-8"))
r, rows = d["report"], d["rows"]
o = r["oracle_over_5"]
means = {k: float(np.mean([x[k] for x in rows])) for k in ("prod", "pick_e1", "pick_e0", "random", "oracle_min")}
same = sum(1 for x in rows if abs(x["pick_e1"] - x["prod"]) < 1e-9)
nconv = sum(1 for x in rows if x["n_converged"] < 5)

ENTRY = f"""`s26/ph_branch.py solutions | relax | report`, artefacts `s26/results/ph_branch_solutions.json`
(126/126, G1 exact 0.0 on every target), `ph_branch_relax.json` (126/126; jobs `ph_branch_relax`
41 cells + `ph_branch_relax2` 85 cells, AMBER, peak RSS 0.27 GB probe / 0.3 GB run, about 65 s
per target), `ph_branch_report.json`. Pre-registered in `s26/PREREG_branch_select.md` (sections
1 to 7, addendum 1 with the gates); tournament item 5 (L51). Per-target cells under
`s26/results/ph_branch_{{sol_,}}cells/`.

THE OPERATOR, native-free. The production projection chooses among FIVE solutions at its
lam = 0.3 rung (the warm start from lam = 0 and the four generic starts alpha / beta / PPII /
extended) by the lowest objective; the choice was reconstructed and equals `I.project` to
0.0 A on 126/126. The arm relaxes all five built chains with the production operator
(`refine_coords(k=10, steps=0)`), reads the converged energy with the restraint off (`e1`) and
emits the BUILT chain of the lowest converged `e1` (ties averaged by `ST.argmin_tied`; 0 ties
occurred; on {nconv} targets fewer than five converged and the pick is among the converged). All 126
targets have at least two distinct solutions (4.77 on average), so the effective n is 126 and
the moved subset equals the whole. BASIS: BUILT CHAIN (`rmsd_arm`) for the primary; the relaxed
chain of the same pick as a secondary. Controls in the operator's space: the production choice
(anchor, the re-projected 3.2126); the exact expectation of a uniformly random pick among the
five; the raw single point `e0` as the declared secondary picker. `ST.fmt` verbatim:

{blocks}

    means, built chain:  production {means['prod']:.4f}   e1 pick {means['pick_e1']:.4f}   e0 pick {means['pick_e0']:.4f}   random pick {means['random']:.4f}   ORACLE min over five {means['oracle_min']:.4f}
    e1 pick equals the production choice on {same} of 126 targets (the 42 exact ties in the first block)
    ORACLE DIAGNOSTIC: the e1 pick is the per-target best on {100 * r['frac_e1_pick_is_oracle']:.1f}% of targets, the production choice on {100 * r['frac_prod_is_oracle']:.1f}%
    ORACLE min over the five (an order statistic): {o['observed_gain']:+.4f}; valid across-target null {o['null_across_targets']:+.4f} (share {o['share_accounted']:.2f});
      split-half transfer {o['split_half']:+.4f} ({100 * o['split_half_frac']:.0f}%); k_eff {o['k_eff']:.2f}; argmin counts over the five columns {[int(c) for c in o['argmin_counts']]}

READING. (1) The falsifier (`PREREG_branch_select.md` section 4) required the relaxed-energy
pick to BEAT the production choice past its MDE with the fold CI excluding zero. It does not:
+0.0055 A, 0.18x MDE, 39W/45L/42T, the fold CI [+0.0004, +0.0106] on the wrong side of zero
anyway. NOT MEASURED, and at SE 0.011 an improvement of 0.03 A or more would have been seen:
underpowered below that, null above it. The idea is CLOSED as an accuracy step. (2) What the
energy DOES do, and it is the interesting half: the converged relaxed energy beats a uniformly
random choice among the same five solutions by -0.102 A [fold -0.127, -0.075], 91W/35L, 5/5
folds, 2.29x MDE, MEASURED. The production objective (CA-RMSD to the cloud plus 0.3 x ramah)
beats the random pick by 0.108 (2.65x MDE). So the converged all-atom energy carries the SAME
discriminating power among the projection's branches as the 2D torsion prior, and adds nothing
to it: on 42 targets it picks the identical solution, on the rest it trades one near-equivalent
branch for another (median difference 0.000). This is the first place in the record where an
AMBER quantity ranks a real discrete choice as well as the structural objective does; it is
also the first place where that choice is among only five candidates that all came from the
same cloud. (3) The raw single point `e0` is WORSE than the objective by +0.092 (1.10x MDE,
Type-M, 39W/58L/29T): the unrelaxed energy, which L23 showed is the builder's side-chain clash,
picks the wrong branch where the relaxed energy does not. Relaxation is what makes the energy a
usable discriminator here, the same operator fact as S20 L6 ("the relaxation is what makes the
AMBER objective defined"). (4) On the relaxed-chain basis the e1 pick is a null against the
production choice relaxed (-0.004, 0.09x MDE). (5) The ORACLE minimum over the five solutions is
-0.19 A below the production choice, but the valid best-of-5 null accounts for 130% of it and
the split-half transfer is 56% (-0.107): the branch degeneracy is real (a per-target-consistent
column exists) and it is worth about 0.1 A to a perfect chooser, which neither the objective
nor the energy is (each finds the per-target best on about 30% of targets, chance being 20%).

POWER. n = 126 with 42 exact ties; SE of the primary 0.011, MDE 0.031: a gain of 0.03 A would
have cleared. No positive result on the falsifier, so no seed replication is due (the arm has
no randomness; the relaxation is deterministic at threads = 1). The measured negative-control
contrast (energy beats random by 0.10) is not a proposal and needs no replication to be quoted
as a diagnostic. DISPOSITION: branch_select CLOSED as an accuracy lever; the record gains one
measured sentence for the report: "the converged force-field energy discriminates among the
projection's own branches exactly as well as the torsion prior already does, and the
unrelaxed energy does not."
"""
out = os.path.join(HERE, "results", "_entry_branch.md")
open(out, "w", encoding="utf-8", newline="\n").write(ENTRY)
title = ("BRANCH SELECT (tournament 5): THE RELAXED ENERGY PICKS THE PROJECTION BRANCH AS WELL AS THE OBJECTIVE "
         "AND NO BETTER (+0.0055 vs production, 0.18x MDE, 42 ties; -0.102 vs a random branch, 2.29x, 5/5); "
         "THE RAW SINGLE POINT PICKS WORSE (+0.092); CLOSED AS AN ACCURACY STEP (2026-09-13, PH)")
p = subprocess.run([sys.executable, os.path.join(HERE, "ph_ledger.py"), "--title", title, "--body", out],
                   capture_output=True, text=True)
print(p.stdout.strip(), p.stderr.strip())
os.remove(out)
